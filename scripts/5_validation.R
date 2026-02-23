###############################################################################
# 5_validation.R
#
# Purpose: Validate occupancy model predictions against the Spanish
#          Biodiversity Atlas using spatial block cross-validation.
#          Computes AUC, TSS, RMSE, Spearman correlation, and Moran's I.
#
# Inputs:  data/processed/{sp}/occ_{sp}_prediction.csv   (from step 4)
#          data/raw/validation/atlas_biodiversidad/aves_spain.shp
#
# Outputs: results/{sp}_validation_cv.csv
#          results/{sp}_validation_summary.txt
#          figs/{sp}_validation_calibration.png
#          figs/{sp}_validation_residuals_map.png
#          figs/{sp}_validation_maps.png
#
# Replaces: 5_validation_otitar.R (now parameterised for all species)
#
# Note:    The atlas shapefile must contain columns with species codes
#          (OTITAR, PTEALC, PTEORI, TETTET) as 0/1 presence data.
###############################################################################

# -- Load packages --
library(here)
library(sf)
library(terra)
library(dplyr)
library(ggplot2)
library(pROC)
library(Metrics)
library(spdep)
library(blockCV)

set.seed(123)  # Reproducible spatial CV

# -- Species to validate --
species_codes <- c("otitar", "ptealc", "pteori", "tettet")

# Map species code to atlas column name (uppercase)
atlas_col_map <- c(
  otitar = "OTITAR",
  ptealc = "PTEALC",
  pteori = "PTEORI",
  tettet = "TETTET"
)

# -- Load atlas (shared across species) --
atlas_path <- here("data", "raw", "validation", "atlas_biodiversidad", "aves_spain.shp")
if (!file.exists(atlas_path)) {
  stop("Validation shapefile not found: ", atlas_path,
       "\nSee data-raw/get_data.R for instructions.")
}
atlas_full <- st_read(atlas_path, quiet = TRUE)
message("Atlas loaded: ", nrow(atlas_full), " polygons")


###############################################################################
# MAIN LOOP
###############################################################################

for (sp in species_codes) {

  message("\n=== Validation: ", sp, " ===")

  atlas_col <- atlas_col_map[sp]
  if (!atlas_col %in% names(atlas_full)) {
    warning("  Column '", atlas_col, "' not found in atlas. Skipping ", sp)
    next
  }

  # -- Load predictions --
  pred_path <- here("data", "processed", sp,
                     paste0("occ_", sp, "_prediction.csv"))
  if (!file.exists(pred_path)) {
    warning("  Prediction file not found: ", pred_path, ". Skipping ", sp)
    next
  }

  pred_occ <- read.csv(pred_path)
  names(pred_occ) <- tolower(names(pred_occ))
  stopifnot(all(c("longitude", "latitude", "occ_prob") %in% names(pred_occ)))

  # -- Rasterize predictions --
  r_pred <- terra::rast(
    data.frame(x = pred_occ$longitude, y = pred_occ$latitude,
               occ_prob = pred_occ$occ_prob),
    type = "xyz"
  )
  crs(r_pred) <- "+proj=longlat"

  # Reproject to atlas CRS
  r_pred_proj <- project(r_pred, crs(atlas_full))

  # Extract median probability per atlas polygon
  atlas <- atlas_full
  ext_prob <- terra::extract(r_pred_proj, atlas, fun = median, na.rm = TRUE)
  atlas$pred_prob_mean <- ext_prob$occ_prob

  atlas_valid <- atlas %>% filter(!is.na(pred_prob_mean))
  message("  Valid atlas polygons: ", nrow(atlas_valid))

  if (nrow(atlas_valid) < 20) {
    warning("  Too few valid polygons for ", sp, ". Skipping.")
    next
  }

  ##############################################################################
  # SPATIAL BLOCK CROSS-VALIDATION
  ##############################################################################
  sf_cv <- atlas_valid

  # cv_spatial replaces the deprecated spatialBlock() in blockCV >= 3.0
  sb <- cv_spatial(
    x = sf_cv,
    column = atlas_col,
    size = 50000,
    k = 5,
    selection = "random",
    iteration = 100,
    progress = FALSE,
    plot = FALSE
  )

  sf_cv$fold <- sb$folds_ids
  results_cv <- list()

  for (i in 1:max(sf_cv$fold)) {
    train <- sf_cv[sf_cv$fold != i, ]
    test  <- sf_cv[sf_cv$fold == i, ]

    # Threshold from training ROC
    roc_train <- pROC::roc(train[[atlas_col]], train$pred_prob_mean, quiet = TRUE)
    th <- pROC::coords(roc_train, "best", ret = "threshold") %>% as.numeric()

    # Binary predictions for test fold
    test$bin <- ifelse(test$pred_prob_mean >= th, 1, 0)

    # Metrics
    auc_i  <- pROC::auc(pROC::roc(test[[atlas_col]], test$pred_prob_mean,
                                     quiet = TRUE))
    rmse_i <- Metrics::rmse(test[[atlas_col]], test$pred_prob_mean)

    cm <- table(factor(test[[atlas_col]], levels = c(0, 1)),
                factor(test$bin, levels = c(0, 1)))
    TP <- cm["1", "1"]; TN <- cm["0", "0"]
    FP <- cm["0", "1"]; FN <- cm["1", "0"]

    sens_i <- TP / (TP + FN)
    spec_i <- TN / (TN + FP)
    TSS_i  <- sens_i + spec_i - 1

    results_cv[[i]] <- data.frame(
      fold = i, threshold = th, AUC = auc_i, RMSE = rmse_i,
      Sens = sens_i, Spec = spec_i, TSS = TSS_i
    )
  }

  df_cv <- do.call(rbind, results_cv)
  write.csv(df_cv, here("results", paste0(sp, "_validation_cv.csv")),
            row.names = FALSE)

  # Optimal threshold
  opt_th <- mean(df_cv$threshold)

  ##############################################################################
  # APPLY OPTIMAL THRESHOLD
  ##############################################################################
  r_pred_bin <- r_pred >= opt_th
  names(r_pred_bin) <- "occ_bin"
  r_pred_bin_proj <- project(r_pred_bin, crs(atlas_full), method = "bilinear")
  ext_bin <- terra::extract(r_pred_bin_proj, atlas, fun = max, na.rm = TRUE)
  atlas$pred_bin_final <- ifelse(
    is.na(atlas$pred_prob_mean), NA,
    ifelse(ext_bin$occ_bin > 0, 1, 0)
  )

  ##############################################################################
  # CONTINUOUS VALIDATION
  ##############################################################################

  # --- Point-biserial correlation on INDIVIDUAL polygons (honest metric) ---
  # NOTE: Spearman/R² must be computed on raw polygon-level data, NOT on
  # aggregated decile bins. Computing on deciles (n=10 points) artificially
  # inflates correlation to ~1.0 because binning removes within-group variance.
  # See docs/validation_diagnosis.md for full explanation.

  spearman_raw <- cor.test(atlas_valid$pred_prob_mean,
                           as.numeric(atlas_valid[[atlas_col]]),
                           method = "spearman")

  # Pseudo-R² from logistic regression (appropriate for binary response)
  glm_calib <- glm(reformulate("pred_prob_mean", response = atlas_col),
                    data = atlas_valid, family = binomial)
  # McFadden's pseudo-R²
  null_dev <- glm_calib$null.deviance
  res_dev  <- glm_calib$deviance
  pseudo_r2 <- 1 - (res_dev / null_dev)

  # --- Per-fold Spearman (spatially honest) ---
  spearman_folds <- sapply(1:max(sf_cv$fold), function(i) {
    test_fold <- sf_cv[sf_cv$fold == i, ]
    suppressWarnings(
      cor(test_fold$pred_prob_mean,
          as.numeric(test_fold[[atlas_col]]),
          method = "spearman", use = "complete.obs")
    )
  })
  spearman_cv_mean <- mean(spearman_folds, na.rm = TRUE)
  spearman_cv_sd   <- sd(spearman_folds, na.rm = TRUE)

  # --- Calibration plot (decile bins — for VISUALISATION only) ---
  atlas_deciles <- atlas_valid %>%
    mutate(prob_bin = ntile(pred_prob_mean, 10)) %>%
    group_by(prob_bin) %>%
    summarise(mean_prob = mean(pred_prob_mean),
              obs_prev  = mean(as.numeric(.data[[atlas_col]])),
              n = n(), .groups = "drop")

  p_calib <- ggplot(atlas_deciles, aes(x = mean_prob, y = obs_prev)) +
    geom_point(aes(size = n)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_smooth(method = "lm", se = TRUE, color = "blue") +
    labs(x = "Mean predicted probability (decile bins)",
         y = "Observed prevalence",
         title = paste0(sp, ": Calibration plot"),
         subtitle = paste0("Spearman \u03c1 (raw polygons) = ",
                           round(spearman_raw$estimate, 3),
                           " | Pseudo-R\u00b2 = ", round(pseudo_r2, 3))) +
    theme_minimal()
  ggsave(here("figs", paste0(sp, "_validation_calibration.png")),
         p_calib, width = 7, height = 6)

  ##############################################################################
  # SPATIAL RESIDUALS + MORAN'S I
  ##############################################################################
  atlas_resid <- atlas_valid %>%
    mutate(resid = pred_prob_mean - as.numeric(.data[[atlas_col]]))

  p_resid <- ggplot(atlas_resid) +
    geom_sf(aes(fill = resid), color = NA) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 0) +
    ggtitle(paste0(sp, ": Residuals (Predicted - Observed)")) +
    theme_minimal()
  ggsave(here("figs", paste0(sp, "_validation_residuals_map.png")),
         p_resid, width = 8, height = 6)

  sf_valid <- atlas_resid %>% filter(!is.na(resid))
  sp_valid <- sf::as_Spatial(sf_valid)
  nb <- poly2nb(sp_valid, queen = TRUE)
  lw <- nb2listw(nb, style = "W", zero.policy = TRUE)
  moran_res <- moran.test(sf_valid$resid, lw, zero.policy = TRUE)

  ##############################################################################
  # PREDICTION vs ATLAS MAPS
  ##############################################################################
  map_pred <- ggplot(atlas) +
    geom_sf(aes(fill = factor(pred_bin_final))) +
    scale_fill_manual(values = c("white", "darkgreen"),
                      name = "Prediction",
                      labels = c("Absence", "Presence"),
                      na.value = "grey90") +
    ggtitle(paste0(sp, ": Model prediction")) + theme_bw()

  map_obs <- ggplot(atlas) +
    geom_sf(aes(fill = factor(.data[[atlas_col]]))) +
    scale_fill_manual(values = c("grey80", "red"),
                      name = "Atlas",
                      labels = c("Absence", "Presence")) +
    ggtitle(paste0(sp, ": Biodiversity Atlas")) + theme_bw()

  ggsave(here("figs", paste0(sp, "_validation_maps.png")),
         gridExtra::grid.arrange(map_pred, map_obs, ncol = 2),
         width = 14, height = 6)

  ##############################################################################
  # SAVE SUMMARY
  ##############################################################################
  sink(here("results", paste0(sp, "_validation_summary.txt")))
  cat("Validation Summary:", sp, "\n")
  cat("Date:", as.character(Sys.time()), "\n\n")

  cat("-- Spatial Block CV (5-fold) --\n")
  print(df_cv)
  cat("\nMean AUC:", round(mean(df_cv$AUC), 3), "\n")
  cat("Mean TSS:", round(mean(df_cv$TSS), 3), "\n")
  cat("Mean RMSE:", round(mean(df_cv$RMSE), 3), "\n")
  cat("Optimal threshold:", round(opt_th, 3), "\n\n")

  cat("-- Continuous Calibration (polygon-level) --\n")
  cat("Spearman rho (raw):", round(spearman_raw$estimate, 3),
      "| p =", signif(spearman_raw$p.value, 3), "\n")
  cat("Spearman rho (CV mean +/- sd):", round(spearman_cv_mean, 3),
      "+/-", round(spearman_cv_sd, 3), "\n")
  cat("McFadden pseudo-R^2:", round(pseudo_r2, 3), "\n\n")

  cat("-- Spatial Autocorrelation (Moran's I on residuals) --\n")
  print(moran_res)
  cat("\nNote: Significant spatial autocorrelation is expected for\n")
  cat("environmental-only models that do not account for dispersal\n")
  cat("limitation or historical biogeographic processes.\n")
  cat("See docs/validation_diagnosis.md for interpretation.\n")
  sink()

  message("  Validation complete for ", sp,
          " | Mean AUC=", round(mean(df_cv$AUC), 3),
          " | Mean TSS=", round(mean(df_cv$TSS), 3),
          " | Spearman(raw)=", round(spearman_raw$estimate, 3))
}

message("\nStep 5 complete: validation done for all species.")
