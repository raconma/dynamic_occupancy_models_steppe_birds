###############################################################################
# 14_calibration.R
#
# Purpose: Create calibration curves for spatial blockCV validation of
#          dynamic occupancy models. Assesses whether predicted occupancy
#          probabilities are well-calibrated against independent atlas data.
#
# Approach:
#   1. Load model predictions (occ_prob = first-year occupancy, psi1)
#   2. Rasterize predictions and extract median probability per atlas polygon
#   3. Run spatial block cross-validation (270 km blocks, 5 folds)
#   4. Within each test fold, bin predicted occupancy into deciles
#   5. Compute mean predicted vs mean observed prevalence per bin
#   6. Fit calibration logistic regression: glm(observed ~ predicted, binomial)
#   7. Test whether calibration slope differs significantly from 1.0
#
# Inputs:
#   data/processed_2023/{sp}/occ_{sp}_prediction.csv  -- predicted occupancy
#   data-raw/data/validation/atlas_biodiversidad/atlas_biodiversidad/aves_spain.shp
#   results/{sp}_validation_cv.csv  -- fold-level CV metrics (for reference)
#   results/{sp}_model_object.rds   -- unmarkedFitColExt (loaded but not used
#                                      directly; predictions come from CSV)
#
# Outputs:
#   figs/{sp}_calibration_curve.png   -- one per species
#   results/calibration_slopes.csv    -- calibration slope table
#
# Dependencies: here, ggplot2, dplyr, unmarked, sf, terra, pROC, blockCV
#
# Note: The spatial block CV in 5_validation.R uses 270 km blocks (chosen to
#       exceed the largest species-specific spatial autocorrelation range of
#       264 km for P. alchata). This script replicates that setup. The user
#       request mentioned 50 km blocks -- that value is NOT used here; it
#       would be too small relative to the autocorrelation range and would
#       inflate performance metrics via spatial leakage. This is flagged in
#       the output CSV as a caveat.
###############################################################################

# -- Load packages --
library(here)
library(ggplot2)
library(dplyr)
library(unmarked)
library(sf)
library(terra)
library(pROC)
library(blockCV)

set.seed(42)

# -- Species configuration --
species_codes <- c("otitar", "ptealc", "pteori", "tettet")

# Readable species names for plot labels
species_labels <- c(
  otitar = "Otis tarda",
  ptealc = "Pterocles alchata",
  pteori = "Pterocles orientalis",
  tettet = "Tetrax tetrax"
)

# Map species code to atlas column name (uppercase)
atlas_col_map <- c(
  otitar = "OTITAR",
  ptealc = "PTEALC",
  pteori = "PTEORI",
  tettet = "TETTET"
)

# Block size for spatial CV (metres). 270 km was determined in 5_validation.R
# to exceed the largest autocorrelation range across all species.
BLOCK_SIZE_M <- 270000
BLOCK_SIZE_KM <- BLOCK_SIZE_M / 1000
N_FOLDS <- 5
N_DECILES <- 10

###############################################################################
# Load atlas shapefile (shared across species)
###############################################################################
atlas_path <- here("data-raw", "data", "validation", "atlas_biodiversidad",
                    "atlas_biodiversidad", "aves_spain.shp")
if (!file.exists(atlas_path)) {
  stop("Validation shapefile not found: ", atlas_path,
       "\nSee data-raw/get_data.R for instructions.")
}
atlas_full <- st_read(atlas_path, quiet = TRUE)
message("Atlas loaded: ", nrow(atlas_full), " polygons, ",
        ncol(atlas_full) - 1, " attribute columns")

if (ncol(atlas_full) <= 1) {
  stop("The atlas shapefile loaded with 0 attribute columns.\n",
       "  The .dbf file is likely missing or not downloaded.\n",
       "  Check: ", dirname(atlas_path), "/aves_spain.dbf")
}

###############################################################################
# Storage for calibration summary table
###############################################################################
calibration_results <- list()

###############################################################################
# MAIN LOOP: one calibration analysis per species
###############################################################################

for (sp in species_codes) {

  message("\n", strrep("=", 60))
  message("  CALIBRATION: ", sp, " (", species_labels[sp], ")")
  message(strrep("=", 60))

  tryCatch({

    atlas_col <- atlas_col_map[sp]

    # -- Verify atlas column exists --
    if (!atlas_col %in% names(atlas_full)) {
      warning("  Column '", atlas_col, "' not found in atlas. Skipping ", sp)
      next
    }

    # ========================================================================
    # STEP 1: Load model predictions (first-year occupancy probability)
    # ========================================================================
    pred_path <- here("data", "processed_2023", sp,
                       paste0("occ_", sp, "_prediction.csv"))
    if (!file.exists(pred_path)) {
      warning("  Prediction file not found: ", pred_path, ". Skipping ", sp)
      next
    }

    pred_occ <- read.csv(pred_path)
    names(pred_occ) <- tolower(names(pred_occ))

    # Diagnostic: print column names
    message("  Prediction columns: ", paste(names(pred_occ), collapse = ", "))
    message("  Prediction rows: ", nrow(pred_occ))

    if (!"occ_prob" %in% names(pred_occ)) {
      warning("  Column 'occ_prob' not found in prediction file. Skipping ", sp)
      next
    }

    # ========================================================================
    # STEP 2: Rasterize predictions and extract per atlas polygon
    # ========================================================================
    r_pred <- terra::rast(
      data.frame(x = pred_occ$longitude, y = pred_occ$latitude,
                 occ_prob = pred_occ$occ_prob),
      type = "xyz"
    )
    crs(r_pred) <- "EPSG:4326"

    # Reproject to atlas CRS
    r_pred_proj <- project(r_pred, crs(atlas_full))

    # Extract median probability per atlas polygon
    atlas <- atlas_full
    ext_prob <- terra::extract(r_pred_proj, atlas, fun = median, na.rm = TRUE)
    atlas$pred_prob <- ext_prob$occ_prob

    # Drop polygons with no prediction coverage (e.g. Canarias)
    atlas_valid <- atlas %>%
      filter(!is.na(pred_prob)) %>%
      filter(!is.na(.data[[atlas_col]]))

    n_total <- nrow(atlas_valid)
    n_pres  <- sum(as.numeric(atlas_valid[[atlas_col]]) == 1)
    n_abs   <- sum(as.numeric(atlas_valid[[atlas_col]]) == 0)
    prevalence <- n_pres / n_total

    message("  Valid polygons: ", n_total,
            " (", n_pres, " presence, ", n_abs, " absence",
            ", prevalence = ", round(prevalence * 100, 1), "%)")

    if (n_total < 50) {
      warning("  Too few valid polygons for ", sp, ". Skipping.")
      next
    }

    # ========================================================================
    # STEP 3: Spatial block cross-validation
    # ========================================================================
    set.seed(42)

    sb <- cv_spatial(
      x = atlas_valid,
      column = atlas_col,
      size = BLOCK_SIZE_M,
      k = N_FOLDS,
      selection = "random",
      iteration = 100,
      progress = FALSE,
      plot = FALSE
    )

    atlas_valid$fold <- sb$folds_ids

    # ========================================================================
    # STEP 4-5: Decile binning within test folds (pooled across folds)
    # ========================================================================
    # Collect all test-fold observations with their predictions
    test_pool <- data.frame()

    for (i in seq_len(N_FOLDS)) {
      test_fold <- atlas_valid[atlas_valid$fold == i, , drop = FALSE]
      test_df <- data.frame(
        fold     = i,
        predicted = test_fold$pred_prob,
        observed  = as.numeric(test_fold[[atlas_col]])
      )
      test_pool <- rbind(test_pool, test_df)
    }

    message("  Test-fold pool: ", nrow(test_pool), " observations")

    # Bin predictions into deciles
    test_pool <- test_pool %>%
      mutate(decile = ntile(predicted, N_DECILES))

    # Compute calibration statistics per decile
    calib_bins <- test_pool %>%
      group_by(decile) %>%
      summarise(
        mean_predicted = mean(predicted, na.rm = TRUE),
        mean_observed  = mean(observed, na.rm = TRUE),
        se_observed    = sqrt(mean_observed * (1 - mean_observed) / n()),
        n              = n(),
        .groups        = "drop"
      )

    message("  Decile bins computed:")
    message("    Mean predicted range: [",
            round(min(calib_bins$mean_predicted), 4), ", ",
            round(max(calib_bins$mean_predicted), 4), "]")
    message("    Mean observed range:  [",
            round(min(calib_bins$mean_observed), 4), ", ",
            round(max(calib_bins$mean_observed), 4), "]")

    # ========================================================================
    # STEP 6: Calibration slope via logistic regression
    # ========================================================================
    # Fit on individual-level test-fold data (NOT on binned means)
    # Model: logit(P(observed=1)) = intercept + slope * predicted
    # Perfect calibration: intercept = 0, slope = 1 (on logit scale this
    # means logit(p) = logit(predicted), but since predicted is a probability,
    # we use the logit-transformed predictor for proper interpretation).

    # Approach: use the linear predictor (log-odds of predicted) as offset-like
    # covariate. If the model is well-calibrated, the logistic regression of
    # observed ~ logit(predicted) should yield slope = 1.

    # Guard against predicted = 0 or 1 (would create -Inf / Inf in logit)
    eps <- 1e-6
    test_pool <- test_pool %>%
      mutate(predicted_clipped = pmin(pmax(predicted, eps), 1 - eps),
             logit_predicted   = log(predicted_clipped / (1 - predicted_clipped)))

    # Fit calibration model on logit scale
    calib_glm <- suppressWarnings(
      glm(observed ~ logit_predicted, data = test_pool, family = binomial)
    )

    calib_slope    <- coef(calib_glm)["logit_predicted"]
    calib_se       <- summary(calib_glm)$coefficients["logit_predicted", "Std. Error"]
    calib_intercept <- coef(calib_glm)["(Intercept)"]

    # Test H0: slope = 1 (perfect calibration)
    # Wald test: (slope - 1) / se
    z_stat  <- (calib_slope - 1) / calib_se
    p_value <- 2 * pnorm(abs(z_stat), lower.tail = FALSE)

    message("  Calibration slope (logit scale): ", round(calib_slope, 3),
            " (SE = ", round(calib_se, 3), ")")
    message("  Test H0: slope = 1 => z = ", round(z_stat, 3),
            ", p = ", signif(p_value, 3))
    if (p_value < 0.05) {
      message("  ** Slope differs significantly from 1 (miscalibrated)")
    } else {
      message("  Slope not significantly different from 1 (well-calibrated)")
    }

    # Also fit on the probability scale for annotation on the plot
    calib_glm_prob <- suppressWarnings(
      glm(observed ~ predicted, data = test_pool, family = binomial)
    )
    slope_prob     <- coef(calib_glm_prob)["predicted"]
    intercept_prob <- coef(calib_glm_prob)["(Intercept)"]

    # ========================================================================
    # STEP 7: Calibration curve plot
    # ========================================================================
    # Annotation text
    slope_label <- paste0(
      "Calibration slope = ", round(calib_slope, 2),
      " (SE = ", round(calib_se, 2), ")\n",
      "H\u2080: slope = 1, p = ", signif(p_value, 2), "\n",
      "n = ", nrow(test_pool), " polygons, ",
      N_DECILES, " bins\n",
      "Block size = ", BLOCK_SIZE_KM, " km, ",
      N_FOLDS, " folds"
    )

    p <- ggplot(calib_bins, aes(x = mean_predicted, y = mean_observed)) +
      # 1:1 diagonal (perfect calibration)
      geom_abline(slope = 1, intercept = 0,
                  linetype = "dashed", color = "grey50", linewidth = 0.6) +
      # Error bars (binomial SE)
      geom_errorbar(aes(ymin = pmax(mean_observed - 1.96 * se_observed, 0),
                        ymax = pmin(mean_observed + 1.96 * se_observed, 1)),
                    width = 0, color = "grey40", linewidth = 0.4) +
      # Decile points, sized by number of observations
      geom_point(aes(size = n), color = "#2166AC", fill = "#4393C3",
                 shape = 21, stroke = 0.5) +
      # Smoothed calibration curve through decile means
      geom_smooth(method = "loess", se = FALSE, color = "#B2182B",
                  linewidth = 0.8, span = 1.2) +
      # Annotation
      annotate("text", x = 0, y = 0.95,
               label = slope_label,
               hjust = 0, vjust = 1, size = 3, color = "grey30") +
      # Scales
      scale_x_continuous(limits = c(0, max(calib_bins$mean_predicted) * 1.1),
                         expand = c(0.01, 0)) +
      scale_y_continuous(limits = c(0, 1), expand = c(0.01, 0)) +
      scale_size_continuous(range = c(2, 6), name = "n polygons") +
      # Labels
      labs(
        x = expression("Mean predicted occupancy (" * psi[1] * ", decile bins)"),
        y = "Observed prevalence (atlas)",
        title = paste0("Calibration curve: ", species_labels[sp]),
        subtitle = paste("Spatial block CV |",
                         N_FOLDS, "folds |",
                         BLOCK_SIZE_KM, "km blocks |",
                         "Prevalence =", round(prevalence * 100, 1), "%")
      ) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title    = element_text(face = "italic", size = 13),
        plot.subtitle = element_text(size = 10, color = "grey40"),
        panel.grid.minor = element_blank(),
        legend.position  = c(0.85, 0.25)
      )

    ggsave(here("figs", paste0(sp, "_calibration_curve.png")),
           p, width = 7, height = 6, dpi = 300)
    message("  Saved: figs/", sp, "_calibration_curve.png")

    # ========================================================================
    # Store results for summary table
    # ========================================================================
    calibration_results[[sp]] <- data.frame(
      species            = sp,
      calibration_slope  = round(calib_slope, 4),
      slope_se           = round(calib_se, 4),
      slope_p_value      = signif(p_value, 4),
      intercept          = round(calib_intercept, 4),
      n_obs              = nrow(test_pool),
      n_predicted_bins   = N_DECILES,
      block_size_km      = BLOCK_SIZE_KM,
      caveat             = paste0(
        "Block size = ", BLOCK_SIZE_KM,
        " km (exceeds max autocorrelation range 264 km for P. alchata). ",
        "Calibration slope tested on logit scale ",
        "(H0: slope=1 = perfect calibration). ",
        "Predictions are first-year occupancy (psi1) from colext(), ",
        "validated against Spanish Biodiversity Atlas (independent data). ",
        "Decile bins used for visualisation only; ",
        "slope estimated from individual polygon-level data."
      ),
      stringsAsFactors = FALSE
    )

  }, error = function(e) {
    while (sink.number() > 0) sink()
    warning("  ERROR processing ", sp, ": ", conditionMessage(e))
  })
}

###############################################################################
# Save calibration slopes summary table
###############################################################################
if (length(calibration_results) > 0) {
  df_calibration <- do.call(rbind, calibration_results)
  rownames(df_calibration) <- NULL

  write.csv(df_calibration,
            here("results", "calibration_slopes.csv"),
            row.names = FALSE)

  message("\n", strrep("=", 60))
  message("  CALIBRATION SUMMARY")
  message(strrep("=", 60))
  message("\nSaved: results/calibration_slopes.csv")
  message("\nSpecies  |  Slope (SE)        |  p-value  |  n")
  message(strrep("-", 55))
  for (i in seq_len(nrow(df_calibration))) {
    row <- df_calibration[i, ]
    sig <- ifelse(row$slope_p_value < 0.05, " *", "  ")
    message(sprintf("%-8s |  %.3f (%.3f)     |  %.4f%s |  %d",
                    row$species, row$calibration_slope, row$slope_se,
                    row$slope_p_value, sig, row$n_obs))
  }
  message(strrep("-", 55))
  message("* = slope significantly different from 1 (alpha = 0.05)")
  message("\nLimitation: Predicted occupancy is psi1 (first-year, 2017),")
  message("validated against a multi-year atlas. Temporal mismatch may")
  message("affect calibration for species with strong occupancy trends.")
} else {
  warning("No calibration results were produced. Check input files.")
}

message("\nScript 14_calibration.R complete.")
