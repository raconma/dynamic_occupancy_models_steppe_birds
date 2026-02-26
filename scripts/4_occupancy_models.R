###############################################################################
# 4_occupancy_models.R
#
# Purpose: Fit dynamic occupancy models (colext) for each steppe bird species.
#          Includes: data formatting, model fitting, GOF, response curves,
#          occupancy maps, and stochastic simulations.
#
# Inputs:  data/processed/{sp}/{sp}_occ_wide_dynamic.csv  (from step 3)
#          data/processed/{sp}/{sp}_scaling_params.rds     (from step 2)
#          data/raw/environmental_data/...                 (for prediction surface)
#          data/raw/topology_data/...                      (for prediction surface)
#          R/model_configs.R                               (species formulas)
#
# Outputs: results/{sp}_model_summary.txt
#          results/{sp}_gof_parboot.rds
#          results/{sp}_gof_mackenzie_bailey.rds
#          results/{sp}_simulation_prevalence.csv
#          results/{sp}_model_object.rds
#          data/processed/{sp}/occ_{sp}_prediction.csv
#          data/processed/{sp}/{sp}_OccuMap.tif
#          figs/{sp}_response_occupancy.png
#          figs/{sp}_response_colonization.png
#          figs/{sp}_response_extinction.png
#          figs/{sp}_response_detection.png
#          figs/{sp}_occupancy_map.png
#          figs/{sp}_prevalence_over_time.png
#
# Replaces: 4_occupancy_{otitar,ptealc,pteori,tettet}.R
###############################################################################

# -- Load packages --
library(here)
library(unmarked)
library(auk)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(raster)
library(terra)
library(AICcmodavg)
library(rmapshaper)

# Resolve namespace conflicts
select <- dplyr::select
filter <- dplyr::filter

# Load model configurations
source(here("R", "model_configs.R"))

theme_set(theme_bw())

# -- Species to process --
species_codes <- c("otitar", "ptealc", "pteori", "tettet")

# -- Constants --
YEARS   <- 2017:2022
T_years <- length(YEARS)
J_reps  <- 10            # max secondary periods per primary period
NSIM_PARBOOT <- 1000     # publication-grade GOF (was 10)
NSIM_MB_GOF  <- 500      # publication-grade MB GOF (was 3)
NSIM_PREV    <- 5000     # stochastic simulations

###############################################################################
# Helper: expand yearly covariates to observation-level dimensions
###############################################################################
expand_matrix <- function(mat, J) {
  expanded_mat <- kronecker(mat, matrix(1, nrow = 1, ncol = J))
  original_colnames <- colnames(mat)
  new_colnames <- as.vector(sapply(original_colnames,
                                    function(name) paste0(name, ".", 1:J)))
  colnames(expanded_mat) <- new_colnames
  return(expanded_mat)
}

###############################################################################
# Helper: apply TRAINING scaling to new data (fixes scaling mismatch bug)
###############################################################################
scale_with_params <- function(x, center, scale_val) {
  (x - center) / scale_val
}

###############################################################################
# Helper: plot response curves for one model component
###############################################################################
plot_response_curves <- function(formula_list, occ_umf, type, data_source,
                                 y_label, color, sp_code) {
  pred_list <- list()

  for (i in seq_along(formula_list)) {
    frm <- formula_list[[i]]
    # Fit univariate model
    new_model <- colext(
      psiformula     = if (type == "psi") frm else ~ 1,
      gammaformula   = if (type == "col") frm else ~ 1,
      epsilonformula = if (type == "ext") frm else ~ 1,
      pformula       = if (type == "det") frm else ~ 1,
      data = occ_umf
    )

    var_name <- gsub("~", "", as.character(frm))[2]
    var_name <- trimws(var_name)
    var_data <- data_source[, names(data_source) == var_name]

    sim_data <- data.frame(seq(min(var_data, na.rm = TRUE),
                                max(var_data, na.rm = TRUE),
                                length.out = 100))
    colnames(sim_data) <- var_name

    pred_df <- predict(new_model, type = type, newdata = sim_data,
                        appendData = TRUE)
    colnames(pred_df)[5] <- "cov_value"
    pred_df$CovariateName <- var_name
    pred_list <- rbind(pred_list, pred_df)
  }

  pred_list$CovariateName <- as.factor(pred_list$CovariateName)

  p <- ggplot(data = pred_list) +
    geom_ribbon(aes(x = cov_value, ymin = lower, ymax = upper),
                fill = "gray80") +
    geom_line(aes(x = cov_value, y = Predicted), color = color) +
    labs(x = "Covariate Values (standardised)", y = y_label) +
    theme(panel.border = element_rect(color = "black", fill = "transparent"),
          panel.background = element_rect(fill = "white")) +
    facet_wrap(. ~ CovariateName, scales = "free", ncol = 3)

  return(p)
}


###############################################################################
# MAIN LOOP: process each species
###############################################################################

for (sp in species_codes) {

  message("\n", strrep("=", 60))
  message("  SPECIES: ", sp)
  message(strrep("=", 60))

  # -- Load config --
  cfg <- get_model_config(sp)

  # -- Load data --
  occ_wide_clean <- read.csv(
    here("data", "processed", sp, paste0(sp, "_occ_wide_dynamic.csv"))
  )
  message("  Data: ", nrow(occ_wide_clean), " sites x ",
          ncol(occ_wide_clean), " columns")

  # -- Extract matrices by column indices --
  # Survey covariates (10 reps x 6 years = 60 columns)
  duration  <- as.matrix(occ_wide_clean[, c(33:42, 104:113, 175:184,
                                             246:255, 317:326, 388:397)])
  effort    <- as.matrix(occ_wide_clean[, c(43:52, 114:123, 185:194,
                                             256:265, 327:336, 398:407)])
  observers <- as.matrix(occ_wide_clean[, c(53:62, 124:133, 195:204,
                                             266:275, 337:346, 408:417)])
  time      <- occ_wide_clean[, c(23:32, 94:103, 165:174,
                                   236:245, 307:316, 378:387)]

  # Detection histories
  detections <- occ_wide_clean[, c(2:11, 74:83, 145:154,
                                    216:225, 287:296, 358:367)]
  detections <- detections %>% mutate(across(where(is.logical), as.integer))
  y.cross <- as.matrix(detections)
  y.cross[is.na(time) != is.na(y.cross)] <- NA

  # Site covariates (static)
  siteCovs <- occ_wide_clean[, c(17:22)]

  # Yearly site covariates
  EVI     <- occ_wide_clean[, c(428:433)]
  Land_Cover_Type_1_Percent_Class_0  <- occ_wide_clean[, c(434:439)]
  Land_Cover_Type_1_Percent_Class_10 <- occ_wide_clean[, c(440:445)]
  Land_Cover_Type_1_Percent_Class_12 <- occ_wide_clean[, c(446:451)]
  Land_Cover_Type_1_Percent_Class_13 <- occ_wide_clean[, c(452:457)]
  Land_Cover_Type_1_Percent_Class_14 <- occ_wide_clean[, c(458:463)]
  Land_Cover_Type_1_Percent_Class_6  <- occ_wide_clean[, c(464:469)]
  Land_Cover_Type_1_Percent_Class_7  <- occ_wide_clean[, c(470:475)]
  NDVI <- occ_wide_clean[, c(482:487)]
  pr   <- occ_wide_clean[, c(488:493)]
  tmmn <- occ_wide_clean[, c(494:499)]
  tmmx <- occ_wide_clean[, c(500:505)]

  # -- Standardise all covariates --
  time      <- scale(time)
  duration  <- scale(duration)
  effort    <- scale(effort)
  observers <- scale(observers)
  siteCovs  <- scale(siteCovs)

  EVI  <- scale(EVI)
  Land_Cover_Type_1_Percent_Class_0  <- scale(Land_Cover_Type_1_Percent_Class_0)
  Land_Cover_Type_1_Percent_Class_6  <- scale(Land_Cover_Type_1_Percent_Class_6)
  Land_Cover_Type_1_Percent_Class_7  <- scale(Land_Cover_Type_1_Percent_Class_7)
  Land_Cover_Type_1_Percent_Class_10 <- scale(Land_Cover_Type_1_Percent_Class_10)
  Land_Cover_Type_1_Percent_Class_12 <- scale(Land_Cover_Type_1_Percent_Class_12)
  Land_Cover_Type_1_Percent_Class_13 <- scale(Land_Cover_Type_1_Percent_Class_13)
  Land_Cover_Type_1_Percent_Class_14 <- scale(Land_Cover_Type_1_Percent_Class_14)
  NDVI <- scale(NDVI)
  pr   <- scale(pr)
  tmmn <- scale(tmmn)
  tmmx <- scale(tmmx)

  # -- Capture training scaling for dynamic covariates --
  # scale() on a matrix gives per-column center/scale (one per year).
  # We store these keyed by "covariate_year" so the simulation section
  # can apply the EXACT SAME scaling as training — not a per-year re-scaling.
  train_dyn_scale <- list()
  .dyn_mats <- list(
    EVI = EVI, NDVI = NDVI, pr = pr, tmmn = tmmn, tmmx = tmmx,
    Land_Cover_Type_1_Percent_Class_0  = Land_Cover_Type_1_Percent_Class_0,
    Land_Cover_Type_1_Percent_Class_6  = Land_Cover_Type_1_Percent_Class_6,
    Land_Cover_Type_1_Percent_Class_7  = Land_Cover_Type_1_Percent_Class_7,
    Land_Cover_Type_1_Percent_Class_10 = Land_Cover_Type_1_Percent_Class_10,
    Land_Cover_Type_1_Percent_Class_12 = Land_Cover_Type_1_Percent_Class_12,
    Land_Cover_Type_1_Percent_Class_13 = Land_Cover_Type_1_Percent_Class_13,
    Land_Cover_Type_1_Percent_Class_14 = Land_Cover_Type_1_Percent_Class_14
  )
  for (.nm in names(.dyn_mats)) {
    .ctrs <- attr(.dyn_mats[[.nm]], "scaled:center")
    .scls <- attr(.dyn_mats[[.nm]], "scaled:scale")
    for (.yi in seq_along(YEARS)) {
      train_dyn_scale[[paste0(.nm, "_", YEARS[.yi])]] <- list(
        center = .ctrs[.yi], scale = .scls[.yi]
      )
    }
  }
  rm(.dyn_mats, .nm, .ctrs, .scls, .yi)

  # -- Build year factor --
  n <- nrow(occ_wide_clean)
  years_df <- data.frame(matrix(rep(YEARS, each = n), n, T_years))
  years_df <- data.frame(lapply(years_df, as.factor))

  # -- Expand covariates to observation-level --
  NDVI_obs        <- expand_matrix(NDVI, J_reps)
  pr_obs          <- expand_matrix(pr, J_reps)
  topo_aspect_obs <- expand_matrix(siteCovs[, "topo_aspect", drop = FALSE],
                                    J_reps * T_years)
  topo_elev_obs   <- expand_matrix(siteCovs[, "topo_elev", drop = FALSE],
                                    J_reps * T_years)

  # -- Build unmarkedMultFrame --
  occ_umf <- unmarkedMultFrame(
    y = y.cross,
    siteCovs = data.frame(siteCovs),
    yearlySiteCovs = list(
      years  = years_df,
      EVI    = EVI,
      Land_Cover_Type_1_Percent_Class_0  = Land_Cover_Type_1_Percent_Class_0,
      Land_Cover_Type_1_Percent_Class_6  = Land_Cover_Type_1_Percent_Class_6,
      Land_Cover_Type_1_Percent_Class_7  = Land_Cover_Type_1_Percent_Class_7,
      Land_Cover_Type_1_Percent_Class_10 = Land_Cover_Type_1_Percent_Class_10,
      Land_Cover_Type_1_Percent_Class_12 = Land_Cover_Type_1_Percent_Class_12,
      Land_Cover_Type_1_Percent_Class_13 = Land_Cover_Type_1_Percent_Class_13,
      Land_Cover_Type_1_Percent_Class_14 = Land_Cover_Type_1_Percent_Class_14,
      NDVI = NDVI, pr = pr, tmmn = tmmn, tmmx = tmmx
    ),
    obsCovs = list(
      duration = duration, effort = effort, observers = observers, time = time,
      NDVI_obs = NDVI_obs, pr_obs = pr_obs,
      topo_aspect_obs = topo_aspect_obs, topo_elev_obs = topo_elev_obs
    ),
    numPrimary = T_years
  )

  ##############################################################################
  # FIT FINAL MODEL
  ##############################################################################
  message("  Fitting colext model...")

  Mod.final <- colext(
    psiformula     = cfg$psi_formula,
    gammaformula   = cfg$gamma_formula,
    epsilonformula = cfg$epsilon_formula,
    pformula       = cfg$p_formula,
    data = occ_umf
  )

  # Save model summary
  sink(here("results", paste0(sp, "_model_summary.txt")))
  cat("Species:", sp, "\n")
  cat("Date:", as.character(Sys.time()), "\n\n")
  summary(Mod.final)
  sink()

  # Save model object
  saveRDS(Mod.final, here("results", paste0(sp, "_model_object.rds")))
  message("  Model saved. AIC = ", round(Mod.final@AIC, 2))

  ##############################################################################
  # GOODNESS OF FIT (publication-grade nsim)
  ##############################################################################
  message("  Running GOF (parboot nsim=", NSIM_PARBOOT, ")...")
  set.seed(42)  # Reproducible GOF
  GOF <- parboot(Mod.final, nsim = NSIM_PARBOOT)
  saveRDS(GOF, here("results", paste0(sp, "_gof_parboot.rds")))

  message("  Running MacKenzie-Bailey GOF (nsim=", NSIM_MB_GOF, ")...")
  set.seed(42)
  occ_gof <- mb.gof.test(Mod.final, nsim = NSIM_MB_GOF, plot.hist = FALSE)
  saveRDS(occ_gof, here("results", paste0(sp, "_gof_mackenzie_bailey.rds")))

  ##############################################################################
  # RESPONSE VARIABLE PLOTS
  ##############################################################################
  message("  Generating response curves...")

  # Occupancy
  occ_formulas <- lapply(cfg$psi_vars, function(v) as.formula(paste0("~", v)))
  p_occ <- plot_response_curves(occ_formulas, occ_umf, "psi",
                                 as.data.frame(occ_umf@siteCovs),
                                 "Probability of Site Occupancy",
                                 "seagreen4", sp)
  ggsave(here("figs", paste0(sp, "_response_occupancy.png")),
         p_occ, width = 10, height = 6)

  # Colonization
  col_formulas <- lapply(cfg$gamma_vars, function(v) as.formula(paste0("~", v)))
  p_col <- plot_response_curves(col_formulas, occ_umf, "col",
                                 as.data.frame(occ_umf@yearlySiteCovs),
                                 "Colonization Probability",
                                 "royalblue3", sp)
  ggsave(here("figs", paste0(sp, "_response_colonization.png")),
         p_col, width = 10, height = 6)

  # Extinction
  ext_formulas <- lapply(cfg$epsilon_vars, function(v) as.formula(paste0("~", v)))
  p_ext <- plot_response_curves(ext_formulas, occ_umf, "ext",
                                 as.data.frame(occ_umf@yearlySiteCovs),
                                 "Extinction Probability",
                                 "royalblue3", sp)
  ggsave(here("figs", paste0(sp, "_response_extinction.png")),
         p_ext, width = 10, height = 6)

  # Detection
  det_var_names <- all.vars(cfg$p_formula)
  det_formulas <- lapply(det_var_names, function(v) as.formula(paste0("~", v)))
  p_det <- plot_response_curves(det_formulas, occ_umf, "det",
                                 as.data.frame(occ_umf@obsCovs),
                                 "Detection Probability",
                                 "royalblue3", sp)
  ggsave(here("figs", paste0(sp, "_response_detection.png")),
         p_det, width = 10, height = 6)

  ##############################################################################
  # OCCUPANCY MAP (spatial prediction)
  ##############################################################################
  message("  Building occupancy map...")

  # Load rasters for prediction surface
  variables <- brick(stack(
    here("data", "raw", "environmental_data", "environmental_data_occ",
         "variables_spain.grd")))
  aspect <- raster(here("data", "raw", "topology_data", "topo_aspect_masked.tif"))
  names(aspect) <- "topo_aspect"
  elev <- raster(here("data", "raw", "topology_data", "topo_elev_masked.tif"))
  names(elev) <- "topo_elev"

  aspect_resampled <- resample(aspect, variables, method = "bilinear")
  elev_resampled   <- resample(elev, variables, method = "bilinear")
  variables <- addLayer(variables, aspect_resampled, elev_resampled)

  # Select prediction variables
  variables.sel <- variables[[cfg$psi_vars]]
  pred_surface <- data.frame(rasterToPoints(variables.sel))
  names(pred_surface)[names(pred_surface) == "x"] <- "longitude"
  names(pred_surface)[names(pred_surface) == "y"] <- "latitude"
  pred_surface <- pred_surface %>% drop_na()

  # FIX: Scale prediction surface using TRAINING scaling parameters
  scaling_params <- readRDS(
    here("data", "processed", sp, paste0(sp, "_scaling_params.rds"))
  )

  pred_surface_std <- pred_surface
  for (v in cfg$psi_vars) {
    if (v %in% names(scaling_params)) {
      pred_surface_std[[v]] <- scale_with_params(
        pred_surface_std[[v]],
        scaling_params[[v]]$center,
        scaling_params[[v]]$scale
      )
    }
  }

  # Predict occupancy
  occ_pred <- predict(Mod.final,
                       newdata = as.data.frame(pred_surface_std[, cfg$psi_vars]),
                       type = "psi")

  pred_occ <- bind_cols(
    pred_surface_std,
    occ_prob = occ_pred$Predicted,
    occ_se   = occ_pred$SE
  ) %>%
    dplyr::select(longitude, latitude, occ_prob, occ_se)

  # Save prediction CSV
  write.csv(pred_occ,
            here("data", "processed", sp, paste0("occ_", sp, "_prediction.csv")),
            row.names = FALSE)

  # Create and save raster
  map_proj <- "+proj=longlat +datum=WGS84 +no_defs"
  r_pred <- pred_occ %>%
    st_as_sf(coords = c("longitude", "latitude"), crs = map_proj) %>%
    st_transform(crs = raster::projection(variables.sel[[1]])) %>%
    rasterize(variables.sel[[1]])
  r_pred <- r_pred[[c("occ_prob", "occ_se")]]

  spain <- ne_countries(country = "spain", scale = "medium", returnclass = "sf")
  spain_crop <- ms_filter_islands(spain, min_area = 100000000000,
                                   drop_null_geometries = TRUE)

  r_pred_proj <- projectRaster(r_pred, crs = map_proj, method = "ngb")
  r_pred_proj_crop <- crop(r_pred_proj, spain_crop)
  r_pred_proj_masked <- mask(r_pred_proj_crop, spain_crop)

  # Save GeoTIFF
  writeRaster(r_pred_proj_masked[[1]],
              here("data", "processed", sp, paste0(sp, "_OccuMap.tif")),
              format = "GTiff", overwrite = TRUE)

  # Plot occupancy map
  r_pred_df <- as.data.frame(r_pred_proj_masked[[1]], xy = TRUE, na.rm = TRUE)
  names(r_pred_df)[3] <- "occ_prob"

  p_map <- ggplot() +
    geom_tile(data = r_pred_df,
              aes(x = x, y = y, fill = occ_prob)) +
    scale_fill_viridis_c(name = "p(occupancy)") +
    labs(x = "Longitude", y = "Latitude",
         title = paste0("Occupancy probability: ", sp)) +
    theme(panel.border = element_rect(color = "black", fill = "transparent"))
  ggsave(here("figs", paste0(sp, "_occupancy_map.png")),
         p_map, width = 10, height = 8)

  ##############################################################################
  # STOCHASTIC SIMULATIONS
  #
  # IMPORTANT: All covariates must be scaled using the TRAINING scaling
  # parameters, not re-scaled from scratch. Three bugs were fixed here:
  #   Bug 1 (CRITICAL): psi covariates were passed UNSCALED to predict().
  #          The model expects scaled inputs (trained on scale(siteCovs)).
  #          Raw values (e.g. grass_cover=45) caused logit saturation → ψ₁≈1.
  #   Bug 2 (HIGH): gamma/epsilon covariates were re-scaled per year with
  #          their own mean/sd instead of the training scaling parameters.
  #   Bug 3 (HIGH): Land Cover classes 6,7,10,12,14 were not scaled in
  #          training but got scaled in simulation → coefficient mismatch.
  #          Fixed above (now all classes are scaled in training too).
  ##############################################################################
  message("  Running stochastic simulations (nsim=", NSIM_PREV, ")...")

  set.seed(123)  # reproducible simulations

  predict_data <- occ_wide_clean %>%
    drop_na(paste0("NDVI_", YEARS))

  S <- nrow(predict_data)
  T_sim <- occ_umf@numPrimary

  # --- Initial occupancy (ψ₁) ---
  # FIX Bug 1: scale psi covariates using TRAINING params (from step 2).
  # Previously passed raw values to predict() → logit saturation → ψ₁ ≈ 1.
  psi_data <- predict_data[, cfg$psi_vars, drop = FALSE]
  for (v in cfg$psi_vars) {
    if (v %in% names(scaling_params)) {
      psi_data[[v]] <- scale_with_params(
        psi_data[[v]], scaling_params[[v]]$center, scaling_params[[v]]$scale
      )
    }
  }
  psi1 <- predict(Mod.final, type = "psi",
                   newdata = as.data.frame(psi_data))$Predicted

  message("    psi1: mean=", round(mean(psi1, na.rm = TRUE), 4),
          " | range=[", round(min(psi1, na.rm = TRUE), 4), ", ",
          round(max(psi1, na.rm = TRUE), 4), "]")

  # --- Colonisation and extinction per year ---
  # FIX Bug 2: scale dynamic covariates using train_dyn_scale (captured from
  # the training scale() calls), not a fresh scale() per year.
  col_mat <- matrix(NA_real_, nrow = S, ncol = T_sim)
  ext_mat <- matrix(NA_real_, nrow = S, ncol = T_sim)

  for (i in 1:T_sim) {
    yr <- YEARS[i]

    # Colonisation
    new_col <- data.frame(row.names = seq_len(S))
    for (v in cfg$gamma_vars) {
      key <- paste0(v, "_", yr)
      raw_val <- predict_data[[key]]
      sc <- train_dyn_scale[[key]]
      if (!is.null(sc)) {
        new_col[[v]] <- scale_with_params(raw_val, sc$center, sc$scale)
      } else {
        new_col[[v]] <- raw_val
        warning("No training scaling for: ", key, " — using raw values")
      }
    }
    col_mat[, i] <- predict(Mod.final, type = "col",
                             newdata = new_col)$Predicted

    # Extinction
    new_ext <- data.frame(row.names = seq_len(S))
    for (v in cfg$epsilon_vars) {
      key <- paste0(v, "_", yr)
      raw_val <- predict_data[[key]]
      sc <- train_dyn_scale[[key]]
      if (!is.null(sc)) {
        new_ext[[v]] <- scale_with_params(raw_val, sc$center, sc$scale)
      } else {
        new_ext[[v]] <- raw_val
        warning("No training scaling for: ", key, " — using raw values")
      }
    }
    ext_mat[, i] <- predict(Mod.final, type = "ext",
                             newdata = new_ext)$Predicted
  }

  message("    col mean: ", paste(round(colMeans(col_mat, na.rm = TRUE), 4),
                                   collapse = ", "))
  message("    ext mean: ", paste(round(colMeans(ext_mat, na.rm = TRUE), 4),
                                   collapse = ", "))

  # --- Stochastic simulations ---
  prevT <- matrix(NA, NSIM_PREV, T_sim)

  for (nn in 1:NSIM_PREV) {
    Z <- matrix(NA, S, T_sim)
    Z[, 1] <- (runif(S) < psi1) * 1
    for (ii in 2:T_sim) {
      tmp <- Z[, ii - 1] * (1 - ext_mat[, ii]) +
        (1 - Z[, ii - 1]) * col_mat[, ii]
      Z[, ii] <- (runif(S) < tmp) * 1
    }
    prevT[nn, ] <- colMeans(Z, na.rm = TRUE)
  }

  # Save simulation results
  sim_results <- data.frame(
    year     = as.character(YEARS),
    mean_occ = colMeans(prevT),
    sd_occ   = apply(prevT, 2, sd),
    mean_col = colMeans(col_mat, na.rm = TRUE),
    sd_col   = apply(col_mat, 2, sd, na.rm = TRUE),
    mean_ext = colMeans(ext_mat, na.rm = TRUE),
    sd_ext   = apply(ext_mat, 2, sd, na.rm = TRUE)
  )
  write.csv(sim_results,
            here("results", paste0(sp, "_simulation_prevalence.csv")),
            row.names = FALSE)

  # Plot prevalence over time
  p_prev <- ggplot(sim_results, aes(x = year, y = mean_occ, group = 1)) +
    geom_errorbar(aes(ymin = mean_occ - sd_occ, ymax = mean_occ + sd_occ),
                  width = 0.1) +
    geom_line() + geom_point() +
    labs(x = "Year", y = "Predicted mean occupancy",
         title = paste0("Occupancy trend: ", sp))
  ggsave(here("figs", paste0(sp, "_prevalence_over_time.png")),
         p_prev, width = 8, height = 5)

  message("  Step 4 complete for ", sp)
}

message("\nStep 4 complete: all occupancy models fitted.")
