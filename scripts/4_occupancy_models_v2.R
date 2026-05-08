###############################################################################
# 4_occupancy_models_v2.R
#
# Purpose: Same as scripts/4_occupancy_models.R (dynamic occupancy fits via
#          colext) but with the steppe-representativeness covariate added to
#          the detection sub-model. This is the "just run it" script for
#          Raul: it auto-runs scripts/3b_add_stepRep.R if needed, applies
#          the same data preparation as v1, and overrides cfg$p_formula
#          inline so R/model_configs.R does not need to be edited.
#
# Inputs:  data/processed_2023/{sp}/{sp}_occ_wide_dynamic.csv  (from step 3,
#                with stepRep_*_<year> columns added by step 3b)
#          data/processed_2023/{sp}/{sp}_scaling_params.rds     (from step 2)
#          data/derived/stepRep_cellyear_{sp}.csv               (build_stepRep.R)
#          data-raw/data/environmental_data/...                 (prediction surface)
#          data-raw/data/topology_data/...                      (prediction surface)
#          R/model_configs.R                                    (species formulas)
#
# Outputs: All outputs identical to v1 but with a "_v2" tag inserted right
#          after the species code, so the v1 fits are NOT overwritten and
#          can be compared side-by-side. Examples:
#            results/{sp}_v2_model_summary.txt
#            results/{sp}_v2_gof_parboot.rds
#            results/{sp}_v2_simulation_prevalence.csv
#            results/{sp}_v2_model_object.rds
#            results/{sp}_v2_train_dyn_scale.rds
#            data/processed_2023/{sp}/occ_{sp}_v2_prediction.csv
#            data/processed_2023/{sp}/{sp}_v2_OccuMap.tif
#            figs/{sp}_v2_response_*.png
#            figs/{sp}_v2_occupancy_map.png
#            figs/{sp}_v2_prevalence_over_time.png
#
# Notes on integration of stepRep (vs scripts/4_occupancy_models.R):
#   1. Auto-runs scripts/3b_add_stepRep.R at the top if stepRep_*_<year>
#      columns are missing from any species' occ_wide_dynamic.csv.
#   2. Extracts stepRep_strict_500m_2017..2023 as an n_sites x T_years
#      matrix, scales it, and pushes it through the same per-year
#      training-scale capture as NDVI / pr / tmmn / tmmx so the
#      attribution scripts (8, 10) can reapply the same transform.
#   3. Builds stepRep_obs via expand_matrix(stepRep, J_reps), exactly
#      analogous to NDVI_obs / pr_obs.
#   4. Adds stepRep to yearlySiteCovs and stepRep_obs to obsCovs.
#   5. Overrides cfg$p_formula inline:
#         cfg$p_formula <- update(cfg$p_formula, ~ . + stepRep_obs)
#      so R/model_configs.R does NOT need to be edited.
#
# Sensitivity variants (stepRep_strict_1km, stepRep_broad_500m,
# stepRep_broad_1km) are available in the same wide CSV after step 3b;
# change the STEPREP_VARIANT constant below to switch.
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
YEARS   <- 2017:2023
T_years <- length(YEARS)
J_reps  <- 10            # max secondary periods per primary period
NSIM_PARBOOT <- 10       # testing; set to 1000 for publication
NSIM_MB_GOF  <- 10       # testing; set to 500 for publication
NSIM_PREV    <- 500      # stochastic simulations (testing; 5000 for publication)

# -- v2-specific constants --
OUT_SUFFIX      <- "_v2"
STEPREP_VARIANT <- "stepRep_strict_500m"  # primary; switch to *_1km / broad_*
                                          # for sensitivity reruns.

###############################################################################
# Auto-run step 3b if stepRep columns are missing in any species
###############################################################################
needs_3b <- FALSE
for (sp in species_codes) {
  fp <- here("data", "processed_2023", sp,
             paste0(sp, "_occ_wide_dynamic.csv"))
  if (!file.exists(fp)) {
    stop("missing input: ", fp,
         "  -- run scripts/1..3 first.")
  }
  hdr <- readLines(fp, n = 1)
  expected <- paste0(STEPREP_VARIANT, "_", YEARS)
  if (!all(sapply(expected, function(x) grepl(paste0("\\b", x, "\\b"), hdr)))) {
    needs_3b <- TRUE; break
  }
}
if (needs_3b) {
  message("[v2] stepRep_*_<year> columns missing in at least one species --",
          " auto-running scripts/3b_add_stepRep.R")
  source(here("scripts", "3b_add_stepRep.R"))
}

###############################################################################
# Helper: extract observation-level columns in correct year-major, rep-minor order
###############################################################################
get_obs_cols <- function(df, prefix, years, reps = 1:10) {
  cols <- c()
  for (yr in years) {
    for (k in reps) {
      cn <- paste0(prefix, ".", k, ".", yr)
      if (cn %in% names(df)) cols <- c(cols, cn)
    }
  }
  cols
}

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
  message("  SPECIES: ", sp, "  (v2: detection ~ ... + stepRep_obs)")
  message(strrep("=", 60))

  # -- Load config and override p_formula to add stepRep_obs --
  cfg <- get_model_config(sp)
  cfg$p_formula <- update(cfg$p_formula, ~ . + stepRep_obs)
  message("  p_formula (v2): ",
          paste(deparse(cfg$p_formula), collapse = " "))

  # -- Load data --
  occ_wide_clean <- read.csv(
    here("data", "processed_2023", sp, paste0(sp, "_occ_wide_dynamic.csv"))
  )
  message("  Data: ", nrow(occ_wide_clean), " sites x ",
          ncol(occ_wide_clean), " columns")

  # -- Extract matrices by column NAMES (robust to any year count) --
  # Survey covariates (10 reps x T_years)
  dur_cols  <- get_obs_cols(occ_wide_clean, "duration_minutes", YEARS)
  eff_cols  <- get_obs_cols(occ_wide_clean, "effort_distance_km", YEARS)
  obs_cols  <- get_obs_cols(occ_wide_clean, "number_observers", YEARS)
  time_cols <- get_obs_cols(occ_wide_clean, "time_observations_started", YEARS)

  duration  <- as.matrix(occ_wide_clean[, dur_cols])
  effort    <- as.matrix(occ_wide_clean[, eff_cols])
  observers <- as.matrix(occ_wide_clean[, obs_cols])
  time      <- occ_wide_clean[, time_cols]

  # Detection histories
  y_cols <- get_obs_cols(occ_wide_clean, "y", YEARS)
  detections <- occ_wide_clean[, y_cols]
  detections <- detections %>% mutate(across(where(is.logical), as.integer))
  y.cross <- as.matrix(detections)
  y.cross[is.na(time) != is.na(y.cross)] <- NA

  stopifnot(ncol(y.cross) == J_reps * T_years)
  message("  Detection matrix: ", nrow(y.cross), " x ", ncol(y.cross))

  # Site covariates (static, by name)
  siteCovs <- occ_wide_clean[, c("bio1", "bio2", "tree_cover", "grass_cover",
                                   "topo_aspect", "topo_elev")]

  # Yearly site covariates (by name pattern)
  EVI  <- occ_wide_clean[, paste0("EVI_", YEARS)]
  Land_Cover_Type_1_Percent_Class_0  <- occ_wide_clean[, paste0("Land_Cover_Type_1_Percent_Class_0_", YEARS)]
  Land_Cover_Type_1_Percent_Class_6  <- occ_wide_clean[, paste0("Land_Cover_Type_1_Percent_Class_6_", YEARS)]
  Land_Cover_Type_1_Percent_Class_7  <- occ_wide_clean[, paste0("Land_Cover_Type_1_Percent_Class_7_", YEARS)]
  Land_Cover_Type_1_Percent_Class_10 <- occ_wide_clean[, paste0("Land_Cover_Type_1_Percent_Class_10_", YEARS)]
  Land_Cover_Type_1_Percent_Class_12 <- occ_wide_clean[, paste0("Land_Cover_Type_1_Percent_Class_12_", YEARS)]
  Land_Cover_Type_1_Percent_Class_13 <- occ_wide_clean[, paste0("Land_Cover_Type_1_Percent_Class_13_", YEARS)]
  Land_Cover_Type_1_Percent_Class_14 <- occ_wide_clean[, paste0("Land_Cover_Type_1_Percent_Class_14_", YEARS)]
  NDVI <- occ_wide_clean[, paste0("NDVI_", YEARS)]
  pr   <- occ_wide_clean[, paste0("pr_", YEARS)]
  tmmn <- occ_wide_clean[, paste0("tmmn_", YEARS)]
  tmmx <- occ_wide_clean[, paste0("tmmx_", YEARS)]

  # NEW IN V2: stepRep yearly site covariate (from step 3b)
  stepRep <- occ_wide_clean[, paste0(STEPREP_VARIANT, "_", YEARS)]

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
  stepRep <- scale(stepRep)                       # NEW IN V2

  # -- Capture training scaling for dynamic covariates --
  # scale() on a matrix gives per-column center/scale (one per year).
  # We store these keyed by "covariate_year" so the simulation section
  # can apply the EXACT SAME scaling as training -- not a per-year re-scaling.
  train_dyn_scale <- list()
  .dyn_mats <- list(
    EVI = EVI, NDVI = NDVI, pr = pr, tmmn = tmmn, tmmx = tmmx,
    Land_Cover_Type_1_Percent_Class_0  = Land_Cover_Type_1_Percent_Class_0,
    Land_Cover_Type_1_Percent_Class_6  = Land_Cover_Type_1_Percent_Class_6,
    Land_Cover_Type_1_Percent_Class_7  = Land_Cover_Type_1_Percent_Class_7,
    Land_Cover_Type_1_Percent_Class_10 = Land_Cover_Type_1_Percent_Class_10,
    Land_Cover_Type_1_Percent_Class_12 = Land_Cover_Type_1_Percent_Class_12,
    Land_Cover_Type_1_Percent_Class_13 = Land_Cover_Type_1_Percent_Class_13,
    Land_Cover_Type_1_Percent_Class_14 = Land_Cover_Type_1_Percent_Class_14,
    stepRep = stepRep                              # NEW IN V2
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

  # Save training scaling parameters for use by attribution scripts (8, 10)
  saveRDS(train_dyn_scale,
          here("results", paste0(sp, OUT_SUFFIX, "_train_dyn_scale.rds")))

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
  stepRep_obs     <- expand_matrix(stepRep, J_reps)   # NEW IN V2

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
      NDVI = NDVI, pr = pr, tmmn = tmmn, tmmx = tmmx,
      stepRep = stepRep                                # NEW IN V2
    ),
    obsCovs = list(
      duration = duration, effort = effort, observers = observers, time = time,
      NDVI_obs = NDVI_obs, pr_obs = pr_obs,
      topo_aspect_obs = topo_aspect_obs, topo_elev_obs = topo_elev_obs,
      stepRep_obs = stepRep_obs                        # NEW IN V2
    ),
    numPrimary = T_years
  )

  ##############################################################################
  # FIT FINAL MODEL
  ##############################################################################
  message("  Fitting colext model (v2 with stepRep_obs in detection)...")

  Mod.final <- colext(
    psiformula     = cfg$psi_formula,
    gammaformula   = cfg$gamma_formula,
    epsilonformula = cfg$epsilon_formula,
    pformula       = cfg$p_formula,
    data = occ_umf
  )

  # Save model summary
  sink(here("results", paste0(sp, OUT_SUFFIX, "_model_summary.txt")))
  cat("Species:", sp, "  (v2: detection includes stepRep_obs)\n")
  cat("Date:", as.character(Sys.time()), "\n")
  cat("p_formula:", paste(deparse(cfg$p_formula), collapse = " "), "\n\n")
  summary(Mod.final)
  sink()

  # Save model object
  saveRDS(Mod.final, here("results", paste0(sp, OUT_SUFFIX, "_model_object.rds")))
  message("  Model saved. AIC = ", round(Mod.final@AIC, 2))

  ##############################################################################
  # GOODNESS OF FIT (publication-grade nsim)
  ##############################################################################
  message("  Running GOF (parboot nsim=", NSIM_PARBOOT, ")...")
  set.seed(42)  # Reproducible GOF
  GOF <- parboot(Mod.final, nsim = NSIM_PARBOOT)
  saveRDS(GOF, here("results", paste0(sp, OUT_SUFFIX, "_gof_parboot.rds")))

  message("  Running MacKenzie-Bailey GOF (nsim=", NSIM_MB_GOF, ")...")
  set.seed(42)
  occ_gof <- mb.gof.test(Mod.final, nsim = NSIM_MB_GOF, plot.hist = FALSE)
  saveRDS(occ_gof, here("results", paste0(sp, OUT_SUFFIX,
                                           "_gof_mackenzie_bailey.rds")))

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
  ggsave(here("figs", paste0(sp, OUT_SUFFIX, "_response_occupancy.png")),
         p_occ, width = 10, height = 6)

  # Colonization
  col_formulas <- lapply(cfg$gamma_vars, function(v) as.formula(paste0("~", v)))
  p_col <- plot_response_curves(col_formulas, occ_umf, "col",
                                 as.data.frame(occ_umf@yearlySiteCovs),
                                 "Colonization Probability",
                                 "royalblue3", sp)
  ggsave(here("figs", paste0(sp, OUT_SUFFIX, "_response_colonization.png")),
         p_col, width = 10, height = 6)

  # Extinction
  ext_formulas <- lapply(cfg$epsilon_vars, function(v) as.formula(paste0("~", v)))
  p_ext <- plot_response_curves(ext_formulas, occ_umf, "ext",
                                 as.data.frame(occ_umf@yearlySiteCovs),
                                 "Extinction Probability",
                                 "royalblue3", sp)
  ggsave(here("figs", paste0(sp, OUT_SUFFIX, "_response_extinction.png")),
         p_ext, width = 10, height = 6)

  # Detection (now includes stepRep_obs in cfg$p_formula)
  det_var_names <- all.vars(cfg$p_formula)
  det_formulas <- lapply(det_var_names, function(v) as.formula(paste0("~", v)))
  p_det <- plot_response_curves(det_formulas, occ_umf, "det",
                                 as.data.frame(occ_umf@obsCovs),
                                 "Detection Probability",
                                 "royalblue3", sp)
  ggsave(here("figs", paste0(sp, OUT_SUFFIX, "_response_detection.png")),
         p_det, width = 10, height = 6)

  ##############################################################################
  # OCCUPANCY MAP (spatial prediction)
  ##############################################################################
  message("  Building occupancy map...")

  # Load rasters for prediction surface
  variables <- brick(stack(
    here("data-raw", "data", "environmental_data", "environmental_data_occ",
         "variables_spain.grd")))
  aspect <- raster(here("data-raw", "data", "topology_data", "topo_aspect_masked.tif"))
  names(aspect) <- "topo_aspect"
  elev <- raster(here("data-raw", "data", "topology_data", "topo_elev_masked.tif"))
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
    here("data", "processed_2023", sp, paste0(sp, "_scaling_params.rds"))
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
            here("data", "processed_2023", sp,
                 paste0("occ_", sp, OUT_SUFFIX, "_prediction.csv")),
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
              here("data", "processed_2023", sp,
                   paste0(sp, OUT_SUFFIX, "_OccuMap.tif")),
              format = "GTiff", overwrite = TRUE)

  # Plot occupancy map
  r_pred_df <- as.data.frame(r_pred_proj_masked[[1]], xy = TRUE, na.rm = TRUE)
  names(r_pred_df)[3] <- "occ_prob"

  p_map <- ggplot() +
    geom_tile(data = r_pred_df,
              aes(x = x, y = y, fill = occ_prob)) +
    scale_fill_viridis_c(name = "p(occupancy)") +
    labs(x = "Longitude", y = "Latitude",
         title = paste0("Occupancy probability: ", sp, " (v2 + stepRep)")) +
    theme(panel.border = element_rect(color = "black", fill = "transparent"))
  ggsave(here("figs", paste0(sp, OUT_SUFFIX, "_occupancy_map.png")),
         p_map, width = 10, height = 8)

  ##############################################################################
  # STOCHASTIC SIMULATIONS
  #
  # IMPORTANT: All covariates must be scaled using the TRAINING scaling
  # parameters, not re-scaled from scratch. Three bugs fixed in v1 (psi
  # covariates unscaled, dynamic covariates re-scaled per year, Land Cover
  # classes not scaled in training) carry over here.
  #
  # Note: simulation predicts col/ext per year, not p; stepRep does not
  # appear in cfg$gamma_vars or cfg$epsilon_vars, so no change here.
  ##############################################################################
  message("  Running stochastic simulations (nsim=", NSIM_PREV, ")...")

  set.seed(123)  # reproducible simulations

  predict_data <- occ_wide_clean %>%
    drop_na(paste0("NDVI_", YEARS))

  S <- nrow(predict_data)
  T_sim <- occ_umf@numPrimary

  # --- Initial occupancy (psi_1) ---
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
        warning("No training scaling for: ", key, " -- using raw values")
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
        warning("No training scaling for: ", key, " -- using raw values")
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
            here("results", paste0(sp, OUT_SUFFIX, "_simulation_prevalence.csv")),
            row.names = FALSE)

  # Plot prevalence over time
  p_prev <- ggplot(sim_results, aes(x = year, y = mean_occ, group = 1)) +
    geom_errorbar(aes(ymin = mean_occ - sd_occ, ymax = mean_occ + sd_occ),
                  width = 0.1) +
    geom_line() + geom_point() +
    labs(x = "Year", y = "Predicted mean occupancy",
         title = paste0("Occupancy trend: ", sp, " (v2 + stepRep)"))
  ggsave(here("figs", paste0(sp, OUT_SUFFIX, "_prevalence_over_time.png")),
         p_prev, width = 8, height = 5)

  message("  Step 4 v2 complete for ", sp)
}

message("\nStep 4 v2 complete: all occupancy models fitted with stepRep_obs",
        " in the detection sub-model.")
message("Outputs tagged with '", OUT_SUFFIX,
        "' to preserve v1 fits for comparison.")
