###############################################################################
# 9_collinearity_ndvi.R
#
# Purpose: Collinearity diagnostics for dynamic covariates, NDVI decomposition
#          into climate vs land-management components, and revised three-way
#          counterfactual attribution analysis.
#
# Inputs:  results/{sp}_model_object.rds               (fitted colext models)
#          data/processed_2023/{sp}/{sp}_occ_wide_dynamic.csv (site-year data)
#          R/model_configs.R                             (model formulas)
#
# Outputs: results/diagnostics/vif_summary.csv
#          results/diagnostics/coefficient_stability.csv
#          results/diagnostics/ndvi_climate_r2.tif
#          results/diagnostics/ndvi_climate_component.rds
#          results/diagnostics/ndvi_residual_component.rds
#          results/diagnostics/ndvi_model_revision_recommendations.csv
#          results/diagnostics/ndvi_revision_summary.txt
#          results/pub_table_attribution_revised.csv
#          figs/diagnostics/collinearity_heatmap_{sp}.png
#          figs/diagnostics/fig_ndvi_climate_r2_map.png
#          figs/fig_attribution_revised_summary_{sp}.png
#          figs/pub_fig_collinearity_diagnostics.png
#
# Structure: Tasks 1-5 matching the analysis specification.
###############################################################################

# -- Load packages --
library(here)
library(unmarked)
library(dplyr)
library(tidyr)
library(ggplot2)
library(terra)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(corrplot)
library(patchwork)

select  <- dplyr::select
filter  <- dplyr::filter
theme_set(theme_bw())
set.seed(42)

# Load model configurations
source(here("R", "model_configs.R"))

###############################################################################
# GLOBAL DEFINITIONS
###############################################################################

SPECIES <- list(
  otitar = "Otis tarda",
  ptealc = "Pterocles alchata",
  pteori = "Pterocles orientalis",
  tettet = "Tetrax tetrax"
)

YEARS   <- 2017:2023
T_YEARS <- length(YEARS)

# Covariate classification
CLIMATE_COVS <- c("NDVI", "pr", "tmmn", "tmmx")
LANDUSE_COVS <- c("Land_Cover_Type_1_Percent_Class_0",
                   "Land_Cover_Type_1_Percent_Class_6",
                   "Land_Cover_Type_1_Percent_Class_7",
                   "Land_Cover_Type_1_Percent_Class_10",
                   "Land_Cover_Type_1_Percent_Class_12",
                   "Land_Cover_Type_1_Percent_Class_13",
                   "Land_Cover_Type_1_Percent_Class_14")

# Helper: apply training scaling
scale_with_params <- function(x, center, scale_val) {
  (x - center) / scale_val
}

# Session tracking
session_log   <- list()
output_files  <- character(0)
track_file    <- function(f) { output_files <<- c(output_files, f); f }

# Spain basemap
spain <- ne_countries(country = "spain", scale = "medium", returnclass = "sf")
spain_crop <- st_crop(spain, xmin = -10, xmax = 5, ymin = 35, ymax = 44)

###############################################################################
# INPUT FILE VERIFICATION
###############################################################################

message("\n", strrep("=", 70))
message("  INPUT FILE VERIFICATION")
message(strrep("=", 70))

all_ok <- TRUE
for (sp in names(SPECIES)) {
  f_mod  <- here("results", paste0(sp, "_model_object.rds"))
  f_data <- here("data", "processed_2023", sp, paste0(sp, "_occ_wide_dynamic.csv"))
  for (f in c(f_mod, f_data)) {
    ok <- file.exists(f)
    message(ifelse(ok, "  [OK] ", "  [MISSING] "), basename(f))
    if (!ok) all_ok <- FALSE
  }
}
f_cfg <- here("R", "model_configs.R")
message(ifelse(file.exists(f_cfg), "  [OK] ", "  [MISSING] "), "model_configs.R")
if (!file.exists(f_cfg)) all_ok <- FALSE

if (!all_ok) stop("Missing required input files. See checklist above.")
message("  All required inputs present.\n")

# Ensure output directories
dir.create(here("results", "diagnostics"), showWarnings = FALSE, recursive = TRUE)
dir.create(here("figs", "diagnostics"), showWarnings = FALSE, recursive = TRUE)

###############################################################################
# HELPER: Load species data (adapted from script 8)
###############################################################################

load_species_data <- function(sp) {

  occ_wide <- read.csv(
    here("data", "processed_2023", sp, paste0(sp, "_occ_wide_dynamic.csv"))
  )

  # Drop sites with missing NDVI (same filter as script 4)
  predict_data <- occ_wide %>% drop_na(paste0("NDVI_", YEARS))

  # Extract raw dynamic covariate matrices and compute training scaling
  all_covs <- unique(c(CLIMATE_COVS, LANDUSE_COVS))
  raw_matrices <- list()
  train_scale  <- list()

  for (cov in all_covs) {
    cols <- paste0(cov, "_", YEARS)
    if (all(cols %in% names(predict_data))) {
      mat <- as.matrix(predict_data[, cols])
      scaled_mat <- scale(mat)
      ctrs <- attr(scaled_mat, "scaled:center")
      scls <- attr(scaled_mat, "scaled:scale")
      raw_matrices[[cov]] <- mat
      for (yi in seq_along(YEARS)) {
        train_scale[[paste0(cov, "_", YEARS[yi])]] <- list(
          center = ctrs[yi], scale = scls[yi]
        )
      }
    }
  }

  # Site coordinates
  coords <- predict_data[, c("latitude", "longitude")]

  # Detection histories (for occupancy state estimation)
  y_cols <- c()
  for (yr in YEARS) {
    y_cols <- c(y_cols, paste0("y.", 1:10, ".", yr))
  }
  y_cols_present <- intersect(y_cols, names(predict_data))
  y_matrix <- as.matrix(predict_data[, y_cols_present])

  list(
    predict_data = predict_data,
    raw_matrices = raw_matrices,
    train_scale  = train_scale,
    coords       = coords,
    S            = nrow(predict_data),
    y_matrix     = y_matrix
  )
}

###############################################################################
# HELPER: Estimate naive occupancy state from detection history
###############################################################################

get_naive_occupancy <- function(sp_data) {
  # Returns S x T matrix: 1 if detected at least once in year t, 0 otherwise
  S <- sp_data$S
  y <- sp_data$y_matrix
  J <- 10  # secondary periods
  occ_state <- matrix(0, nrow = S, ncol = T_YEARS)

  for (t in seq_len(T_YEARS)) {
    col_start <- (t - 1) * J + 1
    col_end   <- t * J
    if (col_end <= ncol(y)) {
      year_obs <- y[, col_start:col_end]
      occ_state[, t] <- apply(year_obs, 1, function(row) {
        if (all(is.na(row))) return(NA)
        as.integer(any(row == 1, na.rm = TRUE))
      })
    }
  }
  occ_state
}

###############################################################################
# HELPER: Build scenario newdata (adapted from script 8)
###############################################################################

build_scenario_newdata <- function(sp_data, vars, scenario_logic, year_idx) {
  # scenario_logic: named list mapping variable name -> "observed" or "frozen"
  yr <- YEARS[year_idx]
  S  <- sp_data$S
  newdata <- data.frame(row.names = seq_len(S))

  for (v in vars) {
    key      <- paste0(v, "_", yr)
    key_2017 <- paste0(v, "_", YEARS[1])

    use_frozen <- isTRUE(scenario_logic[[v]] == "frozen")

    if (use_frozen) {
      raw_val <- sp_data$raw_matrices[[v]][, 1]
      sc      <- sp_data$train_scale[[key_2017]]
    } else {
      raw_val <- sp_data$raw_matrices[[v]][, year_idx]
      sc      <- sp_data$train_scale[[key]]
    }

    if (!is.null(sc)) {
      newdata[[v]] <- scale_with_params(raw_val, sc$center, sc$scale)
    } else {
      newdata[[v]] <- raw_val
    }
  }
  newdata
}

###############################################################################
#                          TASK 1: COLLINEARITY DIAGNOSTICS
###############################################################################

message("\n", strrep("#", 70))
message("# TASK 1: COLLINEARITY DIAGNOSTICS FOR DYNAMIC COVARIATES")
message(strrep("#", 70))

# Storage for results
all_corr_matrices  <- list()
all_vif_results    <- data.frame()
all_corr_flags     <- data.frame()

# The six focal pairs to report explicitly
FOCAL_PAIRS <- list(
  c("NDVI", "pr"), c("NDVI", "tmmx"), c("NDVI", "tmmn"),
  c("pr", "tmmx"), c("pr", "tmmn"), c("tmmx", "tmmn")
)

for (sp in names(SPECIES)) {

  message("\n--- ", sp, " (", SPECIES[[sp]], ") ---")

  cfg     <- get_model_config(sp)
  sp_data <- load_species_data(sp)

  #---------------------------------------------------------------------------
  # 1.1 Pairwise correlations
  #---------------------------------------------------------------------------

  # All dynamic covariates used in ANY submodel for this species
  all_dyn_vars <- unique(c(cfg$gamma_vars, cfg$epsilon_vars))
  message("  Dynamic covariates: ", paste(all_dyn_vars, collapse = ", "))

  # Build long-form covariate matrix: stack all years (S*T rows x ncov cols)
  cov_long <- data.frame(row.names = seq_len(sp_data$S * T_YEARS))
  for (v in all_dyn_vars) {
    if (!is.null(sp_data$raw_matrices[[v]])) {
      cov_long[[v]] <- as.vector(sp_data$raw_matrices[[v]])
    }
  }
  # Remove columns that are all NA or constant
  keep_cols <- sapply(cov_long, function(x) !all(is.na(x)) && sd(x, na.rm = TRUE) > 0)
  cov_long  <- cov_long[, keep_cols, drop = FALSE]

  # Correlation matrix
  corr_mat <- cor(cov_long, use = "pairwise.complete.obs")
  all_corr_matrices[[sp]] <- corr_mat

  # Shorten names for display
  short_names <- gsub("Land_Cover_Type_1_Percent_", "LC", colnames(corr_mat))
  colnames(corr_mat) <- short_names
  rownames(corr_mat) <- short_names

  # Heatmap
  heatmap_path <- here("figs", "diagnostics",
                        paste0("collinearity_heatmap_", sp, ".png"))
  png(heatmap_path, width = 8, height = 7, units = "in", res = 300)
  corrplot(corr_mat, method = "color", type = "lower",
           addCoef.col = "black", number.cex = 0.7, tl.cex = 0.8,
           title = paste0(SPECIES[[sp]], " — Dynamic Covariate Correlations"),
           mar = c(0, 0, 2, 0))
  dev.off()
  track_file(heatmap_path)
  message("  Heatmap saved: ", basename(heatmap_path))

  # Flag collinear pairs
  corr_full <- cor(cov_long, use = "pairwise.complete.obs")
  for (i in seq_len(ncol(corr_full))) {
    for (j in seq_len(i - 1)) {
      r_val <- corr_full[i, j]
      flag  <- ifelse(abs(r_val) > 0.7, "SEVERE",
                      ifelse(abs(r_val) > 0.5, "MODERATE", "OK"))
      all_corr_flags <- rbind(all_corr_flags, data.frame(
        species  = sp,
        var1     = colnames(corr_full)[j],
        var2     = colnames(corr_full)[i],
        r        = round(r_val, 4),
        flag     = flag,
        stringsAsFactors = FALSE
      ))
    }
  }

  # Report focal pairs
  message("  Focal pair correlations:")
  for (pair in FOCAL_PAIRS) {
    if (all(pair %in% colnames(corr_full))) {
      r_val <- corr_full[pair[1], pair[2]]
      flag <- ifelse(abs(r_val) > 0.7, " ** SEVERE **",
                     ifelse(abs(r_val) > 0.5, " * MODERATE *", ""))
      message("    r(", pair[1], ", ", pair[2], ") = ",
              round(r_val, 3), flag)
    }
  }

  #---------------------------------------------------------------------------
  # 1.2 Variance Inflation Factors
  #---------------------------------------------------------------------------

  # Get naive occupancy states for subsetting

  occ_state <- get_naive_occupancy(sp_data)

  for (submodel in c("gamma", "epsilon")) {
    vars <- if (submodel == "gamma") cfg$gamma_vars else cfg$epsilon_vars

    if (length(vars) < 2) {
      message("  VIF: ", submodel, " has < 2 covariates, skipping.")
      for (v in vars) {
        all_vif_results <- rbind(all_vif_results, data.frame(
          species = sp, submodel = submodel, covariate = v,
          VIF = 1.0, flag = "OK", stringsAsFactors = FALSE
        ))
      }
      next
    }

    # Build covariate matrix for relevant transitions
    # Colonisation: unoccupied in year t, any state in t+1 (transitions 0->0, 0->1)
    # Extinction: occupied in year t, any state in t+1 (transitions 1->0, 1->1)
    trans_rows <- c()
    trans_year_idx <- c()
    for (t in seq_len(T_YEARS - 1)) {
      if (submodel == "gamma") {
        # Unoccupied at time t
        eligible <- which(!is.na(occ_state[, t]) & occ_state[, t] == 0)
      } else {
        # Occupied at time t
        eligible <- which(!is.na(occ_state[, t]) & occ_state[, t] == 1)
      }
      if (length(eligible) > 0) {
        trans_rows     <- c(trans_rows, eligible)
        trans_year_idx <- c(trans_year_idx, rep(t, length(eligible)))
      }
    }

    if (length(trans_rows) < 10) {
      message("  VIF: ", submodel, " has < 10 transition rows, skipping.")
      for (v in vars) {
        all_vif_results <- rbind(all_vif_results, data.frame(
          species = sp, submodel = submodel, covariate = v,
          VIF = NA, flag = "INSUFFICIENT_DATA", stringsAsFactors = FALSE
        ))
      }
      next
    }

    # Extract covariates at those transitions
    trans_df <- data.frame(row.names = seq_along(trans_rows))
    for (v in vars) {
      vals <- numeric(length(trans_rows))
      for (k in seq_along(trans_rows)) {
        vals[k] <- sp_data$raw_matrices[[v]][trans_rows[k], trans_year_idx[k]]
      }
      trans_df[[v]] <- vals
    }
    trans_df <- na.omit(trans_df)

    if (nrow(trans_df) < length(vars) + 2) {
      message("  VIF: ", submodel, " insufficient data after NA removal.")
      for (v in vars) {
        all_vif_results <- rbind(all_vif_results, data.frame(
          species = sp, submodel = submodel, covariate = v,
          VIF = NA, flag = "INSUFFICIENT_DATA", stringsAsFactors = FALSE
        ))
      }
      next
    }

    # Compute VIF via R² method
    for (v in vars) {
      other_vars <- setdiff(vars, v)
      fml <- as.formula(paste0("`", v, "` ~ ",
                                paste0("`", other_vars, "`", collapse = " + ")))
      lm_fit <- tryCatch(lm(fml, data = trans_df), error = function(e) NULL)
      if (is.null(lm_fit)) {
        vif_val <- NA; flag <- "FIT_ERROR"
      } else {
        r2 <- summary(lm_fit)$r.squared
        vif_val <- 1 / (1 - r2)
        flag <- ifelse(vif_val > 10, "SEVERE",
                       ifelse(vif_val > 5, "PROBLEMATIC", "OK"))
      }
      all_vif_results <- rbind(all_vif_results, data.frame(
        species = sp, submodel = submodel, covariate = v,
        VIF = round(vif_val, 2), flag = flag, stringsAsFactors = FALSE
      ))
      if (!is.na(vif_val) && vif_val > 5) {
        message("  ** VIF FLAG: ", sp, "/", submodel, "/", v,
                " VIF=", round(vif_val, 2), " (", flag, ")")
      }
    }
    message("  VIF computed for ", sp, "/", submodel, ": ", nrow(trans_df),
            " transition rows, ", length(vars), " covariates.")
  }
}

# Save correlation flags
corr_flags_path <- here("results", "diagnostics", "correlation_flags.csv")
write.csv(all_corr_flags, track_file(corr_flags_path), row.names = FALSE)

# Save VIF results
vif_path <- here("results", "diagnostics", "vif_summary.csv")
write.csv(all_vif_results, track_file(vif_path), row.names = FALSE)
message("\n  VIF summary saved: ", vif_path)

# Focal pairs summary table
focal_summary <- data.frame()
for (sp in names(SPECIES)) {
  corr_full <- all_corr_matrices[[sp]]
  for (pair in FOCAL_PAIRS) {
    if (all(pair %in% colnames(corr_full))) {
      r_val <- corr_full[pair[1], pair[2]]
      focal_summary <- rbind(focal_summary, data.frame(
        species = sp, var1 = pair[1], var2 = pair[2],
        r = round(r_val, 4),
        flag = ifelse(abs(r_val) > 0.7, "SEVERE",
                      ifelse(abs(r_val) > 0.5, "MODERATE", "OK")),
        stringsAsFactors = FALSE
      ))
    }
  }
}
focal_path <- here("results", "diagnostics", "focal_pair_correlations.csv")
write.csv(focal_summary, track_file(focal_path), row.names = FALSE)
message("  Focal pair correlations saved: ", focal_path)

###############################################################################
# 1.3 Coefficient stability test
###############################################################################

message("\n", strrep("-", 70))
message("  TASK 1.3: Coefficient stability test (NDVI removal)")
message(strrep("-", 70))

stability_results <- data.frame()

for (sp in names(SPECIES)) {

  cfg <- get_model_config(sp)
  mod <- readRDS(here("results", paste0(sp, "_model_object.rds")))
  orig_aic <- mod@AIC

  # Check which submodels contain NDVI
  has_ndvi_gamma   <- "NDVI" %in% cfg$gamma_vars
  has_ndvi_epsilon <- "NDVI" %in% cfg$epsilon_vars

  # Check if collinearity was flagged
  sp_vif <- all_vif_results %>% filter(species == sp)
  sp_corr <- all_corr_flags %>% filter(species == sp)
  flagged_gamma <- any(sp_vif$submodel == "gamma" &
                       (sp_vif$flag %in% c("PROBLEMATIC", "SEVERE"))) ||
    any(sp_corr$flag %in% c("MODERATE", "SEVERE") &
        (sp_corr$var1 == "NDVI" | sp_corr$var2 == "NDVI"))
  flagged_epsilon <- any(sp_vif$submodel == "epsilon" &
                         (sp_vif$flag %in% c("PROBLEMATIC", "SEVERE"))) ||
    any(sp_corr$flag %in% c("MODERATE", "SEVERE") &
        (sp_corr$var1 == "NDVI" | sp_corr$var2 == "NDVI"))

  for (submodel in c("gamma", "epsilon")) {
    has_ndvi <- if (submodel == "gamma") has_ndvi_gamma else has_ndvi_epsilon
    flagged  <- if (submodel == "gamma") flagged_gamma else flagged_epsilon

    if (!has_ndvi) {
      message("  ", sp, "/", submodel, ": no NDVI in submodel, skipping.")
      next
    }

    if (!flagged) {
      message("  ", sp, "/", submodel,
              ": no collinearity flagged, running stability test anyway.")
    }

    # Build formula without NDVI
    orig_vars <- if (submodel == "gamma") cfg$gamma_vars else cfg$epsilon_vars
    new_vars  <- setdiff(orig_vars, "NDVI")

    if (length(new_vars) == 0) {
      new_formula <- ~ 1
    } else {
      new_formula <- as.formula(paste0("~ ", paste(new_vars, collapse = " + ")))
    }

    # Refit model with NDVI removed from this submodel only
    message("  Refitting ", sp, " without NDVI in ", submodel, "...")
    refit_args <- list(
      psiformula     = cfg$psi_formula,
      gammaformula   = if (submodel == "gamma") new_formula else cfg$gamma_formula,
      epsilonformula = if (submodel == "epsilon") new_formula else cfg$epsilon_formula,
      pformula       = cfg$p_formula,
      data           = mod@data
    )

    refit_mod <- tryCatch(
      do.call(colext, refit_args),
      error = function(e) {
        message("    REFIT FAILED: ", e$message)
        NULL
      }
    )

    if (is.null(refit_mod)) {
      stability_results <- rbind(stability_results, data.frame(
        species = sp, submodel = submodel, covariate = "ALL",
        original_beta = NA, original_se = NA,
        refit_beta = NA, refit_se = NA,
        beta_shift_se = NA, se_change_pct = NA,
        sign_change = NA, aic_original = orig_aic, aic_refit = NA,
        delta_aic = NA, flag = "REFIT_FAILED", stringsAsFactors = FALSE
      ))
      next
    }

    # Save refitted model
    refit_path <- here("results", "diagnostics",
                       paste0("models_ndvi_excluded_", sp, "_", submodel, ".rds"))
    saveRDS(refit_mod, track_file(refit_path))

    # Compare coefficients for remaining covariates
    orig_coefs <- coef(mod)
    refit_coefs <- coef(refit_mod)
    orig_se <- tryCatch(sqrt(diag(vcov(mod))), error = function(e) rep(NA, length(orig_coefs)))
    refit_se <- tryCatch(sqrt(diag(vcov(refit_mod))), error = function(e) rep(NA, length(refit_coefs)))

    # Match coefficients by name in the relevant submodel
    submodel_prefix <- if (submodel == "gamma") "col" else "ext"

    for (v in new_vars) {
      # Find coefficient name in original model
      orig_name <- paste0(submodel_prefix, "(", v, ")")
      refit_name <- orig_name  # same name in refitted model

      orig_idx  <- match(orig_name, names(orig_coefs))
      refit_idx <- match(refit_name, names(refit_coefs))

      if (!is.na(orig_idx) && !is.na(refit_idx)) {
        ob <- orig_coefs[orig_idx]
        os <- orig_se[orig_idx]
        rb <- refit_coefs[refit_idx]
        rs <- refit_se[refit_idx]

        beta_shift <- abs(rb - ob) / os
        se_change  <- ((rs - os) / os) * 100
        sign_chg   <- sign(ob) != sign(rb) & abs(ob) > 0.01 & abs(rb) > 0.01

        flag <- "OK"
        if (!is.na(sign_chg) && sign_chg) flag <- "SIGN_CHANGE"
        else if (!is.na(beta_shift) && beta_shift > 1) flag <- "BETA_SHIFT>1SE"
        else if (!is.na(se_change) && se_change > 50) flag <- "SE_INCREASE>50%"

        stability_results <- rbind(stability_results, data.frame(
          species = sp, submodel = submodel, covariate = v,
          original_beta = round(ob, 4), original_se = round(os, 4),
          refit_beta = round(rb, 4), refit_se = round(rs, 4),
          beta_shift_se = round(beta_shift, 3),
          se_change_pct = round(se_change, 1),
          sign_change = sign_chg,
          aic_original = round(orig_aic, 2),
          aic_refit = round(refit_mod@AIC, 2),
          delta_aic = round(refit_mod@AIC - orig_aic, 2),
          flag = flag, stringsAsFactors = FALSE
        ))

        message("    ", v, ": beta ", round(ob, 3), " -> ", round(rb, 3),
                " (shift=", round(beta_shift, 2), " SE)",
                ifelse(sign_chg, " *** SIGN CHANGE ***", ""))
      }
    }

    message("    AIC: original=", round(orig_aic, 2),
            " refit=", round(refit_mod@AIC, 2),
            " delta=", round(refit_mod@AIC - orig_aic, 2))
  }
}

stab_path <- here("results", "diagnostics", "coefficient_stability.csv")
write.csv(stability_results, track_file(stab_path), row.names = FALSE)
message("\n  Coefficient stability saved: ", stab_path)

session_log[["task1"]] <- "COMPLETE"


###############################################################################
#                     TASK 2: NDVI DECOMPOSITION
###############################################################################

message("\n", strrep("#", 70))
message("# TASK 2: NDVI DECOMPOSITION — CLIMATE VS LAND-USE SIGNAL")
message(strrep("#", 70))

# We use otitar data as reference grid since all species share similar spatial
# coverage, but we compute decomposition once on the pooled grid.
# Actually, compute per-species since site sets may differ.

all_r2_data <- data.frame()
ndvi_climate_all <- list()
ndvi_residual_all <- list()

for (sp in names(SPECIES)) {

  message("\n--- ", sp, " ---")
  sp_data <- load_species_data(sp)

  ndvi_mat <- sp_data$raw_matrices[["NDVI"]]
  pr_mat   <- sp_data$raw_matrices[["pr"]]
  tmmx_mat <- sp_data$raw_matrices[["tmmx"]]
  tmmn_mat <- sp_data$raw_matrices[["tmmn"]]

  if (is.null(ndvi_mat) || is.null(pr_mat) || is.null(tmmx_mat) || is.null(tmmn_mat)) {
    message("  Missing climate matrices for ", sp, ", skipping NDVI decomposition.")
    next
  }

  S <- sp_data$S
  r2_vals        <- numeric(S)
  ndvi_climate   <- matrix(NA, nrow = S, ncol = T_YEARS)
  ndvi_residual  <- matrix(NA, nrow = S, ncol = T_YEARS)

  for (i in seq_len(S)) {
    df_site <- data.frame(
      ndvi = ndvi_mat[i, ],
      pr   = pr_mat[i, ],
      tmmx = tmmx_mat[i, ],
      tmmn = tmmn_mat[i, ]
    )

    if (any(is.na(df_site))) {
      r2_vals[i]        <- NA
      ndvi_climate[i, ] <- NA
      ndvi_residual[i, ] <- NA
      next
    }

    lm_fit <- tryCatch(lm(ndvi ~ pr + tmmx + tmmn, data = df_site),
                        error = function(e) NULL)

    if (is.null(lm_fit)) {
      r2_vals[i] <- NA
      next
    }

    r2_vals[i]         <- summary(lm_fit)$r.squared
    ndvi_climate[i, ]  <- fitted(lm_fit)
    ndvi_residual[i, ] <- residuals(lm_fit)
  }

  r2_data <- data.frame(
    species = sp,
    lat     = sp_data$coords$latitude,
    lon     = sp_data$coords$longitude,
    r2      = r2_vals,
    stringsAsFactors = FALSE
  )
  all_r2_data <- rbind(all_r2_data, r2_data)

  ndvi_climate_all[[sp]]  <- ndvi_climate
  ndvi_residual_all[[sp]] <- ndvi_residual

  # Summary stats
  r2_clean <- r2_vals[!is.na(r2_vals)]
  message("  R² summary: mean=", round(mean(r2_clean), 3),
          " SD=", round(sd(r2_clean), 3),
          " prop>0.5=", round(mean(r2_clean > 0.5), 3),
          " prop>0.7=", round(mean(r2_clean > 0.7), 3))
}

# Save decomposed components
saveRDS(ndvi_climate_all,
        track_file(here("results", "diagnostics", "ndvi_climate_component.rds")))
saveRDS(ndvi_residual_all,
        track_file(here("results", "diagnostics", "ndvi_residual_component.rds")))

#---------------------------------------------------------------------------
# 2.2 R² map
#---------------------------------------------------------------------------

message("\n  Creating R² map...")

# Create raster from the largest species dataset
r2_for_raster <- all_r2_data %>%
  group_by(lon, lat) %>%
  summarise(r2 = mean(r2, na.rm = TRUE), .groups = "drop")

r2_vect <- vect(r2_for_raster, geom = c("lon", "lat"), crs = "EPSG:4326")
r_template <- rast(ext(r2_vect), res = 0.05, crs = "EPSG:4326")
r2_rast <- rasterize(r2_vect, r_template, field = "r2", fun = mean)
writeRaster(r2_rast,
            track_file(here("results", "diagnostics", "ndvi_climate_r2.tif")),
            overwrite = TRUE)

# Map per species (faceted)
p_r2_map <- ggplot() +
  geom_sf(data = spain_crop, fill = "grey95", colour = "grey50") +
  geom_point(data = all_r2_data %>% filter(!is.na(r2)),
             aes(x = lon, y = lat, colour = r2), size = 0.4, alpha = 0.8) +
  scale_colour_viridis_c(name = expression(R^2), option = "inferno",
                          limits = c(0, 1)) +
  facet_wrap(~ species, ncol = 2,
             labeller = labeller(species = unlist(SPECIES))) +
  coord_sf(xlim = c(-10, 5), ylim = c(35, 44)) +
  labs(title = expression(paste("NDVI ~ climate ", R^2, " per site")),
       subtitle = "High = NDVI tracks climate; Low = land-management signal",
       x = "Longitude", y = "Latitude") +
  theme(strip.text = element_text(face = "italic"))

ggsave(track_file(here("figs", "diagnostics", "fig_ndvi_climate_r2_map.png")),
       p_r2_map, width = 12, height = 10, dpi = 300)
message("  R² map saved.")

session_log[["task2"]] <- "COMPLETE"


###############################################################################
#            TASK 3: REVISED ATTRIBUTION WITH THREE-WAY DECOMPOSITION
###############################################################################

message("\n", strrep("#", 70))
message("# TASK 3: REVISED ATTRIBUTION (THREE-WAY DECOMPOSITION)")
message(strrep("#", 70))

# Revised scenario design:
# S0: all frozen (null baseline)
# S1: climate only (pr, tmmn, tmmx observed; NDVI_climate frozen; NDVI_residual frozen; LC frozen)
# S2: NDVI_climate only (pr, tmmn, tmmx frozen; NDVI_climate observed; NDVI_residual frozen; LC frozen)
# S3: NDVI_landuse only (pr, tmmn, tmmx frozen; NDVI_climate frozen; NDVI_residual observed; LC frozen)
# S4: land use only (pr, tmmn, tmmx frozen; NDVI_climate frozen; NDVI_residual frozen; LC observed)
# S5: combined (all observed)

# For models where NDVI appears, we need to reconstruct it from components:
# NDVI_observed = NDVI_climate + NDVI_residual (additive decomposition)
# For scenarios: NDVI = NDVI_climate(scenario) + NDVI_residual(scenario)

SCENARIOS_REV <- list(
  S0 = list(name = "Null baseline",
            pr = "frozen", tmmn = "frozen", tmmx = "frozen",
            ndvi_clim = "frozen", ndvi_resid = "frozen", landuse = "frozen"),
  S1 = list(name = "Climate only",
            pr = "observed", tmmn = "observed", tmmx = "observed",
            ndvi_clim = "frozen", ndvi_resid = "frozen", landuse = "frozen"),
  S2 = list(name = "NDVI climate",
            pr = "frozen", tmmn = "frozen", tmmx = "frozen",
            ndvi_clim = "observed", ndvi_resid = "frozen", landuse = "frozen"),
  S3 = list(name = "NDVI landuse",
            pr = "frozen", tmmn = "frozen", tmmx = "frozen",
            ndvi_clim = "frozen", ndvi_resid = "observed", landuse = "frozen"),
  S4 = list(name = "Land use only",
            pr = "frozen", tmmn = "frozen", tmmx = "frozen",
            ndvi_clim = "frozen", ndvi_resid = "frozen", landuse = "observed"),
  S5 = list(name = "Combined",
            pr = "observed", tmmn = "observed", tmmx = "observed",
            ndvi_clim = "observed", ndvi_resid = "observed", landuse = "observed")
)

# Helper: build newdata for revised scenarios
build_revised_newdata <- function(sp_data, vars, scenario, year_idx,
                                   ndvi_clim_mat, ndvi_resid_mat) {
  yr <- YEARS[year_idx]
  S  <- sp_data$S
  newdata <- data.frame(row.names = seq_len(S))

  for (v in vars) {
    key      <- paste0(v, "_", yr)
    key_2017 <- paste0(v, "_", YEARS[1])

    if (v == "NDVI") {
      # Reconstruct NDVI from decomposed components
      # NDVI_climate component
      if (scenario$ndvi_clim == "frozen") {
        clim_part <- ndvi_clim_mat[, 1]  # 2017
      } else {
        clim_part <- ndvi_clim_mat[, year_idx]
      }
      # NDVI_residual component
      if (scenario$ndvi_resid == "frozen") {
        resid_part <- ndvi_resid_mat[, 1]  # 2017
      } else {
        resid_part <- ndvi_resid_mat[, year_idx]
      }
      # Reconstructed NDVI
      raw_val <- clim_part + resid_part

      # Scale with appropriate year's params
      # For frozen scenario: use 2017 scaling (since values approximate 2017 level)
      # For observed: use that year's scaling
      if (scenario$ndvi_clim == "frozen" && scenario$ndvi_resid == "frozen") {
        sc <- sp_data$train_scale[[key_2017]]
      } else {
        sc <- sp_data$train_scale[[key]]
      }

    } else if (v %in% c("pr", "tmmn", "tmmx")) {
      use_frozen <- (scenario[[v]] == "frozen")
      if (use_frozen) {
        raw_val <- sp_data$raw_matrices[[v]][, 1]
        sc      <- sp_data$train_scale[[key_2017]]
      } else {
        raw_val <- sp_data$raw_matrices[[v]][, year_idx]
        sc      <- sp_data$train_scale[[key]]
      }

    } else {
      # Land-use covariate
      use_frozen <- (scenario$landuse == "frozen")
      if (use_frozen) {
        raw_val <- sp_data$raw_matrices[[v]][, 1]
        sc      <- sp_data$train_scale[[key_2017]]
      } else {
        raw_val <- sp_data$raw_matrices[[v]][, year_idx]
        sc      <- sp_data$train_scale[[key]]
      }
    }

    if (!is.null(sc)) {
      newdata[[v]] <- scale_with_params(raw_val, sc$center, sc$scale)
    } else {
      newdata[[v]] <- raw_val
    }
  }
  newdata
}

# Main attribution loop
all_rev_attribution <- list()

for (sp in names(SPECIES)) {

  message("\n--- ", sp, " ---")

  mod     <- readRDS(here("results", paste0(sp, "_model_object.rds")))
  cfg     <- get_model_config(sp)
  sp_data <- load_species_data(sp)
  S       <- sp_data$S

  has_ndvi_gamma   <- "NDVI" %in% cfg$gamma_vars
  has_ndvi_epsilon <- "NDVI" %in% cfg$epsilon_vars

  # Get NDVI decomposition for this species
  ndvi_clim  <- ndvi_climate_all[[sp]]
  ndvi_resid <- ndvi_residual_all[[sp]]

  # Determine which scenarios to run
  # If NDVI not in any submodel, S2 and S3 are identical to S0
  scenarios_to_run <- c("S0", "S1", "S4", "S5")
  if (has_ndvi_gamma || has_ndvi_epsilon) {
    scenarios_to_run <- c("S0", "S1", "S2", "S3", "S4", "S5")
  }

  # Predictions for all scenarios
  preds <- list()
  for (sc_name in scenarios_to_run) {
    sc <- SCENARIOS_REV[[sc_name]]
    preds[[sc_name]] <- list()

    for (yi in seq_len(T_YEARS)) {
      # Colonisation predictions
      col_nd <- build_revised_newdata(sp_data, cfg$gamma_vars, sc, yi,
                                       ndvi_clim, ndvi_resid)
      col_pred <- tryCatch(
        predict(mod, type = "col", newdata = col_nd),
        error = function(e) NULL
      )

      # Extinction predictions
      ext_nd <- build_revised_newdata(sp_data, cfg$epsilon_vars, sc, yi,
                                       ndvi_clim, ndvi_resid)
      ext_pred <- tryCatch(
        predict(mod, type = "ext", newdata = ext_nd),
        error = function(e) NULL
      )

      preds[[sc_name]][[as.character(YEARS[yi])]] <- list(
        gamma   = if (!is.null(col_pred)) col_pred$Predicted else rep(NA, S),
        epsilon = if (!is.null(ext_pred)) ext_pred$Predicted else rep(NA, S)
      )
    }
  }

  # Attribution decomposition per site-year
  attr_df <- data.frame()

  for (yi in seq_len(T_YEARS)) {
    yr_char <- as.character(YEARS[yi])

    g0 <- preds[["S0"]][[yr_char]]$gamma
    g1 <- preds[["S1"]][[yr_char]]$gamma
    g4 <- preds[["S4"]][[yr_char]]$gamma
    g5 <- preds[["S5"]][[yr_char]]$gamma

    e0 <- preds[["S0"]][[yr_char]]$epsilon
    e1 <- preds[["S1"]][[yr_char]]$epsilon
    e4 <- preds[["S4"]][[yr_char]]$epsilon
    e5 <- preds[["S5"]][[yr_char]]$epsilon

    # NDVI decomposition (only if NDVI in model)
    if ("S2" %in% scenarios_to_run) {
      g2 <- preds[["S2"]][[yr_char]]$gamma
      g3 <- preds[["S3"]][[yr_char]]$gamma
      e2 <- preds[["S2"]][[yr_char]]$epsilon
      e3 <- preds[["S3"]][[yr_char]]$epsilon
    } else {
      g2 <- g0; g3 <- g0
      e2 <- e0; e3 <- e0
    }

    yr_attr <- data.frame(
      site_id = seq_len(S),
      lat     = sp_data$coords$latitude,
      lon     = sp_data$coords$longitude,
      year    = YEARS[yi],
      # Colonisation attribution
      delta_gamma_climate      = g1 - g0,
      delta_gamma_ndvi_climate = g2 - g0,
      delta_gamma_ndvi_landuse = g3 - g0,
      delta_gamma_landuse      = g4 - g0,
      delta_gamma_total        = g5 - g0,
      # Extinction attribution
      delta_epsilon_climate      = e1 - e0,
      delta_epsilon_ndvi_climate = e2 - e0,
      delta_epsilon_ndvi_landuse = e3 - e0,
      delta_epsilon_landuse      = e4 - e0,
      delta_epsilon_total        = e5 - e0,
      stringsAsFactors = FALSE
    )
    attr_df <- rbind(attr_df, yr_attr)
  }

  all_rev_attribution[[sp]] <- attr_df

  # Summary
  mean_effects <- attr_df %>%
    summarise(across(starts_with("delta_"), ~ mean(.x, na.rm = TRUE)))
  message("  Mean attribution effects:")
  for (nm in names(mean_effects)) {
    message("    ", nm, " = ", format(mean_effects[[nm]], digits = 3, scientific = TRUE))
  }
}

#---------------------------------------------------------------------------
# Revised cross-species summary table
#---------------------------------------------------------------------------

rev_summary <- data.frame()
for (sp in names(SPECIES)) {
  attr_df <- all_rev_attribution[[sp]]

  means <- attr_df %>%
    summarise(across(starts_with("delta_"), ~ mean(.x, na.rm = TRUE)))
  sds <- attr_df %>%
    summarise(across(starts_with("delta_"), ~ sd(.x, na.rm = TRUE)))

  rev_summary <- rbind(rev_summary, data.frame(
    species = sp,
    species_name = SPECIES[[sp]],
    # Colonisation
    delta_gamma_climate_mean      = means$delta_gamma_climate,
    delta_gamma_climate_sd        = sds$delta_gamma_climate,
    delta_gamma_ndvi_climate_mean = means$delta_gamma_ndvi_climate,
    delta_gamma_ndvi_climate_sd   = sds$delta_gamma_ndvi_climate,
    delta_gamma_ndvi_landuse_mean = means$delta_gamma_ndvi_landuse,
    delta_gamma_ndvi_landuse_sd   = sds$delta_gamma_ndvi_landuse,
    delta_gamma_landuse_mean      = means$delta_gamma_landuse,
    delta_gamma_landuse_sd        = sds$delta_gamma_landuse,
    delta_gamma_total_mean        = means$delta_gamma_total,
    delta_gamma_total_sd          = sds$delta_gamma_total,
    # Extinction
    delta_epsilon_climate_mean      = means$delta_epsilon_climate,
    delta_epsilon_climate_sd        = sds$delta_epsilon_climate,
    delta_epsilon_ndvi_climate_mean = means$delta_epsilon_ndvi_climate,
    delta_epsilon_ndvi_climate_sd   = sds$delta_epsilon_ndvi_climate,
    delta_epsilon_ndvi_landuse_mean = means$delta_epsilon_ndvi_landuse,
    delta_epsilon_ndvi_landuse_sd   = sds$delta_epsilon_ndvi_landuse,
    delta_epsilon_landuse_mean      = means$delta_epsilon_landuse,
    delta_epsilon_landuse_sd        = sds$delta_epsilon_landuse,
    delta_epsilon_total_mean        = means$delta_epsilon_total,
    delta_epsilon_total_sd          = sds$delta_epsilon_total,
    stringsAsFactors = FALSE
  ))
}

rev_table_path <- here("results", "pub_table_attribution_revised.csv")
write.csv(rev_summary, track_file(rev_table_path), row.names = FALSE)
message("\n  Revised attribution table saved: ", rev_table_path)

#---------------------------------------------------------------------------
# Revised attribution maps per species
#---------------------------------------------------------------------------

for (sp in names(SPECIES)) {

  attr_df <- all_rev_attribution[[sp]]

  # Mean attribution per site (across years)
  site_means <- attr_df %>%
    group_by(site_id, lat, lon) %>%
    summarise(across(starts_with("delta_"), ~ mean(.x, na.rm = TRUE)),
              .groups = "drop")

  # 4-panel summary: climate, NDVI_climate, NDVI_landuse, land_use
  # for gamma (top) and epsilon (bottom) = 8 panels total, use 4x2
  plot_data <- site_means %>%
    pivot_longer(cols = starts_with("delta_"),
                 names_to = "variable", values_to = "value") %>%
    filter(variable %in% c("delta_gamma_climate", "delta_gamma_ndvi_climate",
                           "delta_gamma_ndvi_landuse", "delta_gamma_landuse",
                           "delta_epsilon_climate", "delta_epsilon_ndvi_climate",
                           "delta_epsilon_ndvi_landuse", "delta_epsilon_landuse"))

  var_labels <- c(
    delta_gamma_climate      = "Col: Climate (pr,tmmn,tmmx)",
    delta_gamma_ndvi_climate = "Col: NDVI-climate",
    delta_gamma_ndvi_landuse = "Col: NDVI-landuse",
    delta_gamma_landuse      = "Col: Land use (LC)",
    delta_epsilon_climate      = "Ext: Climate (pr,tmmn,tmmx)",
    delta_epsilon_ndvi_climate = "Ext: NDVI-climate",
    delta_epsilon_ndvi_landuse = "Ext: NDVI-landuse",
    delta_epsilon_landuse      = "Ext: Land use (LC)"
  )

  plot_data$panel <- factor(var_labels[plot_data$variable],
                            levels = var_labels)

  max_abs <- max(abs(plot_data$value), na.rm = TRUE)
  if (max_abs == 0 || !is.finite(max_abs)) max_abs <- 0.001

  p_rev <- ggplot() +
    geom_sf(data = spain_crop, fill = "grey95", colour = "grey50") +
    geom_point(data = plot_data,
               aes(x = lon, y = lat, colour = value),
               size = 0.3, alpha = 0.8) +
    scale_colour_gradient2(low = "red3", mid = "white", high = "blue3",
                           midpoint = 0, limits = c(-max_abs, max_abs),
                           name = "Effect") +
    facet_wrap(~ panel, ncol = 2) +
    coord_sf(xlim = c(-10, 5), ylim = c(35, 44)) +
    labs(title = paste0(SPECIES[[sp]], ": Revised 3-way Attribution"),
         x = "Longitude", y = "Latitude") +
    theme(strip.text = element_text(size = 8))

  fig_path <- here("figs",
                    paste0("fig_attribution_revised_summary_", sp, ".png"))
  ggsave(track_file(fig_path), p_rev, width = 12, height = 14, dpi = 300)
  message("  Revised attribution figure saved: ", basename(fig_path))
}

session_log[["task3"]] <- "COMPLETE"


###############################################################################
#              TASK 4: MODEL REVISION RECOMMENDATION
###############################################################################

message("\n", strrep("#", 70))
message("# TASK 4: MODEL REVISION RECOMMENDATIONS")
message(strrep("#", 70))

recommendations <- data.frame()
summary_text <- c()

for (sp in names(SPECIES)) {

  cfg <- get_model_config(sp)
  sp_vif   <- all_vif_results %>% filter(species == sp)
  sp_stab  <- stability_results %>% filter(species == sp)

  # R² summary for this species
  r2_sp <- all_r2_data %>% filter(species == sp)
  r2_clean <- r2_sp$r2[!is.na(r2_sp$r2)]
  mean_r2    <- mean(r2_clean)
  prop_r2_05 <- mean(r2_clean > 0.5)
  prop_r2_07 <- mean(r2_clean > 0.7)

  for (submodel in c("gamma", "epsilon")) {
    vars <- if (submodel == "gamma") cfg$gamma_vars else cfg$epsilon_vars

    if (!"NDVI" %in% vars) {
      recommendations <- rbind(recommendations, data.frame(
        species = sp, submodel = submodel,
        has_ndvi = FALSE, recommendation = "N/A",
        vif_ndvi = NA, max_r_ndvi = NA,
        mean_r2_climate = round(mean_r2, 3),
        prop_r2_gt05 = round(prop_r2_05, 3),
        prop_r2_gt07 = round(prop_r2_07, 3),
        sign_change = NA, delta_aic = NA,
        beta_shift_max = NA, se_change_max = NA,
        justification = "NDVI not in submodel",
        stringsAsFactors = FALSE
      ))
      next
    }

    # Get VIF for NDVI in this submodel
    vif_row <- sp_vif %>%
      filter(submodel == !!submodel, covariate == "NDVI")
    vif_ndvi <- if (nrow(vif_row) > 0) vif_row$VIF[1] else NA

    # Max |r| between NDVI and other covariates in this submodel
    sp_corr <- all_corr_flags %>% filter(species == sp)
    ndvi_cors <- sp_corr %>%
      filter((var1 == "NDVI" & var2 %in% vars) |
             (var2 == "NDVI" & var1 %in% vars))
    max_r <- if (nrow(ndvi_cors) > 0) max(abs(ndvi_cors$r)) else 0

    # Stability results
    stab_sub <- sp_stab %>% filter(submodel == !!submodel)
    any_sign_change <- any(stab_sub$sign_change == TRUE, na.rm = TRUE)
    daic <- if (nrow(stab_sub) > 0) stab_sub$delta_aic[1] else NA
    max_beta_shift <- if (nrow(stab_sub) > 0) max(stab_sub$beta_shift_se, na.rm = TRUE) else NA
    max_se_change  <- if (nrow(stab_sub) > 0) max(stab_sub$se_change_pct, na.rm = TRUE) else NA

    # Decision logic
    if (!is.na(vif_ndvi) && vif_ndvi > 5) {
      rec <- "(c) REMOVE"
      just <- paste0("VIF=", round(vif_ndvi, 1), " exceeds threshold of 5")
    } else if (any_sign_change) {
      rec <- "(c) REMOVE"
      just <- "Coefficient sign change when NDVI removed indicates instability"
    } else if (!is.na(max_se_change) && max_se_change > 50) {
      rec <- "(c) REMOVE"
      just <- paste0("SE increase of ", round(max_se_change, 0),
                      "% when collinear variable removed")
    } else if (prop_r2_07 > 0.5) {
      rec <- "(b) REPLACE with decomposed components"
      just <- paste0("R²(NDVI~climate)>0.7 at ", round(prop_r2_07 * 100, 0),
                      "% of sites — NDVI largely redundant with climate")
    } else if (prop_r2_05 > 0.5 && prop_r2_07 <= 0.5) {
      rec <- "(d) RETAIN but reclassify in attribution"
      just <- paste0("Intermediate R² (mean=", round(mean_r2, 2),
                      "); VIF<5; stable coefficients; use decomposed attribution")
    } else {
      rec <- "(a) RETAIN as-is"
      just <- paste0("VIF<5, |r|<", round(max_r, 2),
                      ", R²(NDVI~climate)<0.5 at ",
                      round((1 - prop_r2_05) * 100, 0),
                      "% of sites — substantial non-climate NDVI signal")
    }

    recommendations <- rbind(recommendations, data.frame(
      species = sp, submodel = submodel,
      has_ndvi = TRUE, recommendation = rec,
      vif_ndvi = round(vif_ndvi, 2), max_r_ndvi = round(max_r, 3),
      mean_r2_climate = round(mean_r2, 3),
      prop_r2_gt05 = round(prop_r2_05, 3),
      prop_r2_gt07 = round(prop_r2_07, 3),
      sign_change = any_sign_change,
      delta_aic = round(daic, 2),
      beta_shift_max = round(max_beta_shift, 3),
      se_change_max = round(max_se_change, 1),
      justification = just,
      stringsAsFactors = FALSE
    ))
  }

  # Narrative summary per species
  sp_recs <- recommendations %>% filter(species == sp, has_ndvi == TRUE)
  if (nrow(sp_recs) > 0) {
    txt <- paste0("\n", SPECIES[[sp]], " (", sp, "):\n")
    txt <- paste0(txt, "NDVI~climate R² across sites: mean=",
                  round(mean_r2, 3), ", ", round(prop_r2_05 * 100, 0),
                  "% of sites > 0.5, ", round(prop_r2_07 * 100, 0),
                  "% > 0.7. ")
    for (i in seq_len(nrow(sp_recs))) {
      row <- sp_recs[i, ]
      txt <- paste0(txt, "For the ", row$submodel, " submodel, NDVI ",
                    "co-occurs with ", paste(setdiff(
                      if (row$submodel == "gamma") cfg$gamma_vars else cfg$epsilon_vars,
                      "NDVI"), collapse = ", "), ". ")
      if (!is.na(row$vif_ndvi)) {
        txt <- paste0(txt, "VIF(NDVI)=", row$vif_ndvi, ". ")
      }
      if (!is.na(row$delta_aic)) {
        txt <- paste0(txt, "Removing NDVI changes AIC by ", row$delta_aic,
                      ifelse(abs(row$delta_aic) > 2, " (>2, meaningful)", " (<2, negligible)"),
                      ". ")
      }
      if (!is.na(row$sign_change) && row$sign_change) {
        txt <- paste0(txt, "WARNING: coefficient sign change detected. ")
      }
      txt <- paste0(txt, "Recommendation: ", row$recommendation,
                    " — ", row$justification, ". ")
    }
    summary_text <- c(summary_text, txt)
  } else {
    summary_text <- c(summary_text,
      paste0("\n", SPECIES[[sp]], " (", sp, "):\n",
             "NDVI does not appear in either colonisation or extinction submodel. ",
             "No revision needed. The attribution analysis correctly classifies ",
             "all covariates as either pure climate (pr, tmmn, tmmx) or pure ",
             "land use (LC classes). NDVI~climate R² is reported for reference: ",
             "mean=", round(mean_r2, 3), ".\n"))
  }
}

# Save recommendations
rec_path <- here("results", "diagnostics", "ndvi_model_revision_recommendations.csv")
write.csv(recommendations, track_file(rec_path), row.names = FALSE)

# Save narrative
sum_path <- here("results", "diagnostics", "ndvi_revision_summary.txt")
writeLines(c("NDVI MODEL REVISION RECOMMENDATIONS",
             paste0("Generated: ", Sys.time()),
             strrep("=", 70),
             summary_text),
           track_file(sum_path))
message("\n  Recommendations saved: ", rec_path)
message("  Narrative saved: ", sum_path)

session_log[["task4"]] <- "COMPLETE"


###############################################################################
#     TASK 5: PUBLICATION-READY COLLINEARITY DIAGNOSTIC FIGURE
###############################################################################

message("\n", strrep("#", 70))
message("# TASK 5: PUBLICATION-READY COLLINEARITY DIAGNOSTIC FIGURE")
message(strrep("#", 70))

# Layout: 4 rows (species) x 3 columns (A: corr heatmap, B: VIF bars, C: R² map)

plot_list <- list()

for (idx in seq_along(names(SPECIES))) {
  sp <- names(SPECIES)[idx]
  sp_name <- SPECIES[[sp]]
  cfg <- get_model_config(sp)

  #--- Panel A: Correlation heatmap ---
  corr_mat <- all_corr_matrices[[sp]]
  short_names <- gsub("Land_Cover_Type_1_Percent_", "LC", colnames(corr_mat))
  colnames(corr_mat) <- short_names
  rownames(corr_mat) <- short_names

  # Convert to long format for ggplot
  corr_long <- expand.grid(Var1 = rownames(corr_mat), Var2 = colnames(corr_mat),
                           stringsAsFactors = FALSE)
  corr_long$value <- as.vector(corr_mat)
  # Lower triangle only
  corr_long <- corr_long[match(corr_long$Var1, rownames(corr_mat)) >=
                          match(corr_long$Var2, colnames(corr_mat)), ]

  pA <- ggplot(corr_long, aes(x = Var2, y = Var1, fill = value)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = round(value, 2)), size = 2.2) +
    scale_fill_gradient2(low = "blue3", mid = "white", high = "red3",
                         midpoint = 0, limits = c(-1, 1), name = "r") +
    labs(title = paste0(sp_name), x = NULL, y = NULL) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
          axis.text.y = element_text(size = 6),
          title = element_text(size = 8, face = "italic"),
          legend.key.width = unit(0.3, "cm"),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 7))

  #--- Panel B: VIF bar chart ---
  sp_vif <- all_vif_results %>%
    filter(species == sp, !is.na(VIF)) %>%
    mutate(cov_short = gsub("Land_Cover_Type_1_Percent_", "LC", covariate))

  if (nrow(sp_vif) > 0) {
    pB <- ggplot(sp_vif, aes(x = cov_short, y = VIF, fill = submodel)) +
      geom_col(position = "dodge", width = 0.7) +
      geom_hline(yintercept = 5, linetype = "dashed", colour = "red", linewidth = 0.7) +
      geom_hline(yintercept = 10, linetype = "dotted", colour = "red2", linewidth = 0.5) +
      scale_fill_manual(values = c(gamma = "#2196F3", epsilon = "#FF9800"),
                        name = "Submodel") +
      labs(x = NULL, y = "VIF") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
            axis.text.y = element_text(size = 6),
            legend.text = element_text(size = 6),
            legend.title = element_text(size = 7),
            legend.key.size = unit(0.3, "cm"))
  } else {
    pB <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No VIF data") +
      theme_void()
  }

  #--- Panel C: R² map ---
  r2_sp <- all_r2_data %>% filter(species == sp, !is.na(r2))

  pC <- ggplot() +
    geom_sf(data = spain_crop, fill = "grey95", colour = "grey50") +
    geom_point(data = r2_sp, aes(x = lon, y = lat, colour = r2),
               size = 0.3, alpha = 0.8) +
    scale_colour_viridis_c(name = expression(R^2), option = "inferno",
                            limits = c(0, 1)) +
    coord_sf(xlim = c(-10, 5), ylim = c(35, 44)) +
    labs(x = NULL, y = NULL) +
    theme(axis.text = element_text(size = 5),
          legend.key.width = unit(0.3, "cm"),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 7))

  # Combine row
  plot_list[[idx]] <- (pA | pB | pC) +
    plot_layout(widths = c(1, 0.8, 1))
}

# Stack all 4 species
pub_fig <- plot_list[[1]] / plot_list[[2]] / plot_list[[3]] / plot_list[[4]] +
  plot_annotation(
    title = "Dynamic Covariate Collinearity Diagnostics",
    subtitle = "A: Correlation matrix | B: Variance Inflation Factors | C: NDVI~Climate R²",
    theme = theme(plot.title = element_text(size = 12, face = "bold"),
                  plot.subtitle = element_text(size = 9))
  )

pub_path <- here("figs", "pub_fig_collinearity_diagnostics.png")
ggsave(track_file(pub_path), pub_fig, width = 14, height = 18, dpi = 300)
message("  Publication figure saved: ", pub_path)

session_log[["task5"]] <- "COMPLETE"


###############################################################################
# SESSION SUMMARY
###############################################################################

message("\n", strrep("=", 70))
message("  SESSION SUMMARY")
message(strrep("=", 70))

# Task completion
for (task in names(session_log)) {
  message("  ", task, ": ", session_log[[task]])
}

# Collinearity flags
flagged_spp <- all_corr_flags %>%
  filter(flag %in% c("MODERATE", "SEVERE")) %>%
  distinct(species)
message("\n  Species with collinearity flags: ",
        paste(flagged_spp$species, collapse = ", "))

vif_flagged <- all_vif_results %>%
  filter(flag %in% c("PROBLEMATIC", "SEVERE"))
if (nrow(vif_flagged) > 0) {
  message("  VIF flags:")
  for (i in seq_len(nrow(vif_flagged))) {
    message("    ", vif_flagged$species[i], "/", vif_flagged$submodel[i],
            "/", vif_flagged$covariate[i], ": VIF=", vif_flagged$VIF[i])
  }
} else {
  message("  No VIF flags (all < 5).")
}

# Output files
message("\n  Output files generated (", length(output_files), "):")
for (f in output_files) {
  message("    ", f)
}

message("\n  Script completed: ", Sys.time())
