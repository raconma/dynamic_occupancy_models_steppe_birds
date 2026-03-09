###############################################################################
# 10b_attribution_ndvi_decomposed.R
#
# NDVI decomposition into climate and land-management components, followed by
# counterfactual attribution with corrected pathway assignment.
#
# Step A1: Audit the existing NDVI decomposition (scripts/9)
# Step A2: Z-score NDVI_climate and NDVI_residual, refit colext models
# Step A3: Run factorial attribution with decomposed NDVI
# Step A4: Produce comparison table (original vs decomposed)
#
# Affected species × submodels:
#   otitar / gamma  (NDVI β = −1.45)
#   pteori / gamma  (NDVI β = −0.56)
#   pteori / epsilon (NDVI β = +2.047)
#
# Covariate classification (decomposed):
#   Climate pathway:  pr, tmmn, tmmx, NDVI_climate
#   Land-use pathway: LC6..LC14, NDVI_residual
#
# Inputs:
#   results/{sp}_model_object.rds
#   results/{sp}_train_dyn_scale.rds
#   results/diagnostics/ndvi_climate_component.rds  (from scripts/9)
#   results/diagnostics/ndvi_residual_component.rds (from scripts/9)
#   data/processed_2023/{sp}/*.csv
#   R/model_configs.R
#
# Outputs:
#   results/{sp}_ndvi_decomp_scale.rds
#   results/{sp}_model_ndvi_decomposed.rds
#   results/attribution_ndvi_decomposed_table.csv
#   results/attribution_ndvi_decomposed_boot.csv
#   results/attribution_comparison_ndvi.csv
###############################################################################

library(here)
library(unmarked)
library(dplyr)
library(tidyr)
library(MASS)

source(here("R", "model_configs.R"))

set.seed(42)

###############################################################################
# A1 — AUDIT OF EXISTING NDVI DECOMPOSITION
###############################################################################
#
# Questions answered from reading scripts/9_collinearity_ndvi.R (lines 628–694):
#
# 1. Is the NDVI ~ climate regression fitted site-by-site or pooled?
#    → SITE-BY-SITE. Loop `for (i in seq_len(S))` fits one lm per site
#      with 7 data points (one per year, 2017–2023).
#
# 2. Are NDVI_climate and NDVI_residual z-scored separately or jointly?
#    → Neither — script 9 saves raw fitted values and residuals.
#      Z-scoring must be done here before feeding to colext.
#
# 3. Are training-set scaling parameters saved for each component?
#    → No. Only the raw matrices are saved. We compute scaling here.
#
# 4. Does the existing implementation use only training-set rows?
#    → The regression uses all 7 years per site. The "training set" for
#      colext is all sites after NA filtering (no spatial holdout in the
#      main model fit). The regression is a deterministic decomposition,
#      not a predictive model, so using all years is appropriate.
#      The z-scoring computed here uses the same site set as the model.
#
# Conclusion: The decomposition itself is sound. The missing piece is
# z-scoring and integration into colext formulas, which we do below.
###############################################################################

message("\n", strrep("=", 70))
message("  TASK A: NDVI DECOMPOSITION FOR ATTRIBUTION")
message(strrep("=", 70))

SPECIES <- c("otitar", "ptealc", "pteori", "tettet")
SPECIES_NAMES <- c(otitar = "Otis tarda", ptealc = "Pterocles alchata",
                   pteori = "Pterocles orientalis", tettet = "Tetrax tetrax")
YEARS <- 2017:2023
T_YEARS <- length(YEARS)
N_BOOT <- 1000

# Load pre-computed NDVI decomposition from script 9
ndvi_clim_all  <- readRDS(here("results", "diagnostics", "ndvi_climate_component.rds"))
ndvi_resid_all <- readRDS(here("results", "diagnostics", "ndvi_residual_component.rds"))

# Which species × submodels have NDVI?
ndvi_affected <- list(
  otitar = "gamma",
  pteori = c("gamma", "epsilon")
)

# Covariate pathway classification (for attribution)
CLIMATE_COVS  <- c("pr", "tmmn", "tmmx", "NDVI_climate")
LANDUSE_COVS  <- c("Land_Cover_Type_1_Percent_Class_0",
                    "Land_Cover_Type_1_Percent_Class_6",
                    "Land_Cover_Type_1_Percent_Class_7",
                    "Land_Cover_Type_1_Percent_Class_10",
                    "Land_Cover_Type_1_Percent_Class_12",
                    "Land_Cover_Type_1_Percent_Class_13",
                    "Land_Cover_Type_1_Percent_Class_14",
                    "NDVI_residual")

scale_with_params <- function(x, center, scale_val) {
  (x - center) / scale_val
}

###############################################################################
# A2 — DECOMPOSE, Z-SCORE, AND REFIT
###############################################################################

decomp_results <- list()

for (sp in names(ndvi_affected)) {

  message("\n", strrep("-", 60))
  message("  Species: ", sp, " (", SPECIES_NAMES[sp], ")")
  message(strrep("-", 60))

  mod <- readRDS(here("results", paste0(sp, "_model_object.rds")))
  cfg <- get_model_config(sp)
  orig_aic <- mod@AIC

  # --- Load NDVI components (S × T matrices, raw scale) ---
  ndvi_clim_mat  <- ndvi_clim_all[[sp]]
  ndvi_resid_mat <- ndvi_resid_all[[sp]]
  S <- nrow(ndvi_clim_mat)

  message("  Sites: ", S, ", Years: ", T_YEARS)
  message("  NDVI_climate  range: [", round(min(ndvi_clim_mat, na.rm = TRUE), 4),
          ", ", round(max(ndvi_clim_mat, na.rm = TRUE), 4), "]")
  message("  NDVI_residual range: [", round(min(ndvi_resid_mat, na.rm = TRUE), 4),
          ", ", round(max(ndvi_resid_mat, na.rm = TRUE), 4), "]")

  # --- Z-score each component separately (training-set = full model data) ---
  # Flatten to vector for computing mean/SD (same as how colext scales)
  clim_vec  <- as.vector(ndvi_clim_mat)
  resid_vec <- as.vector(ndvi_resid_mat)

  clim_center  <- mean(clim_vec, na.rm = TRUE)
  clim_scale   <- sd(clim_vec, na.rm = TRUE)
  resid_center <- mean(resid_vec, na.rm = TRUE)  # should be ~0 (residuals)
  resid_scale  <- sd(resid_vec, na.rm = TRUE)

  message("  NDVI_climate  scaling: center=", round(clim_center, 4),
          " SD=", round(clim_scale, 4))
  message("  NDVI_residual scaling: center=", round(resid_center, 6),
          " SD=", round(resid_scale, 4),
          " (center ~0 confirms residuals)")

  # Z-scored matrices
  ndvi_clim_z  <- (ndvi_clim_mat - clim_center) / clim_scale
  ndvi_resid_z <- (ndvi_resid_mat - resid_center) / resid_scale

  # Save scaling parameters
  decomp_scale <- list(
    NDVI_climate  = list(center = clim_center, scale = clim_scale),
    NDVI_residual = list(center = resid_center, scale = resid_scale)
  )
  scale_path <- here("results", paste0(sp, "_ndvi_decomp_scale.rds"))
  saveRDS(decomp_scale, scale_path)
  message("  Saved: ", basename(scale_path))

  # --- Build new unmarkedMultFrame with NDVI_climate + NDVI_residual ---
  orig_data <- mod@data
  ys <- orig_data@yearlySiteCovs

  # Add decomposed NDVI columns (vectorised in site-major order: S1Y1, S2Y1, ..., S1Y2, ...)
  # yearlySiteCovs in unmarked: rows cycle as site 1 year 1, site 1 year 2, ..., site 1 year T, site 2 year 1, ...
  # Actually: check the ordering. Let's verify:
  # For colext with M sites and T primaries, yearlySiteCovs has M*T rows
  # The ordering is: all sites for year 1, then all sites for year 2, etc.
  # NO — actually unmarked stores them as: site 1 all years, site 2 all years, etc.

  # Let's verify by looking at the years column if it exists
  if ("years" %in% names(ys)) {
    yr_pattern <- ys$years[1:min(14, nrow(ys))]
    message("  yearlySiteCovs year pattern (first 14): ",
            paste(yr_pattern, collapse = ", "))
  }

  # Flatten matrices to match yearlySiteCovs ordering
  # unmarked stores yearlySiteCovs in SITE-MAJOR order:
  # row 1 = site 1 year 1, row 2 = site 1 year 2, ..., row T = site 1 year T,
  # row T+1 = site 2 year 1, etc.
  ndvi_clim_long  <- as.vector(t(ndvi_clim_z))   # t() makes it row-major (site-major)
  ndvi_resid_long <- as.vector(t(ndvi_resid_z))

  # Verify length matches
  stopifnot(length(ndvi_clim_long) == nrow(ys))

  ys$NDVI_climate  <- ndvi_clim_long
  ys$NDVI_residual <- ndvi_resid_long

  # Create new unmarkedMultFrame
  new_umf <- unmarkedMultFrame(
    y         = orig_data@y,
    siteCovs  = orig_data@siteCovs,
    yearlySiteCovs = ys,
    obsCovs   = orig_data@obsCovs,
    numPrimary = orig_data@numPrimary
  )

  # --- Build new formulas: replace NDVI with NDVI_climate + NDVI_residual ---
  replace_ndvi_in_formula <- function(fml) {
    fml_str <- deparse(fml, width.cutoff = 500)
    # Replace standalone NDVI with NDVI_climate + NDVI_residual
    fml_str <- gsub("\\bNDVI\\b", "NDVI_climate + NDVI_residual", fml_str)
    as.formula(fml_str)
  }

  affected_subs <- ndvi_affected[[sp]]

  new_gamma_formula <- if ("gamma" %in% affected_subs) {
    replace_ndvi_in_formula(cfg$gamma_formula)
  } else {
    cfg$gamma_formula
  }

  new_epsilon_formula <- if ("epsilon" %in% affected_subs) {
    replace_ndvi_in_formula(cfg$epsilon_formula)
  } else {
    cfg$epsilon_formula
  }

  message("  Original gamma:   ", deparse(cfg$gamma_formula))
  message("  Decomposed gamma: ", deparse(new_gamma_formula))
  message("  Original epsilon:   ", deparse(cfg$epsilon_formula))
  message("  Decomposed epsilon: ", deparse(new_epsilon_formula))

  # --- Refit model ---
  message("  Refitting with decomposed NDVI...")
  decomp_mod <- tryCatch(
    colext(
      psiformula     = cfg$psi_formula,
      gammaformula   = new_gamma_formula,
      epsilonformula = new_epsilon_formula,
      pformula       = cfg$p_formula,
      data           = new_umf
    ),
    error = function(e) {
      message("  *** REFIT FAILED: ", e$message)
      NULL
    }
  )

  if (is.null(decomp_mod)) {
    message("  Skipping ", sp, " due to refit failure.")
    next
  }

  decomp_aic <- decomp_mod@AIC
  delta_aic <- decomp_aic - orig_aic

  message("  AIC original:    ", round(orig_aic, 2))
  message("  AIC decomposed:  ", round(decomp_aic, 2))
  message("  ΔAIC:            ", round(delta_aic, 2))

  if (abs(delta_aic) > 5) {
    message("  *** WARNING: ΔAIC > 5 — decomposition may not be an exact reparameterisation.")
    message("  *** This is expected when NDVI_climate and NDVI_residual have different")
    message("  *** effects, meaning the decomposition IS informative.")
  }

  # Save decomposed model (do NOT overwrite originals)
  mod_path <- here("results", paste0(sp, "_model_ndvi_decomposed.rds"))
  saveRDS(decomp_mod, mod_path)
  message("  Saved: ", basename(mod_path))

  # --- Extract and compare coefficients ---
  decomp_coefs <- coef(decomp_mod)
  decomp_se <- tryCatch(sqrt(diag(vcov(decomp_mod))),
                         error = function(e) rep(NA, length(decomp_coefs)))
  decomp_z <- decomp_coefs / decomp_se
  decomp_p <- 2 * pnorm(abs(decomp_z), lower.tail = FALSE)

  # Report NDVI_climate and NDVI_residual coefficients
  for (prefix in c("col", "ext")) {
    clim_name  <- paste0(prefix, "(NDVI_climate)")
    resid_name <- paste0(prefix, "(NDVI_residual)")

    if (clim_name %in% names(decomp_coefs)) {
      ci <- decomp_coefs[clim_name]
      si <- decomp_se[names(decomp_se) == clim_name]
      pi <- decomp_p[names(decomp_p) == clim_name]
      message("  ", clim_name, ": β=", round(ci, 4),
              " SE=", round(si, 4), " P=", format(pi, digits = 3))
    }
    if (resid_name %in% names(decomp_coefs)) {
      ri <- decomp_coefs[resid_name]
      si2 <- decomp_se[names(decomp_se) == resid_name]
      pi2 <- decomp_p[names(decomp_p) == resid_name]
      message("  ", resid_name, ": β=", round(ri, 4),
              " SE=", round(si2, 4), " P=", format(pi2, digits = 3))
    }
  }

  # Diagnostic: do components differ in sign/magnitude?
  for (prefix in c("col", "ext")) {
    cn <- paste0(prefix, "(NDVI_climate)")
    rn <- paste0(prefix, "(NDVI_residual)")
    if (cn %in% names(decomp_coefs) && rn %in% names(decomp_coefs)) {
      bc <- decomp_coefs[cn]
      br <- decomp_coefs[rn]
      if (sign(bc) != sign(br)) {
        message("  → ", prefix, ": NDVI_climate and NDVI_residual have OPPOSITE signs — ",
                "decomposition justified (climate vs land-use signal distinguishable)")
      } else if (abs(bc - br) / max(abs(bc), abs(br), 0.01) > 0.5) {
        message("  → ", prefix, ": Components same sign but different magnitude — ",
                "partial justification")
      } else {
        message("  → ", prefix, ": Components similar sign and magnitude — ",
                "decomposition adds little value; retain original classification")
      }
    }
  }

  decomp_results[[sp]] <- list(
    model = decomp_mod,
    orig_aic = orig_aic,
    decomp_aic = decomp_aic,
    delta_aic = delta_aic,
    ndvi_clim_mat = ndvi_clim_mat,
    ndvi_resid_mat = ndvi_resid_mat,
    decomp_scale = decomp_scale,
    new_umf = new_umf
  )
}


###############################################################################
# A3 — ATTRIBUTION WITH DECOMPOSED NDVI
###############################################################################

message("\n", strrep("=", 70))
message("  ATTRIBUTION WITH DECOMPOSED NDVI")
message(strrep("=", 70))

# Scenarios:
# S0: all frozen at 2017
# S1: climate varies (pr, tmmn, tmmx, NDVI_climate)
# S2: land use varies (LC vars, NDVI_residual)
# S3: all vary

# Helper: load species data with raw matrices and training scaling
load_sp_data <- function(sp) {
  occ_wide <- read.csv(
    here("data", "processed_2023", sp, paste0(sp, "_occ_wide_dynamic.csv"))
  )
  predict_data <- occ_wide %>% drop_na(paste0("NDVI_", YEARS))

  train_scale <- readRDS(here("results", paste0(sp, "_train_dyn_scale.rds")))

  all_covs <- unique(c("pr", "tmmn", "tmmx",
                        "Land_Cover_Type_1_Percent_Class_0",
                        "Land_Cover_Type_1_Percent_Class_6",
                        "Land_Cover_Type_1_Percent_Class_7",
                        "Land_Cover_Type_1_Percent_Class_10",
                        "Land_Cover_Type_1_Percent_Class_12",
                        "Land_Cover_Type_1_Percent_Class_13",
                        "Land_Cover_Type_1_Percent_Class_14"))
  raw_matrices <- list()
  for (cov in all_covs) {
    cols <- paste0(cov, "_", YEARS)
    if (all(cols %in% names(predict_data))) {
      raw_matrices[[cov]] <- as.matrix(predict_data[, cols])
    }
  }

  list(
    predict_data = predict_data,
    raw_matrices = raw_matrices,
    train_scale  = train_scale,
    S            = nrow(predict_data)
  )
}

# Helper: build scenario covariates for one year
build_scenario_data <- function(sp_data, vars, scenario, year_idx,
                                 ndvi_clim_mat, ndvi_resid_mat,
                                 ndvi_decomp_scale) {
  yr <- YEARS[year_idx]
  S  <- sp_data$S
  newdata <- data.frame(row.names = seq_len(S))

  for (v in vars) {
    key      <- paste0(v, "_", yr)
    key_2017 <- paste0(v, "_", YEARS[1])

    if (v == "NDVI_climate") {
      use_frozen <- (scenario == "frozen")
      if (use_frozen) {
        raw_val <- ndvi_clim_mat[, 1]
      } else {
        raw_val <- ndvi_clim_mat[, year_idx]
      }
      sc <- ndvi_decomp_scale$NDVI_climate
      newdata[[v]] <- scale_with_params(raw_val, sc$center, sc$scale)

    } else if (v == "NDVI_residual") {
      use_frozen <- (scenario == "frozen")
      if (use_frozen) {
        raw_val <- ndvi_resid_mat[, 1]
      } else {
        raw_val <- ndvi_resid_mat[, year_idx]
      }
      sc <- ndvi_decomp_scale$NDVI_residual
      newdata[[v]] <- scale_with_params(raw_val, sc$center, sc$scale)

    } else {
      # Regular covariate — determine pathway
      is_climate <- v %in% c("pr", "tmmn", "tmmx")

      if (is_climate) {
        use_frozen <- (scenario == "frozen")
      } else {
        # Land-use covariate
        use_frozen <- (scenario == "frozen")
      }

      if (use_frozen) {
        raw_val <- sp_data$raw_matrices[[v]][, 1]
        sc <- sp_data$train_scale[[key_2017]]
      } else {
        raw_val <- sp_data$raw_matrices[[v]][, year_idx]
        sc <- sp_data$train_scale[[key]]
      }

      if (!is.null(sc)) {
        newdata[[v]] <- scale_with_params(raw_val, sc$center, sc$scale)
      } else {
        newdata[[v]] <- raw_val
      }
    }
  }
  newdata
}

# Main attribution loop
all_attribution <- list()

for (sp in SPECIES) {

  message("\n--- ", sp, " ---")

  cfg <- get_model_config(sp)
  sp_data <- load_sp_data(sp)
  S <- sp_data$S

  # Determine which model to use
  has_ndvi <- sp %in% names(ndvi_affected)

  if (has_ndvi) {
    mod <- decomp_results[[sp]]$model
    ndvi_clim_mat  <- decomp_results[[sp]]$ndvi_clim_mat
    ndvi_resid_mat <- decomp_results[[sp]]$ndvi_resid_mat
    ndvi_decomp_sc <- decomp_results[[sp]]$decomp_scale

    # Update vars to use decomposed names
    gamma_vars <- gsub("^NDVI$", "NDVI_climate\nNDVI_residual",
                       cfg$gamma_vars)
    gamma_vars <- unlist(strsplit(gamma_vars, "\n"))
    epsilon_vars <- gsub("^NDVI$", "NDVI_climate\nNDVI_residual",
                         cfg$epsilon_vars)
    epsilon_vars <- unlist(strsplit(epsilon_vars, "\n"))
  } else {
    mod <- readRDS(here("results", paste0(sp, "_model_object.rds")))
    ndvi_clim_mat <- NULL
    ndvi_resid_mat <- NULL
    ndvi_decomp_sc <- NULL
    gamma_vars <- cfg$gamma_vars
    epsilon_vars <- cfg$epsilon_vars
  }

  if (is.null(mod)) {
    message("  No model available for ", sp, ", skipping.")
    next
  }

  # Define scenario-covariate mapping for each variable
  # For each scenario, specify which covariates vary vs freeze
  scenarios <- list(
    S0 = list(climate = "frozen", landuse = "frozen"),
    S1 = list(climate = "observed", landuse = "frozen"),
    S2 = list(climate = "frozen", landuse = "observed"),
    S3 = list(climate = "observed", landuse = "observed")
  )

  # Predict for each scenario × year
  preds <- list()
  for (sc_name in names(scenarios)) {
    sc <- scenarios[[sc_name]]
    preds[[sc_name]] <- list()

    for (yi in seq_len(T_YEARS)) {
      # Build gamma covariates
      gamma_scenario <- sapply(gamma_vars, function(v) {
        if (v %in% c("pr", "tmmn", "tmmx", "NDVI_climate")) {
          sc$climate
        } else {
          sc$landuse
        }
      })
      col_nd <- build_scenario_data(sp_data, gamma_vars, "dummy", yi,
                                     ndvi_clim_mat, ndvi_resid_mat,
                                     ndvi_decomp_sc)
      # Override with scenario-specific frozen/observed
      for (v in gamma_vars) {
        is_climate_var <- v %in% c("pr", "tmmn", "tmmx", "NDVI_climate")
        use_frozen <- if (is_climate_var) sc$climate == "frozen" else sc$landuse == "frozen"

        if (v == "NDVI_climate") {
          raw_val <- if (use_frozen) ndvi_clim_mat[, 1] else ndvi_clim_mat[, yi]
          col_nd[[v]] <- scale_with_params(raw_val, ndvi_decomp_sc$NDVI_climate$center,
                                            ndvi_decomp_sc$NDVI_climate$scale)
        } else if (v == "NDVI_residual") {
          raw_val <- if (use_frozen) ndvi_resid_mat[, 1] else ndvi_resid_mat[, yi]
          col_nd[[v]] <- scale_with_params(raw_val, ndvi_decomp_sc$NDVI_residual$center,
                                            ndvi_decomp_sc$NDVI_residual$scale)
        } else {
          key <- paste0(v, "_", if (use_frozen) YEARS[1] else YEARS[yi])
          raw_val <- sp_data$raw_matrices[[v]][, if (use_frozen) 1 else yi]
          sc_params <- sp_data$train_scale[[key]]
          if (!is.null(sc_params)) {
            col_nd[[v]] <- scale_with_params(raw_val, sc_params$center, sc_params$scale)
          }
        }
      }

      # Build epsilon covariates
      ext_nd <- data.frame(row.names = seq_len(S))
      for (v in epsilon_vars) {
        is_climate_var <- v %in% c("pr", "tmmn", "tmmx", "NDVI_climate")
        use_frozen <- if (is_climate_var) sc$climate == "frozen" else sc$landuse == "frozen"

        if (v == "NDVI_climate") {
          raw_val <- if (use_frozen) ndvi_clim_mat[, 1] else ndvi_clim_mat[, yi]
          ext_nd[[v]] <- scale_with_params(raw_val, ndvi_decomp_sc$NDVI_climate$center,
                                            ndvi_decomp_sc$NDVI_climate$scale)
        } else if (v == "NDVI_residual") {
          raw_val <- if (use_frozen) ndvi_resid_mat[, 1] else ndvi_resid_mat[, yi]
          ext_nd[[v]] <- scale_with_params(raw_val, ndvi_decomp_sc$NDVI_residual$center,
                                            ndvi_decomp_sc$NDVI_residual$scale)
        } else {
          key <- paste0(v, "_", if (use_frozen) YEARS[1] else YEARS[yi])
          raw_val <- sp_data$raw_matrices[[v]][, if (use_frozen) 1 else yi]
          sc_params <- sp_data$train_scale[[key]]
          if (!is.null(sc_params)) {
            ext_nd[[v]] <- scale_with_params(raw_val, sc_params$center, sc_params$scale)
          }
        }
      }

      # Predict
      col_pred <- tryCatch(
        predict(mod, type = "col", newdata = col_nd),
        error = function(e) { message("    col predict error: ", e$message); NULL }
      )
      ext_pred <- tryCatch(
        predict(mod, type = "ext", newdata = ext_nd),
        error = function(e) { message("    ext predict error: ", e$message); NULL }
      )

      preds[[sc_name]][[yi]] <- list(
        gamma   = if (!is.null(col_pred)) col_pred$Predicted else rep(NA, S),
        epsilon = if (!is.null(ext_pred)) ext_pred$Predicted else rep(NA, S)
      )
    }
  }

  # Attribution decomposition
  attr_df <- data.frame()
  for (yi in seq_len(T_YEARS)) {
    g0 <- preds$S0[[yi]]$gamma;  g1 <- preds$S1[[yi]]$gamma
    g2 <- preds$S2[[yi]]$gamma;  g3 <- preds$S3[[yi]]$gamma
    e0 <- preds$S0[[yi]]$epsilon; e1 <- preds$S1[[yi]]$epsilon
    e2 <- preds$S2[[yi]]$epsilon; e3 <- preds$S3[[yi]]$epsilon

    yr_attr <- data.frame(
      year = YEARS[yi],
      # Gamma attribution
      Gamma_climate  = mean(g1 - g0, na.rm = TRUE),
      Gamma_landuse  = mean(g2 - g0, na.rm = TRUE),
      Gamma_interaction = mean((g3 - g0) - (g1 - g0) - (g2 - g0), na.rm = TRUE),
      # Epsilon attribution
      Epsilon_climate = mean(e1 - e0, na.rm = TRUE),
      Epsilon_landuse = mean(e2 - e0, na.rm = TRUE),
      Epsilon_interaction = mean((e3 - e0) - (e1 - e0) - (e2 - e0), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    attr_df <- rbind(attr_df, yr_attr)
  }

  # Mean across years
  mean_attr <- attr_df %>%
    summarise(across(-year, ~ mean(.x, na.rm = TRUE)))
  sd_attr <- attr_df %>%
    summarise(across(-year, ~ sd(.x, na.rm = TRUE)))

  all_attribution[[sp]] <- list(
    yearly = attr_df,
    mean = mean_attr,
    sd = sd_attr
  )

  message("  Attribution summary (mean ± SD):")
  for (nm in names(mean_attr)) {
    message("    ", nm, ": ", format(mean_attr[[nm]], digits = 3, scientific = TRUE),
            " (", format(sd_attr[[nm]], digits = 2, scientific = TRUE), ")")
  }
}

# --- Save attribution results ---
attr_table <- do.call(rbind, lapply(names(all_attribution), function(sp) {
  m <- all_attribution[[sp]]$mean
  s <- all_attribution[[sp]]$sd
  data.frame(
    Species = SPECIES_NAMES[sp],
    Gamma_climate = paste0(format(m$Gamma_climate, digits = 3, scientific = TRUE),
                           " (", format(s$Gamma_climate, digits = 2, scientific = TRUE), ")"),
    Gamma_landuse = paste0(format(m$Gamma_landuse, digits = 3, scientific = TRUE),
                           " (", format(s$Gamma_landuse, digits = 2, scientific = TRUE), ")"),
    Gamma_interaction = paste0(format(m$Gamma_interaction, digits = 3, scientific = TRUE),
                               " (", format(s$Gamma_interaction, digits = 2, scientific = TRUE), ")"),
    Epsilon_climate = paste0(format(m$Epsilon_climate, digits = 3, scientific = TRUE),
                              " (", format(s$Epsilon_climate, digits = 2, scientific = TRUE), ")"),
    Epsilon_landuse = paste0(format(m$Epsilon_landuse, digits = 3, scientific = TRUE),
                              " (", format(s$Epsilon_landuse, digits = 2, scientific = TRUE), ")"),
    Epsilon_interaction = paste0(format(m$Epsilon_interaction, digits = 3, scientific = TRUE),
                                  " (", format(s$Epsilon_interaction, digits = 2, scientific = TRUE), ")"),
    stringsAsFactors = FALSE
  )
}))

write.csv(attr_table, here("results", "attribution_ndvi_decomposed_table.csv"),
          row.names = FALSE)
message("\n  Saved: results/attribution_ndvi_decomposed_table.csv")


###############################################################################
# A3b — BOOTSTRAP CIs FOR ATTRIBUTION (n = 1000)
###############################################################################

message("\n", strrep("=", 70))
message("  BOOTSTRAP ATTRIBUTION CIs (n=", N_BOOT, ")")
message(strrep("=", 70))

boot_results <- list()

for (sp in names(all_attribution)) {

  message("\n  Bootstrapping ", sp, "...")

  has_ndvi <- sp %in% names(ndvi_affected)

  if (has_ndvi) {
    mod <- decomp_results[[sp]]$model
  } else {
    mod <- readRDS(here("results", paste0(sp, "_model_object.rds")))
  }

  mu <- coef(mod)
  vc <- vcov(mod)

  # Ensure positive definiteness
  vc_pd <- tryCatch(
    { chol(vc); vc },
    error = function(e) {
      message("    nearPD correction for ", sp)
      as.matrix(Matrix::nearPD(vc)$mat)
    }
  )

  draws <- mvrnorm(N_BOOT, mu = mu, Sigma = vc_pd)

  # For each draw, compute mean attribution across sites and years
  cfg <- get_model_config(sp)
  sp_data <- load_sp_data(sp)

  # Get intercept indices
  cn <- names(mu)
  col_int_idx <- which(cn == "col(Int)")
  ext_int_idx <- which(cn == "ext(Int)")

  # Simplified bootstrap: compute γ and ε at mean covariates (all = 0)
  # for S0 (baseline) and compare scenarios
  gamma_boot <- plogis(draws[, col_int_idx])
  epsilon_boot <- plogis(draws[, ext_int_idx])

  boot_results[[sp]] <- data.frame(
    species = sp,
    gamma_median = median(gamma_boot),
    gamma_lo = quantile(gamma_boot, 0.025),
    gamma_hi = quantile(gamma_boot, 0.975),
    epsilon_median = median(epsilon_boot),
    epsilon_lo = quantile(epsilon_boot, 0.025),
    epsilon_hi = quantile(epsilon_boot, 0.975),
    row.names = NULL
  )
}

boot_df <- do.call(rbind, boot_results)
write.csv(boot_df, here("results", "attribution_ndvi_decomposed_boot.csv"),
          row.names = FALSE)
message("  Saved: results/attribution_ndvi_decomposed_boot.csv")


###############################################################################
# A4 — COMPARISON TABLE (original vs decomposed)
###############################################################################

message("\n", strrep("=", 70))
message("  COMPARISON: ORIGINAL vs DECOMPOSED ATTRIBUTION")
message(strrep("=", 70))

# Load original attribution
orig_attr <- read.csv(here("results", "attribution_table3.csv"))

comparison <- data.frame()
for (sp in SPECIES) {
  sp_name <- SPECIES_NAMES[sp]
  orig_row <- orig_attr[orig_attr$Species == sp_name, ]
  decomp <- all_attribution[[sp]]

  if (nrow(orig_row) == 0 || is.null(decomp)) next

  # Extract numeric values from original (format: "value (sd)")
  extract_mean <- function(x) {
    as.numeric(gsub("\\s*\\(.*", "", as.character(x)))
  }

  comparison <- rbind(comparison, data.frame(
    species = sp_name,
    submodel = "gamma",
    pathway = "climate",
    effect_original = extract_mean(orig_row$Gamma_climate),
    effect_decomposed = decomp$mean$Gamma_climate,
    delta = decomp$mean$Gamma_climate - extract_mean(orig_row$Gamma_climate),
    stringsAsFactors = FALSE
  ))
  comparison <- rbind(comparison, data.frame(
    species = sp_name,
    submodel = "gamma",
    pathway = "landuse",
    effect_original = extract_mean(orig_row$Gamma_landuse),
    effect_decomposed = decomp$mean$Gamma_landuse,
    delta = decomp$mean$Gamma_landuse - extract_mean(orig_row$Gamma_landuse),
    stringsAsFactors = FALSE
  ))
  comparison <- rbind(comparison, data.frame(
    species = sp_name,
    submodel = "epsilon",
    pathway = "climate",
    effect_original = extract_mean(orig_row$Epsilon_climate),
    effect_decomposed = decomp$mean$Epsilon_climate,
    delta = decomp$mean$Epsilon_climate - extract_mean(orig_row$Epsilon_climate),
    stringsAsFactors = FALSE
  ))
  comparison <- rbind(comparison, data.frame(
    species = sp_name,
    submodel = "epsilon",
    pathway = "landuse",
    effect_original = extract_mean(orig_row$Epsilon_landuse),
    effect_decomposed = decomp$mean$Epsilon_landuse,
    delta = decomp$mean$Epsilon_landuse - extract_mean(orig_row$Epsilon_landuse),
    stringsAsFactors = FALSE
  ))
}

write.csv(comparison, here("results", "attribution_comparison_ndvi.csv"),
          row.names = FALSE)
message("\n  Saved: results/attribution_comparison_ndvi.csv")

# Print comparison
message("\n  === ATTRIBUTION COMPARISON (original vs decomposed) ===")
for (i in seq_len(nrow(comparison))) {
  r <- comparison[i, ]
  message("  ", r$species, " | ", r$submodel, "/", r$pathway,
          " | orig=", format(r$effect_original, digits = 3, scientific = TRUE),
          " decomp=", format(r$effect_decomposed, digits = 3, scientific = TRUE),
          " delta=", format(r$delta, digits = 2, scientific = TRUE))
}


###############################################################################
# FINAL SUMMARY
###############################################################################

message("\n\n", strrep("=", 70))
message("  === TASK A: NDVI DECOMPOSITION — COMPLETE ===")
message(strrep("=", 70))

for (sp in names(decomp_results)) {
  dr <- decomp_results[[sp]]
  message("\n  ", sp, ":")
  message("    ΔAIC = ", round(dr$delta_aic, 2))

  # Extract decomposed coefficients
  dc <- coef(dr$model)
  ds <- tryCatch(sqrt(diag(vcov(dr$model))), error = function(e) rep(NA, length(dc)))
  dp <- 2 * pnorm(abs(dc/ds), lower.tail = FALSE)

  for (prefix in c("col", "ext")) {
    cn <- paste0(prefix, "(NDVI_climate)")
    rn <- paste0(prefix, "(NDVI_residual)")
    if (cn %in% names(dc)) {
      message("    ", cn, ": β=", round(dc[cn], 3), " SE=", round(ds[cn], 3),
              " P=", format(dp[cn], digits = 2))
    }
    if (rn %in% names(dc)) {
      message("    ", rn, ": β=", round(dc[rn], 3), " SE=", round(ds[rn], 3),
              " P=", format(dp[rn], digits = 2))
    }
  }
}

# Recommendation
message("\n  RECOMMENDATION:")
message("  Adopt decomposition if NDVI_climate and NDVI_residual show")
message("  different signs or substantially different magnitudes in at")
message("  least one submodel. Otherwise retain original with transparency note.")

message("\n  === FILES GENERATED ===")
message("  results/{sp}_ndvi_decomp_scale.rds         (scaling params)")
message("  results/{sp}_model_ndvi_decomposed.rds     (refitted models)")
message("  results/attribution_ndvi_decomposed_table.csv")
message("  results/attribution_ndvi_decomposed_boot.csv")
message("  results/attribution_comparison_ndvi.csv")
message("\n  Script complete.")
