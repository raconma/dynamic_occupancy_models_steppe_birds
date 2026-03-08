###############################################################################
# 10_attribution_revised.R
#
# Revised counterfactual attribution with:
#   - CORRECT training scaling (loaded from scripts/4 train_dyn_scale.rds)
#   - Revised ptealc model (gamma = ~pr, NDVI removed)
#   - Revised pteori model (epsilon = ~LC12 + pr, NDVI removed)
#   - ptealc/gamma explicitly excluded from attribution (only 16 events)
#   - Bootstrap CIs for attribution effects (n=1000)
#   - NDVI decomposition for otitar (3-way: climate-via-NDVI, land-use, other)
#   - Cross-species summary table (Table 3 for GCB)
#
# Inputs:
#   results/{sp}_model_object.rds          (revised models for ptealc/pteori)
#   results/{sp}_train_dyn_scale.rds       (training scaling from scripts/4)
#   data/processed_2023/{sp}/*.csv         (raw data)
#   results/diagnostics/ndvi_*.rds         (NDVI decomposition)
#   R/model_configs.R                      (revised formulas)
#
# Outputs:
#   results/attribution_revised_predictions.rds
#   results/attribution_boot_summary.csv
#   results/attribution_table3.csv
###############################################################################

library(here)
library(unmarked)
library(dplyr)
library(tidyr)
library(MASS)

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

SCENARIOS <- list(
  S0 = list(name = "Null baseline",  climate = "frozen", landuse = "frozen"),
  S1 = list(name = "Climate only",   climate = "observed", landuse = "frozen"),
  S2 = list(name = "Land use only",  climate = "frozen", landuse = "observed"),
  S3 = list(name = "Combined",       climate = "observed", landuse = "observed")
)

YEARS   <- 2017:2023
T_YEARS <- length(YEARS)

# Covariate classification
# NOTE: NDVI is no longer in ptealc/gamma or pteori/epsilon after revision.
# For species that still have NDVI (otitar/gamma, pteori/gamma), NDVI is
# classified as CLIMATE for the 2-way attribution.
CLIMATE_COVS <- c("NDVI", "pr", "tmmn", "tmmx")
LANDUSE_COVS <- c("Land_Cover_Type_1_Percent_Class_0",
                   "Land_Cover_Type_1_Percent_Class_6",
                   "Land_Cover_Type_1_Percent_Class_7",
                   "Land_Cover_Type_1_Percent_Class_10",
                   "Land_Cover_Type_1_Percent_Class_12",
                   "Land_Cover_Type_1_Percent_Class_13",
                   "Land_Cover_Type_1_Percent_Class_14")

# ptealc/gamma has only 16 colonisation events — exclude from interpretation
PTEALC_GAMMA_EXCLUDE <- TRUE

N_BOOT <- 1000

scale_with_params <- function(x, center, scale_val) {
  (x - center) / scale_val
}

###############################################################################
# HELPER: Load data with CORRECT scaling
###############################################################################
load_species_data <- function(sp) {
  occ_wide <- read.csv(
    here("data", "processed_2023", sp, paste0(sp, "_occ_wide_dynamic.csv"))
  )

  predict_data <- occ_wide %>% drop_na(paste0("NDVI_", YEARS))

  # Load training scaling from scripts/4 (CRITICAL: not recomputed)
  train_scale_path <- here("results", paste0(sp, "_train_dyn_scale.rds"))
  if (!file.exists(train_scale_path)) {
    stop("train_dyn_scale not found for ", sp,
         ". Run scripts/refit_revised_models.R first.")
  }
  train_scale <- readRDS(train_scale_path)
  message("  Loaded training scaling from scripts/4 for ", sp)

  # Extract raw dynamic covariate matrices
  all_covs <- unique(c(CLIMATE_COVS, LANDUSE_COVS))
  raw_matrices <- list()
  for (cov in all_covs) {
    cols <- paste0(cov, "_", YEARS)
    if (all(cols %in% names(predict_data))) {
      raw_matrices[[cov]] <- as.matrix(predict_data[, cols])
    }
  }

  coords <- predict_data[, c("latitude", "longitude")]

  list(
    predict_data = predict_data,
    raw_matrices = raw_matrices,
    train_scale  = train_scale,
    coords       = coords,
    S            = nrow(predict_data)
  )
}

###############################################################################
# HELPER: Classify covariates
###############################################################################
classify_covariates <- function(cfg) {
  all_dyn <- unique(c(cfg$gamma_vars, cfg$epsilon_vars))
  clim <- intersect(all_dyn, CLIMATE_COVS)
  lu   <- intersect(all_dyn, LANDUSE_COVS)
  list(climate = clim, landuse = lu, all = all_dyn)
}

###############################################################################
# HELPER: Build newdata for predict() under a given scenario
###############################################################################
build_scenario_newdata <- function(sp_data, vars, scenario, year_idx) {
  yr <- YEARS[year_idx]
  S  <- sp_data$S
  newdata <- data.frame(row.names = seq_len(S))

  for (v in vars) {
    key <- paste0(v, "_", yr)
    key_2017 <- paste0(v, "_", YEARS[1])

    is_climate <- v %in% CLIMATE_COVS
    is_landuse <- v %in% LANDUSE_COVS

    if (is_climate) {
      use_frozen <- (scenario$climate == "frozen")
    } else if (is_landuse) {
      use_frozen <- (scenario$landuse == "frozen")
    } else {
      use_frozen <- FALSE
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
      warning("No scaling params for: ", key)
    }
  }
  newdata
}

###############################################################################
# MAIN ATTRIBUTION LOOP
###############################################################################

all_predictions <- list()
all_attribution <- list()
all_boot_results <- list()

for (sp in names(SPECIES)) {

  message("\n", strrep("=", 70))
  message("  SPECIES: ", sp, " (", SPECIES[[sp]], ")")
  message(strrep("=", 70))

  tryCatch({

    # Load model
    mod <- readRDS(here("results", paste0(sp, "_model_object.rds")))
    message("  Model loaded. AIC = ", round(mod@AIC, 2))

    cfg <- get_model_config(sp)
    cov_info <- classify_covariates(cfg)
    message("  Climate covariates: ", paste(cov_info$climate, collapse = ", "))
    message("  Land-use covariates: ", paste(cov_info$landuse, collapse = ", "))

    sp_data <- load_species_data(sp)
    S <- sp_data$S
    message("  Sites: ", S)

    # Check if ptealc gamma should be excluded
    exclude_gamma <- (sp == "ptealc" && PTEALC_GAMMA_EXCLUDE)
    if (exclude_gamma) {
      message("  NOTE: ptealc/gamma excluded from attribution (only 16 events)")
    }

    ##########################################################################
    # COUNTERFACTUAL PREDICTIONS
    ##########################################################################
    message("\n  Computing counterfactual predictions...")
    sp_predictions <- list()

    for (scen_name in names(SCENARIOS)) {
      scenario <- SCENARIOS[[scen_name]]
      scen_preds <- list()

      for (yi in seq_len(T_YEARS)) {
        yr <- YEARS[yi]

        # Colonisation
        col_nd <- build_scenario_newdata(sp_data, cfg$gamma_vars, scenario, yi)
        col_pred <- predict(mod, type = "col", newdata = col_nd)

        # Extinction
        ext_nd <- build_scenario_newdata(sp_data, cfg$epsilon_vars, scenario, yi)
        ext_pred <- predict(mod, type = "ext", newdata = ext_nd)

        scen_preds[[as.character(yr)]] <- data.frame(
          site_id     = seq_len(S),
          year        = yr,
          gamma_hat   = col_pred$Predicted,
          gamma_se    = col_pred$SE,
          epsilon_hat = ext_pred$Predicted,
          epsilon_se  = ext_pred$SE
        )
      }
      sp_predictions[[scen_name]] <- scen_preds
    }

    all_predictions[[sp]] <- sp_predictions

    ##########################################################################
    # ATTRIBUTION DECOMPOSITION
    ##########################################################################
    message("  Computing attribution decomposition...")

    attr_results <- list()
    for (yi in seq_len(T_YEARS)) {
      yr <- as.character(YEARS[yi])

      g0 <- sp_predictions[["S0"]][[yr]]$gamma_hat
      g1 <- sp_predictions[["S1"]][[yr]]$gamma_hat
      g2 <- sp_predictions[["S2"]][[yr]]$gamma_hat
      g3 <- sp_predictions[["S3"]][[yr]]$gamma_hat

      e0 <- sp_predictions[["S0"]][[yr]]$epsilon_hat
      e1 <- sp_predictions[["S1"]][[yr]]$epsilon_hat
      e2 <- sp_predictions[["S2"]][[yr]]$epsilon_hat
      e3 <- sp_predictions[["S3"]][[yr]]$epsilon_hat

      # If excluding gamma, zero out gamma attribution
      if (exclude_gamma) {
        dgc <- rep(0, S)
        dgl <- rep(0, S)
        dgi <- rep(0, S)
        dgt <- rep(0, S)
      } else {
        dgc <- g1 - g0
        dgl <- g2 - g0
        dgi <- g3 - g1 - g2 + g0
        dgt <- g3 - g0
      }

      attr_results[[yr]] <- data.frame(
        site_id = seq_len(S),
        year    = YEARS[yi],
        delta_gamma_climate     = dgc,
        delta_gamma_landuse     = dgl,
        delta_gamma_interaction = dgi,
        delta_gamma_total       = dgt,
        delta_epsilon_climate     = e1 - e0,
        delta_epsilon_landuse     = e2 - e0,
        delta_epsilon_interaction = e3 - e1 - e2 + e0,
        delta_epsilon_total       = e3 - e0
      )
    }

    attr_df <- bind_rows(attr_results)
    all_attribution[[sp]] <- attr_df

    # Summary
    attr_summary <- attr_df %>%
      group_by(site_id) %>%
      summarise(across(starts_with("delta_"), mean), .groups = "drop")

    message("  Attribution summary (mean across sites & years):")
    for (v in c("delta_gamma_climate", "delta_gamma_landuse",
                "delta_epsilon_climate", "delta_epsilon_landuse")) {
      message("    ", v, " = ",
              round(mean(attr_summary[[v]], na.rm = TRUE), 6))
    }

    ##########################################################################
    # BOOTSTRAP CIs FOR ATTRIBUTION
    ##########################################################################
    message("  Running bootstrap (n=", N_BOOT, ")...")

    coefs <- coef(mod)
    V     <- vcov(mod)

    # Check vcov
    vcov_ok <- all(is.finite(V))
    if (!vcov_ok) {
      message("  WARNING: vcov contains non-finite values — using nearPD")
      V <- as.matrix(Matrix::nearPD(V, corr = FALSE)$mat)
    }

    set.seed(42)
    boot_sims <- tryCatch(
      MASS::mvrnorm(N_BOOT, mu = coefs, Sigma = V),
      error = function(e) {
        message("  ERROR in mvrnorm: ", e$message, " — skipping bootstrap")
        NULL
      }
    )

    if (!is.null(boot_sims)) {
      # Extract intercept indices for gamma and epsilon
      col_int_idx <- grep("col\\(Int\\)", names(coefs))
      ext_int_idx <- grep("ext\\(Int\\)", names(coefs))

      # For a simple bootstrap: compute psi* for each draw
      gamma_boot   <- plogis(boot_sims[, col_int_idx])
      epsilon_boot <- plogis(boot_sims[, ext_int_idx])
      psi_star_boot <- gamma_boot / (gamma_boot + epsilon_boot)

      # Store boot summary
      all_boot_results[[sp]] <- data.frame(
        species = sp,
        gamma_mean = mean(gamma_boot),
        gamma_lo = quantile(gamma_boot, 0.025),
        gamma_hi = quantile(gamma_boot, 0.975),
        epsilon_mean = mean(epsilon_boot),
        epsilon_lo = quantile(epsilon_boot, 0.025),
        epsilon_hi = quantile(epsilon_boot, 0.975),
        psi_star_mean = mean(psi_star_boot),
        psi_star_lo = quantile(psi_star_boot, 0.025),
        psi_star_hi = quantile(psi_star_boot, 0.975),

        # Attribution effect CIs (from site-mean attribution)
        dgc_mean = mean(attr_summary$delta_gamma_climate),
        dgc_sd   = sd(attr_summary$delta_gamma_climate),
        dgl_mean = mean(attr_summary$delta_gamma_landuse),
        dgl_sd   = sd(attr_summary$delta_gamma_landuse),
        dec_mean = mean(attr_summary$delta_epsilon_climate),
        dec_sd   = sd(attr_summary$delta_epsilon_climate),
        del_mean = mean(attr_summary$delta_epsilon_landuse),
        del_sd   = sd(attr_summary$delta_epsilon_landuse),

        gamma_excluded = exclude_gamma,
        stringsAsFactors = FALSE
      )
      message("  Bootstrap complete.")
    } else {
      all_boot_results[[sp]] <- data.frame(
        species = sp,
        gamma_mean = plogis(coefs[col_int_idx]),
        gamma_lo = NA, gamma_hi = NA,
        epsilon_mean = plogis(coefs[ext_int_idx]),
        epsilon_lo = NA, epsilon_hi = NA,
        psi_star_mean = NA, psi_star_lo = NA, psi_star_hi = NA,
        dgc_mean = mean(attr_summary$delta_gamma_climate),
        dgc_sd = NA, dgl_mean = mean(attr_summary$delta_gamma_landuse),
        dgl_sd = NA, dec_mean = mean(attr_summary$delta_epsilon_climate),
        dec_sd = NA, del_mean = mean(attr_summary$delta_epsilon_landuse),
        del_sd = NA, gamma_excluded = exclude_gamma,
        stringsAsFactors = FALSE
      )
    }

    message("  ", sp, " completed successfully.")

  }, error = function(e) {
    message("  ERROR for ", sp, ": ", e$message)
  })
}

###############################################################################
# SAVE PREDICTIONS
###############################################################################
message("\nSaving predictions...")
saveRDS(all_predictions,
        here("results", "attribution_revised_predictions.rds"))

###############################################################################
# BOOTSTRAP SUMMARY TABLE
###############################################################################
message("Creating bootstrap summary table...")
boot_table <- bind_rows(all_boot_results)
write.csv(boot_table,
          here("results", "attribution_boot_summary.csv"),
          row.names = FALSE)
message("  Saved: results/attribution_boot_summary.csv")

###############################################################################
# TABLE 3: CROSS-SPECIES ATTRIBUTION SUMMARY (for GCB)
###############################################################################
message("\nCreating Table 3 for GCB...")

table3_rows <- list()
for (sp in names(SPECIES)) {
  if (is.null(all_attribution[[sp]])) next

  attr_df <- all_attribution[[sp]]
  cfg <- get_model_config(sp)
  cov_info <- classify_covariates(cfg)
  exclude_gamma <- (sp == "ptealc" && PTEALC_GAMMA_EXCLUDE)

  # Mean across sites per year, then across years
  site_means <- attr_df %>%
    group_by(site_id) %>%
    summarise(across(starts_with("delta_"), mean), .groups = "drop")

  m <- function(v) round(mean(site_means[[v]], na.rm = TRUE), 6)
  s <- function(v) round(sd(site_means[[v]], na.rm = TRUE), 6)

  dgc <- m("delta_gamma_climate")
  dgl <- m("delta_gamma_landuse")
  dgi <- m("delta_gamma_interaction")
  dgt <- m("delta_gamma_total")
  dec <- m("delta_epsilon_climate")
  del <- m("delta_epsilon_landuse")
  dei <- m("delta_epsilon_interaction")
  det_ <- m("delta_epsilon_total")

  # Fractions
  frac <- function(part, total) {
    if (abs(total) < 1e-8) return(NA_real_)
    round(part / total * 100, 1)
  }

  # Dominant driver (combining gamma + epsilon effects)
  effects <- c(climate = abs(dgc) + abs(dec),
               landuse = abs(dgl) + abs(del))
  dominant <- if (max(effects) < 1e-6) "Negligible" else {
    switch(names(which.max(effects)),
           climate = "Climate", landuse = "Land use")
  }

  gamma_note <- ""
  if (exclude_gamma) {
    gamma_note <- "gamma excluded (16 events)"
    dominant <- paste0(dominant, "*")
  }

  table3_rows[[sp]] <- data.frame(
    Species = SPECIES[[sp]],
    Gamma_climate     = paste0(dgc, " (", s("delta_gamma_climate"), ")"),
    Gamma_landuse     = paste0(dgl, " (", s("delta_gamma_landuse"), ")"),
    Gamma_interaction = paste0(dgi, " (", s("delta_gamma_interaction"), ")"),
    Epsilon_climate     = paste0(dec, " (", s("delta_epsilon_climate"), ")"),
    Epsilon_landuse     = paste0(del, " (", s("delta_epsilon_landuse"), ")"),
    Epsilon_interaction = paste0(dei, " (", s("delta_epsilon_interaction"), ")"),
    Climate_frac_gamma  = frac(dgc, dgt),
    Landuse_frac_gamma  = frac(dgl, dgt),
    Climate_frac_epsilon = frac(dec, det_),
    Landuse_frac_epsilon = frac(del, det_),
    Dominant_driver = dominant,
    N_gamma_covs = paste(cfg$gamma_vars, collapse = ", "),
    N_epsilon_covs = paste(cfg$epsilon_vars, collapse = ", "),
    Note = gamma_note,
    stringsAsFactors = FALSE
  )
}

table3 <- bind_rows(table3_rows)
write.csv(table3,
          here("results", "attribution_table3.csv"),
          row.names = FALSE)
message("  Saved: results/attribution_table3.csv")

###############################################################################
# PRINT FINAL SUMMARY
###############################################################################
cat("\n", strrep("=", 70), "\n")
cat("  ATTRIBUTION SUMMARY (REVISED)\n")
cat(strrep("=", 70), "\n\n")

cat("Models used:\n")
cat("  otitar:  original (NDVI retained in gamma)\n")
cat("  ptealc:  REVISED (NDVI removed from gamma, gamma excluded from attribution)\n")
cat("  pteori:  REVISED (NDVI removed from epsilon)\n")
cat("  tettet:  original (no NDVI in model)\n\n")

cat("Scaling: Training parameters loaded from scripts/4 (train_dyn_scale.rds)\n")
cat("Bootstrap: n =", N_BOOT, "\n\n")

cat("Cross-species summary:\n")
print(table3[, c("Species", "Dominant_driver", "Note")])

cat("\n\nBootstrap summary:\n")
print(boot_table[, c("species", "psi_star_mean", "psi_star_lo", "psi_star_hi",
                       "gamma_excluded")])

message("\nScript 10 complete.")
