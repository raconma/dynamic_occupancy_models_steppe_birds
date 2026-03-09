###############################################################################
# refit_revised_models.R
#
# Purpose: Refit ptealc and pteori with revised formulas (NDVI removed),
#          save train_dyn_scale for all species, compute naive vs corrected
#          gamma, and recompute equilibrium occupancy.
#
# Changes from original scripts/4:
#   - ptealc: gamma = ~ pr  (was ~ NDVI + pr)
#   - pteori: epsilon = ~ LC12 + pr  (was ~ LC12 + NDVI + pr)
#
# Also saves train_dyn_scale.rds for ALL 4 species (needed by scripts/8,10).
###############################################################################

library(here)
library(unmarked)
library(dplyr)
library(tidyr)
library(MASS)

source(here("R", "model_configs.R"))

YEARS   <- 2017:2023
T_years <- length(YEARS)
J_reps  <- 10

###############################################################################
# Helper functions (from scripts/4)
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

expand_matrix <- function(mat, J) {
  expanded_mat <- kronecker(mat, matrix(1, nrow = 1, ncol = J))
  original_colnames <- colnames(mat)
  new_colnames <- as.vector(sapply(original_colnames,
                                    function(name) paste0(name, ".", 1:J)))
  colnames(expanded_mat) <- new_colnames
  return(expanded_mat)
}

scale_with_params <- function(x, center, scale_val) {
  (x - center) / scale_val
}

###############################################################################
# Build UMF and scaling for a species (replicates scripts/4 data prep exactly)
###############################################################################
build_umf_and_scale <- function(sp) {
  cfg <- get_model_config(sp)

  occ_wide_clean <- read.csv(
    here("data", "processed_2023", sp, paste0(sp, "_occ_wide_dynamic.csv"))
  )
  message("  Data: ", nrow(occ_wide_clean), " sites x ",
          ncol(occ_wide_clean), " columns")

  # Survey covariates
  dur_cols  <- get_obs_cols(occ_wide_clean, "duration_minutes", YEARS)
  eff_cols  <- get_obs_cols(occ_wide_clean, "effort_distance_km", YEARS)
  obs_cols  <- get_obs_cols(occ_wide_clean, "number_observers", YEARS)
  time_cols <- get_obs_cols(occ_wide_clean, "time_observations_started", YEARS)

  duration  <- as.matrix(occ_wide_clean[, dur_cols])
  effort    <- as.matrix(occ_wide_clean[, eff_cols])
  observers <- as.matrix(occ_wide_clean[, obs_cols])
  time_mat  <- occ_wide_clean[, time_cols]

  # Detection histories
  y_cols <- get_obs_cols(occ_wide_clean, "y", YEARS)
  detections <- occ_wide_clean[, y_cols]
  detections <- detections %>% mutate(across(where(is.logical), as.integer))
  y.cross <- as.matrix(detections)
  y.cross[is.na(time_mat) != is.na(y.cross)] <- NA

  n <- nrow(occ_wide_clean)

  # Site covariates
  siteCovs <- occ_wide_clean[, c("bio1", "bio2", "tree_cover", "grass_cover",
                                   "topo_aspect", "topo_elev")]

  # Yearly covariates
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

  # Standardise ALL covariates (on ALL rows, before any filtering)
  time_mat  <- scale(time_mat)
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

  # Capture training scaling for dynamic covariates
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

  # Save train_dyn_scale
  saveRDS(train_dyn_scale,
          here("results", paste0(sp, "_train_dyn_scale.rds")))
  message("  Saved train_dyn_scale for ", sp)

  # Build year factor
  years_df <- data.frame(matrix(rep(YEARS, each = n), n, T_years))
  years_df <- data.frame(lapply(years_df, as.factor))

  # Expand to observation level
  NDVI_obs        <- expand_matrix(NDVI, J_reps)
  pr_obs          <- expand_matrix(pr, J_reps)
  topo_aspect_obs <- expand_matrix(siteCovs[, "topo_aspect", drop = FALSE],
                                    J_reps * T_years)
  topo_elev_obs   <- expand_matrix(siteCovs[, "topo_elev", drop = FALSE],
                                    J_reps * T_years)

  occ_umf <- unmarkedMultFrame(
    y = y.cross,
    siteCovs = data.frame(siteCovs),
    yearlySiteCovs = list(
      years  = years_df,
      EVI    = EVI, NDVI = NDVI, pr = pr, tmmn = tmmn, tmmx = tmmx,
      Land_Cover_Type_1_Percent_Class_0  = Land_Cover_Type_1_Percent_Class_0,
      Land_Cover_Type_1_Percent_Class_6  = Land_Cover_Type_1_Percent_Class_6,
      Land_Cover_Type_1_Percent_Class_7  = Land_Cover_Type_1_Percent_Class_7,
      Land_Cover_Type_1_Percent_Class_10 = Land_Cover_Type_1_Percent_Class_10,
      Land_Cover_Type_1_Percent_Class_12 = Land_Cover_Type_1_Percent_Class_12,
      Land_Cover_Type_1_Percent_Class_13 = Land_Cover_Type_1_Percent_Class_13,
      Land_Cover_Type_1_Percent_Class_14 = Land_Cover_Type_1_Percent_Class_14
    ),
    obsCovs = list(
      duration = duration, effort = effort, observers = observers, time = time_mat,
      NDVI_obs = NDVI_obs, pr_obs = pr_obs,
      topo_aspect_obs = topo_aspect_obs, topo_elev_obs = topo_elev_obs
    ),
    numPrimary = T_years
  )

  list(umf = occ_umf, cfg = cfg, train_dyn_scale = train_dyn_scale,
       occ_wide_clean = occ_wide_clean)
}

###############################################################################
# PART 1: Save train_dyn_scale for otitar and tettet (no refit needed)
###############################################################################
message("\n=== Saving train_dyn_scale for otitar and tettet ===")

for (sp in c("otitar", "tettet")) {
  message("\n  Processing ", sp, "...")
  sp_data <- build_umf_and_scale(sp)
  message("  Done: train_dyn_scale saved for ", sp)
}

###############################################################################
# PART 2: Refit ptealc with gamma = ~pr (NDVI removed)
###############################################################################
message("\n=== Refitting ptealc ===")

ptealc_data <- build_umf_and_scale("ptealc")

message("  Fitting revised colext model (gamma = ~pr)...")
ptealc_mod_revised <- colext(
  psiformula     = ptealc_data$cfg$psi_formula,
  gammaformula   = ptealc_data$cfg$gamma_formula,  # now ~ pr
  epsilonformula = ptealc_data$cfg$epsilon_formula,
  pformula       = ptealc_data$cfg$p_formula,
  data = ptealc_data$umf
)

message("  ptealc revised AIC = ", round(ptealc_mod_revised@AIC, 2))

# Compare with original
ptealc_orig <- readRDS(here("results", "ptealc_model_object.rds"))
message("  ptealc original AIC = ", round(ptealc_orig@AIC, 2))
message("  Delta AIC = ", round(ptealc_mod_revised@AIC - ptealc_orig@AIC, 2))

# Print coefficients
cat("\n  Revised ptealc coefficients:\n")
print(summary(ptealc_mod_revised))

# Save revised model
saveRDS(ptealc_mod_revised, here("results", "ptealc_model_revised.rds"))
# Also overwrite the main model object (this IS the new canonical model)
saveRDS(ptealc_mod_revised, here("results", "ptealc_model_object.rds"))
message("  ptealc revised model saved.")

# Save summary
sink(here("results", "ptealc_model_summary_revised.txt"))
cat("Species: ptealc (REVISED - NDVI removed from gamma)\n")
cat("Date:", as.character(Sys.time()), "\n")
cat("Original AIC:", round(ptealc_orig@AIC, 2), "\n")
cat("Revised AIC:", round(ptealc_mod_revised@AIC, 2), "\n")
cat("Delta AIC:", round(ptealc_mod_revised@AIC - ptealc_orig@AIC, 2), "\n\n")
summary(ptealc_mod_revised)
sink()

###############################################################################
# PART 3: Refit pteori with epsilon = ~LC12 + pr (NDVI removed)
###############################################################################
message("\n=== Refitting pteori ===")

pteori_data <- build_umf_and_scale("pteori")

message("  Fitting revised colext model (epsilon = ~LC12 + pr)...")
pteori_mod_revised <- colext(
  psiformula     = pteori_data$cfg$psi_formula,
  gammaformula   = pteori_data$cfg$gamma_formula,
  epsilonformula = pteori_data$cfg$epsilon_formula,  # now ~ LC12 + pr
  pformula       = pteori_data$cfg$p_formula,
  data = pteori_data$umf
)

message("  pteori revised AIC = ", round(pteori_mod_revised@AIC, 2))

# Compare with original
pteori_orig <- readRDS(here("results", "pteori_model_object.rds"))
message("  pteori original AIC = ", round(pteori_orig@AIC, 2))
message("  Delta AIC = ", round(pteori_mod_revised@AIC - pteori_orig@AIC, 2))

# Print coefficients
cat("\n  Revised pteori coefficients:\n")
print(summary(pteori_mod_revised))

# Check LC12 coefficient recovery
pteori_ext_coefs <- coef(pteori_mod_revised)
lc12_idx <- grep("ext\\(Land_Cover", names(pteori_ext_coefs))
if (length(lc12_idx) > 0) {
  lc12_coef <- pteori_ext_coefs[lc12_idx]
  message("  LC12 coefficient in epsilon: ", round(lc12_coef, 4))
  if (lc12_coef < -0.09) {
    message("  SUCCESS: LC12 recovered to more negative value than -0.09")
  } else {
    message("  NOTE: LC12 coefficient did not recover as expected")
  }
}

# Save revised model
saveRDS(pteori_mod_revised, here("results", "pteori_epsilon_revised.rds"))
saveRDS(pteori_mod_revised, here("results", "pteori_model_object.rds"))
message("  pteori revised model saved.")

# Save summary
sink(here("results", "pteori_model_summary_revised.txt"))
cat("Species: pteori (REVISED - NDVI removed from epsilon)\n")
cat("Date:", as.character(Sys.time()), "\n")
cat("Original AIC:", round(pteori_orig@AIC, 2), "\n")
cat("Revised AIC:", round(pteori_mod_revised@AIC, 2), "\n")
cat("Delta AIC:", round(pteori_mod_revised@AIC - pteori_orig@AIC, 2), "\n\n")
summary(pteori_mod_revised)
sink()

###############################################################################
# PART 4: Naive vs corrected gamma comparison
###############################################################################
message("\n=== Naive vs corrected gamma comparison ===")

species_codes <- c("otitar", "ptealc", "pteori", "tettet")
naive_results <- list()

for (sp in species_codes) {
  message("\n  Processing ", sp, "...")

  occ_wide <- read.csv(
    here("data", "processed_2023", sp, paste0(sp, "_occ_wide_dynamic.csv"))
  )

  # Build 3D detection array: sites x reps x years
  y_cols <- c()
  for (yr in YEARS) {
    for (k in 1:J_reps) {
      cn <- paste0("y.", k, ".", yr)
      if (cn %in% names(occ_wide)) y_cols <- c(y_cols, cn)
    }
  }

  y_mat <- as.matrix(occ_wide[, y_cols])
  y_3d <- array(y_mat, dim = c(nrow(occ_wide), J_reps, T_years))

  # Per-site per-year: was the species detected at least once?
  detected <- apply(y_3d > 0, c(1, 3), max, na.rm = TRUE)
  # Replace -Inf (all NA) with NA
  detected[!is.finite(detected)] <- NA

  # Naive colonisation: undetected in year t, detected in year t+1
  n_col_events <- 0
  n_col_opportunities <- 0
  for (t in 1:(T_years - 1)) {
    unoccupied_t  <- (detected[, t] == 0)
    occupied_t1   <- (detected[, t + 1] == 1)
    valid <- !is.na(unoccupied_t) & !is.na(occupied_t1)
    n_col_events <- n_col_events + sum(unoccupied_t & occupied_t1 & valid, na.rm = TRUE)
    n_col_opportunities <- n_col_opportunities + sum(unoccupied_t & valid, na.rm = TRUE)
  }

  naive_gamma <- if (n_col_opportunities > 0) n_col_events / n_col_opportunities else NA

  # Naive extinction: detected in year t, undetected in year t+1
  n_ext_events <- 0
  n_ext_opportunities <- 0
  for (t in 1:(T_years - 1)) {
    occupied_t    <- (detected[, t] == 1)
    unoccupied_t1 <- (detected[, t + 1] == 0)
    valid <- !is.na(occupied_t) & !is.na(unoccupied_t1)
    n_ext_events <- n_ext_events + sum(occupied_t & unoccupied_t1 & valid, na.rm = TRUE)
    n_ext_opportunities <- n_ext_opportunities + sum(occupied_t & valid, na.rm = TRUE)
  }

  naive_epsilon <- if (n_ext_opportunities > 0) n_ext_events / n_ext_opportunities else NA

  # Model-corrected values (mean across sites and years)
  mod <- readRDS(here("results", paste0(sp, "_model_object.rds")))
  corrected_gamma   <- mean(predict(mod, type = "col")$Predicted, na.rm = TRUE)
  corrected_epsilon <- mean(predict(mod, type = "ext")$Predicted, na.rm = TRUE)

  naive_results[[sp]] <- data.frame(
    species = sp,
    naive_gamma = round(naive_gamma, 4),
    corrected_gamma = round(corrected_gamma, 4),
    gamma_ratio = round(corrected_gamma / naive_gamma, 2),
    naive_epsilon = round(naive_epsilon, 4),
    corrected_epsilon = round(corrected_epsilon, 4),
    epsilon_ratio = round(corrected_epsilon / naive_epsilon, 2),
    n_col_events = n_col_events,
    n_col_opportunities = n_col_opportunities,
    n_ext_events = n_ext_events,
    n_ext_opportunities = n_ext_opportunities,
    stringsAsFactors = FALSE
  )

  message("  ", sp, ": naive_gamma=", round(naive_gamma, 4),
          " corrected_gamma=", round(corrected_gamma, 4),
          " ratio=", round(corrected_gamma / naive_gamma, 2))
  message("        naive_epsilon=", round(naive_epsilon, 4),
          " corrected_epsilon=", round(corrected_epsilon, 4),
          " ratio=", round(corrected_epsilon / naive_epsilon, 2))
}

naive_table <- bind_rows(naive_results)
write.csv(naive_table,
          here("results", "naive_vs_corrected_gamma.csv"),
          row.names = FALSE)
message("\n  Saved: results/naive_vs_corrected_gamma.csv")
print(naive_table)

###############################################################################
# PART 5: Recompute equilibrium occupancy with revised models
###############################################################################
message("\n=== Recomputing equilibrium occupancy ===")

species_all <- c("otitar", "ptealc", "pteori", "tettet")
species_names <- c("Otis tarda", "Pterocles alchata",
                   "Pterocles orientalis", "Tetrax tetrax")

compute_equilibrium_ci <- function(model, sp_code, n_boot = 5000) {
  coefs <- coef(model)
  V     <- vcov(model)
  vcov_ok <- all(is.finite(V))
  if (!vcov_ok) warning(paste(sp_code, ": vcov contains non-finite values"))

  set.seed(42)
  sims <- MASS::mvrnorm(n_boot, mu = coefs, Sigma = V)
  col_idx <- grep("col\\(Int\\)", names(coefs))
  ext_idx <- grep("ext\\(Int\\)", names(coefs))

  gamma_sims   <- plogis(sims[, col_idx])
  epsilon_sims <- plogis(sims[, ext_idx])
  psi_star          <- gamma_sims / (gamma_sims + epsilon_sims)
  recolonisation_yr <- 1 / gamma_sims
  persistence_yr    <- 1 / epsilon_sims

  data.frame(
    species           = sp_code,
    vcov_ok           = vcov_ok,
    gamma_intercept   = round(coefs[col_idx], 4),
    epsilon_intercept = round(coefs[ext_idx], 4),
    gamma_baseline_pct    = round(plogis(coefs[col_idx]) * 100, 4),
    epsilon_baseline_pct  = round(plogis(coefs[ext_idx]) * 100, 4),
    psi_star_med      = round(median(psi_star) * 100, 4),
    psi_star_lo       = round(quantile(psi_star, 0.025) * 100, 4),
    psi_star_hi       = round(quantile(psi_star, 0.975) * 100, 4),
    recol_yr_med      = round(median(recolonisation_yr)),
    recol_yr_lo       = round(quantile(recolonisation_yr, 0.025)),
    recol_yr_hi       = round(quantile(recolonisation_yr, 0.975)),
    persist_yr_med    = round(median(persistence_yr)),
    persist_yr_lo     = round(quantile(persistence_yr, 0.025)),
    persist_yr_hi     = round(quantile(persistence_yr, 0.975)),
    stringsAsFactors = FALSE
  )
}

results_list <- list()
for (i in seq_along(species_all)) {
  sp <- species_all[i]
  message("\n  ", sp, " (", species_names[i], "):")
  mod <- readRDS(here("results", paste0(sp, "_model_object.rds")))
  message("  AIC = ", round(mod@AIC, 2))

  res <- compute_equilibrium_ci(mod, sp)
  results_list[[sp]] <- res

  message("  psi* = ", res$psi_star_med, "% [",
          res$psi_star_lo, ", ", res$psi_star_hi, "]")

  # Flag ptealc
  if (sp == "ptealc") {
    ci_width <- res$psi_star_hi - res$psi_star_lo
    if (ci_width > 50) {
      message("  WARNING: ptealc psi* CI width = ", round(ci_width, 1),
              "% — uninformative due to epsilon separation")
      message("  EXCLUDE from Abstract numbers or add explicit caveat")
    }
  }
}

eq_table <- do.call(rbind, results_list)
write.csv(eq_table,
          here("results", "equilibrium_occupancy_table.csv"),
          row.names = FALSE)
message("\n  Saved: results/equilibrium_occupancy_table.csv")

cat("\n==========================================\n")
cat("REVISED SUMMARY FOR ABSTRACT:\n")
cat("==========================================\n")
# Exclude ptealc from reporting due to uninformative CI
reportable <- eq_table[eq_table$species != "ptealc", ]
cat("Reportable species (excluding P. alchata):\n")
print(reportable[, c("species", "psi_star_med", "psi_star_lo", "psi_star_hi")])
cat("\nP. alchata (EXCLUDED from headline number):\n")
print(eq_table[eq_table$species == "ptealc",
               c("species", "psi_star_med", "psi_star_lo", "psi_star_hi")])

message("\n=== All refitting tasks complete ===")
