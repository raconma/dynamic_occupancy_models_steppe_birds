###############################################################################
# execute_decisions_v4.R
#
# Executes all decisions from decisions_gcb_v4.md:
#   D1: Revert pteori/epsilon to original (with NDVI), save no-NDVI as sensitivity
#   D2: Resolved by D1 (gamma recovers convergence)
#   D3: No code action (ptealc epsilon kept as-is)
#   D5: Compute extinction debt (current occ − ψ*)
#   D6: Produce reportable ψ* summary
#   D7: Prepare data for isocline plot
#
# Also re-runs attribution with NDVI classified as "climate-adjacent"
###############################################################################

library(here)
library(unmarked)
library(dplyr)
library(tidyr)
library(MASS)
library(Matrix)

source(here("R", "model_configs.R"))

YEARS   <- 2017:2023
T_years <- length(YEARS)
J_reps  <- 10

# Helper functions
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
  colnames(expanded_mat) <- as.vector(sapply(colnames(mat),
                                              function(name) paste0(name, ".", 1:J)))
  return(expanded_mat)
}

scale_with_params <- function(x, center, scale_val) {
  (x - center) / scale_val
}

build_umf_and_scale <- function(sp) {
  cfg <- get_model_config(sp)
  occ_wide_clean <- read.csv(
    here("data", "processed_2023", sp, paste0(sp, "_occ_wide_dynamic.csv"))
  )

  dur_cols  <- get_obs_cols(occ_wide_clean, "duration_minutes", YEARS)
  eff_cols  <- get_obs_cols(occ_wide_clean, "effort_distance_km", YEARS)
  obs_cols  <- get_obs_cols(occ_wide_clean, "number_observers", YEARS)
  time_cols <- get_obs_cols(occ_wide_clean, "time_observations_started", YEARS)

  duration  <- as.matrix(occ_wide_clean[, dur_cols])
  effort    <- as.matrix(occ_wide_clean[, eff_cols])
  observers <- as.matrix(occ_wide_clean[, obs_cols])
  time_mat  <- occ_wide_clean[, time_cols]

  y_cols <- get_obs_cols(occ_wide_clean, "y", YEARS)
  detections <- occ_wide_clean[, y_cols]
  detections <- detections %>% mutate(across(where(is.logical), as.integer))
  y.cross <- as.matrix(detections)
  y.cross[is.na(time_mat) != is.na(y.cross)] <- NA

  n <- nrow(occ_wide_clean)
  siteCovs <- occ_wide_clean[, c("bio1","bio2","tree_cover","grass_cover","topo_aspect","topo_elev")]

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

  time_mat  <- scale(time_mat);  duration  <- scale(duration)
  effort    <- scale(effort);    observers <- scale(observers)
  siteCovs  <- scale(siteCovs)

  EVI  <- scale(EVI)
  Land_Cover_Type_1_Percent_Class_0  <- scale(Land_Cover_Type_1_Percent_Class_0)
  Land_Cover_Type_1_Percent_Class_6  <- scale(Land_Cover_Type_1_Percent_Class_6)
  Land_Cover_Type_1_Percent_Class_7  <- scale(Land_Cover_Type_1_Percent_Class_7)
  Land_Cover_Type_1_Percent_Class_10 <- scale(Land_Cover_Type_1_Percent_Class_10)
  Land_Cover_Type_1_Percent_Class_12 <- scale(Land_Cover_Type_1_Percent_Class_12)
  Land_Cover_Type_1_Percent_Class_13 <- scale(Land_Cover_Type_1_Percent_Class_13)
  Land_Cover_Type_1_Percent_Class_14 <- scale(Land_Cover_Type_1_Percent_Class_14)
  NDVI <- scale(NDVI); pr <- scale(pr); tmmn <- scale(tmmn); tmmx <- scale(tmmx)

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
  saveRDS(train_dyn_scale, here("results", paste0(sp, "_train_dyn_scale.rds")))

  years_df <- data.frame(matrix(rep(YEARS, each = n), n, T_years))
  years_df <- data.frame(lapply(years_df, as.factor))

  NDVI_obs        <- expand_matrix(NDVI, J_reps)
  pr_obs          <- expand_matrix(pr, J_reps)
  topo_aspect_obs <- expand_matrix(siteCovs[, "topo_aspect", drop = FALSE], J_reps * T_years)
  topo_elev_obs   <- expand_matrix(siteCovs[, "topo_elev", drop = FALSE], J_reps * T_years)

  occ_umf <- unmarkedMultFrame(
    y = y.cross, siteCovs = data.frame(siteCovs),
    yearlySiteCovs = list(
      years = years_df, EVI = EVI, NDVI = NDVI, pr = pr, tmmn = tmmn, tmmx = tmmx,
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
# PART 1: Refit pteori with ORIGINAL epsilon (with NDVI)
###############################################################################
message("\n", strrep("=", 70))
message("  DECISION 1: Refit pteori with original epsilon (~LC12 + NDVI + pr)")
message(strrep("=", 70))

pteori_data <- build_umf_and_scale("pteori")

message("  Fitting colext model...")
pteori_mod_original <- colext(
  psiformula     = pteori_data$cfg$psi_formula,
  gammaformula   = pteori_data$cfg$gamma_formula,
  epsilonformula = pteori_data$cfg$epsilon_formula,  # now ~ LC12 + NDVI + pr
  pformula       = pteori_data$cfg$p_formula,
  data = pteori_data$umf
)

message("  AIC (original, with NDVI): ", round(pteori_mod_original@AIC, 2))
cat("\n  Coefficients:\n")
print(summary(pteori_mod_original))

# Save as the canonical model
saveRDS(pteori_mod_original, here("results", "pteori_model_object.rds"))

# Rename previous no-NDVI model as sensitivity
file.copy(
  here("results", "pteori_epsilon_revised.rds"),
  here("results", "pteori_epsilon_sensitivity_noNDVI.rds"),
  overwrite = TRUE
)
message("  Saved no-NDVI model as sensitivity: pteori_epsilon_sensitivity_noNDVI.rds")

# Save summary
sink(here("results", "pteori_model_summary.txt"))
cat("Species: pteori (REVERTED to original epsilon with NDVI)\n")
cat("Date:", as.character(Sys.time()), "\n")
cat("AIC:", round(pteori_mod_original@AIC, 2), "\n\n")
summary(pteori_mod_original)
sink()


###############################################################################
# PART 2: Sensitivity table (with vs without NDVI for pteori/epsilon)
###############################################################################
message("\n", strrep("=", 70))
message("  Creating sensitivity table: pteori epsilon with vs without NDVI")
message(strrep("=", 70))

pteori_noNDVI <- readRDS(here("results", "pteori_epsilon_sensitivity_noNDVI.rds"))

# Extract epsilon coefficients from both
coefs_with <- coef(pteori_mod_original)
coefs_without <- coef(pteori_noNDVI)

ext_idx_with    <- grep("^ext\\(", names(coefs_with))
ext_idx_without <- grep("^ext\\(", names(coefs_without))

se_with <- tryCatch(sqrt(diag(vcov(pteori_mod_original)))[ext_idx_with],
                     warning = function(w) rep(NA, length(ext_idx_with)))
se_without <- tryCatch(sqrt(diag(vcov(pteori_noNDVI)))[ext_idx_without],
                        warning = function(w) rep(NA, length(ext_idx_without)))

sens_with <- data.frame(
  model = "With NDVI (primary)",
  AIC = round(pteori_mod_original@AIC, 2),
  parameter = names(coefs_with)[ext_idx_with],
  estimate = round(coefs_with[ext_idx_with], 4),
  SE = round(se_with, 4),
  stringsAsFactors = FALSE
)

sens_without <- data.frame(
  model = "Without NDVI (sensitivity)",
  AIC = round(pteori_noNDVI@AIC, 2),
  parameter = names(coefs_without)[ext_idx_without],
  estimate = round(coefs_without[ext_idx_without], 4),
  SE = round(se_without, 4),
  stringsAsFactors = FALSE
)

sensitivity_table <- rbind(sens_with, sens_without)
write.csv(sensitivity_table,
          here("results", "pteori_sensitivity_ndvi_epsilon.csv"),
          row.names = FALSE)
message("  Saved: results/pteori_sensitivity_ndvi_epsilon.csv")
print(sensitivity_table)


###############################################################################
# PART 3: Recompute equilibrium with reverted pteori + all species
###############################################################################
message("\n", strrep("=", 70))
message("  Recomputing equilibrium occupancy (all species)")
message(strrep("=", 70))

species_all <- c("otitar", "ptealc", "pteori", "tettet")
species_names <- c("Otis tarda", "Pterocles alchata",
                   "Pterocles orientalis", "Tetrax tetrax")

compute_equilibrium_ci <- function(model, sp_code, n_boot = 5000) {
  coefs <- coef(model)
  V     <- vcov(model)
  vcov_ok <- all(is.finite(V))

  col_idx <- grep("col\\(Int\\)", names(coefs))
  ext_idx <- grep("ext\\(Int\\)", names(coefs))

  if (!vcov_ok || inherits(try(chol(V), silent = TRUE), "try-error")) {
    message("  vcov not PD for ", sp_code, " — using nearPD")
    V <- as.matrix(Matrix::nearPD(V, corr = FALSE)$mat)
    vcov_ok <- FALSE
  }

  set.seed(42)
  sims <- MASS::mvrnorm(n_boot, mu = coefs, Sigma = V)

  gamma_sims   <- plogis(sims[, col_idx])
  epsilon_sims <- plogis(sims[, ext_idx])
  psi_star     <- gamma_sims / (gamma_sims + epsilon_sims)
  recol_yr     <- 1 / gamma_sims
  persist_yr   <- 1 / epsilon_sims

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
    recol_yr_med      = round(median(recol_yr)),
    recol_yr_lo       = round(quantile(recol_yr, 0.025)),
    recol_yr_hi       = round(quantile(recol_yr, 0.975)),
    persist_yr_med    = round(median(persist_yr)),
    persist_yr_lo     = round(quantile(persist_yr, 0.025)),
    persist_yr_hi     = round(quantile(persist_yr, 0.975)),
    stringsAsFactors = FALSE
  )
}

eq_results <- list()
for (i in seq_along(species_all)) {
  sp <- species_all[i]
  message("\n  ", sp, " (", species_names[i], "):")
  mod <- readRDS(here("results", paste0(sp, "_model_object.rds")))
  message("  AIC = ", round(mod@AIC, 2))
  res <- compute_equilibrium_ci(mod, sp)
  eq_results[[sp]] <- res
  message("  psi* = ", res$psi_star_med, "% [",
          res$psi_star_lo, ", ", res$psi_star_hi, "]")
}

eq_table <- do.call(rbind, eq_results)
write.csv(eq_table, here("results", "equilibrium_occupancy_table.csv"),
          row.names = FALSE)
message("\n  Saved: results/equilibrium_occupancy_table.csv")


###############################################################################
# PART 4: Compute extinction debt
###############################################################################
message("\n", strrep("=", 70))
message("  DECISION 5: Computing extinction debt")
message(strrep("=", 70))

# Current occupancy = last-year simulated prevalence from scripts/4
debt_rows <- list()
for (sp in species_all) {
  sim <- read.csv(here("results", paste0(sp, "_simulation_prevalence.csv")))
  current_occ <- sim$mean_occ[nrow(sim)] * 100  # last year, convert to %

  eq <- eq_table[eq_table$species == sp, ]

  debt_rows[[sp]] <- data.frame(
    species = sp,
    current_occupancy_pct = round(current_occ, 4),
    psi_star_pct = eq$psi_star_med,
    psi_star_lo = eq$psi_star_lo,
    psi_star_hi = eq$psi_star_hi,
    extinction_debt_pct = round(current_occ - eq$psi_star_med, 4),
    debt_fraction = round((current_occ - eq$psi_star_med) / current_occ * 100, 1),
    recol_yr_med = eq$recol_yr_med,
    vcov_ok = eq$vcov_ok,
    stringsAsFactors = FALSE
  )

  message("  ", sp, ": current=", round(current_occ, 3), "%, psi*=",
          eq$psi_star_med, "%, debt=",
          round(current_occ - eq$psi_star_med, 3), "% (",
          round((current_occ - eq$psi_star_med)/current_occ*100, 1), "% of current)")
}

debt_table <- do.call(rbind, debt_rows)
write.csv(debt_table, here("results", "extinction_debt_table.csv"),
          row.names = FALSE)
message("\n  Saved: results/extinction_debt_table.csv")
print(debt_table)


###############################################################################
# PART 5: Prepare isocline plot data
###############################################################################
message("\n", strrep("=", 70))
message("  Preparing isocline plot data")
message(strrep("=", 70))

# Compute corrected gamma and epsilon means for each species
# from the REVISED models (ptealc without NDVI in gamma, pteori with NDVI in epsilon)
iso_rows <- list()
for (sp in species_all) {
  mod <- readRDS(here("results", paste0(sp, "_model_object.rds")))

  # Mean corrected gamma and epsilon across sites × years
  gamma_pred <- predict(mod, type = "col")$Predicted
  eps_pred   <- predict(mod, type = "ext")$Predicted

  # Baseline (intercept-only) rates
  coefs <- coef(mod)
  col_int <- grep("col\\(Int\\)", names(coefs))
  ext_int <- grep("ext\\(Int\\)", names(coefs))
  gamma_baseline <- plogis(coefs[col_int])
  eps_baseline   <- plogis(coefs[ext_int])

  # Bootstrap CIs for baseline rates
  V <- vcov(mod)
  vcov_ok <- all(is.finite(V))
  if (!vcov_ok || inherits(try(chol(V), silent = TRUE), "try-error")) {
    V <- as.matrix(Matrix::nearPD(V, corr = FALSE)$mat)
  }
  set.seed(42)
  sims <- MASS::mvrnorm(5000, mu = coefs, Sigma = V)
  g_sims <- plogis(sims[, col_int])
  e_sims <- plogis(sims[, ext_int])

  iso_rows[[sp]] <- data.frame(
    species = sp,
    gamma_mean = mean(gamma_pred, na.rm = TRUE),
    epsilon_mean = mean(eps_pred, na.rm = TRUE),
    gamma_baseline = gamma_baseline,
    epsilon_baseline = eps_baseline,
    gamma_lo = quantile(g_sims, 0.025),
    gamma_hi = quantile(g_sims, 0.975),
    epsilon_lo = quantile(e_sims, 0.025),
    epsilon_hi = quantile(e_sims, 0.975),
    psi_star = gamma_baseline / (gamma_baseline + eps_baseline),
    stringsAsFactors = FALSE
  )

  message("  ", sp, ": gamma_baseline=", round(gamma_baseline, 6),
          " epsilon_baseline=", round(eps_baseline, 4),
          " psi*=", round(gamma_baseline/(gamma_baseline+eps_baseline)*100, 3), "%")
}

iso_data <- do.call(rbind, iso_rows)
write.csv(iso_data, here("results", "isocline_plot_data.csv"), row.names = FALSE)
message("\n  Saved: results/isocline_plot_data.csv")

message("\n=== All decisions executed successfully ===")
