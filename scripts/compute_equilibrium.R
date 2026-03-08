###############################################################################
# compute_equilibrium.R
#
# Compute equilibrium occupancy psi* = gamma/(gamma+epsilon) with bootstrap CI
# for all four steppe bird species.
#
# This uses the intercept-only baseline (at mean covariate values = 0 on
# standardised scale), which gives the spatially-averaged equilibrium.
###############################################################################

library(unmarked)
library(MASS)

species <- c("otitar", "ptealc", "pteori", "tettet")
species_names <- c("Otis tarda", "Pterocles alchata",
                   "Pterocles orientalis", "Tetrax tetrax")

compute_equilibrium_ci <- function(model, sp_code, n_boot = 5000) {
  coefs <- coef(model)
  V     <- vcov(model)

  # Check for NAs or Inf in vcov
  vcov_ok <- all(is.finite(V))
  if (!vcov_ok) {
    warning(paste(sp_code, ": vcov contains non-finite values"))
  }

  set.seed(42)
  sims <- MASS::mvrnorm(n_boot, mu = coefs, Sigma = V)

  # Extract gamma and epsilon intercept indices
  col_idx <- grep("col\\(Int\\)", names(coefs))
  ext_idx <- grep("ext\\(Int\\)", names(coefs))

  cat("  col intercept name:", names(coefs)[col_idx], "=", coefs[col_idx], "\n")
  cat("  ext intercept name:", names(coefs)[ext_idx], "=", coefs[ext_idx], "\n")

  # Baseline gamma and epsilon at mean covariates
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
for (i in seq_along(species)) {
  sp <- species[i]
  cat("\n==========================================\n")
  cat("SPECIES:", sp, "(", species_names[i], ")\n")

  mod_path <- paste0("results/", sp, "_model_object.rds")
  if (!file.exists(mod_path)) {
    cat("  MODEL NOT FOUND:", mod_path, "\n")
    next
  }

  mod <- readRDS(mod_path)
  cat("  AIC:", mod@AIC, "\n")

  res <- compute_equilibrium_ci(mod, sp)
  results_list[[sp]] <- res

  cat("  gamma baseline:", res$gamma_baseline_pct, "%\n")
  cat("  epsilon baseline:", res$epsilon_baseline_pct, "%\n")
  cat("  psi* median:", res$psi_star_med, "% [",
      res$psi_star_lo, ",", res$psi_star_hi, "]\n")
  cat("  Recolonisation time:", res$recol_yr_med, "yr [",
      res$recol_yr_lo, ",", res$recol_yr_hi, "]\n")
  cat("  Persistence time:", res$persist_yr_med, "yr [",
      res$persist_yr_lo, ",", res$persist_yr_hi, "]\n")
}

eq_table <- do.call(rbind, results_list)
write.csv(eq_table, "results/equilibrium_occupancy_table.csv", row.names = FALSE)
cat("\n\nTable saved to results/equilibrium_occupancy_table.csv\n")
cat("\n==========================================\n")
cat("SUMMARY FOR ABSTRACT:\n")
cat("==========================================\n")
print(eq_table[, c("species", "psi_star_med", "psi_star_lo", "psi_star_hi",
                    "recol_yr_med", "recol_yr_lo", "recol_yr_hi",
                    "persist_yr_med", "persist_yr_lo", "persist_yr_hi")])
