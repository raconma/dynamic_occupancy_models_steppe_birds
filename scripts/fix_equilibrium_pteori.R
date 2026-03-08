###############################################################################
# fix_equilibrium_pteori.R
#
# Compute equilibrium for pteori with non-PD vcov fix, then produce
# the final combined equilibrium table.
###############################################################################

library(here)
library(unmarked)
library(MASS)
library(Matrix)

compute_equilibrium_ci <- function(model, sp_code, n_boot = 5000) {
  coefs <- coef(model)
  V     <- vcov(model)
  vcov_ok <- all(is.finite(V))

  col_idx <- grep("col\\(Int\\)", names(coefs))
  ext_idx <- grep("ext\\(Int\\)", names(coefs))

  cat("  ", sp_code, ": col intercept =", coefs[col_idx],
      " ext intercept =", coefs[ext_idx], "\n")

  if (!vcov_ok || inherits(try(chol(V), silent = TRUE), "try-error")) {
    # Use nearPD to get approximate positive definite matrix
    message("  vcov not PD for ", sp_code, " — using nearPD approximation")
    V_pd <- as.matrix(Matrix::nearPD(V, corr = FALSE)$mat)
    vcov_ok <- FALSE
  } else {
    V_pd <- V
  }

  set.seed(42)
  sims <- MASS::mvrnorm(n_boot, mu = coefs, Sigma = V_pd)

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

species_all <- c("otitar", "ptealc", "pteori", "tettet")

results_list <- list()
for (sp in species_all) {
  message("\n  Processing ", sp, "...")
  mod <- readRDS(here("results", paste0(sp, "_model_object.rds")))
  message("  AIC = ", round(mod@AIC, 2))
  res <- compute_equilibrium_ci(mod, sp)
  results_list[[sp]] <- res
  message("  psi* = ", res$psi_star_med, "% [",
          res$psi_star_lo, ", ", res$psi_star_hi, "]")
}

eq_table <- do.call(rbind, results_list)
write.csv(eq_table,
          here("results", "equilibrium_occupancy_table.csv"),
          row.names = FALSE)

cat("\n==========================================\n")
cat("REVISED EQUILIBRIUM TABLE:\n")
cat("==========================================\n")
print(eq_table[, c("species", "vcov_ok",
                    "gamma_baseline_pct", "epsilon_baseline_pct",
                    "psi_star_med", "psi_star_lo", "psi_star_hi")])

cat("\n\nFOR ABSTRACT (excluding ptealc if CI > 50%):\n")
for (i in 1:nrow(eq_table)) {
  sp <- eq_table$species[i]
  ci_width <- eq_table$psi_star_hi[i] - eq_table$psi_star_lo[i]
  flag <- if (ci_width > 50) " [EXCLUDE - uninformative CI]" else ""
  cat(sprintf("  %s: psi* = %.2f%% [%.2f%%, %.2f%%]%s\n",
              sp, eq_table$psi_star_med[i],
              eq_table$psi_star_lo[i], eq_table$psi_star_hi[i], flag))
}

message("\nDone.")
