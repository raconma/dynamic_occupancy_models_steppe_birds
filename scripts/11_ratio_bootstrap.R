###############################################################################
# 11_ratio_bootstrap.R
#
# Compute the extinction-to-colonisation ratio (epsilon/gamma) with bootstrap
# confidence intervals for 4 steppe bird species.
#
# At mean covariate values (= 0 on the z-scored scale), gamma and epsilon
# reduce to plogis(intercept) because all other covariate terms vanish.
# We draw 5000 coefficient vectors from the multivariate normal defined by
# the fitted coefficients and their vcov matrix, then compute the ratio for
# each draw.
#
# Inputs:
#   results/{sp}_model_object.rds   (unmarkedFitColExt objects)
#
# Outputs:
#   results/ratio_bootstrap.csv     (summary table)
#   results/bootstrap_draws_5000.rds (raw draws for reuse by scripts 12/13)
###############################################################################

library(here)
library(unmarked)
library(MASS)

set.seed(42)

###############################################################################
# SPECIES DEFINITIONS
###############################################################################

species <- c("otitar", "ptealc", "pteori", "tettet")
species_names <- c("Otis tarda", "Pterocles alchata",
                   "Pterocles orientalis", "Tetrax tetrax")

N_BOOT <- 5000

###############################################################################
# MAIN LOOP
###############################################################################

results_list   <- list()
draws_list     <- list()

for (i in seq_along(species)) {

  sp <- species[i]
  cat("\n===", species_names[i], "(", sp, ") ===\n")

  # --- Load model ----------------------------------------------------------
  mod_path <- here("results", paste0(sp, "_model_object.rds"))
  mod      <- readRDS(mod_path)

  coefs <- coef(mod)
  V     <- vcov(mod)

  # --- Identify intercept indices ------------------------------------------
  col_idx <- grep("^col\\(Int\\)$", names(coefs))
  ext_idx <- grep("^ext\\(Int\\)$", names(coefs))

  cat("  col intercept:", names(coefs)[col_idx], "=",
      round(coefs[col_idx], 4), "\n")
  cat("  ext intercept:", names(coefs)[ext_idx], "=",
      round(coefs[ext_idx], 4), "\n")

  # --- Ensure vcov is positive definite ------------------------------------
  V_pd <- tryCatch({
    chol(V)
    V
  }, error = function(e) {
    cat("  vcov not positive definite; applying nearPD correction\n")
    as.matrix(Matrix::nearPD(V, corr = FALSE)$mat)
  })

  # --- Draw coefficient vectors --------------------------------------------
  sims <- MASS::mvrnorm(n = N_BOOT, mu = coefs, Sigma = V_pd)

  # --- Compute gamma and epsilon at mean covariates (= intercept only) -----
  gamma_draws   <- plogis(sims[, col_idx])
  epsilon_draws <- plogis(sims[, ext_idx])
  ratio_draws   <- epsilon_draws / gamma_draws

  # --- Store draws for reuse -----------------------------------------------
  draws_list[[sp]] <- list(
    gamma_draws   = gamma_draws,
    epsilon_draws = epsilon_draws,
    ratio_draws   = ratio_draws
  )

  # --- Summarise -----------------------------------------------------------
  ratio_q   <- quantile(ratio_draws, probs = c(0.025, 0.5, 0.975))
  gamma_q   <- quantile(gamma_draws, probs = c(0.025, 0.5, 0.975))
  epsilon_q <- quantile(epsilon_draws, probs = c(0.025, 0.5, 0.975))

  pct_gt100   <- 100 * mean(ratio_draws > 100)
  pct_gt1000  <- 100 * mean(ratio_draws > 1000)
  pct_gt10000 <- 100 * mean(ratio_draws > 10000)

  cat("  gamma   median:", round(gamma_q[2], 4),
      " [", round(gamma_q[1], 4), ",", round(gamma_q[3], 4), "]\n")
  cat("  epsilon median:", round(epsilon_q[2], 4),
      " [", round(epsilon_q[1], 4), ",", round(epsilon_q[3], 4), "]\n")
  cat("  ratio   median:", round(ratio_q[2], 1),
      " [", round(ratio_q[1], 1), ",", round(ratio_q[3], 1), "]\n")
  cat("  P(ratio > 100):", round(pct_gt100, 1), "%\n")
  cat("  P(ratio > 1000):", round(pct_gt1000, 1), "%\n")
  cat("  P(ratio > 10000):", round(pct_gt10000, 1), "%\n")

  results_list[[i]] <- data.frame(
    species        = sp,
    ratio_median   = round(ratio_q[2], 2),
    ratio_lo       = round(ratio_q[1], 2),
    ratio_hi       = round(ratio_q[3], 2),
    pct_gt100      = round(pct_gt100, 2),
    pct_gt1000     = round(pct_gt1000, 2),
    pct_gt10000    = round(pct_gt10000, 2),
    gamma_median   = round(gamma_q[2], 6),
    gamma_lo       = round(gamma_q[1], 6),
    gamma_hi       = round(gamma_q[3], 6),
    epsilon_median = round(epsilon_q[2], 6),
    epsilon_lo     = round(epsilon_q[1], 6),
    epsilon_hi     = round(epsilon_q[3], 6),
    stringsAsFactors = FALSE
  )
}

###############################################################################
# COMBINE AND SAVE
###############################################################################

results_df <- do.call(rbind, results_list)
rownames(results_df) <- NULL

# Save summary CSV
csv_path <- here("results", "ratio_bootstrap.csv")
write.csv(results_df, csv_path, row.names = FALSE)
cat("\nSaved:", csv_path, "\n")

# Save raw draws for reuse by scripts 12 and 13
draws_path <- here("results", "bootstrap_draws_5000.rds")
saveRDS(draws_list, draws_path)
cat("Saved:", draws_path, "\n")

###############################################################################
# PRINT FORMATTED TABLE
###############################################################################

cat("\n")
cat("=======================================================================\n")
cat("  epsilon/gamma RATIO BOOTSTRAP RESULTS (n = 5000, seed = 42)\n")
cat("=======================================================================\n")
cat(sprintf("  %-8s %10s %22s %8s %8s %8s\n",
            "Species", "Median", "95% CI", ">100", ">1000", ">10000"))
cat("-----------------------------------------------------------------------\n")
for (j in seq_len(nrow(results_df))) {
  r <- results_df[j, ]
  ci_str <- sprintf("[%.1f, %.1f]", r$ratio_lo, r$ratio_hi)
  cat(sprintf("  %-8s %10.1f %22s %7.1f%% %7.1f%% %7.1f%%\n",
              r$species, r$ratio_median, ci_str,
              r$pct_gt100, r$pct_gt1000, r$pct_gt10000))
}
cat("=======================================================================\n")
cat("\nDone.\n")
