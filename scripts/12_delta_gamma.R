###############################################################################
# 12_delta_gamma.R
#
# What increase in colonisation would achieve a viable equilibrium occupancy?
#
# For each species, we compute the multiplicative increase in baseline gamma
# needed to reach target equilibrium occupancies (psi* = 0.05, 0.10).
#
# At equilibrium: psi* = gamma / (gamma + epsilon)
# Rearranging:    gamma_required = psi*_target * epsilon / (1 - psi*_target)
# Multiplier:     gamma_required / gamma_current
#
# Inputs:
#   results/bootstrap_draws_5000.rds   (list by species code, each with
#                                       gamma_draws, epsilon_draws, ratio_draws)
#
# Outputs:
#   results/delta_gamma.csv
###############################################################################

library(here)

set.seed(42)

###############################################################################
# 1. Load bootstrap draws
###############################################################################

draws_path <- here("results", "bootstrap_draws_5000.rds")
stopifnot(file.exists(draws_path))
boot_draws <- readRDS(draws_path)

species_codes <- c("otitar", "ptealc", "pteori", "tettet")
species_names <- c(
  otitar = "Otis tarda",
  ptealc = "Pterocles alchata",
  pteori = "Pterocles orientalis",
  tettet = "Tetrax tetrax"
)

psi_targets <- c(0.05, 0.10)

###############################################################################
# 2. Compute required gamma and multiplier for each species x target
###############################################################################

results_list <- list()

for (sp in species_codes) {

  sp_draws <- boot_draws[[sp]]
  gamma_current <- sp_draws$gamma_draws
  epsilon       <- sp_draws$epsilon_draws

  for (psi_t in psi_targets) {

    # Required colonisation rate to achieve target equilibrium
    gamma_required <- psi_t * epsilon / (1 - psi_t)

    # Multiplicative increase needed
    multiplier <- gamma_required / gamma_current

    results_list[[paste0(sp, "_", psi_t)]] <- data.frame(
      species              = sp,
      psi_target           = psi_t,
      gamma_required_median = median(gamma_required),
      gamma_required_lo    = quantile(gamma_required, 0.025, names = FALSE),
      gamma_required_hi    = quantile(gamma_required, 0.975, names = FALSE),
      multiplier_median    = median(multiplier),
      multiplier_lo        = quantile(multiplier, 0.025, names = FALSE),
      multiplier_hi        = quantile(multiplier, 0.975, names = FALSE),
      stringsAsFactors     = FALSE
    )
  }
}

delta_gamma <- do.call(rbind, results_list)
rownames(delta_gamma) <- NULL

###############################################################################
# 3. Save output
###############################################################################

out_path <- here("results", "delta_gamma.csv")
write.csv(delta_gamma, out_path, row.names = FALSE)
cat("Saved:", out_path, "\n\n")

###############################################################################
# 4. Print formatted table to console
###############################################################################

cat("==========================================================================\n")
cat("  Delta-gamma analysis: colonisation multiplier for viable equilibrium\n")
cat("==========================================================================\n\n")

for (sp in species_codes) {
  cat(sprintf("--- %s (%s) ---\n", sp, species_names[sp]))

  sp_rows <- delta_gamma[delta_gamma$species == sp, ]

  for (i in seq_len(nrow(sp_rows))) {
    r <- sp_rows[i, ]
    cat(sprintf(
      "  psi* target = %.0f%%:  gamma_req = %.6f [%.6f, %.6f]  |  multiplier = %.1fx [%.1f, %.1f]\n",
      r$psi_target * 100,
      r$gamma_required_median, r$gamma_required_lo, r$gamma_required_hi,
      r$multiplier_median, r$multiplier_lo, r$multiplier_hi
    ))
  }
  cat("\n")
}

cat("==========================================================================\n")
cat("Done.\n")
