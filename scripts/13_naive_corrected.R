###############################################################################
# 13_naive_corrected.R
#
# Compare naive vs detection-corrected transition rates (gamma, epsilon)
# with bootstrap confidence intervals.
#
# Naive estimates treat sites as "occupied" if detected at least once in a
# primary period, ignoring imperfect detection. Corrected estimates come
# from the dynamic occupancy model bootstrap draws (at mean covariate values).
#
# Inputs:
#   data/processed_2023/{sp}/{sp}_occ_wide_dynamic.csv
#   results/bootstrap_draws_5000.rds   (from script 11)
#
# Outputs:
#   results/naive_vs_corrected_full.csv
###############################################################################

library(here)
library(dplyr)
library(tidyr)

set.seed(42)

###############################################################################
# GLOBAL DEFINITIONS
###############################################################################

SPECIES <- c("otitar", "ptealc", "pteori", "tettet")
YEARS   <- 2017:2023
N_YEARS <- length(YEARS)

###############################################################################
# FUNCTION: compute naive gamma and epsilon from raw detection histories
###############################################################################

compute_naive_rates <- function(sp) {

  # -- Load data and apply same NDVI filter as model fitting -----------------
  occ_wide <- read.csv(
    here("data", "processed_2023", sp, paste0(sp, "_occ_wide_dynamic.csv"))
  ) %>%
    drop_na(paste0("NDVI_", YEARS))

  # -- Build naive occupancy state matrix ------------------------------------
  # A site is "occupied" in year t if detected at least once (any y.j.t == 1)
  # A site is "surveyed" in year t if at least one detection column is non-NA

  n_sites <- nrow(occ_wide)
  occ_state <- matrix(NA, nrow = n_sites, ncol = N_YEARS)
  surveyed  <- matrix(FALSE, nrow = n_sites, ncol = N_YEARS)

  for (t in seq_along(YEARS)) {
    yr <- YEARS[t]
    det_cols <- paste0("y.", 1:10, ".", yr)
    det_mat  <- occ_wide[, det_cols]

    # Site is surveyed if at least one secondary occasion is non-NA
    surveyed[, t] <- rowSums(!is.na(det_mat)) > 0

    # Site is naively "occupied" if detected at least once
    occ_state[, t] <- ifelse(
      surveyed[, t],
      as.integer(rowSums(det_mat, na.rm = TRUE) > 0),
      NA_integer_
    )
  }

  # -- Count colonisation and extinction events ------------------------------
  n_col_events        <- 0
  n_col_opportunities <- 0
  n_ext_events        <- 0
  n_ext_opportunities <- 0

  for (t in 1:(N_YEARS - 1)) {
    for (i in 1:n_sites) {
      # Need both years to have valid state
      if (is.na(occ_state[i, t]) || is.na(occ_state[i, t + 1])) next

      if (occ_state[i, t] == 0) {
        # Unoccupied in t -> colonisation opportunity
        n_col_opportunities <- n_col_opportunities + 1
        if (occ_state[i, t + 1] == 1) {
          n_col_events <- n_col_events + 1
        }
      } else {
        # Occupied in t -> extinction opportunity
        n_ext_opportunities <- n_ext_opportunities + 1
        if (occ_state[i, t + 1] == 0) {
          n_ext_events <- n_ext_events + 1
        }
      }
    }
  }

  naive_gamma   <- n_col_events / n_col_opportunities
  naive_epsilon <- n_ext_events / n_ext_opportunities

  list(
    naive_gamma          = naive_gamma,
    naive_epsilon        = naive_epsilon,
    n_col_events         = n_col_events,
    n_col_opportunities  = n_col_opportunities,
    n_ext_events         = n_ext_events,
    n_ext_opportunities  = n_ext_opportunities
  )
}

###############################################################################
# MAIN: loop over species and build comparison table
###############################################################################

# Load bootstrap draws (created by script 11)
boot_draws <- readRDS(here("results", "bootstrap_draws_5000.rds"))

results_list <- list()

for (sp in SPECIES) {

  cat("Processing:", sp, "\n")

  # -- Naive rates -----------------------------------------------------------
  naive <- compute_naive_rates(sp)

  # -- Corrected rates from bootstrap draws ----------------------------------
  gamma_draws   <- boot_draws[[sp]]$gamma_draws
  epsilon_draws <- boot_draws[[sp]]$epsilon_draws

  corrected_gamma_median <- median(gamma_draws)
  corrected_gamma_lo     <- quantile(gamma_draws, 0.025)
  corrected_gamma_hi     <- quantile(gamma_draws, 0.975)

  corrected_epsilon_median <- median(epsilon_draws)
  corrected_epsilon_lo     <- quantile(epsilon_draws, 0.025)
  corrected_epsilon_hi     <- quantile(epsilon_draws, 0.975)

  # -- Ratios ----------------------------------------------------------------
  ratio_gamma   <- naive$naive_gamma   / corrected_gamma_median
  ratio_epsilon <- naive$naive_epsilon / corrected_epsilon_median

  results_list[[sp]] <- tibble(
    species                  = sp,
    naive_gamma              = naive$naive_gamma,
    corrected_gamma_median   = corrected_gamma_median,
    corrected_gamma_lo       = corrected_gamma_lo,
    corrected_gamma_hi       = corrected_gamma_hi,
    ratio_gamma              = ratio_gamma,
    naive_epsilon            = naive$naive_epsilon,
    corrected_epsilon_median = corrected_epsilon_median,
    corrected_epsilon_lo     = corrected_epsilon_lo,
    corrected_epsilon_hi     = corrected_epsilon_hi,
    ratio_epsilon            = ratio_epsilon,
    n_col_events             = naive$n_col_events,
    n_col_opportunities      = naive$n_col_opportunities,
    n_ext_events             = naive$n_ext_events,
    n_ext_opportunities      = naive$n_ext_opportunities
  )
}

###############################################################################
# COMBINE AND SAVE
###############################################################################

results_df <- bind_rows(results_list)

write.csv(
  results_df,
  here("results", "naive_vs_corrected_full.csv"),
  row.names = FALSE
)

###############################################################################
# PRINT FORMATTED TABLE
###############################################################################

cat("\n")
cat("=============================================================\n")
cat("  Naive vs Detection-Corrected Transition Rates\n")
cat("=============================================================\n\n")

for (i in seq_len(nrow(results_df))) {
  r <- results_df[i, ]
  cat(sprintf("--- %s ---\n", r$species))
  cat(sprintf("  Colonisation (gamma):\n"))
  cat(sprintf("    Naive:     %.4f  (%d events / %d opportunities)\n",
              r$naive_gamma, r$n_col_events, r$n_col_opportunities))
  cat(sprintf("    Corrected: %.4f  [%.4f, %.4f]\n",
              r$corrected_gamma_median, r$corrected_gamma_lo, r$corrected_gamma_hi))
  cat(sprintf("    Ratio (naive/corrected): %.2f\n", r$ratio_gamma))
  cat(sprintf("  Extinction (epsilon):\n"))
  cat(sprintf("    Naive:     %.4f  (%d events / %d opportunities)\n",
              r$naive_epsilon, r$n_ext_events, r$n_ext_opportunities))
  cat(sprintf("    Corrected: %.4f  [%.4f, %.4f]\n",
              r$corrected_epsilon_median, r$corrected_epsilon_lo, r$corrected_epsilon_hi))
  cat(sprintf("    Ratio (naive/corrected): %.2f\n\n", r$ratio_epsilon))
}

cat("Output saved to: results/naive_vs_corrected_full.csv\n")
