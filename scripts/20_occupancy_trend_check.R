###############################################################################
# 20_occupancy_trend_check.R
#
# Is the occupancy decline projected by the colext models supported by the
# detection data, or is it the model relaxing from psi1 towards its own
# equilibrium?
#
# For each year t, the model implies an expected naive occupancy among the
# cells surveyed that year:
#   E[naive_t] = mean_i( psi_it * (1 - prod_j (1 - p_ijt)) )
# where psi_it is the projected (unconditional) occupancy and p_ijt the
# fitted detection probability of each realised visit. If the projected
# decline were real, observed naive occupancy should fall in step with
# E[naive_t]. Predictive intervals come from simulating latent occupancy and
# detection at the MLE (sampling variation only; parameter uncertainty is not
# propagated). The smoothed occupancy (conditional on the data) is reported
# alongside the projected trajectory.
#
# Inputs:
#   results/{sp}_model_object.rds
#
# Outputs (new files):
#   results/population_averaged_rates/trend_check.csv
#   results/population_averaged_rates/trend_check.png
###############################################################################

library(here)
library(unmarked)
library(ggplot2)

set.seed(42)

species <- c("otitar", "ptealc", "pteori", "tettet")
species_names <- c(otitar = "O. tarda", ptealc = "P. alchata",
                   pteori = "P. orientalis", tettet = "T. tetrax")
YEARS <- 2017:2023
N_SIM <- 2000
out_dir <- here("results", "population_averaged_rates")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

rows <- list()

for (sp in species) {
  message("=== ", sp, " ===")
  mod <- readRDS(here("results", paste0(sp, "_model_object.rds")))
  umf <- mod@data
  n_t <- umf@numPrimary
  keep <- setdiff(seq_len(numSites(umf)), mod@sitesRemoved)
  M <- length(keep)

  y <- getY(umf)[keep, , drop = FALSE]
  J <- ncol(y) / n_t

  # Detection probability per observation, site-major in the column order
  # of y (predict() already drops the removed sites)
  p_vec <- predict(mod, type = "det")$Predicted
  stopifnot(length(p_vec) == M * ncol(y))
  p_mat <- matrix(p_vec, M, ncol(y), byrow = TRUE)
  p_mat[is.na(y)] <- NA

  psi_proj <- t(mod@projected[2, , ])      # M x T, dim 2 = occupied
  psi_smooth <- t(mod@smoothed[2, , ])

  for (t in 1:n_t) {
    cols <- ((t - 1) * J + 1):(t * J)
    yt <- y[, cols, drop = FALSE]
    pt <- p_mat[, cols, drop = FALSE]
    surveyed <- rowSums(!is.na(yt)) > 0
    # P(at least one detection | occupied), over realised visits only
    p_star <- 1 - apply(1 - pt[surveyed, , drop = FALSE], 1, prod,
                        na.rm = TRUE)
    psi_s <- psi_proj[surveyed, t]

    obs_naive <- mean(rowSums(yt[surveyed, , drop = FALSE], na.rm = TRUE) > 0)
    exp_naive <- mean(psi_s * p_star)
    sims <- replicate(N_SIM, {
      z <- rbinom(length(psi_s), 1, psi_s)
      mean(z * rbinom(length(psi_s), 1, p_star))
    })

    rows[[paste(sp, t)]] <- data.frame(
      species = sp, year = YEARS[t], n_surveyed = sum(surveyed),
      obs_naive = obs_naive, exp_naive = exp_naive,
      exp_lo = quantile(sims, 0.025), exp_hi = quantile(sims, 0.975),
      psi_projected = mean(psi_proj[, t]),
      psi_smoothed  = mean(psi_smooth[, t]),
      row.names = NULL
    )
  }
}

res <- do.call(rbind, rows)
write.csv(res, file.path(out_dir, "trend_check.csv"), row.names = FALSE)

# Change from first to last year: observed vs model-expected naive occupancy
chg <- do.call(rbind, lapply(split(res, res$species), function(d) {
  data.frame(species = d$species[1],
             obs_change_pct = 100 * (d$obs_naive[n_t] / d$obs_naive[1] - 1),
             exp_change_pct = 100 * (d$exp_naive[n_t] / d$exp_naive[1] - 1),
             proj_change_pct = 100 * (d$psi_projected[n_t] /
                                      d$psi_projected[1] - 1),
             smooth_change_pct = 100 * (d$psi_smoothed[n_t] /
                                        d$psi_smoothed[1] - 1),
             years_obs_outside_pi = sum(d$obs_naive < d$exp_lo |
                                        d$obs_naive > d$exp_hi))
}))
options(width = 200)
print(res, digits = 3)
print(chg, digits = 3)

p <- ggplot(res, aes(x = year)) +
  geom_ribbon(aes(ymin = exp_lo, ymax = exp_hi), fill = "grey85") +
  geom_line(aes(y = exp_naive, linetype = "Model-expected")) +
  geom_point(aes(y = obs_naive, shape = "Observed")) +
  geom_line(aes(y = obs_naive, linetype = "Observed")) +
  facet_wrap(~ species, scales = "free_y",
             labeller = labeller(species = species_names)) +
  scale_linetype_manual(values = c("Model-expected" = "dashed",
                                   "Observed" = "solid"), name = NULL) +
  scale_shape_manual(values = c("Observed" = 16), guide = "none") +
  labs(x = NULL, y = "Naive occupancy (surveyed cells)") +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "italic"))
ggsave(file.path(out_dir, "trend_check.png"), p, width = 7, height = 5.5,
       dpi = 300)

sessionInfo()
