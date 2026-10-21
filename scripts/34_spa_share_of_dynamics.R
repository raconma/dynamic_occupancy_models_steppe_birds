###############################################################################
# 34_spa_share_of_dynamics.R
#
# Exploratory check: do SPAs protect a static snapshot of the range rather
# than its dynamics? Among cells with >= 1 detection, and on anchored recent
# transitions (2018->19 ... 2022->23, surveyed in both years), we compute the
# share falling in SPA cells (>= 50 % SPA cover) of
#   - expected occupied cell-years            (where the species is)
#   - expected unoccupied cell-years          (where it could colonise)
#   - expected extinctions and colonisations  (where the dynamics happen)
# and the same shares for the naive (raw) events. Under static protection the
# SPA share of colonisations should be lower than the SPA share of occupancy.
#
# Uncertainty: 1,000 coefficient draws combined with resampling of recorded
# cells. All four species are identifiable with the pooled histories
# (script 39).
#
# Inputs:  results/{sp}_model_object.rds, colext inputs (colext_input_path()),
#          results/core_spa/cell_habitat.csv (script 23)
# Output:  results/core_spa/spa_share_of_dynamics.csv
###############################################################################

library(here)
library(unmarked)
library(MASS)

source(here("R", "colext_utils.R"))
set.seed(42)

species <- c("otitar", "tettet", "pteori", "ptealc")
N_DRAWS <- 1000
hab <- read.csv(here("results", "core_spa", "cell_habitat.csv"))
q <- function(x, p) as.numeric(quantile(x, p, na.rm = TRUE))

rows <- list()
for (sp in species) {
  message("=== ", sp)
  mod <- readRDS(colext_model_path(sp))
  D <- colext_design(mod)
  n_t <- D$n_t
  occ <- read.csv(colext_input_path(sp))[D$keep, ]
  h <- hab[match(occ$cells, hab$cells), ]
  stopifnot(!anyNA(h$cells))
  naive <- sapply(1:n_t, function(t) {
    yt <- D$y[, ((t - 1) * D$J + 1):(t * D$J), drop = FALSE]
    ifelse(rowSums(!is.na(yt)) > 0, as.integer(rowSums(yt, na.rm = TRUE) > 0), NA_integer_)
  })
  surv <- !is.na(naive)
  both <- surv[, 1:(n_t - 1)] & surv[, 2:n_t]
  det <- which(rowSums(naive, na.rm = TRUE) > 0)
  recent <- 2:(n_t - 1)
  in_spa <- h$spa_frac[det] >= 0.5

  # Naive events on anchored recent transitions, per recorded cell
  n0 <- naive[det, recent, drop = FALSE]; n1 <- naive[det, recent + 1, drop = FALSE]
  a <- both[det, recent, drop = FALSE]
  obs <- cbind(occ = rowSums(n0 == 1 & a, na.rm = TRUE),
               unocc = rowSums(n0 == 0 & a, na.rm = TRUE),
               ext = rowSums(n0 == 1 & n1 == 0 & a, na.rm = TRUE),
               col = rowSums(n0 == 0 & n1 == 1 & a, na.rm = TRUE))

  draws <- draw_coefs(mod, N_DRAWS)
  point_only <- attr(draws, "point_only")
  share <- matrix(NA_real_, N_DRAWS, 4, dimnames = list(NULL, c("occ", "unocc", "ext", "col")))
  totals <- share
  for (k in seq_len(N_DRAWS)) {
    fb <- forward_backward(colext_rates(draws[k, ], D), D)
    z <- fb$z1[det, 1:(n_t - 1), drop = FALSE]
    bb <- both[det, , drop = FALSE]
    ev <- cbind(occ = rowSums((z * bb)[, recent, drop = FALSE]),
                unocc = rowSums(((1 - z) * bb)[, recent, drop = FALSE]),
                ext = rowSums((fb$ext[det, , drop = FALSE] * bb)[, recent, drop = FALSE]),
                col = rowSums((fb$col[det, , drop = FALSE] * bb)[, recent, drop = FALSE]))
    b <- sample(seq_along(det), replace = TRUE)
    tot <- colSums(ev[b, , drop = FALSE])
    share[k, ] <- colSums(ev[b[in_spa[b]], , drop = FALSE]) / tot
    totals[k, ] <- colSums(ev)
  }
  for (m in colnames(share)) {
    s <- share[, m]
    rows[[paste(sp, m)]] <- data.frame(
      species = sp, point_only = point_only, quantity = m,
      n_recorded_cells = length(det), n_spa_cells = sum(in_spa),
      expected_total_med = median(totals[, m]),
      spa_share_med = median(s, na.rm = TRUE), spa_share_lo = q(s, .025),
      spa_share_hi = q(s, .975),
      # Difference from the SPA share of occupied cell-years, same draws
      diff_vs_occ_med = median(s - share[, "occ"], na.rm = TRUE),
      diff_vs_occ_lo = q(s - share[, "occ"], .025), diff_vs_occ_hi = q(s - share[, "occ"], .975),
      naive_total = sum(obs[, m]), naive_spa = sum(obs[in_spa, m]),
      naive_spa_share = sum(obs[in_spa, m]) / sum(obs[, m]))
  }
}
out <- do.call(rbind, rows); rownames(out) <- NULL
write.csv(out, here("results", "core_spa", "spa_share_of_dynamics.csv"), row.names = FALSE)
print(out, digits = 3)

sessionInfo()
