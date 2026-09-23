###############################################################################
# 19_population_averaged_rates.R
#
# Sensitivity of the headline demographic summaries to how gamma and epsilon
# are summarised across sites.
#
# Scripts 11/compute_equilibrium report gamma and epsilon as plogis(intercept),
# i.e. the rate of a hypothetical cell with all dynamic covariates at their
# z-scored mean. Because the logit link is non-linear and some coefficients are
# large (|beta| up to ~9), the average rate across real cells can differ from
# the rate at the average cell by orders of magnitude (Jensen's inequality).
# The naive rates in script 13, in contrast, are averages over observed
# site-years, so the naive-vs-corrected contrast mixes two different summaries.
#
# This script computes, for each species and for three sets of cells
# (all modelled cells, cells with any detection, top decile of psi1):
#   (a) intercept-based rates (reproduces Table 1)
#   (b) unweighted site-year means of gamma and epsilon
#   (c) occupancy-weighted rates: gamma averaged over expected unoccupied
#       site-years, epsilon over expected occupied site-years, and their
#       ratio (rate asymmetry; the quantity comparable to Table 1 eps/gamma)
#   (d) expected extinctions per expected colonisation (flux ratio). This is
#       an event ratio, which depends on current occupancy: flux ratio =
#       rate ratio x occupied / unoccupied site-years. Not comparable to (c).
#       Note that at equilibrium eps/gamma = (1 - psi)/psi, so for a rare
#       species a large rate ratio is expected even without any imbalance;
#       the flux ratio (> 1 = net loss) is the diagnostic of disequilibrium.
#       flux_first / flux_last isolate the first (2017-18) and last
#       (2022-23) transitions, since the first can reflect psi1 relaxing
#       towards the dynamic equilibrium rather than a real trend.
#   (e) population equilibrium occupancy, mean_i gamma_i / (gamma_i + eps_i)
#   (f) model-projected occupancy in the first and last year, and the
#       extinction debt (last-year projected occupancy minus psistar_pop).
#       Supersedes results/extinction_debt_table.csv, whose "current
#       occupancy" came from the script-4 simulation with double-scaled psi1.
# with uncertainty from 1,000 multivariate-normal coefficient draws, plus the
# naive equivalents computed from the detection histories.
#
# Inputs:
#   results/{sp}_model_object.rds   (unmarkedFitColExt objects)
#
# Outputs (new files; nothing existing is overwritten):
#   results/population_averaged_rates/summary_point.csv
#   results/population_averaged_rates/summary_boot.csv
#   results/population_averaged_rates/naive_summary.csv
###############################################################################

library(here)
library(unmarked)
library(MASS)

set.seed(42)

species <- c("otitar", "ptealc", "pteori", "tettet")
N_DRAWS <- 1000
out_dir <- here("results", "population_averaged_rates")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

###############################################################################
# HELPERS
###############################################################################

# Design matrix that keeps NA rows (model.matrix drops them by default)
design <- function(formula, data) {
  mf <- model.frame(formula, data, na.action = na.pass)
  X <- model.matrix(formula, mf)
  X_full <- matrix(NA_real_, nrow(data), ncol(X),
                   dimnames = list(NULL, colnames(X)))
  X_full[as.integer(rownames(X)), ] <- X
  X_full
}

# Projected occupancy psi[i, t] from psi1 and transition matrices (M x T-1)
project_occ <- function(psi1, gam, eps) {
  n_t <- ncol(gam) + 1
  psi <- matrix(NA_real_, length(psi1), n_t)
  psi[, 1] <- psi1
  for (t in 1:(n_t - 1)) {
    psi[, t + 1] <- psi[, t] * (1 - eps[, t]) + (1 - psi[, t]) * gam[, t]
  }
  psi
}

# All summaries for one coefficient vector and one set of cells
summarise_rates <- function(psi1, gam, eps, b_col0, b_ext0, cells) {
  psi <- project_occ(psi1, gam, eps)
  n_t <- ncol(psi)
  g <- gam[cells, , drop = FALSE]
  e <- eps[cells, , drop = FALSE]
  p <- psi[cells, 1:(n_t - 1), drop = FALSE]

  exp_col <- sum((1 - p) * g, na.rm = TRUE)
  exp_ext <- sum(p * e, na.rm = TRUE)
  g_site <- rowMeans(g, na.rm = TRUE)
  e_site <- rowMeans(e, na.rm = TRUE)

  c(
    gamma_int     = plogis(b_col0),
    eps_int       = plogis(b_ext0),
    ratio_int     = plogis(b_ext0) / plogis(b_col0),
    psistar_int   = plogis(b_col0) / (plogis(b_col0) + plogis(b_ext0)),
    gamma_mean    = mean(g, na.rm = TRUE),
    gamma_median  = median(g, na.rm = TRUE),
    eps_mean      = mean(e, na.rm = TRUE),
    ratio_mean    = mean(e, na.rm = TRUE) / mean(g, na.rm = TRUE),
    gamma_w       = exp_col / sum(1 - p, na.rm = TRUE),
    eps_w         = exp_ext / sum(p, na.rm = TRUE),
    ratio_w       = (exp_ext / sum(p, na.rm = TRUE)) /
                    (exp_col / sum(1 - p, na.rm = TRUE)),
    flux_ratio    = exp_ext / exp_col,
    flux_first    = sum(p[, 1] * e[, 1], na.rm = TRUE) /
                    sum((1 - p[, 1]) * g[, 1], na.rm = TRUE),
    flux_last     = sum(p[, n_t - 1] * e[, n_t - 1], na.rm = TRUE) /
                    sum((1 - p[, n_t - 1]) * g[, n_t - 1], na.rm = TRUE),
    exp_col_yr    = exp_col / (n_t - 1),
    exp_ext_yr    = exp_ext / (n_t - 1),
    psistar_pop   = mean(g_site / (g_site + e_site), na.rm = TRUE),
    occ_first     = mean(psi[cells, 1], na.rm = TRUE),
    occ_last      = mean(psi[cells, n_t], na.rm = TRUE),
    debt          = mean(psi[cells, n_t], na.rm = TRUE) -
                    mean(g_site / (g_site + e_site), na.rm = TRUE)
  )
}

###############################################################################
# MAIN LOOP
###############################################################################

point_rows <- list()
boot_rows  <- list()
naive_rows <- list()

for (sp in species) {
  message("\n=== ", sp, " ===")
  mod <- readRDS(here("results", paste0(sp, "_model_object.rds")))
  umf <- mod@data
  n_t <- umf@numPrimary

  # Drop sites removed by unmarked during fitting (missing covariates)
  keep <- setdiff(seq_len(numSites(umf)), mod@sitesRemoved)
  M <- length(keep)
  sc  <- siteCovs(umf)[keep, , drop = FALSE]
  row_keep <- as.vector(t(outer((keep - 1) * n_t, 1:n_t, "+")))
  ysc <- yearlySiteCovs(umf)[row_keep, , drop = FALSE]
  ysc <- cbind(ysc, sc[rep(seq_len(M), each = n_t), , drop = FALSE])
  rownames(ysc) <- NULL
  rownames(sc)  <- NULL

  fl <- mod@formlist
  X_psi <- design(fl$psi, sc)
  X_col <- design(fl$col, ysc)
  X_ext <- design(fl$ext, ysc)

  # Transitions t -> t+1 use yearly covariates of year t (years 1..T-1)
  trans_rows <- which(rep(1:n_t, M) < n_t)
  X_col <- X_col[trans_rows, , drop = FALSE]
  X_ext <- X_ext[trans_rows, , drop = FALSE]

  cf <- coef(mod)
  V  <- vcov(mod)
  i_psi <- grep("^psi\\(", names(cf))
  i_col <- grep("^col\\(", names(cf))
  i_ext <- grep("^ext\\(", names(cf))
  stopifnot(length(i_psi) == ncol(X_psi), length(i_col) == ncol(X_col),
            length(i_ext) == ncol(X_ext))

  rates_from <- function(b) {
    list(
      psi1 = as.vector(plogis(X_psi %*% b[i_psi])),
      gam  = matrix(plogis(X_col %*% b[i_col]), M, n_t - 1, byrow = TRUE),
      eps  = matrix(plogis(X_ext %*% b[i_ext]), M, n_t - 1, byrow = TRUE)
    )
  }

  r0 <- rates_from(cf)

  # --- Checks against unmarked's own computations ---------------------------
  pg <- predict(mod, type = "col")$Predicted
  pg <- matrix(pg, M, n_t, byrow = TRUE)[, 1:(n_t - 1)]
  message("  max |gamma - predict()|: ",
          signif(max(abs(pg - r0$gam), na.rm = TRUE), 3))
  proj_own <- project_occ(r0$psi1, r0$gam, r0$eps)
  proj_unm <- t(mod@projected[2, , ])  # dim 2 = occupied
  message("  max |projected - unmarked@projected|: ",
          signif(max(abs(proj_own - proj_unm), na.rm = TRUE), 3))

  # --- Naive states from detection histories --------------------------------
  y <- getY(umf)[keep, , drop = FALSE]
  J <- ncol(y) / n_t
  naive <- sapply(1:n_t, function(t) {
    yt <- y[, ((t - 1) * J + 1):(t * J), drop = FALSE]
    ifelse(rowSums(!is.na(yt)) > 0, as.integer(rowSums(yt, na.rm = TRUE) > 0),
           NA_integer_)
  })
  ever_det <- which(rowSums(naive, na.rm = TRUE) > 0)

  z0 <- naive[, 1:(n_t - 1)]
  z1 <- naive[, 2:n_t]
  ok <- !is.na(z0) & !is.na(z1)
  n_col_opp <- sum(ok & z0 == 0); n_col_ev <- sum(ok & z0 == 0 & z1 == 1)
  n_ext_opp <- sum(ok & z0 == 1); n_ext_ev <- sum(ok & z0 == 1 & z1 == 0)
  g_n <- n_col_ev / n_col_opp
  e_n <- n_ext_ev / n_ext_opp
  naive_occ <- colMeans(naive, na.rm = TRUE)
  naive_rows[[sp]] <- data.frame(
    species = sp, n_col_ev, n_col_opp, n_ext_ev, n_ext_opp,
    gamma_naive = g_n, eps_naive = e_n, ratio_naive = e_n / g_n,
    psistar_naive = g_n / (g_n + e_n),
    naive_occ_first = naive_occ[1], naive_occ_last = naive_occ[n_t],
    n_cells = M, n_ever_detected = length(ever_det)
  )

  # --- Cell sets -------------------------------------------------------------
  cell_sets <- list(
    all       = seq_len(M),
    ever_det  = ever_det,
    top10_psi = which(r0$psi1 >= quantile(r0$psi1, 0.9, na.rm = TRUE))
  )

  # --- Point estimates and coefficient draws --------------------------------
  V_pd <- tryCatch({ chol(V); V },
                   error = function(e) as.matrix(Matrix::nearPD(V)$mat))
  draws <- MASS::mvrnorm(N_DRAWS, mu = cf, Sigma = V_pd)
  r_draws <- lapply(seq_len(N_DRAWS), function(k) rates_from(draws[k, ]))

  for (set_name in names(cell_sets)) {
    cells <- cell_sets[[set_name]]
    pt <- summarise_rates(r0$psi1, r0$gam, r0$eps,
                          cf["col(Int)"], cf["ext(Int)"], cells)
    point_rows[[paste(sp, set_name)]] <- data.frame(
      species = sp, cells = set_name, n = length(cells), t(pt))

    bt <- t(sapply(seq_len(N_DRAWS), function(k) {
      summarise_rates(r_draws[[k]]$psi1, r_draws[[k]]$gam, r_draws[[k]]$eps,
                      draws[k, "col(Int)"], draws[k, "ext(Int)"], cells)
    }))
    qs <- apply(bt, 2, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
    boot_rows[[paste(sp, set_name)]] <- data.frame(
      species = sp, cells = set_name,
      metric = colnames(qs), lo = qs[1, ], med = qs[2, ], hi = qs[3, ],
      p_flux_gt1 = mean(bt[, "flux_ratio"] > 1),
      p_flux_last_gt1 = mean(bt[, "flux_last"] > 1),
      p_ratio_w_gt10 = mean(bt[, "ratio_w"] > 10),
      p_debt_gt0 = mean(bt[, "debt"] > 0),
      row.names = NULL)
  }
}

point_df <- do.call(rbind, point_rows); rownames(point_df) <- NULL
boot_df  <- do.call(rbind, boot_rows);  rownames(boot_df)  <- NULL
naive_df <- do.call(rbind, naive_rows); rownames(naive_df) <- NULL

write.csv(point_df, file.path(out_dir, "summary_point.csv"), row.names = FALSE)
write.csv(boot_df,  file.path(out_dir, "summary_boot.csv"),  row.names = FALSE)
write.csv(naive_df, file.path(out_dir, "naive_summary.csv"), row.names = FALSE)

options(width = 200)
print(point_df, digits = 3)
print(naive_df, digits = 3)

sessionInfo()
