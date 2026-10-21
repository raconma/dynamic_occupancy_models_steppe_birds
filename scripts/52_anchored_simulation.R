###############################################################################
# 52_anchored_simulation.R
#
# Does the anchored-transition summary (script 26) create losses? Robustness
# plan B (docs/robustness_plan_BiolCons.md). Great bustard detection histories
# are simulated from the fitted primary model under two truths and the
# primary model is refitted to each replicate:
#   equilibrium  every cell at its stationary state: gamma and epsilon fixed
#                at the cell's mean fitted rates over the transitions and
#                initial occupancy gamma / (gamma + epsilon), so expected
#                extinctions equal expected colonisations in every cell
#   fitted       the fitted initial occupancy, colonisation and extinction
# Observation: the observed survey pattern (which visits exist) and the fitted
# detection probabilities of those visits; in addition, after a simulated
# record in year t the whole of year t + 1 is deleted with probability
# 1 - c_rec / c_none, where c_rec and c_none are the observed proportions of
# recorded cells surveyed again after a year with and without a record
# (Table S11), so that survey continuation depends on the outcome as in the
# data.
# For each replicate: anchored ratio of expected extinctions to colonisations
# on 2018-2023 transitions (posterior means at the refitted coefficients, as
# in script 26), the naive ratio, and the true ratio on the same transitions.
#
# Usage:   Rscript scripts/52_anchored_simulation.R   (env N_SIM, N_CORES)
# Inputs:  results/otitar_model_object.rds, results/manuscript_tables/
#          tableS_missingness.csv
# Outputs: results/robustness/anchored_simulation.csv (one row per replicate)
#          results/robustness/anchored_simulation_summary.csv
###############################################################################

library(here)
library(unmarked)
library(parallel)

source(here("R", "colext_utils.R"))
set.seed(42)

SP <- "otitar"
N_SIM <- as.integer(Sys.getenv("N_SIM", "100"))
N_CORES <- as.integer(Sys.getenv("N_CORES", "6"))
RECENT <- 2:6                                   # transitions 2018->19 ... 2022->23
out_dir <- here("results", "robustness"); dir.create(out_dir, showWarnings = FALSE)

mod <- readRDS(here("results", paste0(SP, "_model_object.rds")))
fl <- mod@formlist
D <- colext_design(mod)
M <- D$M; n_t <- D$n_t; J <- D$J
r_fit <- colext_rates(coef(mod), D)
mask <- !is.na(D$y)                              # observed visits

miss <- read.csv(here("results", "manuscript_tables", "tableS_missingness.csv"))
c_rec <- miss$p_surveyed_next[miss$species == SP & miss$state_t == "recorded"]
c_none <- miss$p_surveyed_next[miss$species == SP & miss$state_t == "not recorded"]
p_drop <- max(0, 1 - c_rec / c_none)
message("deletion probability after a record: ", round(p_drop, 3))

truths <- list(fitted = r_fit, equilibrium = within(r_fit, {
  g <- rowMeans(gam[, 1:(n_t - 1)]); e <- rowMeans(eps[, 1:(n_t - 1)])
  psi1 <- g / (g + e)
  gam <- matrix(g, M, n_t); eps <- matrix(e, M, n_t)
  rm(g, e)
}))

ratio <- function(e, c) if (c > 0) e / c else NA_real_

one_rep <- function(i, truth) {
  set.seed(1000 * match(truth, names(truths)) + i)
  r <- truths[[truth]]
  z <- matrix(0L, M, n_t)
  z[, 1] <- rbinom(M, 1, r$psi1)
  for (t in 1:(n_t - 1)) z[, t + 1] <- rbinom(M, 1, ifelse(z[, t] == 1, 1 - r$eps[, t], r$gam[, t]))
  y <- matrix(NA_integer_, M, n_t * J)
  y[mask] <- rbinom(sum(mask), 1, (r$p * z[, rep(1:n_t, each = J)])[mask])
  # outcome-dependent survey continuation
  for (t in 1:(n_t - 1)) {
    cols_t <- ((t - 1) * J + 1):(t * J); cols_n <- (t * J + 1):((t + 1) * J)
    rec <- rowSums(y[, cols_t, drop = FALSE] == 1, na.rm = TRUE) > 0
    drop <- rec & runif(M) < p_drop
    y[drop, cols_n] <- NA
  }
  # refit the primary model
  umf <- mod@data
  yy <- getY(umf); yy[D$keep, ] <- y; umf@y <- yy
  fit <- tryCatch(colext(psiformula = fl$psi, gammaformula = fl$col, epsilonformula = fl$ext,
                         pformula = fl$det, data = umf, starts = coef(mod), control = list(maxit = 3000)),
                  error = function(e) NULL)
  if (is.null(fit)) return(data.frame(truth = truth, rep = i, converged = FALSE))
  Ds <- D; Ds$y <- y
  fb <- forward_backward(colext_rates(coef(fit), Ds), Ds)
  surv <- sapply(1:n_t, function(t) rowSums(!is.na(y[, ((t - 1) * J + 1):(t * J), drop = FALSE])) > 0)
  recd <- sapply(1:n_t, function(t) rowSums(y[, ((t - 1) * J + 1):(t * J), drop = FALSE] == 1, na.rm = TRUE) > 0)
  det <- rowSums(recd) > 0
  w <- surv[, 1:(n_t - 1)] & surv[, 2:n_t] & det
  w[, -RECENT] <- FALSE
  n0 <- recd[, 1:(n_t - 1)]; n1 <- recd[, 2:n_t]
  t_ext <- z[, 1:(n_t - 1)] == 1 & z[, 2:n_t] == 0; t_col <- z[, 1:(n_t - 1)] == 0 & z[, 2:n_t] == 1
  data.frame(truth = truth, rep = i, converged = fit@opt$convergence == 0,
             n_anchored = sum(w),
             est_ratio = ratio(sum(fb$ext[w]), sum(fb$col[w])),
             naive_ratio = ratio(sum(w & n0 & !n1), sum(w & !n0 & n1)),
             true_ratio_anchored = ratio(sum(t_ext[w]), sum(t_col[w])),
             true_ext_anchored = sum(t_ext[w]), true_col_anchored = sum(t_col[w]))
}

t0 <- Sys.time()
res <- do.call(rbind, lapply(names(truths), function(tr) {
  message("=== ", tr)
  do.call(rbind, mclapply(seq_len(N_SIM), one_rep, truth = tr, mc.cores = N_CORES))
}))
message("simulated in ", round(as.numeric(Sys.time() - t0, units = "mins"), 1), " min")
write.csv(res, file.path(out_dir, "anchored_simulation.csv"), row.names = FALSE)

q <- function(x, p) as.numeric(quantile(x, p, na.rm = TRUE))
summ <- do.call(rbind, lapply(split(res[res$converged %in% TRUE, ], res$truth[res$converged %in% TRUE]), function(d)
  data.frame(truth = d$truth[1], n_ok = nrow(d),
             est_med = median(d$est_ratio, na.rm = TRUE), est_lo = q(d$est_ratio, .025), est_hi = q(d$est_ratio, .975),
             true_anch_med = median(d$true_ratio_anchored, na.rm = TRUE),
             naive_med = median(d$naive_ratio, na.rm = TRUE),
             rel_bias_med = median(d$est_ratio / d$true_ratio_anchored, na.rm = TRUE),
             p_est_gt1 = mean(d$est_ratio > 1, na.rm = TRUE),
             p_est_ge_obs = mean(d$est_ratio >= 1.87, na.rm = TRUE))))
write.csv(summ, file.path(out_dir, "anchored_simulation_summary.csv"), row.names = FALSE)
print(summ, digits = 3)

sessionInfo()
