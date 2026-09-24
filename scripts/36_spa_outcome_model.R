###############################################################################
# 36_spa_outcome_model.R
#
# Two checks on the SPA extinction contrast of script 23 (primary model,
# anchored recent transitions, cells with >= 1 detection):
#
# (1) Outcome-model (regression-adjustment) estimator. Nearest-neighbour
#     matching with replacement reuses few controls (18-25 unique controls for
#     40-66 SPA cells), and resampling can misstate the variance of such
#     estimators (Abadie & Imbens 2008). Here, for every coefficient draw, a
#     quasi-Poisson GLM of expected extinctions with an offset for expected
#     occupied cell-years is fitted to all SPA (>= 50 % cover) and control
#     (< 10 %) cells:
#         E ~ SPA + pseudo-steppe + logit(psi1) + checklists per cell-year
#     No control is reused, and the GLM standard errors describe the
#     between-cell variability. The SPA effect is reported as a rate ratio and
#     as the average difference in extinction rate for SPA cells (predicted
#     with and without SPA status), drawing the GLM coefficients from their
#     sampling distribution within each coefficient draw.
#
# (2) Matching without initial occupancy. psi1 is predicted from static
#     environmental covariates only, but a reviewer may see it as potentially
#     affected by protection; matching on pseudo-steppe alone and on
#     pseudo-steppe plus survey effort removes it.
#
# Inputs:  results/{sp}_model_object.rds, colext inputs,
#          results/core_spa/cell_habitat.csv (script 23)
# Output:  results/core_spa/spa_outcome_model.csv
###############################################################################

library(here)
library(unmarked)
library(MASS)

source(here("R", "colext_utils.R"))
source(here("R", "spa_matching.R"))
set.seed(42)

species <- c("otitar", "tettet", "pteori")
N_DRAWS <- 1000
hab <- read.csv(here("results", "core_spa", "cell_habitat.csv"))
q <- function(x, p) as.numeric(quantile(x, p, na.rm = TRUE))
summ <- function(x) { x <- x[is.finite(x)]; c(med = median(x), lo = q(x, .025), hi = q(x, .975)) }

rows <- list()
for (sp in species) {
  message("=== ", sp)
  s <- spa_setup(sp, "primary", hab)
  m <- spa_matches(s, specs = c("base", "steppe", "steppe_effort"))
  cells <- c(s$tr, s$ctl)
  X <- data.frame(spa = rep(c(1, 0), c(length(s$tr), length(s$ctl))),
                  steppe = s$h$steppe_strict[cells], lpsi = s$lpsi[cells],
                  visits = s$visits[cells])
  draws <- draw_coefs(s$mod, N_DRAWS)
  point_only <- attr(draws, "point_only")
  glm_rr <- glm_diff <- rep(NA_real_, N_DRAWS); n_fail <- 0
  match_d <- match_r <- matrix(NA_real_, N_DRAWS, length(m), dimnames = list(NULL, names(m)))
  rec <- s$recent
  for (k in seq_len(N_DRAWS)) {
    fb <- forward_backward(colext_rates(draws[k, ], s$D), s$D)
    E <- rowSums((fb$ext * s$both)[, rec])
    O <- rowSums((fb$z1[, 1:(s$n_t - 1)] * s$both)[, rec])
    # (1) Outcome model on cells with expected occupancy
    d <- cbind(X, E = E[cells], O = O[cells])
    d <- d[d$O > 1e-6, ]
    fit <- tryCatch(suppressWarnings(glm(E ~ spa + steppe + lpsi + visits, offset = log(O),
                                         family = quasipoisson, data = d)),
                    error = function(e) NULL)
    if (is.null(fit) || !fit$converged || anyNA(coef(fit))) {
      n_fail <- n_fail + 1
    } else {
      b <- MASS::mvrnorm(1, coef(fit), vcov(fit))
      glm_rr[k] <- exp(b[["spa"]])
      tt <- d[d$spa == 1, ]
      mm <- model.matrix(~ spa + steppe + lpsi + visits, tt)
      mm0 <- mm; mm0[, "spa"] <- 0
      glm_diff[k] <- (sum(exp(mm %*% b) * tt$O) - sum(exp(mm0 %*% b) * tt$O)) / sum(tt$O)
    }
    # (2) Matching without psi1 (base kept as reference, same resample)
    bb <- sample(seq_along(s$tr), replace = TRUE)
    for (nm in names(m)) {
      t_i <- s$tr[bb]; c_i <- m[[nm]][bb]
      r_t <- sum(E[t_i]) / sum(O[t_i]); r_c <- sum(E[c_i]) / sum(O[c_i])
      match_d[k, nm] <- r_t - r_c; match_r[k, nm] <- r_t / r_c
    }
  }
  add <- function(method, spec, dd, rr, n_unique = NA, smd_psi = NA) {
    a <- summ(dd); r <- summ(rr)
    rows[[paste(sp, method, spec)]] <<- data.frame(
      species = sp, point_only = point_only, method = method, spec = spec,
      n_spa = length(s$tr), n_controls = length(s$ctl), n_unique_controls = n_unique,
      smd_after_lpsi = smd_psi, n_failed = sum(!is.finite(dd)),
      diff_med = a[["med"]], diff_lo = a[["lo"]], diff_hi = a[["hi"]],
      ratio_med = r[["med"]], ratio_lo = r[["lo"]], ratio_hi = r[["hi"]])
  }
  add("outcome_model", "steppe + psi1 + effort", glm_diff, glm_rr)
  for (nm in names(m)) {
    add("matching", nm, match_d[, nm], match_r[, nm], length(unique(m[[nm]])),
        smd(s$lpsi[s$tr], s$lpsi[m[[nm]]]))
  }
  message("GLM failures: ", n_fail)
}
res <- do.call(rbind, rows); rownames(res) <- NULL
write.csv(res, here("results", "core_spa", "spa_outcome_model.csv"), row.names = FALSE)
print(res, digits = 3)

sessionInfo()
