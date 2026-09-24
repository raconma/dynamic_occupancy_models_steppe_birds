###############################################################################
# spa_matching.R
#
# Per-species set-up of the SPA contrast shared by scripts 33 onwards: naive
# detection histories, anchored transitions, recorded cells, initial
# occupancy and survey effort, and the nearest-neighbour Mahalanobis matching
# of SPA cells (>= 50 % cover) to control cells (< 10 %).
#
# Script 23 holds its own copy of this logic (it produced the published
# tables and is kept unchanged); keep the two in sync. Scripts 33 onwards
# should use these helpers.
#
#   spa_setup()  model design, naive histories, anchored transitions, cell
#                covariates and SPA / control indices for one species
#   match_nn()   index of the nearest control for each treated unit
#   spa_matches() matched control for every SPA cell under each specification
###############################################################################

match_nn <- function(Xt, Xc) {
  S_inv <- solve(cov(rbind(Xt, Xc)))
  apply(Xt, 1, function(xt) which.min(mahalanobis(Xc, xt, S_inv, inverted = TRUE)))
}

smd <- function(a, b) (mean(a) - mean(b)) / sqrt((var(a) + var(b)) / 2)

# `hab` is results/core_spa/cell_habitat.csv (script 23)
spa_setup <- function(sp, variant, hab) {
  mod <- readRDS(colext_model_path(sp, variant))
  D <- colext_design(mod)
  n_t <- D$n_t
  occ <- read.csv(colext_input_path(sp, variant))[D$keep, ]
  h <- hab[match(occ$cells, hab$cells), ]
  stopifnot(!anyNA(h$cells))
  naive <- sapply(1:n_t, function(t) {
    yt <- D$y[, ((t - 1) * D$J + 1):(t * D$J), drop = FALSE]
    ifelse(rowSums(!is.na(yt)) > 0, as.integer(rowSums(yt, na.rm = TRUE) > 0), NA_integer_)
  })
  surv <- !is.na(naive)
  both <- surv[, 1:(n_t - 1)] & surv[, 2:n_t]
  det <- rowSums(naive, na.rm = TRUE) > 0
  psi1 <- colext_rates(coef(mod), D)$psi1
  list(mod = mod, D = D, n_t = n_t, h = h, naive = naive, both = both,
       det = det, recent = 2:(n_t - 1),
       lpsi = qlogis(pmin(pmax(psi1, 1e-6), 1 - 1e-6)),
       visits = rowSums(!is.na(D$y)) / pmax(rowSums(surv), 1),
       tr = which(det & h$spa_frac >= 0.5),
       ctl = which(det & h$spa_frac < 0.1))
}

# Matched control for each SPA cell; base, coords and effort as in script 23
spa_matches <- function(s, specs = c("base", "effort"), ctl = s$ctl) {
  covs <- list(base = cbind(s$h$steppe_strict, s$lpsi),
               coords = cbind(s$h$steppe_strict, s$lpsi, s$h$x / 1e5, s$h$y / 1e5),
               effort = cbind(s$h$steppe_strict, s$lpsi, s$visits),
               # without initial occupancy (sensitivity to matching on psi1)
               steppe = cbind(s$h$steppe_strict),
               steppe_effort = cbind(s$h$steppe_strict, s$visits))
  lapply(setNames(specs, specs), function(nm) {
    X <- covs[[nm]]
    ctl[match_nn(X[s$tr, , drop = FALSE], X[ctl, , drop = FALSE])]
  })
}
