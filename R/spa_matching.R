###############################################################################
# spa_matching.R
#
# Per-species set-up of the SPA contrast shared by scripts 33 onwards: naive
# detection histories, anchored transitions, cells with anchored transitions, initial
# occupancy and survey effort, and the nearest-neighbour Mahalanobis matching
# of SPA cells (>= 50 % cover) to control cells (< 10 %).
#
# Script 23 holds its own copy of this logic (with the same anchored-cell
# restriction, block resampling and summaries); keep the two in sync.
#
#   spa_setup()  model design, naive histories, anchored transitions, cell
#                covariates and SPA / control indices for one species
#   match_nn()   index of the nearest control for each treated unit
#   spa_matches() matched control for every SPA cell under each specification
#   block_sample() positions drawn by resampling spatial blocks
#   summ_diff(), summ_ratio() draw summaries (as in script 23)
###############################################################################

block_sample <- function(blk) {
  ub <- unique(blk)
  unlist(lapply(sample(ub, length(ub), replace = TRUE), function(u) which(blk == u)))
}

summ_diff <- function(d) {
  ok <- is.finite(d)
  qd <- function(p) as.numeric(quantile(d[ok], p))
  c(n_undefined = sum(!ok), med = median(d[ok]), lo = qd(.025), hi = qd(.975),
    p_pos = mean(d[ok] > 0))
}

# Ratios keep infinite draws; NA when > 5 % of draws are 0/0
summ_ratio <- function(r) {
  nan <- is.na(r)
  if (mean(nan) > 0.05) return(c(n_undefined = sum(nan), n_inf = sum(is.infinite(r)),
                                 med = NA, lo = NA, hi = NA))
  rr <- r[!nan]
  qq <- function(p) as.numeric(quantile(rr, p, type = 1))
  c(n_undefined = sum(nan), n_inf = sum(is.infinite(rr)), med = qq(.5), lo = qq(.025),
    hi = qq(.975))
}

match_nn <- function(Xt, Xc) {
  S_inv <- solve(cov(rbind(Xt, Xc)))
  apply(Xt, 1, function(xt) which.min(mahalanobis(Xc, xt, S_inv, inverted = TRUE)))
}

smd <- function(a, b) (mean(a) - mean(b)) / sqrt((var(a) + var(b)) / 2)

# `hab` is results/core_spa/cell_habitat.csv (script 23)
# Block side for spatial resampling (m); BLOCK_KM environment variable for
# the block-size sensitivity of script 23
BLOCK_M <- as.numeric(Sys.getenv("BLOCK_KM", "50")) * 1000

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
  recent <- 2:(n_t - 1)
  # Only cells with >= 1 anchored recent transition inform the contrasts
  anch <- det & rowSums(both[, recent, drop = FALSE]) > 0
  list(mod = mod, D = D, n_t = n_t, h = h, naive = naive, both = both,
       det = det, anch = anch, recent = recent,
       blk = paste(floor(h$x / BLOCK_M), floor(h$y / BLOCK_M)),   # 50 x 50 km blocks by default
       lpsi = qlogis(pmin(pmax(psi1, 1e-6), 1 - 1e-6)),
       visits = rowSums(!is.na(D$y)) / pmax(rowSums(surv), 1),
       tr = which(anch & h$spa_frac >= 0.5),
       ctl = which(anch & h$spa_frac < 0.1))
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
