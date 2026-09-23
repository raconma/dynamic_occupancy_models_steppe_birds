###############################################################################
# colext_utils.R
#
# Helpers to recompute site-year quantities from fitted unmarked colext models
# for arbitrary coefficient vectors (e.g. multivariate-normal draws), used by
# scripts 19 and 21.
#
#   colext_design()   design matrices for psi1, gamma, epsilon and p
#   colext_rates()    psi1, gamma, epsilon (M x T; column T = transition from
#                     the last year, used only for one-step-ahead risk) and p
#   project_occ()     unconditional (projected) occupancy trajectory
#   forward_backward() posterior state and transition probabilities given the
#                     detection histories (two-state hidden Markov model)
###############################################################################

# Path to a fitted model: "primary" = results/{sp}_model_object.rds; any
# other variant = results/robustness/{sp}_{variant}_model_object.rds (28)
colext_model_path <- function(sp, variant = "primary") {
  if (identical(variant, "primary")) {
    here::here("results", paste0(sp, "_model_object.rds"))
  } else {
    here::here("results", "robustness", paste0(sp, "_", variant, "_model_object.rds"))
  }
}

# Site-level input table (cell ids, detection histories) matching a model's
# rows: the stepRep fits used the filtered inputs in reports/script_4_v2_inputs
colext_input_path <- function(sp, variant = "primary") {
  if (startsWith(variant, "steprep")) {
    here::here("reports", "script_4_v2_inputs", paste0(sp, "_occ_wide_dynamic.csv"))
  } else {
    here::here("data", "processed_2023", sp, paste0(sp, "_occ_wide_dynamic.csv"))
  }
}

# Multivariate-normal coefficient draws. When the Hessian gives non-finite
# or non-positive variances (parameters at the optimiser boundary) the point
# estimate is repeated and attr(, "point_only") is TRUE, so callers can flag
# intervals as unavailable rather than report an arbitrary spread.
draw_coefs <- function(mod, n) {
  cf <- coef(mod)
  V <- suppressWarnings(vcov(mod))
  if (!all(is.finite(V)) || any(diag(V) <= 0)) {
    d <- matrix(cf, n, length(cf), byrow = TRUE, dimnames = list(NULL, names(cf)))
    attr(d, "point_only") <- TRUE
    return(d)
  }
  V_pd <- tryCatch({ chol(V); V },
                   error = function(e) as.matrix(Matrix::nearPD(V)$mat))
  d <- MASS::mvrnorm(n, mu = cf, Sigma = V_pd)
  attr(d, "point_only") <- FALSE
  d
}

# Design matrix that keeps NA rows (model.matrix drops them by default)
design <- function(formula, data) {
  mf <- model.frame(formula, data, na.action = na.pass)
  X <- model.matrix(formula, mf)
  X_full <- matrix(NA_real_, nrow(data), ncol(X),
                   dimnames = list(NULL, colnames(X)))
  X_full[as.integer(rownames(X)), ] <- X
  X_full
}

colext_design <- function(mod) {
  umf <- mod@data
  n_t <- umf@numPrimary
  keep <- setdiff(seq_len(numSites(umf)), mod@sitesRemoved)
  M <- length(keep)
  y <- getY(umf)[keep, , drop = FALSE]
  J <- ncol(y) / n_t

  sc <- siteCovs(umf)[keep, , drop = FALSE]
  rownames(sc) <- NULL
  row_yr <- as.vector(t(outer((keep - 1) * n_t, 1:n_t, "+")))
  ysc <- yearlySiteCovs(umf)[row_yr, , drop = FALSE]
  ysc <- cbind(ysc, sc[rep(seq_len(M), each = n_t), , drop = FALSE])
  rownames(ysc) <- NULL

  # Observation covariates with site and yearly covariates available too
  row_obs <- as.vector(t(outer((keep - 1) * J * n_t, 1:(J * n_t), "+")))
  oc <- obsCovs(umf)[row_obs, , drop = FALSE]
  oc <- cbind(oc, ysc[rep(seq_len(M * n_t), each = J), , drop = FALSE])
  rownames(oc) <- NULL

  fl <- mod@formlist
  cf <- coef(mod)
  idx <- list(psi = grep("^psi\\(", names(cf)), col = grep("^col\\(", names(cf)),
              ext = grep("^ext\\(", names(cf)), det = grep("^p\\(", names(cf)))
  D <- list(M = M, n_t = n_t, J = J, keep = keep, y = y, idx = idx,
            X_psi = design(fl$psi, sc), X_col = design(fl$col, ysc),
            X_ext = design(fl$ext, ysc), X_det = design(fl$det, oc))
  stopifnot(length(idx$psi) == ncol(D$X_psi), length(idx$col) == ncol(D$X_col),
            length(idx$ext) == ncol(D$X_ext), length(idx$det) == ncol(D$X_det))
  D
}

# Rates for coefficient vector b. gam/eps are M x T: column t is the
# transition t -> t+1 (yearly covariates of year t); column T is only used
# for one-step-ahead projections beyond the study period.
colext_rates <- function(b, D) {
  with(D, list(
    psi1 = as.vector(plogis(X_psi %*% b[idx$psi])),
    gam  = matrix(plogis(X_col %*% b[idx$col]), M, n_t, byrow = TRUE),
    eps  = matrix(plogis(X_ext %*% b[idx$ext]), M, n_t, byrow = TRUE),
    p    = matrix(plogis(X_det %*% b[idx$det]), M, n_t * J, byrow = TRUE)
  ))
}

# Projected occupancy psi[i, t] from psi1 and transitions (M x >= T-1)
project_occ <- function(psi1, gam, eps, n_t = ncol(gam) + 1) {
  psi <- matrix(NA_real_, length(psi1), n_t)
  psi[, 1] <- psi1
  for (t in 1:(n_t - 1)) {
    psi[, t + 1] <- psi[, t] * (1 - eps[, t]) + (1 - psi[, t]) * gam[, t]
  }
  psi
}

# Posterior quantities of the two-state HMM given detection histories.
# Returns z1 (P(z_it = 1 | y), M x T) and ext/col (P of an extinction /
# colonisation event between t and t+1 given y, M x T-1).
forward_backward <- function(r, D) {
  M <- D$M; n_t <- D$n_t; J <- D$J
  y <- D$y
  e1 <- e0 <- matrix(1, M, n_t)
  for (t in 1:n_t) {
    cols <- ((t - 1) * J + 1):(t * J)
    yt <- y[, cols, drop = FALSE]
    pt <- r$p[, cols, drop = FALSE]
    lik1 <- ifelse(is.na(yt), 1, ifelse(yt == 1, pt, 1 - pt))
    e1[, t] <- apply(lik1, 1, prod)
    e0[, t] <- as.numeric(rowSums(yt, na.rm = TRUE) == 0)
  }
  g <- r$gam; e <- r$eps
  A0 <- A1 <- B0 <- B1 <- matrix(NA_real_, M, n_t)
  a0 <- (1 - r$psi1) * e0[, 1]; a1 <- r$psi1 * e1[, 1]
  s <- a0 + a1; A0[, 1] <- a0 / s; A1[, 1] <- a1 / s
  for (t in 1:(n_t - 1)) {
    a0 <- (A0[, t] * (1 - g[, t]) + A1[, t] * e[, t]) * e0[, t + 1]
    a1 <- (A0[, t] * g[, t] + A1[, t] * (1 - e[, t])) * e1[, t + 1]
    s <- a0 + a1; A0[, t + 1] <- a0 / s; A1[, t + 1] <- a1 / s
  }
  B0[, n_t] <- B1[, n_t] <- 1
  for (t in (n_t - 1):1) {
    b0 <- (1 - g[, t]) * e0[, t + 1] * B0[, t + 1] + g[, t] * e1[, t + 1] * B1[, t + 1]
    b1 <- e[, t] * e0[, t + 1] * B0[, t + 1] + (1 - e[, t]) * e1[, t + 1] * B1[, t + 1]
    s <- b0 + b1; B0[, t] <- b0 / s; B1[, t] <- b1 / s
  }
  z1 <- A1 * B1 / (A0 * B0 + A1 * B1)
  ext <- col <- matrix(NA_real_, M, n_t - 1)
  for (t in 1:(n_t - 1)) {
    x00 <- A0[, t] * (1 - g[, t]) * e0[, t + 1] * B0[, t + 1]
    x01 <- A0[, t] * g[, t] * e1[, t + 1] * B1[, t + 1]
    x10 <- A1[, t] * e[, t] * e0[, t + 1] * B0[, t + 1]
    x11 <- A1[, t] * (1 - e[, t]) * e1[, t + 1] * B1[, t + 1]
    s <- x00 + x01 + x10 + x11
    col[, t] <- x01 / s
    ext[, t] <- x10 / s
  }
  list(z1 = z1, ext = ext, col = col)
}
