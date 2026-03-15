###############################################################################
# simulation_robustness_utils.R
#
# Helper functions for the simulation-based robustness analysis.
# Tests whether colext() can recover known colonisation rates under
# the same sampling design as the empirical study.
#
# Sourced by: scripts/17_simulation_robustness.R
###############################################################################


# --- simulate_dynamic_occ -------------------------------------------------
#
# Generate a simulated detection-history dataset from a constant-parameter
# dynamic occupancy process (HMM), returning an unmarkedMultFrame ready
# for colext().
#
# Arguments:
#   N_sites    integer, number of sites
#   T_years    integer, number of primary periods
#   J_visits   integer, max secondary periods per primary
#   psi1       numeric [0,1], initial occupancy probability
#   gamma      numeric [0,1], colonisation probability
#   epsilon    numeric [0,1], extinction probability
#   p          numeric [0,1], per-visit detection probability
#   na_pattern optional N_sites x T_years integer matrix: actual visit
#              counts per site-year. Visits beyond this count are set to NA.
#
# Returns: unmarkedMultFrame (no covariates, numPrimary = T_years)

simulate_dynamic_occ <- function(N_sites, T_years, J_visits,
                                 psi1, gamma, epsilon, p,
                                 na_pattern = NULL) {

  # --- Latent state process (HMM) ---
  Z <- matrix(NA_integer_, nrow = N_sites, ncol = T_years)
  Z[, 1] <- rbinom(N_sites, 1, psi1)

  for (t in 2:T_years) {
    # occupied sites: persist with prob 1 - epsilon
    # unoccupied sites: colonised with prob gamma
    Z[, t] <- Z[, t - 1] * rbinom(N_sites, 1, 1 - epsilon) +
              (1 - Z[, t - 1]) * rbinom(N_sites, 1, gamma)
  }

  # --- Observation process ---
  # y matrix: N_sites x (T_years * J_visits), year-major order
  y <- matrix(NA_integer_, nrow = N_sites, ncol = T_years * J_visits)

  for (t in seq_len(T_years)) {
    cols <- ((t - 1) * J_visits + 1):(t * J_visits)
    for (j in seq_len(J_visits)) {
      # detection conditional on true occupancy
      y[, cols[j]] <- Z[, t] * rbinom(N_sites, 1, p)
    }
  }

  # --- Apply effort heterogeneity (NA pattern) ---
  if (!is.null(na_pattern)) {
    for (t in seq_len(T_years)) {
      cols <- ((t - 1) * J_visits + 1):(t * J_visits)
      for (j in seq_len(J_visits)) {
        # set to NA where the real data had fewer visits
        missing <- na_pattern[, t] < j
        y[missing, cols[j]] <- NA_integer_
      }
    }
  }

  # --- Build unmarkedMultFrame ---
  unmarked::unmarkedMultFrame(
    y = y,
    numPrimary = T_years
  )
}


# --- extract_effort_pattern -----------------------------------------------
#
# Extract the actual number of non-NA visits per site-year from the
# empirical detection histories.
#
# Arguments:
#   sp_code   character, species code (otitar, ptealc, pteori, tettet)
#   data_dir  path to processed data directory
#   years     integer vector of study years
#
# Returns: N_sites x T_years integer matrix (visit counts per site-year)

extract_effort_pattern <- function(sp_code,
                                   data_dir = here::here("data", "processed_2023"),
                                   years = 2017:2023) {

  csv_path <- file.path(data_dir, sp_code,
                         paste0(sp_code, "_occ_wide_dynamic.csv"))
  dat <- read.csv(csv_path, stringsAsFactors = FALSE)

  # Drop sites with missing NDVI in any year (matching model fitting pipeline)
  ndvi_cols <- paste0("NDVI_", years)
  ndvi_cols <- ndvi_cols[ndvi_cols %in% names(dat)]
  if (length(ndvi_cols) > 0) {
    complete <- complete.cases(dat[, ndvi_cols, drop = FALSE])
    dat <- dat[complete, ]
  }

  T_years <- length(years)
  N_sites <- nrow(dat)
  effort <- matrix(0L, nrow = N_sites, ncol = T_years)

  for (t in seq_along(years)) {
    yr <- years[t]
    # detection columns for this year: y.1.YYYY, y.2.YYYY, ..., y.10.YYYY
    y_cols <- paste0("y.", 1:10, ".", yr)
    y_cols <- y_cols[y_cols %in% names(dat)]
    if (length(y_cols) > 0) {
      y_mat <- as.matrix(dat[, y_cols, drop = FALSE])
      effort[, t] <- rowSums(!is.na(y_mat))
    }
  }

  effort
}


# --- fit_intercept_colext -------------------------------------------------
#
# Fit an intercept-only colext model with error handling and retries.
#
# Arguments:
#   umf        unmarkedMultFrame
#   max_tries  number of fitting attempts (jittered starts on retry)
#
# Returns: list with converged, parameter estimates, CIs, AIC, warnings

fit_intercept_colext <- function(umf, max_tries = 3) {

  result <- list(
    converged    = FALSE,
    psi1_hat     = NA_real_,
    gamma_hat    = NA_real_,
    epsilon_hat  = NA_real_,
    p_hat        = NA_real_,
    gamma_se     = NA_real_,
    epsilon_se   = NA_real_,
    gamma_ci_lo  = NA_real_,
    gamma_ci_hi  = NA_real_,
    epsilon_ci_lo = NA_real_,
    epsilon_ci_hi = NA_real_,
    aic          = NA_real_,
    boundary     = FALSE,
    warnings     = character(0)
  )

  fit <- NULL

  for (attempt in seq_len(max_tries)) {

    # starting values: jitter on retries
    starts <- if (attempt == 1) NULL else rnorm(4, 0, 1)

    fit <- tryCatch({
      w <- character(0)
      mod <- withCallingHandlers(
        unmarked::colext(
          psiformula     = ~ 1,
          gammaformula   = ~ 1,
          epsilonformula = ~ 1,
          pformula       = ~ 1,
          data   = umf,
          starts = starts,
          se     = TRUE
        ),
        warning = function(cond) {
          w <<- c(w, conditionMessage(cond))
          invokeRestart("muffleWarning")
        }
      )
      list(mod = mod, warnings = w)
    },
    error = function(e) {
      list(mod = NULL, warnings = conditionMessage(e))
    })

    if (!is.null(fit$mod)) break
  }

  if (is.null(fit$mod)) {
    result$warnings <- fit$warnings
    return(result)
  }

  # extract coefficients
  coefs <- coef(fit$mod)
  result$converged <- TRUE
  result$warnings  <- fit$warnings
  result$aic       <- fit$mod@AIC

  result$psi1_hat    <- plogis(coefs["psi(Int)"])
  result$gamma_hat   <- plogis(coefs["col(Int)"])
  result$epsilon_hat <- plogis(coefs["ext(Int)"])
  result$p_hat       <- plogis(coefs["p(Int)"])

  # boundary detection: gamma intercept pushed far negative
  if (coefs["col(Int)"] < -15) result$boundary <- TRUE

  # standard errors and Wald CIs on probability scale
  se_ok <- tryCatch({
    ses <- sqrt(diag(vcov(fit$mod)))
    names(ses) <- names(coefs)

    col_se <- ses["col(Int)"]
    ext_se <- ses["ext(Int)"]

    result$gamma_se   <- col_se
    result$epsilon_se <- ext_se

    # 95% Wald CI on logit scale, then transform
    col_logit <- coefs["col(Int)"]
    ext_logit <- coefs["ext(Int)"]

    result$gamma_ci_lo   <- plogis(col_logit - 1.96 * col_se)
    result$gamma_ci_hi   <- plogis(col_logit + 1.96 * col_se)
    result$epsilon_ci_lo <- plogis(ext_logit - 1.96 * ext_se)
    result$epsilon_ci_hi <- plogis(ext_logit + 1.96 * ext_se)

    TRUE
  }, error = function(e) FALSE)

  result
}


# --- classify_regime ------------------------------------------------------
#
# Classify the demographic regime based on the epsilon/gamma ratio.
#
# Returns: character, one of "trap", "asymmetric", "mild_decline", "growth"

classify_regime <- function(gamma, epsilon) {
  if (is.na(gamma) || is.na(epsilon) || gamma <= 0) return("trap")
  ratio <- epsilon / gamma
  if (ratio > 100)  return("trap")
  if (ratio > 10)   return("asymmetric")
  if (ratio > 1)    return("mild_decline")
  return("growth")
}


# --- compute_sim_metrics --------------------------------------------------
#
# Compute evaluation metrics for a single simulation replicate.
#
# Arguments:
#   true_params  list(psi1, gamma, epsilon, p)
#   est_params   list(psi1, gamma, epsilon, p) -- estimated
#   est_ci       list(gamma_ci_lo, gamma_ci_hi, epsilon_ci_lo, epsilon_ci_hi)
#
# Returns: one-row data.frame of metrics

compute_sim_metrics <- function(true_params, est_params, est_ci) {

  bias_gamma   <- est_params$gamma - true_params$gamma
  bias_epsilon <- est_params$epsilon - true_params$epsilon
  bias_p       <- est_params$p - true_params$p

  rel_bias_gamma <- if (true_params$gamma > 1e-12) {
    (est_params$gamma - true_params$gamma) / true_params$gamma
  } else {
    NA_real_
  }
  rel_bias_epsilon <- (est_params$epsilon - true_params$epsilon) / true_params$epsilon

  # CI coverage
  ci_covers_gamma <- !is.na(est_ci$gamma_ci_lo) &&
    est_ci$gamma_ci_lo <= true_params$gamma &&
    est_ci$gamma_ci_hi >= true_params$gamma

  ci_covers_epsilon <- !is.na(est_ci$epsilon_ci_lo) &&
    est_ci$epsilon_ci_lo <= true_params$epsilon &&
    est_ci$epsilon_ci_hi >= true_params$epsilon

  # Ratio and regime
  ratio_hat  <- est_params$epsilon / max(est_params$gamma, 1e-20)
  ratio_true <- true_params$epsilon / max(true_params$gamma, 1e-20)

  regime_hat  <- classify_regime(est_params$gamma, est_params$epsilon)
  regime_true <- classify_regime(true_params$gamma, true_params$epsilon)

  data.frame(
    bias_gamma      = bias_gamma,
    bias_epsilon     = bias_epsilon,
    bias_p           = bias_p,
    rel_bias_gamma   = rel_bias_gamma,
    rel_bias_epsilon = rel_bias_epsilon,
    ci_covers_gamma  = ci_covers_gamma,
    ci_covers_epsilon = ci_covers_epsilon,
    ratio_hat        = ratio_hat,
    ratio_true       = ratio_true,
    regime_hat       = regime_hat,
    regime_true      = regime_true,
    regime_correct   = (regime_hat == regime_true),
    underest_10x     = est_params$gamma < true_params$gamma / 10,
    underest_100x    = est_params$gamma < true_params$gamma / 100,
    stringsAsFactors = FALSE
  )
}
