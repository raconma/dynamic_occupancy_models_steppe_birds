###############################################################################
# 6_spatial_occupancy_test.R
#
# Purpose: Test spatio-temporal occupancy model with explicit spatial structure
#          (Gaussian Process via NNGP) for Otis tarda, as a proof of concept.
#
# Approach: spOccupancy::stPGOcc — multi-season spatial occupancy with
#           Nearest Neighbor Gaussian Process (NNGP). This extends the current
#           colext() model by adding a spatial random effect, which should
#           absorb the residual spatial autocorrelation (Moran's I ~ 0.38).
#
# Why spOccupancy and not INLA:
#   - INLA requires aggregating detections within site-years (no obs-level
#     detection covariates like effort, duration, observers). That would
#     fundamentally change our detection model.
#   - stPGOcc handles the full observation-level detection structure.
#   - NNGP with m=5 neighbours scales to ~3000 sites in hours, not days.
#
# Important conceptual note:
#   stPGOcc is a "phenomenological" spatial model — it estimates year-specific
#   occupancy probability with a spatial random effect, but does NOT explicitly
#   parameterise colonisation/extinction as Markov transitions (unlike colext).
#   The comparison is therefore between:
#     (a) colext: dynamic (gamma/epsilon) but no spatial structure
#     (b) stPGOcc: spatial structure but phenomenological temporal dynamics
#   Both are valid and complementary. If colonisation/extinction parameters are
#   central to the paper, stPGOcc serves as a sensitivity check, not a replacement.
#
# Strategy (3 phases in one script):
#   Phase A: tPGOcc  — non-spatial multi-season baseline (fast, ~minutes)
#   Phase B: spPGOcc — single-year spatial proof of concept (~1-2 hours)
#   Phase C: stPGOcc — full spatio-temporal model (~6-24 hours)
#   Each phase is self-contained; you can run only Phase A+B to assess
#   feasibility before committing to Phase C.
#
# Inputs:  data-raw/data/otitar/otitar_occ_wide_dynamic.csv
# Outputs: results/otitar_spatial_*
#          figs/otitar_spatial_*
#
# Requirements: install.packages("spOccupancy")
###############################################################################

library(here)
library(sf)
library(dplyr)

# spOccupancy must be installed separately (not in renv yet)
if (!requireNamespace("spOccupancy", quietly = TRUE)) {
  stop("Install spOccupancy first: install.packages('spOccupancy')\n",
       "Then re-run this script.")
}
library(spOccupancy)

set.seed(42)
t_start_global <- Sys.time()

cat("====================================================================\n")
cat(" Spatial Occupancy Test — Otis tarda (otitar)\n")
cat("====================================================================\n\n")

###############################################################################
# 0. LOAD AND RESHAPE DATA
#    spOccupancy expects:
#      y:        3D array [J sites x T years x K replicates]
#      occ.covs: list of vectors (J) or matrices (J x T)
#      det.covs: list of vectors (J), matrices (J x T), or arrays (J x T x K)
#      coords:   J x 2 matrix in PROJECTED CRS (not lon/lat)
###############################################################################

cat("Loading data...\n")
d <- read.csv(here("data-raw", "data", "otitar", "otitar_occ_wide_dynamic.csv"))
cat("  Raw data:", nrow(d), "sites x", ncol(d), "columns\n")

years <- 2017:2022
n_years <- length(years)
n_reps  <- 10
J <- nrow(d)

# --- Detection history: y[J, T, K] ---
y_array <- array(NA, dim = c(J, n_years, n_reps))
for (t in seq_along(years)) {
  for (k in 1:n_reps) {
    col_name <- paste0("y.", k, ".", years[t])
    if (col_name %in% names(d)) {
      y_array[, t, k] <- d[[col_name]]
    }
  }
}
cat("  Detection array: ", paste(dim(y_array), collapse=" x "), "\n")
cat("  Total observations:", sum(!is.na(y_array)), "\n")
cat("  Detections:", sum(y_array == 1, na.rm = TRUE), "\n")
cat("  Naive occupancy:", round(mean(apply(y_array, 1, function(x) any(x == 1, na.rm = TRUE))), 3), "\n")

# --- Coordinates: project to ETRS89 / UTM 30N ---
coords_lonlat <- cbind(d$longitude, d$latitude)
pts_sf <- st_as_sf(data.frame(x = d$longitude, y = d$latitude),
                   coords = c("x", "y"), crs = 4326)
pts_utm <- st_transform(pts_sf, 25830)  # ETRS89 / UTM zone 30N
coords <- st_coordinates(pts_utm)
cat("  Coordinates reprojected to ETRS89/UTM30N\n")
cat("  X range:", round(range(coords[,1])/1000), "km\n")
cat("  Y range:", round(range(coords[,2])/1000), "km\n")

# Maximum inter-site distance (needed for phi prior)
# Use a subsample for speed
set.seed(1)
idx_sub <- sample(J, min(500, J))
max_dist <- max(dist(coords[idx_sub, ]))
cat("  Approx max distance:", round(max_dist/1000), "km\n\n")

# --- Occupancy covariates ---
# Static (site-level): bio1, bio2, tree_cover, grass_cover, topo_elev
# Already scaled in the dataset

occ_covs <- list(
  bio1       = d$bio1,
  bio2       = d$bio2,
  tree_cover = d$tree_cover,
  grass_cover = d$grass_cover,
  topo_elev  = d$topo_elev
)

# Dynamic (site x year): NDVI, tmmx
# These are stored as NDVI_2017, NDVI_2018, ... in the data
for (var_name in c("NDVI", "tmmx")) {
  mat <- matrix(NA, nrow = J, ncol = n_years)
  for (t in seq_along(years)) {
    col <- paste0(var_name, "_", years[t])
    if (col %in% names(d)) {
      vals <- d[[col]]
      # Scale if not already scaled
      if (abs(mean(vals, na.rm = TRUE)) > 10 || sd(vals, na.rm = TRUE) > 10) {
        vals <- as.numeric(scale(vals))
      }
      mat[, t] <- vals
    }
  }
  occ_covs[[var_name]] <- mat
}

cat("Occupancy covariates:\n")
for (nm in names(occ_covs)) {
  obj <- occ_covs[[nm]]
  if (is.matrix(obj)) {
    cat("  ", nm, ": matrix", nrow(obj), "x", ncol(obj),
        "| range:", round(range(obj, na.rm=TRUE), 2), "\n")
  } else {
    cat("  ", nm, ": vector length", length(obj),
        "| range:", round(range(obj, na.rm=TRUE), 2), "\n")
  }
}

# --- Detection covariates ---
# Observation-level (site x year x replicate): duration, observers

det_covs <- list()

for (var_base in c("duration_minutes", "number_observers")) {
  arr <- array(NA, dim = c(J, n_years, n_reps))
  for (t in seq_along(years)) {
    for (k in 1:n_reps) {
      col <- paste0(var_base, ".", k, ".", years[t])
      if (col %in% names(d)) {
        arr[, t, k] <- d[[col]]
      }
    }
  }
  # Scale
  vals_flat <- as.numeric(arr)
  mu <- mean(vals_flat, na.rm = TRUE)
  sigma <- sd(vals_flat, na.rm = TRUE)
  if (sigma > 0) arr <- (arr - mu) / sigma

  short_name <- ifelse(var_base == "duration_minutes", "duration", "observers")
  det_covs[[short_name]] <- arr
}

cat("\nDetection covariates:\n")
for (nm in names(det_covs)) {
  obj <- det_covs[[nm]]
  cat("  ", nm, ": array", paste(dim(obj), collapse=" x "),
      "| non-NA:", sum(!is.na(obj)), "\n")
}

# --- Assemble data list ---
data_list <- list(
  y = y_array,
  occ.covs = occ_covs,
  det.covs = det_covs,
  coords = coords
)

cat("\n")


###############################################################################
# PHASE A: tPGOcc — NON-SPATIAL MULTI-SEASON BASELINE
#
# This is a Bayesian multi-season occupancy model WITHOUT spatial random effects.
# Fast to fit (~minutes). Serves as:
#   1. Baseline for comparison with spatial model
#   2. Validates that the data formatting is correct
#   3. Provides initial parameter estimates
###############################################################################

cat("====================================================================\n")
cat(" Phase A: tPGOcc (non-spatial multi-season baseline)\n")
cat("====================================================================\n\n")

t_start_A <- Sys.time()

# Priors
priors_nonspatial <- list(
  beta.normal  = list(mean = 0, var = 2.72),   # occupancy coefficients
  alpha.normal = list(mean = 0, var = 2.72)     # detection coefficients
)

# MCMC settings (conservative for exploratory run)
n_batch_A  <- 400
batch_len  <- 25
n_samples_A <- n_batch_A * batch_len  # 10,000 per chain
n_burn_A   <- 5000
n_thin_A   <- 5
n_chains_A <- 3

cat("Fitting tPGOcc...\n")
cat("  MCMC: ", n_samples_A, "samples/chain x", n_chains_A, "chains\n")
cat("  Burn-in:", n_burn_A, "| Thin:", n_thin_A, "\n")
cat("  Post-thinning samples per chain:", (n_samples_A - n_burn_A) / n_thin_A, "\n\n")

fit_tpgocc <- tPGOcc(
  occ.formula  = ~ bio1 + bio2 + tree_cover + grass_cover + topo_elev + NDVI + tmmx,
  det.formula  = ~ duration + observers,
  data         = data_list,
  priors       = priors_nonspatial,
  n.batch      = n_batch_A,
  batch.length = batch_len,
  n.burn       = n_burn_A,
  n.thin       = n_thin_A,
  n.chains     = n_chains_A,
  ar1          = TRUE,
  verbose      = TRUE,
  n.report     = 100
)

t_end_A <- Sys.time()
cat("\n  Phase A time:", round(difftime(t_end_A, t_start_A, units = "mins"), 1), "min\n\n")

# --- Summary ---
cat("--- tPGOcc Summary ---\n")
summary(fit_tpgocc)

# Save
saveRDS(fit_tpgocc, here("results", "otitar_spatial_tPGOcc.rds"))

# WAIC for model comparison
waic_A <- waicOcc(fit_tpgocc)
cat("\ntPGOcc WAIC:", waic_A$WAIC, "\n\n")


###############################################################################
# PHASE B: spPGOcc — SINGLE-YEAR SPATIAL PROOF OF CONCEPT
#
# Pick the year with best data coverage (2021: 1218 sites) and fit a
# single-season spatial occupancy model with NNGP. This tests:
#   1. Whether NNGP machinery works with our data/coords
#   2. What the spatial range parameter (phi) looks like
#   3. Approximate wall-clock time per MCMC iteration
###############################################################################

cat("====================================================================\n")
cat(" Phase B: spPGOcc (single-year spatial, year = 2021)\n")
cat("====================================================================\n\n")

t_start_B <- Sys.time()

# Select year 2021 (index 5 in 2017:2022)
yr_idx <- 5
yr <- years[yr_idx]

# Extract single-year data
y_single <- y_array[, yr_idx, ]   # J x K matrix

# Keep only sites with at least 1 observation in this year
has_data <- rowSums(!is.na(y_single)) > 0
cat("  Year", yr, ":", sum(has_data), "sites with data\n")

y_sub <- y_single[has_data, ]
coords_sub <- coords[has_data, ]

# Occupancy covariates (static only for single-year)
occ_covs_sub <- data.frame(
  bio1        = d$bio1[has_data],
  bio2        = d$bio2[has_data],
  tree_cover  = d$tree_cover[has_data],
  grass_cover = d$grass_cover[has_data],
  topo_elev   = d$topo_elev[has_data]
)

# Add NDVI for this specific year
ndvi_col <- paste0("NDVI_", yr)
if (ndvi_col %in% names(d)) {
  ndvi_vals <- d[[ndvi_col]][has_data]
  if (abs(mean(ndvi_vals, na.rm = TRUE)) > 10) ndvi_vals <- as.numeric(scale(ndvi_vals))
  occ_covs_sub$NDVI <- ndvi_vals
}

# Detection covariates (observation-level)
det_covs_sub <- list(
  duration  = det_covs$duration[has_data, yr_idx, ],
  observers = det_covs$observers[has_data, yr_idx, ]
)

data_list_single <- list(
  y = y_sub,
  occ.covs = occ_covs_sub,
  det.covs = det_covs_sub,
  coords = coords_sub
)

# NNGP spatial model
# Phi prior: 3/max_dist to 3/min_dist covers spatial ranges from
# the study extent down to ~local neighbourhood
max_dist_sub <- max(dist(coords_sub[sample(nrow(coords_sub), min(300, nrow(coords_sub))), ]))
phi_low  <- 3 / max_dist_sub         # long-range correlation
phi_high <- 3 / (max_dist_sub * 0.01) # ~1% of extent

priors_spatial <- list(
  beta.normal  = list(mean = 0, var = 2.72),
  alpha.normal = list(mean = 0, var = 2.72),
  sigma.sq.ig  = c(2, 1),
  phi.unif     = c(phi_low, phi_high)
)

# Inits — beta: intercept + 6 occ covs; alpha: intercept + 2 det covs
inits_spatial <- list(
  beta = rep(0, ncol(occ_covs_sub) + 1),
  alpha = rep(0, length(det_covs_sub) + 1),
  sigma.sq = 1,
  phi = (phi_low + phi_high) / 2,
  z = ifelse(rowSums(y_sub, na.rm = TRUE) > 0, 1, 0)
)

n_batch_B  <- 400
n_burn_B   <- 5000
n_thin_B   <- 5
n_chains_B <- 3

cat("  Fitting spPGOcc with NNGP (n.neighbors=5)...\n")
cat("  MCMC:", n_batch_B * batch_len, "samples/chain x", n_chains_B, "chains\n")

fit_sppgocc <- spPGOcc(
  occ.formula  = ~ bio1 + bio2 + tree_cover + grass_cover + topo_elev + NDVI,
  det.formula  = ~ duration + observers,
  data         = data_list_single,
  priors       = priors_spatial,
  inits        = inits_spatial,
  cov.model    = "exponential",
  NNGP         = TRUE,
  n.neighbors  = 5,
  n.batch      = n_batch_B,
  batch.length = batch_len,
  n.burn       = n_burn_B,
  n.thin       = n_thin_B,
  n.chains     = n_chains_B,
  n.omp.threads = parallel::detectCores() - 1,
  verbose      = TRUE,
  n.report     = 100
)

t_end_B <- Sys.time()
time_B <- difftime(t_end_B, t_start_B, units = "mins")
cat("\n  Phase B time:", round(time_B, 1), "min\n\n")

cat("--- spPGOcc Summary ---\n")
summary(fit_sppgocc)

# Effective spatial range: distance at which correlation drops to ~5%
# For exponential covariance: effective range = 3 / phi
phi_post <- fit_sppgocc$beta.samples  # need to extract phi
if (!is.null(fit_sppgocc$theta.samples)) {
  # theta.samples contains sigma.sq and phi
  phi_samples <- fit_sppgocc$theta.samples[, "phi"]
  eff_range_km <- 3 / quantile(phi_samples, c(0.025, 0.5, 0.975)) / 1000
  cat("\nEffective spatial range (km):\n")
  cat("  2.5%:", round(eff_range_km[1], 1), "km\n")
  cat("  50%: ", round(eff_range_km[2], 1), "km\n")
  cat("  97.5%:", round(eff_range_km[3], 1), "km\n\n")
}

# WAIC
waic_B <- waicOcc(fit_sppgocc)
cat("spPGOcc WAIC (single year):", waic_B$WAIC, "\n\n")

saveRDS(fit_sppgocc, here("results", "otitar_spatial_spPGOcc_2021.rds"))


###############################################################################
# PHASE C: stPGOcc — FULL SPATIO-TEMPORAL MODEL
#
# This is the main model. Runs for hours.
# Comment out this section if you only want the proof of concept (A+B).
#
# Set RUN_PHASE_C <- TRUE to enable. Default is FALSE for safety.
###############################################################################

RUN_PHASE_C <- FALSE   # <-- Change to TRUE when ready

if (RUN_PHASE_C) {

  cat("====================================================================\n")
  cat(" Phase C: stPGOcc (full spatio-temporal, 6 years)\n")
  cat("====================================================================\n\n")

  t_start_C <- Sys.time()

  # Use phi prior informed by Phase B
  if (exists("phi_samples")) {
    phi_median <- median(phi_samples)
    phi_prior <- c(phi_median * 0.1, phi_median * 10)
  } else {
    phi_prior <- c(phi_low, phi_high)
  }

  priors_stpgocc <- list(
    beta.normal    = list(mean = 0, var = 2.72),
    alpha.normal   = list(mean = 0, var = 2.72),
    sigma.sq.ig    = c(2, 1),
    phi.unif       = phi_prior
  )

  # Production MCMC settings — increase for publication
  n_batch_C  <- 800      # 20,000 samples per chain
  n_burn_C   <- 10000
  n_thin_C   <- 20
  n_chains_C <- 3

  cat("  Fitting stPGOcc (NNGP, n.neighbors=5, ar1=TRUE)...\n")
  cat("  MCMC:", n_batch_C * batch_len, "samples/chain x", n_chains_C, "chains\n")
  cat("  This will take several hours. Go get coffee.\n\n")

  fit_stpgocc <- stPGOcc(
    occ.formula  = ~ bio1 + bio2 + tree_cover + grass_cover + topo_elev + NDVI + tmmx,
    det.formula  = ~ duration + observers,
    data         = data_list,
    priors       = priors_stpgocc,
    cov.model    = "exponential",
    NNGP         = TRUE,
    n.neighbors  = 5,
    n.batch      = n_batch_C,
    batch.length = batch_len,
    n.burn       = n_burn_C,
    n.thin       = n_thin_C,
    n.chains     = n_chains_C,
    ar1          = TRUE,
    n.omp.threads = parallel::detectCores() - 1,
    verbose      = TRUE,
    n.report     = 50
  )

  t_end_C <- Sys.time()
  time_C <- difftime(t_end_C, t_start_C, units = "hours")
  cat("\n  Phase C time:", round(time_C, 1), "hours\n\n")

  cat("--- stPGOcc Summary ---\n")
  summary(fit_stpgocc)

  # Spatial range
  if (!is.null(fit_stpgocc$theta.samples)) {
    phi_st <- fit_stpgocc$theta.samples[, "phi"]
    eff_range_st <- 3 / quantile(phi_st, c(0.025, 0.5, 0.975)) / 1000
    cat("\nEffective spatial range (km):\n")
    cat("  2.5%:", round(eff_range_st[1], 1), "km\n")
    cat("  50%: ", round(eff_range_st[2], 1), "km\n")
    cat("  97.5%:", round(eff_range_st[3], 1), "km\n")
  }

  # WAIC
  waic_C <- waicOcc(fit_stpgocc)
  cat("\nstPGOcc WAIC:", waic_C$WAIC, "\n")

  saveRDS(fit_stpgocc, here("results", "otitar_spatial_stPGOcc.rds"))

} else {
  cat("====================================================================\n")
  cat(" Phase C: SKIPPED (set RUN_PHASE_C <- TRUE to run)\n")
  cat("====================================================================\n\n")
}


###############################################################################
# DIAGNOSTICS
#
# Design rationale for Moran's I comparison:
#   The valid comparison is tPGOcc (Phase A) vs stPGOcc (Phase C), because
#   both use the SAME J=3052 sites over the SAME T=6 years. Comparing against
#   spPGOcc (Phase B) would be invalid: different site set (1218 vs 3052),
#   different temporal extent (1 year vs 6), and different neighbour matrix.
#
#   If Phase C was not run (RUN_PHASE_C = FALSE), we report only Phase A
#   Moran's I and flag that the spatial comparison requires Phase C.
#
# What we call "residuals" here:
#   resid_j = E[z_j | data] - y_naive_j
#   where E[z_j | data] is the posterior mean of the latent occupancy state
#   (averaged over years), and y_naive_j = 1 if species was ever detected at
#   site j across all years, 0 otherwise. This is a standard approximation
#   in occupancy modelling (see Doser et al. 2022, Ecological Monographs).
#   It is NOT a classical Pearson/deviance residual, but is the appropriate
#   metric for assessing spatial structure in latent-state models.
###############################################################################

cat("====================================================================\n")
cat(" Diagnostics\n")
cat("====================================================================\n\n")

sp_code <- "otitar"   # for output filenames

library(spdep)
library(ggplot2)


# =========================================================================
# 1. MCMC CONVERGENCE: Rhat and ESS
# =========================================================================

cat("--- 1. MCMC Convergence ---\n\n")

# Helper: extract Rhat and ESS from an spOccupancy summary object
extract_mcmc_diag <- function(fit, model_name) {
  s <- summary(fit)

  rows <- list()

  # Occupancy coefficients (beta)
  if (!is.null(s$beta.samples)) {
    b <- s$beta.samples
    for (i in seq_len(nrow(b))) {
      rows[[length(rows) + 1]] <- data.frame(
        Model = model_name, Parameter = paste0("beta[", rownames(b)[i], "]"),
        Mean = b[i, "Mean"], Rhat = b[i, "Rhat"], ESS = b[i, "ESS"],
        stringsAsFactors = FALSE
      )
    }
  }

  # Detection coefficients (alpha)
  if (!is.null(s$alpha.samples)) {
    a <- s$alpha.samples
    for (i in seq_len(nrow(a))) {
      rows[[length(rows) + 1]] <- data.frame(
        Model = model_name, Parameter = paste0("alpha[", rownames(a)[i], "]"),
        Mean = a[i, "Mean"], Rhat = a[i, "Rhat"], ESS = a[i, "ESS"],
        stringsAsFactors = FALSE
      )
    }
  }

  # Spatial / covariance parameters (theta)
  if (!is.null(s$theta.samples)) {
    th <- s$theta.samples
    for (i in seq_len(nrow(th))) {
      rows[[length(rows) + 1]] <- data.frame(
        Model = model_name, Parameter = paste0("theta[", rownames(th)[i], "]"),
        Mean = th[i, "Mean"], Rhat = th[i, "Rhat"], ESS = th[i, "ESS"],
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(rows) == 0) return(NULL)
  do.call(rbind, rows)
}

diag_frames <- list()

if (exists("fit_tpgocc")) {
  diag_A <- extract_mcmc_diag(fit_tpgocc, "tPGOcc")
  if (!is.null(diag_A)) {
    diag_frames[["A"]] <- diag_A
    cat("  tPGOcc: max Rhat =", round(max(diag_A$Rhat, na.rm = TRUE), 3),
        "| min ESS =", round(min(diag_A$ESS, na.rm = TRUE), 0), "\n")
    n_bad <- sum(diag_A$Rhat > 1.1, na.rm = TRUE)
    if (n_bad > 0) cat("  WARNING:", n_bad, "parameters with Rhat > 1.1\n")
  }
}

if (exists("fit_sppgocc")) {
  diag_B <- extract_mcmc_diag(fit_sppgocc, "spPGOcc")
  if (!is.null(diag_B)) {
    diag_frames[["B"]] <- diag_B
    cat("  spPGOcc: max Rhat =", round(max(diag_B$Rhat, na.rm = TRUE), 3),
        "| min ESS =", round(min(diag_B$ESS, na.rm = TRUE), 0), "\n")
  }
}

if (RUN_PHASE_C && exists("fit_stpgocc")) {
  diag_C <- extract_mcmc_diag(fit_stpgocc, "stPGOcc")
  if (!is.null(diag_C)) {
    diag_frames[["C"]] <- diag_C
    cat("  stPGOcc: max Rhat =", round(max(diag_C$Rhat, na.rm = TRUE), 3),
        "| min ESS =", round(min(diag_C$ESS, na.rm = TRUE), 0), "\n")
  }
}

# Save combined MCMC diagnostics
if (length(diag_frames) > 0) {
  mcmc_diag_df <- do.call(rbind, diag_frames)
  mcmc_path <- here("results", paste0(sp_code, "_mcmc_diagnostics.csv"))
  write.csv(mcmc_diag_df, mcmc_path, row.names = FALSE)
  cat("\n  Saved:", mcmc_path, "\n")
}
cat("\n")


# =========================================================================
# 2. MORAN'S I ON RESIDUALS
#    CRITICAL: both models compared on the SAME sites / SAME neighbour matrix.
#    tPGOcc vs stPGOcc (both J=3052, T=6). NOT vs spPGOcc (different J/T).
# =========================================================================

cat("--- 2. Moran's I on residuals ---\n\n")

# Build a SINGLE neighbour matrix for the full dataset (J=3052).
# This is reused for both tPGOcc and stPGOcc to ensure a fair comparison.
knn_full <- knearneigh(coords, k = 8)
nb_full  <- knn2nb(knn_full)
lw_full  <- nb2listw(nb_full, style = "W")

# Naive detection per site (ever detected across all years)
y_naive <- apply(y_array, 1, function(x) as.numeric(any(x == 1, na.rm = TRUE)))

# --- tPGOcc residuals (non-spatial baseline) ---
moran_A <- NULL
resid_A <- NULL

if (exists("fit_tpgocc")) {
  # z.samples: [n_samples x J x T]
  # Posterior mean of z per site, averaged over years
  z_mean_A <- apply(fit_tpgocc$z.samples, c(2, 3), mean)  # [J x T]
  psi_mean_A <- rowMeans(z_mean_A, na.rm = TRUE)           # [J]

  resid_A <- psi_mean_A - y_naive

  moran_A <- moran.test(resid_A, lw_full, zero.policy = TRUE)
  cat("  tPGOcc (non-spatial, J=", J, ", T=", n_years, "):\n", sep = "")
  cat("    Moran's I =", round(moran_A$estimate["Moran I statistic"], 4),
      "| p =", signif(moran_A$p.value, 3), "\n\n")
}

# --- stPGOcc residuals (spatio-temporal model) ---
moran_C <- NULL
resid_C <- NULL
w_spatial <- NULL

if (RUN_PHASE_C && exists("fit_stpgocc")) {
  z_mean_C <- apply(fit_stpgocc$z.samples, c(2, 3), mean)  # [J x T]
  psi_mean_C <- rowMeans(z_mean_C, na.rm = TRUE)            # [J]

  resid_C <- psi_mean_C - y_naive

  moran_C <- moran.test(resid_C, lw_full, zero.policy = TRUE)
  cat("  stPGOcc (spatial, J=", J, ", T=", n_years, "):\n", sep = "")
  cat("    Moran's I =", round(moran_C$estimate["Moran I statistic"], 4),
      "| p =", signif(moran_C$p.value, 3), "\n\n")

  # Spatial random effect (w) — posterior mean per site
  if (!is.null(fit_stpgocc$w.samples)) {
    w_spatial <- apply(fit_stpgocc$w.samples, 2, mean)  # [J]
  }
} else if (!RUN_PHASE_C) {
  cat("  stPGOcc: NOT RUN (Phase C disabled).\n")
  cat("  -> Set RUN_PHASE_C <- TRUE to enable the spatial comparison.\n")
  cat("  -> Without Phase C, Moran's I comparison is not possible.\n\n")
}

# Reduction summary
if (!is.null(moran_A) && !is.null(moran_C)) {
  I_A <- moran_A$estimate["Moran I statistic"]
  I_C <- moran_C$estimate["Moran I statistic"]
  reduction_pct <- (I_A - I_C) / I_A * 100
  cat("  --> Moran's I reduction: ", round(I_A, 4), " -> ", round(I_C, 4),
      " (", round(reduction_pct, 1), "%)\n\n", sep = "")
}


# =========================================================================
# 3. SPATIAL RANGE (from Phase C if available, else Phase B)
# =========================================================================

cat("--- 3. Spatial range ---\n\n")

eff_range_median <- NA
eff_range_lo     <- NA
eff_range_hi     <- NA

if (RUN_PHASE_C && exists("fit_stpgocc") && !is.null(fit_stpgocc$theta.samples)) {
  phi_post <- fit_stpgocc$theta.samples[, "phi"]
  # Exponential covariance: effective range (correlation ~ 0.05) = 3 / phi
  eff_range <- 3 / phi_post / 1000  # convert m to km
  eff_range_lo     <- quantile(eff_range, 0.025)
  eff_range_median <- quantile(eff_range, 0.500)
  eff_range_hi     <- quantile(eff_range, 0.975)
  cat("  stPGOcc effective spatial range:\n")
  cat("    Median:", round(eff_range_median, 1), "km\n")
  cat("    95% CI: [", round(eff_range_lo, 1), ",", round(eff_range_hi, 1), "] km\n\n")
} else if (exists("fit_sppgocc") && !is.null(fit_sppgocc$theta.samples)) {
  phi_post <- fit_sppgocc$theta.samples[, "phi"]
  eff_range <- 3 / phi_post / 1000
  eff_range_lo     <- quantile(eff_range, 0.025)
  eff_range_median <- quantile(eff_range, 0.500)
  eff_range_hi     <- quantile(eff_range, 0.975)
  cat("  spPGOcc effective spatial range (Phase B only, single year):\n")
  cat("    Median:", round(eff_range_median, 1), "km\n")
  cat("    95% CI: [", round(eff_range_lo, 1), ",", round(eff_range_hi, 1), "] km\n")
  cat("    (Note: from single-year model. Run Phase C for multi-year estimate.)\n\n")
} else {
  cat("  Spatial range: not available (no spatial model fitted).\n\n")
}


# =========================================================================
# 4. MODEL COMPARISON TABLE -> results/{sp}_model_comparison.csv
# =========================================================================

cat("--- 4. Model comparison table ---\n\n")

comp_rows <- list()

# tPGOcc
if (exists("waic_A")) {
  comp_rows[["tPGOcc"]] <- data.frame(
    Model = "tPGOcc",
    WAIC = round(waic_A$WAIC, 1),
    Moran_I = if (!is.null(moran_A)) round(moran_A$estimate["Moran I statistic"], 4) else NA,
    Moran_p = if (!is.null(moran_A)) signif(moran_A$p.value, 3) else NA,
    Spatial_range_km = NA,
    Runtime_min = round(as.numeric(difftime(t_end_A, t_start_A, units = "mins")), 1),
    stringsAsFactors = FALSE
  )
}

# spPGOcc (Phase B — informational only, not directly comparable)
if (exists("waic_B")) {
  range_B <- NA
  if (exists("fit_sppgocc") && !is.null(fit_sppgocc$theta.samples)) {
    range_B <- round(median(3 / fit_sppgocc$theta.samples[, "phi"] / 1000), 1)
  }
  comp_rows[["spPGOcc"]] <- data.frame(
    Model = "spPGOcc (2021 only)",
    WAIC = round(waic_B$WAIC, 1),
    Moran_I = NA,   # not comparable (different site set)
    Moran_p = NA,
    Spatial_range_km = range_B,
    Runtime_min = round(as.numeric(difftime(t_end_B, t_start_B, units = "mins")), 1),
    stringsAsFactors = FALSE
  )
}

# stPGOcc (Phase C)
if (RUN_PHASE_C && exists("waic_C")) {
  comp_rows[["stPGOcc"]] <- data.frame(
    Model = "stPGOcc",
    WAIC = round(waic_C$WAIC, 1),
    Moran_I = if (!is.null(moran_C)) round(moran_C$estimate["Moran I statistic"], 4) else NA,
    Moran_p = if (!is.null(moran_C)) signif(moran_C$p.value, 3) else NA,
    Spatial_range_km = round(eff_range_median, 1),
    Runtime_min = round(as.numeric(difftime(t_end_C, t_start_C, units = "mins")), 1),
    stringsAsFactors = FALSE
  )
}

if (length(comp_rows) > 0) {
  comp_df <- do.call(rbind, comp_rows)
  rownames(comp_df) <- NULL
  comp_path <- here("results", paste0(sp_code, "_model_comparison.csv"))
  write.csv(comp_df, comp_path, row.names = FALSE)
  cat("  Model comparison table:\n")
  print(comp_df, row.names = FALSE)
  cat("\n  Saved:", comp_path, "\n")

  # WAIC difference
  if (exists("waic_A") && RUN_PHASE_C && exists("waic_C")) {
    delta_waic <- waic_A$WAIC - waic_C$WAIC
    cat("\n  DELTA WAIC (tPGOcc - stPGOcc) =", round(delta_waic, 1))
    if (delta_waic > 10) {
      cat("  -> STRONG support for spatial model\n")
    } else if (delta_waic > 2) {
      cat("  -> Moderate support for spatial model\n")
    } else {
      cat("  -> Weak or no support for spatial model\n")
    }
  }
}
cat("\n")


# =========================================================================
# 5. DIAGNOSTIC FIGURE -> figs/{sp}_spatial_diagnostics.png
#    Panel 1: Non-spatial residuals map (tPGOcc)
#    Panel 2: Spatial random effect w (stPGOcc), or Phase B w if C not run
# =========================================================================

cat("--- 5. Diagnostic figure ---\n\n")

# We need at least tPGOcc residuals to make the figure
if (!is.null(resid_A)) {

  # Data frame with coordinates + residuals for plotting
  diag_plot_df <- data.frame(
    lon = coords_lonlat[, 1],
    lat = coords_lonlat[, 2],
    resid_nonspatial = resid_A
  )

  # Panel 1: non-spatial residuals
  p1 <- ggplot(diag_plot_df, aes(x = lon, y = lat, colour = resid_nonspatial)) +
    geom_point(size = 0.4, alpha = 0.7) +
    scale_colour_gradient2(low = "blue", mid = "white", high = "red",
                           midpoint = 0, name = "Residual") +
    labs(title = "tPGOcc residuals (non-spatial)",
         subtitle = paste0("Moran's I = ",
                           if (!is.null(moran_A))
                             round(moran_A$estimate["Moran I statistic"], 3)
                           else "NA")) +
    coord_quickmap() + theme_minimal() +
    theme(legend.position = "bottom")

  # Panel 2: spatial random effect w (stPGOcc preferred, else spPGOcc)
  p2 <- NULL
  if (!is.null(w_spatial) && length(w_spatial) == J) {
    # stPGOcc w (full multi-season, same J sites)
    diag_plot_df$w <- w_spatial
    p2 <- ggplot(diag_plot_df, aes(x = lon, y = lat, colour = w)) +
      geom_point(size = 0.4, alpha = 0.7) +
      scale_colour_gradient2(low = "blue", mid = "white", high = "red",
                             midpoint = 0, name = "w") +
      labs(title = "stPGOcc spatial random effect (w)",
           subtitle = paste0("Effective range: ",
                             round(eff_range_median, 0), " km [",
                             round(eff_range_lo, 0), "-",
                             round(eff_range_hi, 0), "]")) +
      coord_quickmap() + theme_minimal() +
      theme(legend.position = "bottom")

  } else if (exists("fit_sppgocc") && !is.null(fit_sppgocc$w.samples)) {
    # Fallback: spPGOcc w (Phase B, subset of sites)
    w_B <- apply(fit_sppgocc$w.samples, 2, mean)
    w_df <- data.frame(
      lon = coords_lonlat[has_data, 1],
      lat = coords_lonlat[has_data, 2],
      w = w_B
    )
    p2 <- ggplot(w_df, aes(x = lon, y = lat, colour = w)) +
      geom_point(size = 0.5, alpha = 0.7) +
      scale_colour_gradient2(low = "blue", mid = "white", high = "red",
                             midpoint = 0, name = "w") +
      labs(title = "spPGOcc spatial effect (2021 only)",
           subtitle = paste0("N = ", sum(has_data),
                             " sites | Range: ", round(eff_range_median, 0), " km")) +
      coord_quickmap() + theme_minimal() +
      theme(legend.position = "bottom")

  } else {
    # No spatial model available — second panel shows stPGOcc residuals or placeholder
    if (!is.null(resid_C)) {
      diag_plot_df$resid_spatial <- resid_C
      p2 <- ggplot(diag_plot_df, aes(x = lon, y = lat, colour = resid_spatial)) +
        geom_point(size = 0.4, alpha = 0.7) +
        scale_colour_gradient2(low = "blue", mid = "white", high = "red",
                               midpoint = 0, name = "Residual") +
        labs(title = "stPGOcc residuals (spatial)",
             subtitle = paste0("Moran's I = ",
                               if (!is.null(moran_C))
                                 round(moran_C$estimate["Moran I statistic"], 3)
                               else "NA")) +
        coord_quickmap() + theme_minimal() +
        theme(legend.position = "bottom")
    } else {
      p2 <- ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "Phase C not run\nSet RUN_PHASE_C <- TRUE",
                 size = 5, hjust = 0.5) +
        theme_void() +
        labs(title = "Spatial effect not available")
    }
  }

  fig_path <- here("figs", paste0(sp_code, "_spatial_diagnostics.png"))
  ggsave(fig_path,
         gridExtra::grid.arrange(p1, p2, ncol = 2),
         width = 14, height = 6, dpi = 150)
  cat("  Saved:", fig_path, "\n\n")

} else {
  cat("  Cannot generate figure: tPGOcc residuals not available.\n\n")
}


# =========================================================================
# 6. PAPER-READY SUMMARY
# =========================================================================

t_end_global <- Sys.time()
total_time <- difftime(t_end_global, t_start_global, units = "mins")

cat("====================================================================\n")
cat(" PAPER-READY SUMMARY\n")
cat("====================================================================\n\n")

cat("Species:", sp_code, "\n")
cat("Sites:", J, "| Years:", n_years, "| Replicates:", n_reps, "\n")
cat("Total elapsed time:", round(total_time, 1), "min\n\n")

cat("--- Timing ---\n")
if (exists("t_end_A")) cat("  tPGOcc  (non-spatial): ", round(difftime(t_end_A, t_start_A, units = "mins"), 1), "min\n")
if (exists("t_end_B")) cat("  spPGOcc (spatial, 2021): ", round(difftime(t_end_B, t_start_B, units = "mins"), 1), "min\n")
if (RUN_PHASE_C && exists("t_end_C")) cat("  stPGOcc (spatio-temporal): ", round(as.numeric(time_C) * 60, 1), "min (", round(as.numeric(time_C), 1), "h)\n")

cat("\n--- Model selection (WAIC) ---\n")
if (exists("waic_A")) cat("  tPGOcc WAIC: ", round(waic_A$WAIC, 1), "\n")
if (RUN_PHASE_C && exists("waic_C")) {
  cat("  stPGOcc WAIC:", round(waic_C$WAIC, 1), "\n")
  delta <- waic_A$WAIC - waic_C$WAIC
  cat("  DELTA WAIC:  ", round(delta, 1),
      if (delta > 10) " (strong)" else if (delta > 2) " (moderate)" else " (weak)", "\n")
} else {
  cat("  stPGOcc WAIC: not available (Phase C not run)\n")
}

cat("\n--- Spatial autocorrelation (Moran's I) ---\n")
if (!is.null(moran_A)) {
  cat("  tPGOcc  (no spatial effect): I =", round(moran_A$estimate["Moran I statistic"], 4),
      ", p =", signif(moran_A$p.value, 3), "\n")
}
if (!is.null(moran_C)) {
  cat("  stPGOcc (with spatial effect): I =", round(moran_C$estimate["Moran I statistic"], 4),
      ", p =", signif(moran_C$p.value, 3), "\n")
  cat("  Reduction:", round(reduction_pct, 1), "%\n")
} else {
  cat("  stPGOcc: Phase C not run — comparison pending.\n")
}

cat("\n--- Spatial range ---\n")
if (!is.na(eff_range_median)) {
  cat("  Effective range:", round(eff_range_median, 1),
      "km [95% CI:", round(eff_range_lo, 1), "-", round(eff_range_hi, 1), "km]\n")
  if (eff_range_median > 50) {
    cat("  Interpretation: broad-scale spatial structure, consistent with\n")
    cat("  biogeographic / macroecological processes.\n")
  } else if (eff_range_median > 20) {
    cat("  Interpretation: mesoscale spatial structure, consistent with\n")
    cat("  landscape-level habitat connectivity and metapopulation dynamics.\n")
  } else {
    cat("  Interpretation: fine-scale spatial structure, consistent with\n")
    cat("  local dispersal limitation and conspecific attraction.\n")
  }
} else {
  cat("  Not available.\n")
}

cat("\n--- Output files ---\n")
cat("  results/", sp_code, "_mcmc_diagnostics.csv\n", sep = "")
cat("  results/", sp_code, "_model_comparison.csv\n", sep = "")
cat("  figs/", sp_code, "_spatial_diagnostics.png\n", sep = "")
cat("  results/", sp_code, "_spatial_tPGOcc.rds\n", sep = "")
cat("  results/", sp_code, "_spatial_spPGOcc_2021.rds\n", sep = "")
if (RUN_PHASE_C) cat("  results/", sp_code, "_spatial_stPGOcc.rds\n", sep = "")

cat("\n")
cat("====================================================================\n")
cat(" DECISION FRAMEWORK: Scale to other species?\n")
cat("====================================================================\n")
cat("
Evaluate based on the results above:

1. SPATIAL RANGE:
   > 50 km  -> biogeographic-scale structure (broad patterns)
   20-50 km -> landscape-scale connectivity (metapopulations)
   < 20 km  -> local dispersal limitation (fine-scale)

2. MORAN'S I REDUCTION:
   > 50% reduction -> spatial effect highly effective; include in paper
   20-50%          -> moderate; report as supplementary analysis
   < 20%           -> spatial structure minimal; stick with colext

3. WAIC IMPROVEMENT (tPGOcc vs stPGOcc):
   delta > 10  -> strong support for spatial model
   delta 2-10  -> moderate support
   delta < 2   -> no meaningful improvement

4. COST-BENEFIT:
   Phase C < 12h -> viable for all 4 species
   Phase C 12-48h -> consider subsampling to ~1000 sites
   Phase C > 48h  -> not viable at full scale

ALTERNATIVES if stPGOcc is too slow:
  a) Subsample to ~1000 sites (stratified by region)
  b) tPGOcc with ar1=TRUE captures temporal but not spatial autocorrelation
  c) Add spatial eigenvectors (MEMs) as covariates in colext() —
     computationally free approximation to spatial structure
  d) Trend surface: add poly(x, 2) + poly(y, 2) to colext psi formula —
     captures broad spatial gradients without GP overhead
")
