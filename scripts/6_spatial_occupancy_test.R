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
###############################################################################

cat("====================================================================\n")
cat(" Diagnostics\n")
cat("====================================================================\n\n")

# --- 1. Convergence: Rhat and ESS ---
cat("--- Convergence (Phase A: tPGOcc) ---\n")
if (exists("fit_tpgocc")) {
  # Check Rhat from summary
  s <- summary(fit_tpgocc)
  cat("  (See summary above for Rhat and ESS columns)\n\n")
}

# --- 2. Moran's I on spatial residuals ---
cat("--- Moran's I on residuals ---\n")

# For the non-spatial model (tPGOcc), extract posterior mean of z
# and compute residuals against detection history
if (exists("fit_tpgocc")) {
  library(spdep)

  # Posterior mean occupancy probability per site (averaged over years)
  z_mean <- apply(fit_tpgocc$z.samples, c(2, 3), mean)  # [J x T]
  psi_mean <- rowMeans(z_mean, na.rm = TRUE)  # average over years

  # Naive detection per site
  y_detected <- apply(y_array, 1, function(x) as.numeric(any(x == 1, na.rm = TRUE)))

  resid_nonspatial <- psi_mean - y_detected

  # Build neighbours (k-nearest for point data)
  knn <- knearneigh(coords, k = 8)
  nb <- knn2nb(knn)
  lw <- nb2listw(nb, style = "W")

  moran_nonspatial <- moran.test(resid_nonspatial, lw, zero.policy = TRUE)
  cat("\n  tPGOcc (non-spatial) residual Moran's I:\n")
  print(moran_nonspatial)
  cat("  Interpretation: Moran's I > 0 with p < 0.05 means spatial autocorrelation\n")
  cat("  remains unaccounted for.\n\n")
}

# For the spatial model (spPGOcc, single year), if available
if (exists("fit_sppgocc")) {
  z_sp <- apply(fit_sppgocc$z.samples, 2, mean)
  y_det_sub <- apply(y_sub, 1, function(x) as.numeric(any(x == 1, na.rm = TRUE)))
  resid_spatial <- z_sp - y_det_sub

  knn_sub <- knearneigh(coords_sub, k = 8)
  nb_sub <- knn2nb(knn_sub)
  lw_sub <- nb2listw(nb_sub, style = "W")

  moran_spatial <- moran.test(resid_spatial, lw_sub, zero.policy = TRUE)
  cat("  spPGOcc (spatial, 2021) residual Moran's I:\n")
  print(moran_spatial)
  cat("\n")
}

# --- 3. Model comparison: WAIC ---
cat("--- WAIC Comparison ---\n")
cat("  tPGOcc  (non-spatial, multi-season): ", if (exists("waic_A")) round(waic_A$WAIC, 1) else "not run", "\n")
cat("  spPGOcc (spatial, single-year 2021): ", if (exists("waic_B")) round(waic_B$WAIC, 1) else "not run", "\n")
if (RUN_PHASE_C && exists("waic_C")) {
  cat("  stPGOcc (spatial, multi-season):     ", round(waic_C$WAIC, 1), "\n")
}
cat("\n  Lower WAIC = better fit. WAIC difference > 10 is strong.\n\n")


# --- 4. AUC/TSS comparison with colext ---
# Load colext predictions for comparison (if available)
pred_path <- here("data-raw", "data", "otitar", "occ_otitar_prediction.csv")
atlas_path <- here("data", "raw", "validation", "atlas_biodiversidad", "aves_spain.shp")

if (file.exists(pred_path) && file.exists(atlas_path)) {
  cat("--- AUC/TSS vs colext (atlas validation) ---\n")

  library(terra)
  library(pROC)
  library(Metrics)

  atlas <- st_read(atlas_path, quiet = TRUE)

  # colext predictions
  pred_colext <- read.csv(pred_path)
  names(pred_colext) <- tolower(names(pred_colext))
  r_colext <- rast(data.frame(x = pred_colext$longitude, y = pred_colext$latitude,
                               occ_prob = pred_colext$occ_prob), type = "xyz")
  crs(r_colext) <- "+proj=longlat"
  r_colext_proj <- project(r_colext, crs(atlas))
  ext_colext <- terra::extract(r_colext_proj, atlas, fun = median, na.rm = TRUE)
  atlas$pred_colext <- ext_colext$occ_prob

  # spPGOcc predictions (single year, only for sites with data)
  # Map posterior z back to spatial locations and rasterize
  if (exists("fit_sppgocc")) {
    z_sp_mean <- apply(fit_sppgocc$z.samples, 2, mean)
    pred_sp_df <- data.frame(
      x = coords_lonlat[has_data, 1],
      y = coords_lonlat[has_data, 2],
      occ_prob = z_sp_mean
    )
    r_sp <- rast(pred_sp_df, type = "xyz")
    crs(r_sp) <- "+proj=longlat"
    r_sp_proj <- project(r_sp, crs(atlas))
    ext_sp <- terra::extract(r_sp_proj, atlas, fun = median, na.rm = TRUE)
    atlas$pred_spatial <- ext_sp$occ_prob
  }

  # Compute metrics on valid polygons
  atlas_valid <- atlas %>%
    filter(!is.na(pred_colext))

  if (nrow(atlas_valid) > 50) {
    # colext AUC
    roc_colext <- pROC::roc(atlas_valid$OTITAR, atlas_valid$pred_colext, quiet = TRUE)
    auc_colext <- pROC::auc(roc_colext)
    th_colext <- pROC::coords(roc_colext, "best", ret = "threshold") %>% as.numeric()
    bin_colext <- ifelse(atlas_valid$pred_colext >= th_colext, 1, 0)
    cm_colext <- table(factor(atlas_valid$OTITAR, 0:1), factor(bin_colext, 0:1))
    TSS_colext <- cm_colext["1","1"]/(cm_colext["1","1"]+cm_colext["1","0"]) +
                  cm_colext["0","0"]/(cm_colext["0","0"]+cm_colext["0","1"]) - 1

    cat("  colext:   AUC =", round(auc_colext, 3), "| TSS =", round(TSS_colext, 3), "\n")

    # spPGOcc AUC (if available)
    if ("pred_spatial" %in% names(atlas_valid)) {
      av_sp <- atlas_valid %>% filter(!is.na(pred_spatial))
      if (nrow(av_sp) > 50) {
        roc_sp <- pROC::roc(av_sp$OTITAR, av_sp$pred_spatial, quiet = TRUE)
        auc_sp <- pROC::auc(roc_sp)
        th_sp <- pROC::coords(roc_sp, "best", ret = "threshold") %>% as.numeric()
        bin_sp <- ifelse(av_sp$pred_spatial >= th_sp, 1, 0)
        cm_sp <- table(factor(av_sp$OTITAR, 0:1), factor(bin_sp, 0:1))
        TSS_sp <- cm_sp["1","1"]/(cm_sp["1","1"]+cm_sp["1","0"]) +
                  cm_sp["0","0"]/(cm_sp["0","0"]+cm_sp["0","1"]) - 1
        cat("  spPGOcc: AUC =", round(auc_sp, 3), "| TSS =", round(TSS_sp, 3), "\n")
      }
    }
  }
  cat("\n")
} else {
  cat("  (Prediction/atlas files not available for AUC/TSS comparison)\n\n")
}


###############################################################################
# SUMMARY AND RECOMMENDATION
###############################################################################

t_end_global <- Sys.time()
total_time <- difftime(t_end_global, t_start_global, units = "mins")

cat("====================================================================\n")
cat(" SUMMARY\n")
cat("====================================================================\n\n")

cat("Total elapsed time:", round(total_time, 1), "min\n\n")

cat("Timing per phase:\n")
if (exists("t_end_A")) cat("  Phase A (tPGOcc, non-spatial):  ", round(difftime(t_end_A, t_start_A, units="mins"), 1), "min\n")
if (exists("t_end_B")) cat("  Phase B (spPGOcc, spatial 2021):", round(difftime(t_end_B, t_start_B, units="mins"), 1), "min\n")
if (RUN_PHASE_C && exists("t_end_C")) cat("  Phase C (stPGOcc, full ST):     ", round(as.numeric(time_C), 1), "hours\n")

cat("\nModel comparison:\n")
if (exists("waic_A")) cat("  tPGOcc WAIC:  ", round(waic_A$WAIC, 1), "\n")
if (exists("waic_B")) cat("  spPGOcc WAIC: ", round(waic_B$WAIC, 1), "(single year, not directly comparable)\n")
if (RUN_PHASE_C && exists("waic_C")) cat("  stPGOcc WAIC: ", round(waic_C$WAIC, 1), "\n")

cat("\nMoran's I comparison:\n")
if (exists("moran_nonspatial")) cat("  Non-spatial residuals: I =", round(moran_nonspatial$estimate[1], 3), "p =", signif(moran_nonspatial$p.value, 3), "\n")
if (exists("moran_spatial"))    cat("  Spatial residuals:     I =", round(moran_spatial$estimate[1], 3), "p =", signif(moran_spatial$p.value, 3), "\n")

cat("\n")
cat("====================================================================\n")
cat(" DECISION FRAMEWORK: Is it worth scaling to other species?\n")
cat("====================================================================\n")
cat("
Questions to answer from the results above:

1. SPATIAL RANGE: If the effective range is > 50 km, the spatial random
   effect captures broad-scale patterns (biogeography). If < 20 km, it
   captures local clustering (dispersal). Both are ecologically meaningful.

2. MORAN'S I REDUCTION: If the spatial model reduces Moran's I from
   ~0.38 to < 0.10, the spatial effect substantially improves residual
   independence. This justifies the computational cost.

3. WAIC IMPROVEMENT: If stPGOcc WAIC is >10 lower than tPGOcc, the spatial
   structure provides meaningful predictive improvement.

4. AUC/TSS: If AUC improves by > 0.02 or TSS by > 0.05, there is practical
   predictive gain for conservation mapping.

5. COST-BENEFIT: If Phase B (single-year) takes < 2 hours and reduces
   Moran's I, it is worth running Phase C. If Phase C takes < 12 hours
   on your hardware, it is viable for all 4 species.

ALTERNATIVES if stPGOcc is too slow:
  a) Subsample to ~1000 sites (stratified by region) — cuts time ~9x
  b) Use tPGOcc with ar1=TRUE (already done in Phase A) — captures
     temporal autocorrelation without spatial structure
  c) Use a CAR (Conditional AutoRegressive) model via INLA but with
     simplified detection — faster but loses observation-level covariates
  d) Add geographic coordinates as covariates in colext() (trend surface) —
     crude but computationally free. Try: psi ~ bio1 + ... + poly(x, 2) + poly(y, 2)
  e) Use spatial eigenvectors (MEMs) from spdep::ME() as covariates in
     colext() — captures spatial patterns without full GP overhead
")
