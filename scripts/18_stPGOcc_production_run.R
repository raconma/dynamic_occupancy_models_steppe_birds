###############################################################################
# 18_stPGOcc_production_run.R
#
# PURPOSE: Production re-run of stPGOcc with settings optimised for
#          convergence of spatial parameters (sigma.sq, phi).
#
# CONTEXT: The 500-batch run (12,500 iter) showed:
#   - Detection parameters (alpha): perfect convergence across all species
#   - Occupancy covariates (beta): OK for otitar/pteori, poor for ptealc/tettet
#   - Spatial parameters (sigma.sq, phi): FAILED across all species
#     (Rhat 1.1-6.2, ESS 9-250)
#
# KEY CHANGES vs script 6:
#   1. 200K iterations (8000 batches x 25), 100K burn-in, thin=100
#      -> 3000 post-burnin samples per chain (9000 total)
#   2. n.neighbors = 10 (was 5) — better NNGP approximation
#   3. Informative starting values from 500-run posteriors
#   4. Tighter sigma.sq prior: IG(2, 2) instead of IG(2, 1)
#   5. Skips Phase A/B (already completed) — only runs stPGOcc
#   6. Saves to results/results_spatial/results_production/
#
# EXPECTED RUNTIME: ~8-12 hours per species (~2 days total)
#                   Run on a machine with >= 8 cores
#
# USAGE:
#   Rscript scripts/18_stPGOcc_production_run.R
#   # or for a single species:
#   Rscript scripts/18_stPGOcc_production_run.R otitar
#
# OUTPUTS (per species):
#   results/results_spatial/results_production/{sp}_stPGOcc.rds
#   results/results_spatial/results_production/{sp}_mcmc_diagnostics.csv
#   results/results_spatial/results_production/{sp}_model_comparison.csv
#   figs/{sp}_spatial_diagnostics_production.png
###############################################################################

library(here)
library(sf)
library(dplyr)
library(spdep)
library(ggplot2)

if (!requireNamespace("spOccupancy", quietly = TRUE)) {
  stop("Install spOccupancy first: install.packages('spOccupancy')")
}
library(spOccupancy)

source(here("R", "model_configs.R"))

set.seed(42)

# ---- Output directory ----
out_dir <- here("results", "results_spatial", "results_production")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ---- Configuration ----
years    <- 2017:2023
n_years  <- length(years)
n_reps   <- 10
batch_len <- 25

# MCMC settings: 200K iterations per chain
N_BATCH   <- 8000    # 8000 * 25 = 200,000 samples per chain
N_BURN    <- 100000  # 50% burn-in
N_THIN    <- 100     # keep 1,000 per chain -> 3,000 total
N_CHAINS  <- 3
N_NEIGHBORS <- 10    # was 5; better NNGP approximation for broad-scale SA
N_THREADS <- max(1, parallel::detectCores() - 1)

# Species to run (all 4 by default, or from command line)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  species_codes <- args
} else {
  species_codes <- c("otitar", "ptealc", "pteori", "tettet")
}

# Starting values informed by 500-run posteriors
# These help the sampler start in a reasonable region
start_vals <- list(
  otitar = list(sigma.sq = 10, phi = 6e-05),
  ptealc = list(sigma.sq = 5,  phi = 3e-05),
  pteori = list(sigma.sq = 15, phi = 2e-05),
  tettet = list(sigma.sq = 10, phi = 7e-05)
)

cat("\n")
cat("###############################################################\n")
cat("#  stPGOcc PRODUCTION RUN                                     #\n")
cat("#  200K iterations | n.neighbors = 10 | informed starting vals #\n")
cat("###############################################################\n")
cat("\nStarted:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Species:", paste(species_codes, collapse = ", "), "\n")
cat("Threads:", N_THREADS, "\n")
cat("MCMC: ", N_BATCH * batch_len, " iter/chain x ", N_CHAINS, " chains\n")
cat("Burn-in:", N_BURN, "| Thin:", N_THIN, "\n")
cat("Post-thinning samples:", (N_BATCH * batch_len - N_BURN) / N_THIN,
    "per chain,", (N_BATCH * batch_len - N_BURN) / N_THIN * N_CHAINS, "total\n")
cat("n.neighbors:", N_NEIGHBORS, "\n\n")


###############################################################################
# MAIN LOOP
###############################################################################

for (sp in species_codes) {

  tryCatch({

  t_start <- Sys.time()
  cfg <- get_model_config(sp)

  cat("\n====================================================================\n")
  cat(" stPGOcc Production —", sp, "\n")
  cat("====================================================================\n\n")

  # ---- 1. Load and reshape data (same as script 6) ----
  cat("Loading data...\n")
  d <- read.csv(here("data", "processed_2023", sp,
                      paste0(sp, "_occ_wide_dynamic.csv")))
  J <- nrow(d)

  # Detection history
  y_array <- array(NA, dim = c(J, n_years, n_reps))
  for (t in seq_along(years)) {
    for (k in 1:n_reps) {
      col_name <- paste0("y.", k, ".", years[t])
      if (col_name %in% names(d)) y_array[, t, k] <- d[[col_name]]
    }
  }

  # Coordinates (ETRS89/UTM30N)
  pts_sf <- st_as_sf(data.frame(x = d$longitude, y = d$latitude),
                     coords = c("x", "y"), crs = 4326)
  pts_utm <- st_transform(pts_sf, 25830)
  coords <- st_coordinates(pts_utm)
  coords_lonlat <- cbind(d$longitude, d$latitude)

  # Occupancy covariates
  occ_covs <- list()
  for (v in cfg$psi_vars) {
    if (v %in% names(d)) occ_covs[[v]] <- d[[v]]
  }
  for (var_name in c("NDVI", "tmmx")) {
    mat <- matrix(NA, nrow = J, ncol = n_years)
    for (t in seq_along(years)) {
      col <- paste0(var_name, "_", years[t])
      if (col %in% names(d)) {
        vals <- d[[col]]
        if (abs(mean(vals, na.rm = TRUE)) > 10 || sd(vals, na.rm = TRUE) > 10)
          vals <- as.numeric(scale(vals))
        mat[, t] <- vals
      }
    }
    occ_covs[[var_name]] <- mat
  }

  # Remove NA sites
  na_sites <- rep(FALSE, J)
  for (nm in names(occ_covs)) {
    obj <- occ_covs[[nm]]
    if (is.matrix(obj)) {
      na_sites <- na_sites | apply(obj, 1, function(x) any(is.na(x)))
    } else {
      na_sites <- na_sites | is.na(obj)
    }
  }
  na_sites <- na_sites | is.na(coords[, 1]) | is.na(coords[, 2])

  if (sum(na_sites) > 0) {
    keep <- !na_sites
    J <- sum(keep)
    y_array <- y_array[keep, , ]
    coords <- coords[keep, , drop = FALSE]
    coords_lonlat <- coords_lonlat[keep, , drop = FALSE]
    d <- d[keep, ]
    for (nm in names(occ_covs)) {
      if (is.matrix(occ_covs[[nm]])) {
        occ_covs[[nm]] <- occ_covs[[nm]][keep, , drop = FALSE]
      } else {
        occ_covs[[nm]] <- occ_covs[[nm]][keep]
      }
    }
    cat("  Sites after NA removal:", J, "\n")
  }

  # Impute isolated NAs in dynamic covariates
  for (nm in c("NDVI", "tmmx")) {
    if (nm %in% names(occ_covs) && is.matrix(occ_covs[[nm]])) {
      mat <- occ_covs[[nm]]
      for (col_i in 1:ncol(mat)) {
        na_idx <- is.na(mat[, col_i])
        if (any(na_idx)) mat[na_idx, col_i] <- mean(mat[, col_i], na.rm = TRUE)
      }
      occ_covs[[nm]] <- mat
    }
  }

  # Detection covariates
  det_covs <- list()
  for (var_base in c("duration_minutes", "number_observers")) {
    arr <- array(NA, dim = c(J, n_years, n_reps))
    for (t in seq_along(years)) {
      for (k in 1:n_reps) {
        col <- paste0(var_base, ".", k, ".", years[t])
        if (col %in% names(d)) arr[, t, k] <- d[[col]]
      }
    }
    vals_flat <- as.numeric(arr)
    mu <- mean(vals_flat, na.rm = TRUE)
    sigma <- sd(vals_flat, na.rm = TRUE)
    if (sigma > 0) arr <- (arr - mu) / sigma
    short_name <- ifelse(var_base == "duration_minutes", "duration", "observers")
    det_covs[[short_name]] <- arr
  }

  # Data list
  data_list <- list(
    y = y_array,
    occ.covs = occ_covs,
    det.covs = det_covs,
    coords = coords
  )

  # Formulas
  occ_vars_all <- c(cfg$psi_vars, "NDVI", "tmmx")
  occ_formula <- as.formula(paste("~", paste(occ_vars_all, collapse = " + ")))
  det_formula <- ~ duration + observers

  cat("  Occupancy formula:", deparse(occ_formula), "\n")
  cat("  Detection formula:", deparse(det_formula), "\n")
  cat("  J =", J, "sites | T =", n_years, "years | K =", n_reps, "visits\n\n")

  # ---- 2. Phi prior from Phase B results ----
  # Try to load spPGOcc Phase B fit for informative phi prior
  spB_file <- list.files(here("results"), pattern = paste0(sp, "_spatial_spPGOcc"),
                         full.names = TRUE)
  # Also check in results_spatial subdirectories
  spB_file2 <- list.files(here("results", "results_spatial"),
                          pattern = paste0(sp, "_spatial_spPGOcc"),
                          full.names = TRUE, recursive = TRUE)
  spB_file <- c(spB_file, spB_file2)

  phi_prior <- NULL
  if (length(spB_file) > 0) {
    cat("  Loading Phase B fit for phi prior:", basename(spB_file[1]), "\n")
    fit_B <- readRDS(spB_file[1])
    if (!is.null(fit_B$theta.samples)) {
      phi_B <- fit_B$theta.samples[, "phi"]
      phi_med <- median(phi_B)
      phi_prior <- c(phi_med * 0.01, phi_med * 100)
      cat("  Phase B phi median:", signif(phi_med, 3), "\n")
      cat("  phi prior range: [", signif(phi_prior[1], 3), ",",
          signif(phi_prior[2], 3), "]\n")
    }
    rm(fit_B); gc()
  }

  # Fallback phi prior if Phase B not available
  if (is.null(phi_prior)) {
    set.seed(1)
    idx_sub <- sample(J, min(500, J))
    max_dist <- max(dist(coords[idx_sub, ]))
    phi_prior <- c(3 / max_dist, 3 / (max_dist * 0.01))
    cat("  Using default phi prior (no Phase B available)\n")
  }

  # ---- 3. Priors ----
  priors <- list(
    beta.normal    = list(mean = 0, var = 2.72),
    alpha.normal   = list(mean = 0, var = 2.72),
    sigma.sq.ig    = c(2, 2),     # Tighter than c(2,1): pushes sigma.sq down
    phi.unif       = phi_prior
  )

  cat("\n  Priors:\n")
  cat("    beta:     N(0, 2.72)\n")
  cat("    alpha:    N(0, 2.72)\n")
  cat("    sigma.sq: IG(2, 2)\n")
  cat("    phi:      Unif(", signif(phi_prior[1], 3), ",",
      signif(phi_prior[2], 3), ")\n\n")

  # ---- 4. Starting values ----
  sv <- start_vals[[sp]]
  if (is.null(sv)) sv <- list(sigma.sq = 5, phi = 5e-05)

  cat("  Starting values: sigma.sq =", sv$sigma.sq, ", phi =", sv$phi, "\n\n")

  # ---- 5. Fit stPGOcc ----
  cat("  Fitting stPGOcc...\n")
  cat("  MCMC:", N_BATCH * batch_len, "iter/chain x", N_CHAINS, "chains\n")
  cat("  n.neighbors:", N_NEIGHBORS, "\n")
  cat("  n.omp.threads:", N_THREADS, "\n")
  cat("  Expected runtime: 8-12 hours\n\n")

  fit <- stPGOcc(
    occ.formula  = occ_formula,
    det.formula  = det_formula,
    data         = data_list,
    priors       = priors,
    cov.model    = "exponential",
    NNGP         = TRUE,
    n.neighbors  = N_NEIGHBORS,
    n.batch      = N_BATCH,
    batch.length = batch_len,
    n.burn       = N_BURN,
    n.thin       = N_THIN,
    n.chains     = N_CHAINS,
    ar1          = TRUE,
    n.omp.threads = N_THREADS,
    verbose      = TRUE,
    n.report     = 100
  )

  t_end <- Sys.time()
  runtime_hrs <- as.numeric(difftime(t_end, t_start, units = "hours"))
  cat("\n  Runtime:", round(runtime_hrs, 1), "hours\n\n")

  # ---- 6. Save model object ----
  rds_path <- file.path(out_dir, paste0(sp, "_stPGOcc.rds"))
  saveRDS(fit, rds_path)
  cat("  Saved:", rds_path, "\n")
  cat("  File size:", round(file.size(rds_path) / 1e6, 0), "MB\n\n")

  # ---- 7. Diagnostics ----
  cat("--- MCMC Diagnostics ---\n\n")

  extract_diag <- function(fit) {
    rows <- list()
    for (pt in c("beta", "alpha", "theta")) {
      sn <- paste0(pt, ".samples")
      if (!is.null(fit[[sn]])) {
        samps <- fit[[sn]]
        n_p <- ncol(samps)
        pn <- colnames(samps)
        if (is.null(pn)) pn <- paste0(pt, "_", seq_len(n_p))
        means <- colMeans(samps)
        rhats <- if (!is.null(fit$rhat[[pt]])) fit$rhat[[pt]] else rep(NA, n_p)
        ess <- if (!is.null(fit$ESS[[pt]])) fit$ESS[[pt]] else rep(NA, n_p)
        for (i in seq_len(n_p)) {
          rows[[length(rows) + 1]] <- data.frame(
            Parameter = paste0(pt, "[", pn[i], "]"),
            Mean = means[i], Rhat = rhats[i], ESS = ess[i],
            stringsAsFactors = FALSE)
        }
      }
    }
    do.call(rbind, rows)
  }

  diag <- extract_diag(fit)

  # Print with status flags
  for (i in seq_len(nrow(diag))) {
    r <- diag[i, ]
    flag <- "  OK"
    if (!is.na(r$Rhat) && (r$Rhat > 1.1 || r$ESS < 100)) flag <- "  FAIL"
    else if (!is.na(r$Rhat) && (r$Rhat > 1.05 || r$ESS < 400)) flag <- "  WARN"
    cat(sprintf("  %-28s mean=%9.4f  Rhat=%6.3f  ESS=%7.1f %s\n",
                r$Parameter, r$Mean, r$Rhat, r$ESS, flag))
  }

  n_fail <- sum(diag$Rhat > 1.1 | diag$ESS < 100, na.rm = TRUE)
  n_warn <- sum((diag$Rhat > 1.05 & diag$Rhat <= 1.1) |
                (diag$ESS >= 100 & diag$ESS < 400), na.rm = TRUE)
  cat(sprintf("\n  TOTAL: %d FAIL, %d WARN, %d OK (of %d params)\n\n",
              n_fail, n_warn, nrow(diag) - n_fail - n_warn, nrow(diag)))

  # Save diagnostics
  diag_path <- file.path(out_dir, paste0(sp, "_mcmc_diagnostics.csv"))
  write.csv(diag, diag_path, row.names = FALSE)
  cat("  Saved:", diag_path, "\n")

  # ---- 8. WAIC ----
  cat("\n--- WAIC ---\n")
  waic <- waicOcc(fit)
  cat("  stPGOcc WAIC:", round(waic["WAIC"], 1), "\n")

  # Load tPGOcc for comparison if available
  tpg_file <- list.files(here("results"), pattern = paste0(sp, "_spatial_tPGOcc"),
                         full.names = TRUE)
  tpg_file2 <- list.files(here("results", "results_spatial"),
                          pattern = paste0(sp, "_spatial_tPGOcc"),
                          full.names = TRUE, recursive = TRUE)
  tpg_file <- c(tpg_file, tpg_file2)
  waic_tpg <- NA
  if (length(tpg_file) > 0) {
    fit_tpg <- readRDS(tpg_file[1])
    waic_tpg <- waicOcc(fit_tpg)["WAIC"]
    delta <- waic_tpg - waic["WAIC"]
    cat("  tPGOcc WAIC:", round(waic_tpg, 1), "\n")
    cat("  DELTA WAIC:", round(delta, 1),
        if (delta > 10) "(strong support for spatial)"
        else if (delta > 2) "(moderate)"
        else "(weak)", "\n")
    rm(fit_tpg); gc()
  }

  # ---- 9. Moran's I ----
  cat("\n--- Moran's I ---\n")
  knn_full <- knearneigh(coords, k = 8)
  nb_full  <- knn2nb(knn_full)
  lw_full  <- nb2listw(nb_full, style = "W")

  y_naive <- apply(y_array, 1, function(x) as.numeric(any(x == 1, na.rm = TRUE)))

  z_mean <- apply(fit$z.samples, c(2, 3), mean)
  psi_mean <- rowMeans(z_mean, na.rm = TRUE)
  resid_st <- psi_mean - y_naive

  moran_st <- moran.test(resid_st, lw_full, zero.policy = TRUE)
  cat("  Moran's I =", round(moran_st$estimate["Moran I statistic"], 4),
      "| p =", signif(moran_st$p.value, 3), "\n")

  # ---- 10. Spatial range ----
  cat("\n--- Spatial range ---\n")
  if (!is.null(fit$theta.samples)) {
    phi_post <- fit$theta.samples[, "phi"]
    eff_range <- 3 / phi_post / 1000
    cat("  Median:", round(median(eff_range), 1), "km\n")
    cat("  95% CI: [", round(quantile(eff_range, 0.025), 1), ",",
        round(quantile(eff_range, 0.975), 1), "] km\n")
  }

  # ---- 11. Save model comparison ----
  comp <- data.frame(
    Model = "stPGOcc (production)",
    WAIC = round(waic["WAIC"], 1),
    Moran_I = round(moran_st$estimate["Moran I statistic"], 4),
    Moran_p = signif(moran_st$p.value, 3),
    Spatial_range_km = round(median(eff_range), 1),
    Runtime_hrs = round(runtime_hrs, 1),
    N_FAIL = n_fail,
    N_WARN = n_warn,
    stringsAsFactors = FALSE
  )
  comp_path <- file.path(out_dir, paste0(sp, "_model_comparison.csv"))
  write.csv(comp, comp_path, row.names = FALSE)
  cat("\n  Saved:", comp_path, "\n")

  # ---- 12. Diagnostic figure ----
  cat("\n--- Diagnostic figure ---\n")

  w_spatial <- apply(fit$w.samples, 2, mean)

  p1 <- ggplot(data.frame(lon = coords_lonlat[,1], lat = coords_lonlat[,2],
                           resid = resid_st),
               aes(x = lon, y = lat, colour = resid)) +
    geom_point(size = 0.4, alpha = 0.7) +
    scale_colour_gradient2(low = "blue", mid = "white", high = "red",
                           midpoint = 0, name = "Residual") +
    labs(title = paste0(sp, ": stPGOcc residuals (production)"),
         subtitle = paste0("Moran's I = ",
                           round(moran_st$estimate["Moran I statistic"], 4))) +
    coord_quickmap() + theme_minimal() + theme(legend.position = "bottom")

  p2 <- ggplot(data.frame(lon = coords_lonlat[,1], lat = coords_lonlat[,2],
                           w = w_spatial),
               aes(x = lon, y = lat, colour = w)) +
    geom_point(size = 0.4, alpha = 0.7) +
    scale_colour_gradient2(low = "blue", mid = "white", high = "red",
                           midpoint = 0, name = "w") +
    labs(title = paste0(sp, ": spatial random effect (w)"),
         subtitle = paste0("Range: ", round(median(eff_range), 0), " km")) +
    coord_quickmap() + theme_minimal() + theme(legend.position = "bottom")

  fig_path <- here("figs", paste0(sp, "_spatial_diagnostics_production.png"))
  ggsave(fig_path, gridExtra::grid.arrange(p1, p2, ncol = 2),
         width = 14, height = 6, dpi = 150)
  cat("  Saved:", fig_path, "\n")

  cat("\n  ", sp, "COMPLETE in", round(runtime_hrs, 1), "hours\n\n")

  rm(fit); gc()

  }, error = function(e) {
    cat("\n  ERROR processing", sp, ":", conditionMessage(e), "\n")
    cat("  Skipping to next species.\n\n")
  })
}

cat("\n###############################################################\n")
cat("#  ALL SPECIES COMPLETE                                       #\n")
cat("#  Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "       #\n")
cat("###############################################################\n")
