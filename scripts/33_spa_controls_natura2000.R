###############################################################################
# 33_spa_controls_natura2000.R
#
# Are the "unprotected" control cells of the SPA contrast (script 23) really
# outside the protected-area network? Controls are cells with < 10 % SPA
# cover, but they can lie in Natura 2000 sites designated under the Habitats
# Directive (SCIs/SACs). This script
#   (1) computes, for every model cell, the fraction inside any Natura 2000
#       site and inside Habitats Directive sites (SCI + SAC);
#   (2) reproduces the script 23 matching (primary model; base and effort
#       specifications) and reports how many controls lie in Natura 2000;
#   (3) repeats the matched contrast with the control pool restricted to
#       cells with < 10 % cover of any Natura 2000 site, using the same
#       coefficient draws and cell resamples as the original pool, so the two
#       contrasts differ only in the controls.
#
# The matching and the contrast repeat the logic of script 23 (sections 2 and
# its matched contrast) rather than refactoring it, so that script 23 and its
# published outputs stay unchanged. Keep both in sync if either changes.
#
# Inputs:  results/{sp}_model_object.rds, colext inputs (colext_input_path()),
#          results/core_spa/cell_habitat.csv (script 23), CORINE 2018 100 m
#          raster (grid template), Natura 2000 (R/spatial_layers.R)
# Outputs: results/core_spa/cell_natura2000.csv
#          results/core_spa/spa_controls_natura2000.csv
#          results/core_spa/spa_matched_outside_n2k.csv
###############################################################################

library(here)
library(unmarked)
library(MASS)
library(sf)
library(terra)

source(here("R", "colext_utils.R"))
source(here("R", "spatial_layers.R"))
set.seed(42)

species <- c("otitar", "tettet", "pteori")   # extinction identifiable
N_DRAWS <- 1000
N2K_MAX <- 0.1                               # control cells: < 10 % Natura 2000
hab_dir <- here("results", "core_spa")
q <- function(x, p) as.numeric(quantile(x, p, na.rm = TRUE))

###############################################################################
# 1. NATURA 2000 COVER PER CELL (cached; same raster approach as script 23)
###############################################################################

hab <- read.csv(file.path(hab_dir, "cell_habitat.csv"))
n2k_path <- file.path(hab_dir, "cell_natura2000.csv")
if (!file.exists(n2k_path)) {
  grid <- wc_grid()
  xy_c <- xyFromCell(grid, hab$cells)
  grid <- crop(grid, ext(min(xy_c[, 1]) - 0.1, max(xy_c[, 1]) + 0.1,
                         min(xy_c[, 2]) - 0.1, max(xy_c[, 2]) + 0.1))
  clc <- rast(here("data-raw", "data", "corine", "Results",
                   "u2018_clc2018_v2020_20u1_raster100m", "DATA",
                   "U2018_CLC2018_V2020_20u1.tif"))
  clc <- crop(clc, project(ext(grid), crs(grid), crs(clc)) + 2000)
  cover <- function(sites) {
    r <- rasterize(vect(st_transform(st_as_sf(st_union(sites)), crs(clc))),
                   clc, field = 1, background = 0)
    terra::extract(project(r, grid, method = "average"),
                   xyFromCell(wc_grid(), hab$cells))[, 1]
  }
  n2k <- load_n2k_sites()
  hd <- n2k[!grepl("specialProtectionArea", n2k$desig), ]
  write.csv(data.frame(cells = hab$cells, n2k_frac = cover(n2k), hd_frac = cover(hd)),
            n2k_path, row.names = FALSE)
}
n2k_cov <- read.csv(n2k_path)
hab$n2k_frac <- n2k_cov$n2k_frac[match(hab$cells, n2k_cov$cells)]
hab$hd_frac <- n2k_cov$hd_frac[match(hab$cells, n2k_cov$cells)]
stopifnot(!anyNA(hab$n2k_frac))

###############################################################################
# 2-3. MATCHING, CONTROL COVER AND CONTRAST WITH CONTROLS OUTSIDE NATURA 2000
###############################################################################

# As in script 23
match_nn <- function(Xt, Xc) {
  S_inv <- solve(cov(rbind(Xt, Xc)))
  apply(Xt, 1, function(xt) which.min(mahalanobis(Xc, xt, S_inv, inverted = TRUE)))
}
smd <- function(a, b) (mean(a) - mean(b)) / sqrt((var(a) + var(b)) / 2)

desc_rows <- list(); contrast_rows <- list()
for (sp in species) {
  message("=== ", sp)
  mod <- readRDS(colext_model_path(sp))
  D <- colext_design(mod)
  n_t <- D$n_t
  occ <- read.csv(colext_input_path(sp))[D$keep, ]
  h <- hab[match(occ$cells, hab$cells), ]
  stopifnot(!anyNA(h$cells))
  naive <- sapply(1:n_t, function(t) {
    yt <- D$y[, ((t - 1) * D$J + 1):(t * D$J), drop = FALSE]
    ifelse(rowSums(!is.na(yt)) > 0, as.integer(rowSums(yt, na.rm = TRUE) > 0), NA_integer_)
  })
  surv <- !is.na(naive)
  both <- surv[, 1:(n_t - 1)] & surv[, 2:n_t]
  det <- rowSums(naive, na.rm = TRUE) > 0
  recent <- 2:(n_t - 1)
  psi1 <- colext_rates(coef(mod), D)$psi1
  lpsi <- qlogis(pmin(pmax(psi1, 1e-6), 1 - 1e-6))
  visits <- rowSums(!is.na(D$y)) / pmax(rowSums(surv), 1)

  tr <- which(det & h$spa_frac >= 0.5)
  pools <- list(all = which(det & h$spa_frac < 0.1),
                outside_n2k = which(det & h$spa_frac < 0.1 & h$n2k_frac < N2K_MAX))
  covs <- list(base = cbind(h$steppe_strict, lpsi),
               effort = cbind(h$steppe_strict, lpsi, visits))
  matches <- list()
  for (pl in names(pools)) for (nm in names(covs)) {
    ctl <- pools[[pl]]; X <- covs[[nm]]
    m <- ctl[match_nn(X[tr, , drop = FALSE], X[ctl, , drop = FALSE])]
    matches[[paste(pl, nm)]] <- m
    desc_rows[[paste(sp, pl, nm)]] <- data.frame(
      species = sp, pool = pl, spec = nm, n_treated = length(tr),
      n_pool = length(ctl), n_pool_in_n2k = sum(h$n2k_frac[ctl] >= N2K_MAX),
      n_pool_in_hd = sum(h$hd_frac[ctl] >= N2K_MAX),
      n_unique_controls = length(unique(m)),
      n_unique_controls_in_n2k = sum(h$n2k_frac[unique(m)] >= N2K_MAX),
      share_pairs_control_in_n2k = mean(h$n2k_frac[m] >= N2K_MAX),
      mean_n2k_frac_controls = mean(h$n2k_frac[m]),
      smd_after_steppe = smd(h$steppe_strict[tr], h$steppe_strict[m]),
      smd_after_lpsi = smd(lpsi[tr], lpsi[m]),
      smd_after_visits = smd(visits[tr], visits[m]))
  }

  # Contrast on anchored recent transitions; the same draws and the same
  # resample of treated cells serve both control pools
  draws <- draw_coefs(mod, N_DRAWS)
  point_only <- attr(draws, "point_only")
  diff_d <- ratio_d <- matrix(NA_real_, N_DRAWS, length(matches),
                              dimnames = list(NULL, names(matches)))
  for (k in seq_len(N_DRAWS)) {
    fb <- forward_backward(colext_rates(draws[k, ], D), D)
    E <- rowSums((fb$ext * both)[, recent])
    O <- rowSums((fb$z1[, 1:(n_t - 1)] * both)[, recent])
    b <- sample(seq_along(tr), replace = TRUE)
    for (mm in names(matches)) {
      t_i <- tr[b]; c_i <- matches[[mm]][b]
      r_t <- sum(E[t_i]) / sum(O[t_i]); r_c <- sum(E[c_i]) / sum(O[c_i])
      diff_d[k, mm] <- r_t - r_c; ratio_d[k, mm] <- r_t / r_c
    }
  }
  for (mm in names(matches)) {
    d <- diff_d[, mm]; r <- ratio_d[, mm]
    d <- d[is.finite(d)]; r <- r[is.finite(r)]
    contrast_rows[[paste(sp, mm)]] <- data.frame(
      species = sp, pool = sub(" .*", "", mm), spec = sub(".* ", "", mm),
      point_only = point_only, n_undefined = N_DRAWS - length(d),
      diff_med = median(d), diff_lo = q(d, .025), diff_hi = q(d, .975),
      p_spa_higher = mean(d > 0),
      ratio_med = median(r), ratio_lo = q(r, .025), ratio_hi = q(r, .975))
  }
}

desc_df <- do.call(rbind, desc_rows); rownames(desc_df) <- NULL
contrast_df <- do.call(rbind, contrast_rows); rownames(contrast_df) <- NULL
write.csv(desc_df, file.path(hab_dir, "spa_controls_natura2000.csv"), row.names = FALSE)
write.csv(contrast_df, file.path(hab_dir, "spa_matched_outside_n2k.csv"), row.names = FALSE)
print(desc_df, digits = 3)
print(contrast_df, digits = 3)

sessionInfo()
