###############################################################################
# 23_core_spa_analysis.R
#
# Are local extinctions concentrated in core pseudo-steppe habitat, and do
# Special Protection Areas (SPAs) retain occupied cells better than
# comparable unprotected cells?
#
# Core habitat is defined independently of the occupancy model, from the
# CORINE 2018 pseudo-steppe mask used for stepRep (strict: CLC 211, 231, 321,
# 333; broad: + 242, 243, 244): core = cell fraction of strict mask >= 0.5,
# with sensitivity at 0.3 / 0.7 and for the broad mask. SPA cover is the
# fraction of each model cell (WorldClim 2.5' cell) inside the SPA union.
#
# Outcomes use posterior (data-conditioned) events of the colext models over
# the recent transitions (2018->19 ... 2022->23), in cells with >= 1
# detection, and the data-anchored subset of transitions surveyed in both
# years ("both"; script 21):
#   extinction rate   = sum(extinctions) / sum(P(occupied at t))
#   net change        = sum(colonisations) - sum(extinctions)
#
# SPA effect: detected cells mostly inside SPAs (cover >= 0.5) are matched
# with replacement to detected cells essentially outside (cover < 0.1) by
# nearest Mahalanobis distance on strict steppe fraction and logit(psi1)
# (sensitivity: + UTM coordinates). The contrast is the difference in
# extinction rate (SPA - matched controls). Uncertainty: 1,000 draws that
# combine coefficient uncertainty (multivariate normal) with resampling of
# treated cells (and their matched controls).
#
# Inputs:
#   results/{sp}_model_object.rds
#   data/processed_2023/{sp}/{sp}_occ_wide_dynamic.csv
#   CORINE 2018 100 m raster (see build_stepRep.R), SPAs (R/spatial_layers.R)
#
# Outputs (new):
#   results/core_spa/cell_habitat.csv
#   results/core_spa/class_rates.csv
#   results/core_spa/spa_matched.csv
#   results/core_spa/spa_balance.csv
#   results/core_spa/fig_class_rates.png
###############################################################################

library(here)
library(unmarked)
library(MASS)
library(sf)
library(terra)
library(ggplot2)

source(here("R", "colext_utils.R"))
source(here("R", "spatial_layers.R"))
set.seed(42)

species <- c("otitar", "ptealc", "pteori", "tettet")
species_names <- c(otitar = "O. tarda", ptealc = "P. alchata",
                   pteori = "P. orientalis", tettet = "T. tetrax")
N_DRAWS <- 1000
CORE_THRESHOLDS <- c(0.3, 0.5, 0.7)
out_dir <- here("results", "core_spa")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
q <- function(x, p) as.numeric(quantile(x, p, na.rm = TRUE))

###############################################################################
# 1. CELL HABITAT AND SPA COVER (cached)
###############################################################################

occ_list <- lapply(species, function(sp) {
  read.csv(here("data", "processed_2023", sp,
                paste0(sp, "_occ_wide_dynamic.csv")))
})
names(occ_list) <- species
all_cells <- sort(unique(unlist(lapply(occ_list, `[[`, "cells"))))

hab_path <- file.path(out_dir, "cell_habitat.csv")
if (!file.exists(hab_path)) {
  message("Computing CORINE pseudo-steppe and SPA cover for ",
          length(all_cells), " cells ...")
  # Aggregating the 100 m layers onto the WorldClim grid with
  # project(method = "average") gives the cover fraction of every cell in
  # one pass (per-polygon extraction took > 35 min).
  grid <- wc_grid()
  xy_c <- xyFromCell(grid, all_cells)
  grid <- crop(grid, ext(min(xy_c[, 1]) - 0.1, max(xy_c[, 1]) + 0.1,
                         min(xy_c[, 2]) - 0.1, max(xy_c[, 2]) + 0.1))
  clc <- rast(here("data-raw", "data", "corine", "Results",
                   "u2018_clc2018_v2020_20u1_raster100m", "DATA",
                   "U2018_CLC2018_V2020_20u1.tif"))
  clc <- crop(clc, project(ext(grid), crs(grid), crs(clc)) + 2000)
  # Raster values 1..44 are CLC classes (value 48 = NoData); lookup as in
  # build_stepRep.R: 211 = 12, 231 = 18, 242 = 20, 243 = 21, 244 = 22,
  # 321 = 26, 333 = 32
  as_mask <- function(vals) {
    classify(clc, rcl = cbind(0:48, as.integer(0:48 %in% vals)), others = NA)
  }
  to_grid <- function(m) project(m, grid, method = "average")
  strict <- to_grid(as_mask(c(12L, 18L, 26L, 32L)))
  broad  <- to_grid(as_mask(c(12L, 18L, 26L, 32L, 20L, 21L, 22L)))
  spa_r  <- rasterize(vect(st_transform(st_as_sf(load_spa_union()), crs(clc))),
                      clc, field = 1, background = 0)
  spa_g  <- to_grid(spa_r)

  cell_xy <- xyFromCell(wc_grid(), all_cells)
  getv <- function(r) terra::extract(r, cell_xy)[, 1]
  polys <- wc_cell_polygons(all_cells)
  xy <- st_coordinates(st_centroid(polys))
  write.csv(data.frame(cells = all_cells, x = xy[, 1], y = xy[, 2],
                       area_km2 = as.numeric(st_area(polys)) / 1e6,
                       steppe_strict = getv(strict), steppe_broad = getv(broad),
                       spa_frac = getv(spa_g)),
            hab_path, row.names = FALSE)
}
hab <- read.csv(hab_path)
message("Cell area (km2): ", paste(round(range(hab$area_km2), 1), collapse = "-"))

###############################################################################
# 2. POSTERIOR EVENTS, CLASS RATES AND MATCHED SPA CONTRAST
###############################################################################

# Nearest-neighbour Mahalanobis matching with replacement; returns the index
# of the matched control for each treated unit
match_nn <- function(Xt, Xc) {
  S_inv <- solve(cov(rbind(Xt, Xc)))
  apply(Xt, 1, function(xt) {
    d <- mahalanobis(Xc, xt, S_inv, inverted = TRUE)
    which.min(d)
  })
}
smd <- function(a, b) (mean(a) - mean(b)) / sqrt((var(a) + var(b)) / 2)

rate_rows <- list(); match_rows <- list(); bal_rows <- list()

for (sp in species) {
  message("\n=== ", sp, " ===")
  mod <- readRDS(here("results", paste0(sp, "_model_object.rds")))
  D <- colext_design(mod)
  n_t <- D$n_t
  occ <- occ_list[[sp]][D$keep, ]
  y_csv <- sapply(occ[, grep("^y\\.[0-9]+\\.[0-9]{4}$", names(occ))],
                  function(v) as.numeric(as.logical(v)))
  stopifnot(all(y_csv[!is.na(D$y)] == D$y[!is.na(D$y)]))
  h <- hab[match(occ$cells, hab$cells), ]
  stopifnot(!anyNA(h$cells))

  naive <- sapply(1:n_t, function(t) {
    yt <- D$y[, ((t - 1) * D$J + 1):(t * D$J), drop = FALSE]
    ifelse(rowSums(!is.na(yt)) > 0, as.integer(rowSums(yt, na.rm = TRUE) > 0),
           NA_integer_)
  })
  surv <- !is.na(naive)
  both <- surv[, 1:(n_t - 1)] & surv[, 2:n_t]
  det <- rowSums(naive, na.rm = TRUE) > 0
  recent <- 2:(n_t - 1)

  cf <- coef(mod)
  psi1 <- colext_rates(cf, D)$psi1

  # Cell classes (detected cells only)
  classes <- list(det_all = det)
  for (th in CORE_THRESHOLDS) {
    classes[[paste0("core_strict", th)]] <- det & h$steppe_strict >= th
    classes[[paste0("rest_strict", th)]] <- det & h$steppe_strict < th
  }
  classes$core_broad0.5 <- det & h$steppe_broad >= 0.5
  classes$rest_broad0.5 <- det & h$steppe_broad < 0.5
  classes$spa_in  <- det & h$spa_frac >= 0.5
  classes$spa_out <- det & h$spa_frac < 0.1
  classes$core_spa_in  <- classes$core_strict0.5 & h$spa_frac >= 0.5
  classes$core_spa_out <- classes$core_strict0.5 & h$spa_frac < 0.1

  # Matching (fixed across draws)
  tr <- which(classes$spa_in); ctl <- which(classes$spa_out)
  covs <- list(
    base   = cbind(h$steppe_strict, qlogis(pmin(pmax(psi1, 1e-6), 1 - 1e-6))),
    coords = cbind(h$steppe_strict, qlogis(pmin(pmax(psi1, 1e-6), 1 - 1e-6)),
                   h$x / 1e5, h$y / 1e5))
  matches <- lapply(covs, function(X) ctl[match_nn(X[tr, , drop = FALSE],
                                                   X[ctl, , drop = FALSE])])
  for (nm in names(covs)) {
    X <- covs[[nm]]
    for (j in seq_len(ncol(X))) {
      bal_rows[[paste(sp, nm, j)]] <- data.frame(
        species = sp, spec = nm,
        covariate = c("steppe_strict", "logit_psi1", "x", "y")[j],
        smd_before = smd(X[tr, j], X[ctl, j]),
        smd_after = smd(X[tr, j], X[matches[[nm]], j]))
    }
  }

  # Draws
  V <- vcov(mod)
  V_pd <- tryCatch({ chol(V); V },
                   error = function(e) as.matrix(Matrix::nearPD(V)$mat))
  draws <- MASS::mvrnorm(N_DRAWS, mu = cf, Sigma = V_pd)
  per <- c("rec", "both")
  rate_d <- array(NA_real_, c(N_DRAWS, length(classes), 2, 2),
                  dimnames = list(NULL, names(classes), per, c("ext_rate", "net")))
  match_d <- array(NA_real_, c(N_DRAWS, length(covs), 2),
                   dimnames = list(NULL, names(covs), per))

  for (k in seq_len(N_DRAWS)) {
    fb <- forward_backward(colext_rates(draws[k, ], D), D)
    occ_exp <- fb$z1[, 1:(n_t - 1)]
    cell_sum <- list(
      rec = list(E = rowSums(fb$ext[, recent]), O = rowSums(occ_exp[, recent]),
                 C = rowSums(fb$col[, recent])),
      both = list(E = rowSums((fb$ext * both)[, recent]),
                  O = rowSums((occ_exp * both)[, recent]),
                  C = rowSums((fb$col * both)[, recent])))
    for (p in per) {
      cs <- cell_sum[[p]]
      for (j in names(classes)) {
        w <- classes[[j]]
        rate_d[k, j, p, ] <- c(sum(cs$E[w]) / sum(cs$O[w]),
                               sum(cs$C[w]) - sum(cs$E[w]))
      }
      # Bootstrap treated cells (with their matched controls)
      b <- sample(seq_along(tr), replace = TRUE)
      for (nm in names(covs)) {
        t_i <- tr[b]; c_i <- matches[[nm]][b]
        match_d[k, nm, p] <- sum(cs$E[t_i]) / sum(cs$O[t_i]) -
                             sum(cs$E[c_i]) / sum(cs$O[c_i])
      }
    }
  }

  for (j in names(classes)) for (p in per) {
    e <- rate_d[, j, p, "ext_rate"]; n <- rate_d[, j, p, "net"]
    rate_rows[[paste(sp, j, p)]] <- data.frame(
      species = sp, class = j, period = p, n_cells = sum(classes[[j]]),
      ext_rate_med = median(e), ext_rate_lo = q(e, .025), ext_rate_hi = q(e, .975),
      net_med = median(n), net_lo = q(n, .025), net_hi = q(n, .975),
      p_net_loss = mean(n < 0))
  }
  for (nm in names(covs)) for (p in per) {
    d <- match_d[, nm, p]
    match_rows[[paste(sp, nm, p)]] <- data.frame(
      species = sp, spec = nm, period = p, n_treated = length(tr),
      n_control_pool = length(ctl), n_unique_controls = length(unique(matches[[nm]])),
      # Resamples whose treated or control cells carry no expected occupancy
      # give an undefined rate ratio (0/0) and are dropped
      n_undefined = sum(!is.finite(d)),
      diff_ext_rate_med = median(d[is.finite(d)]), diff_lo = q(d[is.finite(d)], .025),
      diff_hi = q(d[is.finite(d)], .975), p_spa_higher = mean(d[is.finite(d)] > 0))
  }
}

rate_df  <- do.call(rbind, rate_rows);  rownames(rate_df)  <- NULL
match_df <- do.call(rbind, match_rows); rownames(match_df) <- NULL
bal_df   <- do.call(rbind, bal_rows);   rownames(bal_df)   <- NULL
write.csv(rate_df,  file.path(out_dir, "class_rates.csv"), row.names = FALSE)
write.csv(match_df, file.path(out_dir, "spa_matched.csv"), row.names = FALSE)
write.csv(bal_df,   file.path(out_dir, "spa_balance.csv"), row.names = FALSE)
options(width = 220)
print(rate_df, digits = 3)
print(match_df, digits = 3)
print(bal_df, digits = 2)

###############################################################################
# 3. FIGURE: extinction rate by class (data-anchored transitions)
###############################################################################

show <- c(core_strict0.5 = "Core (>= 50% pseudo-steppe)",
          rest_strict0.5 = "Rest",
          spa_in = "Mostly inside SPA", spa_out = "Outside SPA",
          core_spa_in = "Core, inside SPA", core_spa_out = "Core, outside SPA")
fd <- rate_df[rate_df$class %in% names(show), ]
fd$class <- factor(show[fd$class], levels = rev(show))
fd$species <- factor(species_names[fd$species], levels = species_names)
fd$period <- factor(fd$period, levels = c("both", "rec"),
                    labels = c("Transitions surveyed in both years",
                               "All recent transitions"))
p <- ggplot(fd, aes(x = ext_rate_med, y = class, colour = period)) +
  geom_pointrange(aes(xmin = ext_rate_lo, xmax = ext_rate_hi),
                  position = position_dodge(width = 0.5), size = 0.25) +
  facet_wrap(~ species, scales = "free_x") +
  labs(x = "Posterior local extinction rate per occupied cell-year, 2018-2023",
       y = NULL, colour = NULL,
       caption = "Cells with >= 1 eBird detection. Bars: 95% intervals over 1,000 coefficient draws.") +
  theme_bw(base_size = 9) +
  theme(strip.text = element_text(face = "italic"), legend.position = "bottom")
ggsave(file.path(out_dir, "fig_class_rates.png"), p, width = 8, height = 6, dpi = 300)

sessionInfo()
