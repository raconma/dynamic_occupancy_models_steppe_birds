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
# (sensitivity: + UTM coordinates; + checklists per surveyed year). The contrast is the difference in
# extinction rate (SPA - matched controls). Uncertainty: 1,000 draws that
# combine coefficient uncertainty (multivariate normal) with resampling of
# treated cells (and their matched controls).
#
# Inputs:
#   results/{sp}_model_object.rds
#   data/processed_2023/{sp}/{sp}_occ_wide_dynamic.csv (stepRep variants:
#   reports/script_4_v2_inputs/, see colext_input_path())
#   CORINE 2018 100 m raster (see build_stepRep.R), SPAs (R/spatial_layers.R)
#
# Outputs (new):
#   results/core_spa/cell_habitat.csv
#   results/core_spa/class_rates.csv
#   results/core_spa/spa_matched.csv  (SPA minus matched difference and SPA / matched ratio)
#   results/core_spa/spa_balance.csv   (matching covariates and survey effort)
#   results/core_spa/core_contrast.csv (marginal minus core extinction rate and
#                                       ratio, with cell resampling per draw)
#   results/core_spa/spa_threshold.csv (SPA contrast under alternative SPA
#                                       cover thresholds, anchored transitions)
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

# Model variant: "primary" (default) or a robustness refit from script 28,
# e.g. Rscript scripts/23_core_spa_analysis.R eps_null
variant <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(variant)) variant <- "primary"

species <- c("otitar", "ptealc", "pteori", "tettet")
species_names <- c(otitar = "O. tarda", ptealc = "P. alchata",
                   pteori = "P. orientalis", tettet = "T. tetrax")
N_DRAWS <- 1000
CORE_THRESHOLDS <- c(0.3, 0.5, 0.7)
hab_dir <- here("results", "core_spa")
out_dir <- if (variant == "primary") hab_dir else file.path(hab_dir, variant)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
q <- function(x, p) as.numeric(quantile(x, p, na.rm = TRUE))

###############################################################################
# 1. CELL HABITAT AND SPA COVER (cached)
###############################################################################

occ_list <- lapply(species, function(sp) read.csv(colext_input_path(sp, variant)))
names(occ_list) <- species
all_cells <- sort(unique(unlist(lapply(occ_list, `[[`, "cells"))))

hab_path <- file.path(hab_dir, "cell_habitat.csv")
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
contrast_store <- list()   # per-draw cell sums for the core-vs-marginal contrast

for (sp in species) {
  message("\n=== ", sp, " ===")
  mod <- readRDS(colext_model_path(sp, variant))
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

  # Survey effort per cell: checklists per surveyed year and number of
  # anchored transitions in the recent period
  effort <- cbind(visits_per_year = rowSums(!is.na(D$y)) / pmax(rowSums(surv), 1),
                  anchored_transitions = rowSums(both[, recent, drop = FALSE]))

  # Matching (fixed across draws). "effort" adds checklists per surveyed year,
  # because SPA cells are surveyed less often and fewer visits give less
  # evidence of absence
  tr <- which(classes$spa_in); ctl <- which(classes$spa_out)
  lpsi <- qlogis(pmin(pmax(psi1, 1e-6), 1 - 1e-6))
  covs <- list(
    base   = cbind(h$steppe_strict, lpsi),
    coords = cbind(h$steppe_strict, lpsi, h$x / 1e5, h$y / 1e5),
    effort = cbind(h$steppe_strict, lpsi, effort[, "visits_per_year"]))
  cov_labels <- list(base = c("steppe_strict", "logit_psi1"),
                     coords = c("steppe_strict", "logit_psi1", "x", "y"),
                     effort = c("steppe_strict", "logit_psi1", "visits_per_year"))
  matches <- lapply(covs, function(X) ctl[match_nn(X[tr, , drop = FALSE],
                                                   X[ctl, , drop = FALSE])])
  # Balance on the matching covariates and on survey effort
  for (nm in names(covs)) {
    extra <- setdiff(colnames(effort), cov_labels[[nm]])
    X <- cbind(covs[[nm]], effort[, extra, drop = FALSE])
    labs <- c(cov_labels[[nm]], extra)
    for (j in seq_len(ncol(X))) {
      bal_rows[[paste(sp, nm, j)]] <- data.frame(
        species = sp, spec = nm, covariate = labs[j],
        matched_on = j <= ncol(covs[[nm]]),
        smd_before = smd(X[tr, j], X[ctl, j]),
        smd_after = smd(X[tr, j], X[matches[[nm]], j]))
    }
  }

  # Draws
  draws <- draw_coefs(mod, N_DRAWS)
  point_only <- attr(draws, "point_only")
  per <- c("rec", "both")
  rate_d <- array(NA_real_, c(N_DRAWS, length(classes), 2, 2),
                  dimnames = list(NULL, names(classes), per, c("ext_rate", "net")))
  match_d <- array(NA_real_, c(N_DRAWS, length(covs), 2),
                   dimnames = list(NULL, names(covs), per))
  match_r <- match_d   # SPA / matched-control ratio, same resamples as match_d
  det_idx <- which(det)
  E_store <- O_store <- array(NA_real_, c(N_DRAWS, length(det_idx), 2),
                              dimnames = list(NULL, NULL, per))

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
      E_store[k, , p] <- cs$E[det_idx]; O_store[k, , p] <- cs$O[det_idx]
      for (j in names(classes)) {
        w <- classes[[j]]
        rate_d[k, j, p, ] <- c(sum(cs$E[w]) / sum(cs$O[w]),
                               sum(cs$C[w]) - sum(cs$E[w]))
      }
      # Bootstrap treated cells (with their matched controls)
      b <- sample(seq_along(tr), replace = TRUE)
      for (nm in names(covs)) {
        t_i <- tr[b]; c_i <- matches[[nm]][b]
        r_t <- sum(cs$E[t_i]) / sum(cs$O[t_i]); r_c <- sum(cs$E[c_i]) / sum(cs$O[c_i])
        match_d[k, nm, p] <- r_t - r_c
        match_r[k, nm, p] <- r_t / r_c
      }
    }
  }

  pairs <- c(setNames(paste0("strict", CORE_THRESHOLDS), paste0("strict", CORE_THRESHOLDS)),
             broad0.5 = "broad0.5")
  # Observed (naive) extinctions on recent anchored transitions, per cell
  obs_ext <- rowSums(naive[, recent] == 1 & naive[, recent + 1] == 0 &
                       both[, recent], na.rm = TRUE)
  contrast_store[[sp]] <- list(
    E = E_store, O = O_store, point_only = point_only,
    cell = data.frame(spa_frac = h$spa_frac, steppe_strict = h$steppe_strict,
                      lpsi = lpsi, visits_per_year = effort[, "visits_per_year"],
                      obs_ext = obs_ext)[det_idx, ],
    groups = lapply(pairs, function(g) list(
      core = classes[[paste0("core_", g)]][det_idx],
      rest = classes[[paste0("rest_", g)]][det_idx])))

  for (j in names(classes)) for (p in per) {
    e <- rate_d[, j, p, "ext_rate"]; n <- rate_d[, j, p, "net"]
    rate_rows[[paste(sp, j, p)]] <- data.frame(
      species = sp, variant = variant, point_only = point_only,
      class = j, period = p, n_cells = sum(classes[[j]]),
      ext_rate_med = median(e), ext_rate_lo = q(e, .025), ext_rate_hi = q(e, .975),
      net_med = median(n), net_lo = q(n, .025), net_hi = q(n, .975),
      p_net_loss = mean(n < 0))
  }
  for (nm in names(covs)) for (p in per) {
    d <- match_d[, nm, p]
    match_rows[[paste(sp, nm, p)]] <- data.frame(
      species = sp, variant = variant, point_only = point_only,
      spec = nm, period = p, n_treated = length(tr),
      n_control_pool = length(ctl), n_unique_controls = length(unique(matches[[nm]])),
      # Resamples whose treated or control cells carry no expected occupancy
      # give an undefined rate ratio (0/0) and are dropped
      n_undefined = sum(!is.finite(d)),
      diff_ext_rate_med = median(d[is.finite(d)]), diff_lo = q(d[is.finite(d)], .025),
      diff_hi = q(d[is.finite(d)], .975), p_spa_higher = mean(d[is.finite(d)] > 0),
      n_ratio_undefined = sum(!is.finite(match_r[, nm, p])),
      ratio_med = median(match_r[is.finite(match_r[, nm, p]), nm, p]),
      ratio_lo = q(match_r[is.finite(match_r[, nm, p]), nm, p], .025),
      ratio_hi = q(match_r[is.finite(match_r[, nm, p]), nm, p], .975))
  }
}

rate_df  <- do.call(rbind, rate_rows);  rownames(rate_df)  <- NULL
match_df <- do.call(rbind, match_rows); rownames(match_df) <- NULL
bal_df   <- do.call(rbind, bal_rows);   rownames(bal_df)   <- NULL
write.csv(rate_df,  file.path(out_dir, "class_rates.csv"), row.names = FALSE)
write.csv(match_df, file.path(out_dir, "spa_matched.csv"), row.names = FALSE)
write.csv(bal_df,   file.path(out_dir, "spa_balance.csv"), row.names = FALSE)

# Core vs marginal contrast with its own uncertainty: per coefficient draw,
# resample cells with replacement within each habitat class. Run after all
# species so that the random stream of the outputs above is unchanged.
contrast_rows <- list()
for (sp in names(contrast_store)) {
  st <- contrast_store[[sp]]
  for (g in names(st$groups)) for (p in dimnames(st$E)[[3]]) {
    ci <- which(st$groups[[g]]$core); ri <- which(st$groups[[g]]$rest)
    if (length(ci) < 2 || length(ri) < 2) next
    d <- r <- numeric(N_DRAWS)
    for (k in seq_len(N_DRAWS)) {
      bc <- sample(ci, replace = TRUE); br <- sample(ri, replace = TRUE)
      rc <- sum(st$E[k, bc, p]) / sum(st$O[k, bc, p])
      rr <- sum(st$E[k, br, p]) / sum(st$O[k, br, p])
      d[k] <- rr - rc; r[k] <- rr / rc
    }
    ok <- is.finite(d) & is.finite(r)
    contrast_rows[[paste(sp, g, p)]] <- data.frame(
      species = sp, variant = variant, point_only = st$point_only,
      definition = g, period = p, n_core = length(ci), n_marginal = length(ri),
      n_undefined = sum(!ok),
      diff_med = median(d[ok]), diff_lo = q(d[ok], .025), diff_hi = q(d[ok], .975),
      ratio_med = median(r[ok]), ratio_lo = q(r[ok], .025), ratio_hi = q(r[ok], .975),
      p_marginal_higher = mean(d[ok] > 0))
  }
}
contrast_df <- do.call(rbind, contrast_rows); rownames(contrast_df) <- NULL
write.csv(contrast_df, file.path(out_dir, "core_contrast.csv"), row.names = FALSE)

# Sensitivity of the SPA contrast to the SPA cover thresholds (anchored
# transitions, 2018-2023). The main analysis compares cells >= 50 % inside
# SPAs with cells < 10 % inside and leaves out intermediate cover; here the
# treated/control split is varied. Matching and resampling as above. Run
# after the other outputs so that their random stream is unchanged.
SPA_DEFS <- list(in50_out10 = c(0.5, 0.1), in50_out50 = c(0.5, 0.5),
                 in30_out30 = c(0.3, 0.3), in10_out10 = c(0.1, 0.1))
thr_rows <- list()
for (sp in names(contrast_store)) {
  st <- contrast_store[[sp]]; cl <- st$cell
  E <- st$E[, , "both"]; O <- st$O[, , "both"]
  rate_of <- function(k, i) sum(E[k, i]) / sum(O[k, i])
  for (dn in names(SPA_DEFS)) {
    thr <- SPA_DEFS[[dn]]
    tr <- which(cl$spa_frac >= thr[1]); ctl <- which(cl$spa_frac < thr[2])
    X_all <- list(base = cbind(cl$steppe_strict, cl$lpsi),
                  effort = cbind(cl$steppe_strict, cl$lpsi, cl$visits_per_year))
    raw_in <- sapply(seq_len(N_DRAWS), rate_of, i = tr)
    raw_out <- sapply(seq_len(N_DRAWS), rate_of, i = ctl)
    for (nm in names(X_all)) {
      X <- X_all[[nm]]
      m <- ctl[match_nn(X[tr, , drop = FALSE], X[ctl, , drop = FALSE])]
      Xb <- X_all$effort   # balance always on habitat, psi1 and effort
      smd_after <- sapply(seq_len(ncol(Xb)), function(j) smd(Xb[tr, j], Xb[m, j]))
      d <- r <- numeric(N_DRAWS)
      for (k in seq_len(N_DRAWS)) {
        b <- sample(seq_along(tr), replace = TRUE)
        r_t <- rate_of(k, tr[b]); r_c <- rate_of(k, m[b])
        d[k] <- r_t - r_c; r[k] <- r_t / r_c
      }
      ok <- is.finite(d); okr <- is.finite(r)
      thr_rows[[paste(sp, dn, nm)]] <- data.frame(
        species = sp, variant = variant, point_only = st$point_only,
        definition = dn, spa_in_min = thr[1], spa_out_max = thr[2], spec = nm,
        n_treated = length(tr), n_control_pool = length(ctl),
        n_unique_controls = length(unique(m)), n_excluded = nrow(cl) - length(tr) - length(ctl),
        obs_ext_treated = sum(cl$obs_ext[tr]), obs_ext_control = sum(cl$obs_ext[ctl]),
        rate_in_med = median(raw_in), rate_in_lo = q(raw_in, .025), rate_in_hi = q(raw_in, .975),
        rate_out_med = median(raw_out), rate_out_lo = q(raw_out, .025), rate_out_hi = q(raw_out, .975),
        n_undefined = sum(!ok), diff_med = median(d[ok]), diff_lo = q(d[ok], .025),
        diff_hi = q(d[ok], .975), p_spa_higher = mean(d[ok] > 0),
        n_ratio_undefined = sum(!okr), ratio_med = median(r[okr]),
        ratio_lo = q(r[okr], .025), ratio_hi = q(r[okr], .975),
        max_abs_smd_after = max(abs(smd_after)))
    }
  }
}
thr_df <- do.call(rbind, thr_rows); rownames(thr_df) <- NULL
write.csv(thr_df, file.path(out_dir, "spa_threshold.csv"), row.names = FALSE)

options(width = 220)
print(contrast_df, digits = 3)
print(thr_df, digits = 3)
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
