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
# SPA effect: cells with >= 1 detection and >= 1 anchored recent transition
# (the only cells that inform the anchored contrasts) mostly inside SPAs
# (cover >= 0.5) are matched with replacement to such cells essentially
# outside (cover < 0.1) by
# nearest Mahalanobis distance on strict steppe fraction and logit(psi1)
# (sensitivity: + UTM coordinates; + checklists per surveyed year). The contrast is the difference in
# extinction rate (SPA - matched controls). Uncertainty: 1,000 draws that
# combine coefficient uncertainty (multivariate normal) with resampling of
# 50 x 50 km spatial blocks of treated cells (with their matched controls),
# because neighbouring cells share birds and management; resampling of
# single cells is kept as a secondary summary (*_cell columns). Ratio
# summaries keep infinite draws (controls without expected extinctions) and
# are set to NA when more than 5 % of draws are undefined (0/0).
#
# Cell cover (pseudo-steppe, SPA) is computed exactly per cell polygon
# (R/spatial_layers.R: cover_by_raster, cover_by_polygons); an earlier
# version used terra::project(method = "average"), which shifted the
# averaging window and misstated cover by up to ~0.1.
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

# Model variant: "primary" (default) or a robustness refit (scripts 28 and 41),
# e.g. Rscript scripts/23_core_spa_analysis.R trend
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
  message("Computing exact CORINE pseudo-steppe and SPA cover for ",
          length(all_cells), " cells ...")
  polys <- wc_cell_polygons(all_cells)
  clc <- rast(here("data-raw", "data", "corine", "Results",
                   "u2018_clc2018_v2020_20u1_raster100m", "DATA",
                   "U2018_CLC2018_V2020_20u1.tif"))
  # Raster values 1..44 are CLC classes (value 48 = NoData, treated as
  # missing); lookup as in build_stepRep.R: 211 = 12, 231 = 18, 242 = 20,
  # 243 = 21, 244 = 22, 321 = 26, 333 = 32
  cv <- cover_by_raster(polys, clc,
                        list(strict = c(12L, 18L, 26L, 32L),
                             broad = c(12L, 18L, 26L, 32L, 20L, 21L, 22L)),
                        na_values = 48L)
  xy <- st_coordinates(st_centroid(polys))
  write.csv(data.frame(cells = all_cells, x = xy[, 1], y = xy[, 2],
                       area_km2 = as.numeric(st_area(polys)) / 1e6,
                       steppe_strict = cv$strict, steppe_broad = cv$broad,
                       spa_frac = cover_by_polygons(polys, load_spa_sites())),
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
# 50 x 50 km spatial blocks (UTM 30N) for cluster resampling
block_of <- function(x, y) paste(floor(x / 5e4), floor(y / 5e4))
# Positions of `ids` drawn by resampling their blocks with replacement
block_sample <- function(blk) {
  ub <- unique(blk)
  unlist(lapply(sample(ub, length(ub), replace = TRUE), function(u) which(blk == u)))
}
# Summaries: differences drop undefined draws (no expected occupancy);
# ratios keep infinite draws and are NA when > 5 % of draws are 0/0
summ_diff <- function(d) {
  ok <- is.finite(d)
  c(n_undefined = sum(!ok), med = median(d[ok]), lo = q(d[ok], .025),
    hi = q(d[ok], .975), p_pos = mean(d[ok] > 0))
}
summ_ratio <- function(r) {
  nan <- is.na(r)
  if (mean(nan) > 0.05) return(c(n_undefined = sum(nan), n_inf = sum(is.infinite(r)),
                                 med = NA, lo = NA, hi = NA))
  rr <- r[!nan]
  qq <- function(p) as.numeric(quantile(rr, p, type = 1))
  c(n_undefined = sum(nan), n_inf = sum(is.infinite(rr)), med = qq(.5), lo = qq(.025),
    hi = qq(.975))
}

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
  # The pooled-history refit (script 38) has its own histories on the same cells
  if (variant != "pooled") stopifnot(all(y_csv[!is.na(D$y)] == D$y[!is.na(D$y)]))
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
  # Cells that inform the anchored contrasts: >= 1 transition surveyed in
  # both years within the recent period
  anch <- det & rowSums(both[, recent, drop = FALSE]) > 0
  blk <- block_of(h$x, h$y)

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
  tr <- which(classes$spa_in & anch); ctl <- which(classes$spa_out & anch)
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
  match_d_cell <- match_d   # single-cell resampling (secondary)
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
      # Resample spatial blocks of treated cells (with their matched
      # controls); single-cell resampling kept as a secondary summary
      bb <- block_sample(blk[tr]); bc <- sample(seq_along(tr), replace = TRUE)
      for (nm in names(covs)) {
        rate <- function(i) sum(cs$E[i]) / sum(cs$O[i])
        r_t <- rate(tr[bb]); r_c <- rate(matches[[nm]][bb])
        match_d[k, nm, p] <- r_t - r_c
        match_r[k, nm, p] <- r_t / r_c
        match_d_cell[k, nm, p] <- rate(tr[bc]) - rate(matches[[nm]][bc])
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
                      obs_ext = obs_ext, anch = anch, blk = blk)[det_idx, ],
    groups = lapply(pairs, function(g) list(
      core = (classes[[paste0("core_", g)]] & anch)[det_idx],
      rest = (classes[[paste0("rest_", g)]] & anch)[det_idx])))

  for (j in names(classes)) for (p in per) {
    e <- rate_d[, j, p, "ext_rate"]; n <- rate_d[, j, p, "net"]
    rate_rows[[paste(sp, j, p)]] <- data.frame(
      species = sp, variant = variant, point_only = point_only,
      class = j, period = p, n_cells = sum(classes[[j]]),
      n_anchored = sum(classes[[j]] & anch),
      ext_rate_med = median(e), ext_rate_lo = q(e, .025), ext_rate_hi = q(e, .975),
      net_med = median(n), net_lo = q(n, .025), net_hi = q(n, .975),
      p_net_loss = mean(n < 0))
  }
  for (nm in names(covs)) for (p in per) {
    sd <- summ_diff(match_d[, nm, p]); sr <- summ_ratio(match_r[, nm, p])
    sc <- summ_diff(match_d_cell[, nm, p])
    match_rows[[paste(sp, nm, p)]] <- data.frame(
      species = sp, variant = variant, point_only = point_only,
      spec = nm, period = p, n_treated = length(tr),
      n_treated_recorded = sum(classes$spa_in), n_control_pool = length(ctl),
      n_unique_controls = length(unique(matches[[nm]])),
      n_blocks_treated = length(unique(blk[tr])),
      n_undefined = sd[["n_undefined"]],
      diff_ext_rate_med = sd[["med"]], diff_lo = sd[["lo"]], diff_hi = sd[["hi"]],
      p_spa_higher = sd[["p_pos"]],
      n_ratio_undefined = sr[["n_undefined"]], n_ratio_inf = sr[["n_inf"]],
      ratio_med = sr[["med"]], ratio_lo = sr[["lo"]], ratio_hi = sr[["hi"]],
      diff_cell_med = sc[["med"]], diff_cell_lo = sc[["lo"]], diff_cell_hi = sc[["hi"]])
  }
}

rate_df  <- do.call(rbind, rate_rows);  rownames(rate_df)  <- NULL
match_df <- do.call(rbind, match_rows); rownames(match_df) <- NULL
bal_df   <- do.call(rbind, bal_rows);   rownames(bal_df)   <- NULL
write.csv(rate_df,  file.path(out_dir, "class_rates.csv"), row.names = FALSE)
write.csv(match_df, file.path(out_dir, "spa_matched.csv"), row.names = FALSE)
write.csv(bal_df,   file.path(out_dir, "spa_balance.csv"), row.names = FALSE)

# Core vs marginal contrast with its own uncertainty: per coefficient draw,
# resample 50 x 50 km blocks of the anchored recorded cells (both habitat
# classes together); resampling cells within each class is kept as a
# secondary summary. Run after all species so that the random stream of the
# outputs above is unchanged.
contrast_rows <- list()
for (sp in names(contrast_store)) {
  st <- contrast_store[[sp]]
  for (g in names(st$groups)) for (p in dimnames(st$E)[[3]]) {
    ci <- which(st$groups[[g]]$core); ri <- which(st$groups[[g]]$rest)
    if (length(ci) < 2 || length(ri) < 2) next
    cr <- c(ci, ri); is_core <- rep(c(TRUE, FALSE), c(length(ci), length(ri)))
    d <- r <- dc <- numeric(N_DRAWS)
    for (k in seq_len(N_DRAWS)) {
      rate <- function(i) sum(st$E[k, i, p]) / sum(st$O[k, i, p])
      b <- block_sample(st$cell$blk[cr])
      rc <- rate(cr[b][is_core[b]]); rr <- rate(cr[b][!is_core[b]])
      d[k] <- rr - rc; r[k] <- rr / rc
      dc[k] <- rate(sample(ri, replace = TRUE)) - rate(sample(ci, replace = TRUE))
    }
    sd <- summ_diff(d); sr <- summ_ratio(r); sc <- summ_diff(dc)
    contrast_rows[[paste(sp, g, p)]] <- data.frame(
      species = sp, variant = variant, point_only = st$point_only,
      definition = g, period = p, n_core = length(ci), n_marginal = length(ri),
      n_undefined = sd[["n_undefined"]],
      diff_med = sd[["med"]], diff_lo = sd[["lo"]], diff_hi = sd[["hi"]],
      n_ratio_undefined = sr[["n_undefined"]], n_ratio_inf = sr[["n_inf"]],
      ratio_med = sr[["med"]], ratio_lo = sr[["lo"]], ratio_hi = sr[["hi"]],
      p_marginal_higher = sd[["p_pos"]],
      diff_cell_med = sc[["med"]], diff_cell_lo = sc[["lo"]], diff_cell_hi = sc[["hi"]])
  }
}
contrast_df <- do.call(rbind, contrast_rows); rownames(contrast_df) <- NULL
write.csv(contrast_df, file.path(out_dir, "core_contrast.csv"), row.names = FALSE)

# Sensitivity of the SPA contrast to the SPA cover thresholds (anchored
# transitions, 2018-2023). The main analysis compares cells >= 50 % inside
# SPAs with cells < 10 % inside and leaves out intermediate cover; here the
# treated/control split is varied, among cells with anchored transitions.
# Matching and block resampling as above. Run
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
    tr <- which(cl$anch & cl$spa_frac >= thr[1]); ctl <- which(cl$anch & cl$spa_frac < thr[2])
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
        b <- block_sample(cl$blk[tr])
        r_t <- rate_of(k, tr[b]); r_c <- rate_of(k, m[b])
        d[k] <- r_t - r_c; r[k] <- r_t / r_c
      }
      sd <- summ_diff(d); sr <- summ_ratio(r)
      thr_rows[[paste(sp, dn, nm)]] <- data.frame(
        species = sp, variant = variant, point_only = st$point_only,
        definition = dn, spa_in_min = thr[1], spa_out_max = thr[2], spec = nm,
        n_treated = length(tr), n_control_pool = length(ctl),
        n_unique_controls = length(unique(m)),
        n_excluded = sum(cl$anch) - length(tr) - length(ctl),
        obs_ext_treated = sum(cl$obs_ext[tr]), obs_ext_control = sum(cl$obs_ext[ctl]),
        rate_in_med = median(raw_in), rate_in_lo = q(raw_in, .025), rate_in_hi = q(raw_in, .975),
        rate_out_med = median(raw_out), rate_out_lo = q(raw_out, .025), rate_out_hi = q(raw_out, .975),
        n_undefined = sd[["n_undefined"]], diff_med = sd[["med"]], diff_lo = sd[["lo"]],
        diff_hi = sd[["hi"]], p_spa_higher = sd[["p_pos"]],
        n_ratio_undefined = sr[["n_undefined"]], n_ratio_inf = sr[["n_inf"]],
        ratio_med = sr[["med"]], ratio_lo = sr[["lo"]], ratio_hi = sr[["hi"]],
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
