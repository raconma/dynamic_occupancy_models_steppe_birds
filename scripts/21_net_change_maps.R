###############################################################################
# 21_net_change_maps.R
#
# Where are local extinctions and colonisations happening?
#
# For each 5x5 km cell, the posterior expected number of local extinctions
# and colonisations over 2017-2023 is computed from the fitted colext model
# conditional on the cell's detection history (two-state HMM forward-
# backward; R/colext_utils.R). Net change = colonisations - extinctions.
# The "recent" period drops the 2017->18 transition, where projected change
# is dominated by psi1 relaxing towards the dynamic equilibrium (script 20).
# Where a cell-year is unsurveyed the posterior follows the model's prior
# dynamics, so events are also reported for (i) cells with at least one
# detection and (ii) transitions with surveys in both years ("both"), the
# subset anchored in data. Most events in never-detected cells are the decay
# of occupancy the model assumes but the data never showed.
# One-step-ahead risk = P(occupied in 2023 | data) x epsilon(2023 covariates).
# Uncertainty: 1,000 multivariate-normal coefficient draws.
#
# Validation / interpretation layers (all independent of eBird):
#   core vs rest   top 10% of cells by psi1 vs remaining cells
#   Atlas II       10x10 km presence 1998-2002 (aves_spain.shp, SEO98 grid)
#   SPA            Natura 2000 Special Protection Areas (MITECO 2023);
#                  a cell is "in SPA" if its centroid falls inside an SPA
#   stPGOcc        spatio-temporal model event rates (script 7 rasters)
# plus the location of observed (naive) colonisation/extinction events.
#
# Inputs:
#   results/{sp}_model_object.rds
#   data/processed_2023/{sp}/{sp}_occ_wide_dynamic.csv  (cell coordinates;
#       same rows and order as the model sites, checked below)
#   data-raw/data/validation/atlas_biodiversidad/atlas_biodiversidad/aves_spain.shp
#   data-raw/data/threats/rn2000/PS_NATURA2000_2023_PB.shp
#   results/stpgocc_{col,ext}_hotspots_{sp}.tif
#
# Outputs (new):
#   results/net_change_maps/{sp}_cells.csv
#   results/net_change_maps/summary_by_class.csv
#   results/net_change_maps/validation.csv
#   results/net_change_maps/share_ext_never_detected.csv  (per-draw share of
#       whole-grid extinctions in never-recorded cells, 2017-23 and 2018-23)
#   results/net_change_maps/fig_net_change_recent.png
#   results/net_change_maps/fig_p_net_loss_recent.png
#   results/net_change_maps/fig_risk_2024.png
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
atlas_col <- c(otitar = "OTITAR", ptealc = "PTEALC", pteori = "PTEORI",
               tettet = "TETTET")
N_DRAWS <- 1000
out_dir <- here("results", "net_change_maps")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

q <- function(x, p) as.numeric(quantile(x, p, na.rm = TRUE))

###############################################################################
# SHARED SPATIAL LAYERS
###############################################################################

atlas2 <- load_atlas2()
spa <- load_spa_union()
spa_plot <- st_simplify(spa, dTolerance = 500)

spain <- NULL
if (requireNamespace("rnaturalearth", quietly = TRUE)) {
  spain <- st_transform(rnaturalearth::ne_countries(
    scale = 50, country = "Spain", returnclass = "sf"), UTM30)
}

###############################################################################
# MAIN LOOP
###############################################################################

class_rows <- list()
valid_rows <- list()
share_rows <- list()
cell_tabs  <- list()

for (sp in species) {
  message("\n=== ", sp, " ===")
  mod <- readRDS(here("results", paste0(sp, "_model_object.rds")))
  D <- colext_design(mod)
  n_t <- D$n_t

  # --- Cell coordinates (verify row order against the model's y) -----------
  occ <- read.csv(here("data", "processed_2023", sp,
                       paste0(sp, "_occ_wide_dynamic.csv")))[D$keep, ]
  y_csv <- sapply(occ[, grep("^y\\.[0-9]+\\.[0-9]{4}$", names(occ))],
                  function(v) as.numeric(as.logical(v)))
  ok <- !is.na(D$y)
  stopifnot(all(y_csv[ok] == D$y[ok]), identical(is.na(y_csv), is.na(D$y)))
  pts <- st_transform(st_as_sf(occ[, c("cells", "longitude", "latitude")],
                               coords = c("longitude", "latitude"), crs = 4326),
                      UTM30)
  xy <- st_coordinates(pts)

  in_spa <- lengths(st_intersects(pts, spa)) > 0
  a2 <- st_join(pts, atlas2[, atlas_col[sp]], join = st_within)
  atlas_pres <- as.numeric(st_drop_geometry(a2)[[atlas_col[sp]]]) > 0
  atlas_pres[is.na(atlas_pres)] <- FALSE

  stp <- lapply(c(col = "col", ext = "ext"), function(k) {
    f <- here("results", sprintf("stpgocc_%s_hotspots_%s.tif", k, sp))
    if (!file.exists(f)) return(rep(NA_real_, D$M))
    r <- rast(f)
    v <- st_transform(pts, crs(r))
    terra::extract(r, vect(v))[, 2]
  })

  # --- Point estimate --------------------------------------------------------
  cf <- coef(mod)
  r0 <- colext_rates(cf, D)
  core <- r0$psi1 >= quantile(r0$psi1, 0.9)
  recent <- 2:(n_t - 1)
  prior_psi <- project_occ(r0$psi1, r0$gam, r0$eps, n_t)
  prior_col <- rowSums((1 - prior_psi[, 1:(n_t - 1)]) * r0$gam[, 1:(n_t - 1)])
  prior_ext <- rowSums(prior_psi[, 1:(n_t - 1)] * r0$eps[, 1:(n_t - 1)])

  # Naive (observed) events per cell
  naive <- sapply(1:n_t, function(t) {
    yt <- D$y[, ((t - 1) * D$J + 1):(t * D$J), drop = FALSE]
    ifelse(rowSums(!is.na(yt)) > 0, as.integer(rowSums(yt, na.rm = TRUE) > 0),
           NA_integer_)
  })
  z0 <- naive[, 1:(n_t - 1)]; z1n <- naive[, 2:n_t]
  obs_col <- rowSums(z0 == 0 & z1n == 1, na.rm = TRUE)
  obs_ext <- rowSums(z0 == 1 & z1n == 0, na.rm = TRUE)
  surv <- !is.na(naive)
  both <- surv[, 1:(n_t - 1)] & surv[, 2:n_t]
  det <- rowSums(naive, na.rm = TRUE) > 0

  # --- Coefficient draws -----------------------------------------------------
  draws <- draw_coefs(mod, N_DRAWS)
  # Point-only fits (non-finite variances) give zero-width intervals; flagged
  # in the outputs rather than stopping (P. alchata, extinction not identifiable)
  point_only <- attr(draws, "point_only")
  if (point_only) warning(sp, ": coefficients fixed at point estimates; intervals not valid")

  ext_full <- col_full <- net_rec <- net_both <- risk <- z1_last <-
    matrix(NA_real_, D$M, N_DRAWS)
  classes <- list(all = rep(TRUE, D$M), never_det = !det, det = det,
                  det_core = det & core, det_rest = det & !core,
                  det_spa_in = det & in_spa, det_spa_out = det & !in_spa,
                  det_atlas2_in = det & atlas_pres,
                  det_atlas2_out = det & !atlas_pres,
                  core = core, rest = !core, spa_in = in_spa, spa_out = !in_spa)
  tot <- array(NA_real_, c(N_DRAWS, length(classes), 6),
               dimnames = list(NULL, names(classes),
                               c("ext_full", "col_full", "ext_rec", "col_rec",
                                 "ext_both", "col_both")))

  for (k in seq_len(N_DRAWS)) {
    r <- colext_rates(draws[k, ], D)
    fb <- forward_backward(r, D)
    ext_full[, k] <- rowSums(fb$ext)
    col_full[, k] <- rowSums(fb$col)
    e_rec <- rowSums(fb$ext[, recent, drop = FALSE])
    c_rec <- rowSums(fb$col[, recent, drop = FALSE])
    net_rec[, k] <- c_rec - e_rec
    e_both <- rowSums((fb$ext * both)[, recent, drop = FALSE])
    c_both <- rowSums((fb$col * both)[, recent, drop = FALSE])
    net_both[, k] <- c_both - e_both
    z1_last[, k] <- fb$z1[, n_t]
    risk[, k] <- fb$z1[, n_t] * r$eps[, n_t]
    for (j in seq_along(classes)) {
      w <- classes[[j]]
      tot[k, j, ] <- c(sum(ext_full[w, k]), sum(col_full[w, k]),
                       sum(e_rec[w]), sum(c_rec[w]),
                       sum(e_both[w]), sum(c_both[w]))
    }
  }

  # --- Per-cell table --------------------------------------------------------
  cell_tabs[[sp]] <- data.frame(
    species = sp, cells = occ$cells, x = xy[, 1], y = xy[, 2],
    psi1 = r0$psi1, core = core, in_spa = in_spa, atlas2 = atlas_pres,
    ever_detected = det, n_years_surveyed = rowSums(surv),
    z1_2023 = rowMeans(z1_last),
    ext_full = rowMeans(ext_full), col_full = rowMeans(col_full),
    net_recent = rowMeans(net_rec),
    net_recent_lo = apply(net_rec, 1, q, 0.025),
    net_recent_hi = apply(net_rec, 1, q, 0.975),
    p_net_loss_recent = rowMeans(net_rec < 0),
    net_both = rowMeans(net_both),
    p_net_loss_both = rowMeans(net_both < 0),
    risk_2024 = rowMeans(risk),
    risk_2024_lo = apply(risk, 1, q, 0.025),
    risk_2024_hi = apply(risk, 1, q, 0.975),
    prior_col = prior_col, prior_ext = prior_ext,
    obs_col = obs_col, obs_ext = obs_ext,
    stpgocc_col = stp$col, stpgocc_ext = stp$ext
  )
  write.csv(cell_tabs[[sp]], file.path(out_dir, paste0(sp, "_cells.csv")),
            row.names = FALSE)

  # --- Totals by class -------------------------------------------------------
  for (j in names(classes)) {
    for (per in c("full", "rec", "both")) {
      e <- tot[, j, paste0("ext_", per)]; c_ <- tot[, j, paste0("col_", per)]
      class_rows[[paste(sp, j, per)]] <- data.frame(
        species = sp, class = j, period = per, n_cells = sum(classes[[j]]),
        ext_med = median(e), ext_lo = q(e, .025), ext_hi = q(e, .975),
        col_med = median(c_), col_lo = q(c_, .025), col_hi = q(c_, .975),
        net_med = median(c_ - e), net_lo = q(c_ - e, .025), net_hi = q(c_ - e, .975),
        p_net_loss = mean(c_ - e < 0))
    }
  }

  # --- Share of whole-grid extinctions in never-recorded cells, per draw -----
  # "full" = 2017-2023 (includes the 2017->18 transition, mostly relaxation of
  # psi1); "rec" = 2018-2023, the period analysed in the paper
  for (per in c("full", "rec")) {
    sh <- tot[, "never_det", paste0("ext_", per)] / tot[, "all", paste0("ext_", per)]
    share_rows[[paste(sp, per)]] <- data.frame(
      species = sp, period = per, point_only = point_only, share_med = median(sh),
      share_lo = q(sh, .025), share_hi = q(sh, .975))
  }

  # --- Validation summaries --------------------------------------------------
  ct <- cell_tabs[[sp]]
  share <- function(v, w) sum(v[w]) / sum(v)
  valid_rows[[sp]] <- data.frame(
    species = sp,
    # Share of posterior events in cells where the species was never detected
    share_ext_never_detected = share(ct$ext_full, !det),
    share_col_never_detected = share(ct$col_full, !det),
    # Where does the model expect vs the data support colonisations?
    share_col_rest_prior = share(ct$prior_col, !core),
    share_col_rest_posterior = share(ct$col_full, !core),
    share_col_rest_observed = share(ct$obs_col, !core),
    share_ext_rest_posterior = share(ct$ext_full, !core),
    share_ext_rest_observed = share(ct$obs_ext, !core),
    # Historical range (Atlas II, 1998-2002)
    share_ext_in_atlas2 = share(ct$ext_full, atlas_pres),
    share_col_outside_atlas2 = share(ct$col_full, !atlas_pres),
    share_obs_col_outside_atlas2 = share(ct$obs_col, !atlas_pres),
    share_cells_atlas2 = mean(atlas_pres),
    # Concordance with the spatio-temporal model (cells ever detected)
    n_stpgocc = sum(!is.na(ct$stpgocc_ext)),
    rho_ext_stpgocc_det = suppressWarnings(cor(ct$ext_full[det], ct$stpgocc_ext[det],
                                               method = "spearman", use = "complete.obs")),
    rho_col_stpgocc_det = suppressWarnings(cor(ct$col_full[det], ct$stpgocc_col[det],
                                               method = "spearman", use = "complete.obs")),
    rho_ext_stpgocc_all = suppressWarnings(cor(ct$ext_full, ct$stpgocc_ext,
                                               method = "spearman", use = "complete.obs"))
  )
}

class_df <- do.call(rbind, class_rows); rownames(class_df) <- NULL
valid_df <- do.call(rbind, valid_rows); rownames(valid_df) <- NULL
write.csv(class_df, file.path(out_dir, "summary_by_class.csv"), row.names = FALSE)
write.csv(valid_df, file.path(out_dir, "validation.csv"), row.names = FALSE)
write.csv(do.call(rbind, share_rows), file.path(out_dir, "share_ext_never_detected.csv"),
          row.names = FALSE)
options(width = 220)
print(class_df, digits = 3)
print(valid_df, digits = 3)

###############################################################################
# MAPS
###############################################################################

cells_all <- do.call(rbind, cell_tabs)
cells_all$species <- factor(species_names[cells_all$species],
                            levels = species_names)

# Light grey tiles = all modelled (surveyed) cells, for coverage
base_map <- function(fill_var, data, fill_scale, title) {
  p <- ggplot()
  if (!is.null(spain)) p <- p + geom_sf(data = spain, fill = "grey97",
                                        colour = "grey60", linewidth = 0.2)
  p + geom_tile(data = cells_all, aes(x = x, y = y), fill = "grey88",
                width = 5000, height = 5000) +
    geom_tile(data = data, aes(x = x, y = y, fill = .data[[fill_var]]),
                width = 5000, height = 5000) +
    geom_sf(data = spa_plot, fill = NA, colour = "grey25", linewidth = 0.15) +
    fill_scale +
    coord_sf(crs = UTM30, xlim = c(-60000, 1130000), ylim = c(3980000, 4880000),
             expand = FALSE) +
    facet_wrap(~ species) +
    labs(title = title, x = NULL, y = NULL) +
    theme_bw(base_size = 9) +
    theme(strip.text = element_text(face = "italic"),
          axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank(), legend.position = "bottom")
}

detected <- cells_all[cells_all$ever_detected, ]
lim <- q(abs(detected$net_recent), 0.99)
p_net <- base_map("net_recent", detected,
                  scale_fill_gradient2(low = "#B2182B", mid = "white",
                                       high = "#2166AC", midpoint = 0,
                                       limits = c(-lim, lim),
                                       oob = scales::squish,
                                       name = "Expected colonisations - extinctions, 2018-2023"),
                  "Posterior net change in cells with >= 1 detection (grey: surveyed, never detected; outlines: SPAs)")
ggsave(file.path(out_dir, "fig_net_change_recent.png"), p_net,
       width = 8, height = 7, dpi = 300)

p_loss <- base_map("p_net_loss_recent", detected,
                   scale_fill_viridis_c(option = "rocket", direction = -1,
                                        limits = c(0, 1),
                                        name = "P(net loss 2018-2023)"),
                   "Probability of net loss 2018-2023 (cells with >= 1 detection)")
ggsave(file.path(out_dir, "fig_p_net_loss_recent.png"), p_loss,
       width = 8, height = 7, dpi = 300)

plaus <- cells_all[cells_all$z1_2023 > 0.05 | cells_all$ever_detected, ]
p_risk <- base_map("risk_2024", plaus,
                   scale_fill_viridis_c(option = "magma", direction = -1,
                                        trans = "sqrt",
                                        name = "P(occupied 2023) x extinction"),
                   "One-step-ahead extinction risk (cells with >= 1 detection or P(occupied 2023) > 0.05)")
ggsave(file.path(out_dir, "fig_risk_2024.png"), p_risk,
       width = 8, height = 7, dpi = 300)

sessionInfo()
