###############################################################################
# 7_spatial_results.R
#
# Purpose: Post-process Bayesian spatial occupancy model results (spOccupancy)
#          and generate publication-ready figures and tables for all 4 species.
#
# Inputs:
#   results/results_spatial/{sp}_spatial_tPGOcc.rds
#   results/results_spatial/{sp}_spatial_stPGOcc.rds
#   results/results_spatial/{sp}_spatial_spPGOcc_2023.rds
#   results/results_spatial/{sp}_mcmc_diagnostics.csv
#   results/results_spatial/{sp}_model_comparison.csv
#   data/processed_2023/{sp}/{sp}_occ_wide_dynamic.csv
#   results/pub_table2_coefficients_psi.csv
#   R/model_configs.R
#
# Outputs (figs/):
#   pub_fig_spatial_w_maps.png              - Spatial random effect maps
#   pub_fig_spatial_forest_betas.png        - Coefficient forest plots
#   pub_fig_spatial_range.png               - Spatial range comparison
#   pub_fig_spatial_waic.png                - WAIC model comparison
#   pub_fig_spatial_moran.png               - Moran's I before/after
#   pub_fig_spatial_colext_comparison.png   - colext vs spOccupancy betas
#   pub_fig_spatial_occupancy_trends.png    - Occupancy trends 2017-2023
#   pub_fig_spatial_col_ext_rates.png       - Derived colonization/extinction
#   pub_fig_spatial_col_ext_maps.png        - Col/ext hotspot maps
#   pub_fig_spatial_delta_psi_maps.png      - Occupancy change maps
#
# Outputs (results/):
#   pub_spatial_table1_model_comparison.csv
#   pub_spatial_table2_stpgocc_coefficients.csv
#
# Requirements: spOccupancy, sf, ggplot2, dplyr, tidyr, gridExtra,
#               rnaturalearth, rnaturalearthdata, spdep
###############################################################################

library(here)
library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(rnaturalearth)
library(spdep)

if (!requireNamespace("spOccupancy", quietly = TRUE)) {
  stop("Install spOccupancy first: install.packages('spOccupancy')")
}
library(spOccupancy)

source(here("R", "model_configs.R"))

set.seed(42)

cat("\n")
cat(strrep("=", 70), "\n")
cat("  SPATIAL OCCUPANCY RESULTS — POST-PROCESSING & VISUALIZATION\n")
cat(strrep("=", 70), "\n\n")

###############################################################################
# 1. CONFIGURATION
###############################################################################

species_codes <- c("otitar", "ptealc", "pteori", "tettet")
sp_latin <- c(otitar = "Otis tarda", ptealc = "Pterocles alchata",
              pteori = "Pterocles orientalis", tettet = "Tetrax tetrax")
sp_short <- c(otitar = "O. tarda", ptealc = "P. alchata",
              pteori = "P. orientalis", tettet = "T. tetrax")
sp_colors <- c("O. tarda" = "#1B9E77", "P. alchata" = "#D95F02",
               "P. orientalis" = "#7570B3", "T. tetrax" = "#E7298A")
panel_labels <- c("otitar" = "A", "ptealc" = "B",
                   "pteori" = "C", "tettet" = "D")

years <- 2017:2023
n_years <- length(years)
n_reps  <- 10

figs_dir    <- here("figs")
results_dir <- here("results")
spatial_dir <- here("results", "results_spatial")
if (!dir.exists(figs_dir)) dir.create(figs_dir, recursive = TRUE)

# Publication theme (matches publication_figures_models.R)
theme_pub <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(face = "italic"),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    axis.title = element_text(size = 10)
  )

# Spain outline for maps
spain <- ne_countries(scale = "medium", country = "Spain", returnclass = "sf")


###############################################################################
# 2. HELPER FUNCTIONS
###############################################################################

#' Extract posterior summary from spOccupancy model
extract_posterior <- function(fit, param_type) {
  samps <- fit[[paste0(param_type, ".samples")]]
  if (is.null(samps)) return(NULL)
  if (inherits(samps, "mcmc.list")) samps <- do.call(rbind, samps)
  if (!is.matrix(samps)) samps <- matrix(samps, ncol = 1)
  pnames <- colnames(samps)
  if (is.null(pnames)) pnames <- paste0(param_type, "_", seq_len(ncol(samps)))
  data.frame(
    Parameter = pnames,
    Mean  = colMeans(samps),
    SD    = apply(samps, 2, sd),
    Q2.5  = apply(samps, 2, quantile, 0.025),
    Q50   = apply(samps, 2, quantile, 0.50),
    Q97.5 = apply(samps, 2, quantile, 0.975),
    Rhat  = if (!is.null(fit$rhat[[param_type]])) fit$rhat[[param_type]] else NA,
    ESS   = if (!is.null(fit$ESS[[param_type]])) fit$ESS[[param_type]] else NA,
    row.names = NULL, stringsAsFactors = FALSE
  )
}

#' Compute effective spatial range from phi posterior
compute_eff_range <- function(fit) {
  if (is.null(fit$theta.samples)) return(list(median = NA, lo = NA, hi = NA))
  samps <- fit$theta.samples
  if (inherits(samps, "mcmc.list")) samps <- do.call(rbind, samps)
  if (!"phi" %in% colnames(samps)) return(list(median = NA, lo = NA, hi = NA))
  phi <- samps[, "phi"]
  range_km <- 3 / phi / 1000
  list(
    median = median(range_km),
    lo = quantile(range_km, 0.025, names = FALSE),
    hi = quantile(range_km, 0.975, names = FALSE)
  )
}

#' Load species coordinates (replicating NA-filtering from script 6)
load_species_data <- function(sp) {
  cfg <- get_model_config(sp)
  d <- read.csv(here("data", "processed_2023", sp,
                      paste0(sp, "_occ_wide_dynamic.csv")))
  J <- nrow(d)

  # Build detection array
  y_array <- array(NA, dim = c(J, n_years, n_reps))
  for (t in seq_along(years)) {
    for (k in 1:n_reps) {
      col_name <- paste0("y.", k, ".", years[t])
      if (col_name %in% names(d)) y_array[, t, k] <- d[[col_name]]
    }
  }

  # Occupancy covariates (replicate script 6 logic exactly)
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
        if (abs(mean(vals, na.rm = TRUE)) > 10 || sd(vals, na.rm = TRUE) > 10) {
          vals <- as.numeric(scale(vals))
        }
        mat[, t] <- vals
      }
    }
    occ_covs[[var_name]] <- mat
  }

  # NA filtering (matches script 6 lines 156-186)
  na_sites <- rep(FALSE, J)
  for (nm in names(occ_covs)) {
    obj <- occ_covs[[nm]]
    if (is.matrix(obj)) {
      na_sites <- na_sites | apply(obj, 1, function(x) any(is.na(x)))
    } else {
      na_sites <- na_sites | is.na(obj)
    }
  }
  coords_lonlat <- cbind(d$longitude, d$latitude)
  na_sites <- na_sites | is.na(coords_lonlat[, 1]) | is.na(coords_lonlat[, 2])

  keep <- !na_sites
  list(
    coords_lonlat = coords_lonlat[keep, ],
    y_array = y_array[keep, , ],
    J = sum(keep),
    d = d[keep, ]
  )
}


###############################################################################
# 3. DATA LOADING
###############################################################################

cat("Loading model objects and data...\n\n")

fits_st <- list()   # stPGOcc (spatio-temporal)
fits_t  <- list()   # tPGOcc (non-spatial baseline)
fits_sp <- list()   # spPGOcc (single-year spatial)
sp_data <- list()   # coordinates and detection data
comp_df <- list()   # model comparison CSVs
diag_df <- list()   # MCMC diagnostics CSVs

for (sp in species_codes) {
  cat("  ", sp, "...")

  fits_st[[sp]] <- readRDS(file.path(spatial_dir,
                                      paste0(sp, "_spatial_stPGOcc.rds")))
  fits_t[[sp]]  <- readRDS(file.path(spatial_dir,
                                      paste0(sp, "_spatial_tPGOcc.rds")))
  fits_sp[[sp]] <- readRDS(file.path(spatial_dir,
                                      paste0(sp, "_spatial_spPGOcc_2023.rds")))

  sp_data[[sp]] <- load_species_data(sp)

  comp_df[[sp]] <- read.csv(file.path(spatial_dir,
                                       paste0(sp, "_model_comparison.csv")))
  diag_df[[sp]] <- read.csv(file.path(spatial_dir,
                                       paste0(sp, "_mcmc_diagnostics.csv")))
  cat(" OK\n")
}

# Load colext coefficients for comparison
colext_psi <- read.csv(here("results", "pub_table2_coefficients_psi.csv"))

cat("\nAll data loaded successfully.\n\n")


###############################################################################
# 4. TABLES
###############################################################################

cat("Generating tables...\n\n")

# --- Table 1: Model comparison summary ---
table1_rows <- list()
for (sp in species_codes) {
  comp <- comp_df[[sp]]
  t_row <- comp[comp$Model == "tPGOcc", ]
  st_row <- comp[grepl("stPGOcc", comp$Model), ]
  range_info <- compute_eff_range(fits_st[[sp]])

  table1_rows[[sp]] <- data.frame(
    Species = sp_latin[sp],
    tPGOcc_WAIC = t_row$WAIC,
    stPGOcc_WAIC = st_row$WAIC,
    Delta_WAIC = round(t_row$WAIC - st_row$WAIC, 1),
    Moran_I_tPGOcc = t_row$Moran_I,
    Moran_I_stPGOcc = st_row$Moran_I,
    Spatial_Range_km = round(range_info$median, 1),
    Range_CI_low = round(range_info$lo, 1),
    Range_CI_high = round(range_info$hi, 1),
    stringsAsFactors = FALSE
  )
}
table1 <- do.call(rbind, table1_rows)
rownames(table1) <- NULL
table1_path <- file.path(results_dir, "pub_spatial_table1_model_comparison.csv")
write.csv(table1, table1_path, row.names = FALSE)
cat("  Saved:", table1_path, "\n")

# --- Table 2: stPGOcc coefficient summary ---
table2_rows <- list()
for (sp in species_codes) {
  for (ptype in c("beta", "alpha", "theta")) {
    post <- extract_posterior(fits_st[[sp]], ptype)
    if (!is.null(post)) {
      post$Species <- sp_latin[sp]
      post$Type <- ptype
      table2_rows[[paste(sp, ptype)]] <- post
    }
  }
}
table2 <- do.call(rbind, table2_rows)
rownames(table2) <- NULL
table2 <- table2[, c("Species", "Type", "Parameter", "Mean", "SD",
                       "Q2.5", "Q50", "Q97.5", "Rhat", "ESS")]
table2_path <- file.path(results_dir,
                          "pub_spatial_table2_stpgocc_coefficients.csv")
write.csv(table2, table2_path, row.names = FALSE)
cat("  Saved:", table2_path, "\n\n")


###############################################################################
# 5. FIGURES
###############################################################################

cat("Generating figures...\n\n")

# =========================================================================
# FIGURE A: Spatial random effect maps (w)
# =========================================================================
cat("  Figure A: Spatial random effect maps...\n")

w_plots <- list()
for (sp in species_codes) {
  fit <- fits_st[[sp]]
  dat <- sp_data[[sp]]

  # Extract posterior mean of w
  w_samps <- fit$w.samples
  if (inherits(w_samps, "mcmc.list")) w_samps <- do.call(rbind, w_samps)
  w_mean <- colMeans(w_samps)

  # Handle dimension mismatch: w may be for J sites
  n_w <- length(w_mean)
  n_coords <- nrow(dat$coords_lonlat)
  if (n_w != n_coords) {
    cat("    WARNING:", sp, "- w has", n_w, "values but", n_coords,
        "coords. Using min.\n")
    n_use <- min(n_w, n_coords)
    w_mean <- w_mean[1:n_use]
    coords_use <- dat$coords_lonlat[1:n_use, ]
  } else {
    coords_use <- dat$coords_lonlat
  }

  range_info <- compute_eff_range(fit)

  w_df <- data.frame(
    lon = coords_use[, 1], lat = coords_use[, 2], w = w_mean
  )

  # Symmetric color limits
  w_lim <- max(abs(quantile(w_mean, c(0.01, 0.99))))

  w_plots[[sp]] <- ggplot() +
    geom_sf(data = spain, fill = "grey95", colour = "grey60", linewidth = 0.3) +
    geom_point(data = w_df, aes(x = lon, y = lat, colour = w),
               size = 0.3, alpha = 0.7) +
    scale_colour_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                           midpoint = 0, limits = c(-w_lim, w_lim),
                           name = "Spatial\neffect (w)") +
    coord_sf(xlim = c(-10, 5), ylim = c(35, 44)) +
    labs(title = paste0(panel_labels[sp], ") ", sp_short[sp]),
         subtitle = paste0("Effective range: ",
                           round(range_info$median, 0), " km [",
                           round(range_info$lo, 0), "\u2013",
                           round(range_info$hi, 0), "]")) +
    theme_pub +
    theme(axis.title = element_blank(),
          plot.title = element_text(face = "italic"))
}

fig_a <- gridExtra::grid.arrange(
  grobs = w_plots, ncol = 2,
  top = grid::textGrob("Spatial random effect (stPGOcc)",
                        gp = grid::gpar(fontsize = 14, fontface = "bold"))
)
ggsave(file.path(figs_dir, "pub_fig_spatial_w_maps.png"),
       fig_a, width = 12, height = 11, dpi = 300)
cat("    Saved: pub_fig_spatial_w_maps.png\n")


# =========================================================================
# FIGURE B: Coefficient forest plots (tPGOcc vs stPGOcc)
# =========================================================================
cat("  Figure B: Coefficient forest plots...\n")

beta_all <- list()
for (sp in species_codes) {
  for (model_name in c("tPGOcc", "stPGOcc")) {
    fit <- if (model_name == "tPGOcc") fits_t[[sp]] else fits_st[[sp]]
    post <- extract_posterior(fit, "beta")
    if (!is.null(post)) {
      post$Species <- sp_short[sp]
      post$Model <- model_name
      beta_all[[paste(sp, model_name)]] <- post
    }
  }
}
beta_df <- do.call(rbind, beta_all)
rownames(beta_df) <- NULL

# Remove intercept
beta_df <- beta_df[beta_df$Parameter != "(Intercept)", ]

fig_b <- ggplot(beta_df, aes(x = Mean, y = Parameter,
                              colour = Species, shape = Model)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_pointrange(aes(xmin = Q2.5, xmax = Q97.5),
                  position = position_dodge(width = 0.6), size = 0.4) +
  scale_colour_manual(values = sp_colors) +
  scale_shape_manual(values = c("tPGOcc" = 1, "stPGOcc" = 16)) +
  labs(title = "Occupancy coefficients: tPGOcc vs stPGOcc",
       subtitle = "Posterior mean with 95% credible interval",
       x = "Coefficient estimate", y = NULL) +
  theme_pub +
  theme(legend.box = "vertical")

ggsave(file.path(figs_dir, "pub_fig_spatial_forest_betas.png"),
       fig_b, width = 10, height = 8, dpi = 300)
cat("    Saved: pub_fig_spatial_forest_betas.png\n")


# =========================================================================
# FIGURE C: Spatial range comparison
# =========================================================================
cat("  Figure C: Spatial range comparison...\n")

range_data <- data.frame(
  Species = character(), Median = numeric(),
  Lo = numeric(), Hi = numeric(), stringsAsFactors = FALSE
)
for (sp in species_codes) {
  r <- compute_eff_range(fits_st[[sp]])
  range_data <- rbind(range_data, data.frame(
    Species = sp_short[sp], Median = r$median,
    Lo = r$lo, Hi = r$hi, stringsAsFactors = FALSE
  ))
}

# Order by range
range_data$Species <- factor(range_data$Species,
                              levels = range_data$Species[order(range_data$Median)])

fig_c <- ggplot(range_data, aes(x = Species, y = Median, colour = Species)) +
  geom_hline(yintercept = c(50, 100), linetype = "dotted", colour = "grey70") +
  geom_pointrange(aes(ymin = Lo, ymax = Hi), size = 1, linewidth = 0.8) +
  scale_colour_manual(values = sp_colors) +
  labs(title = "Effective spatial range (stPGOcc)",
       subtitle = paste("Dashed lines: 50 km (local/landscape threshold),",
                        "100 km (landscape/regional threshold)"),
       x = NULL, y = "Effective range (km)") +
  theme_pub +
  theme(legend.position = "none",
        axis.text.x = element_text(face = "italic"))

ggsave(file.path(figs_dir, "pub_fig_spatial_range.png"),
       fig_c, width = 8, height = 5, dpi = 300)
cat("    Saved: pub_fig_spatial_range.png\n")


# =========================================================================
# FIGURE D: WAIC comparison
# =========================================================================
cat("  Figure D: WAIC comparison...\n")

waic_data <- data.frame(
  Species = character(), Model = character(),
  WAIC = numeric(), stringsAsFactors = FALSE
)
for (sp in species_codes) {
  comp <- comp_df[[sp]]
  waic_data <- rbind(waic_data,
    data.frame(Species = sp_short[sp], Model = "tPGOcc",
               WAIC = comp$WAIC[comp$Model == "tPGOcc"],
               stringsAsFactors = FALSE),
    data.frame(Species = sp_short[sp], Model = "stPGOcc",
               WAIC = comp$WAIC[grepl("stPGOcc", comp$Model)],
               stringsAsFactors = FALSE)
  )
}

# Compute delta for annotation
delta_waic <- waic_data %>%
  pivot_wider(names_from = Model, values_from = WAIC) %>%
  mutate(Delta = round(tPGOcc - stPGOcc, 0),
         x_pos = (tPGOcc + stPGOcc) / 2)

fig_d <- ggplot(waic_data, aes(x = WAIC, y = Species, colour = Model)) +
  geom_line(aes(group = Species), colour = "grey70", linewidth = 0.5) +
  geom_point(size = 3) +
  geom_text(data = delta_waic,
            aes(x = x_pos, y = Species,
                label = paste0("\u0394 = ", Delta)),
            colour = "grey30", size = 3, vjust = -1, inherit.aes = FALSE) +
  scale_colour_manual(values = c("tPGOcc" = "grey50", "stPGOcc" = "#E6550D")) +
  labs(title = "Model comparison: WAIC",
       subtitle = "Lower WAIC = better predictive performance",
       x = "WAIC", y = NULL) +
  theme_pub +
  theme(axis.text.y = element_text(face = "italic"))

ggsave(file.path(figs_dir, "pub_fig_spatial_waic.png"),
       fig_d, width = 8, height = 5, dpi = 300)
cat("    Saved: pub_fig_spatial_waic.png\n")


# =========================================================================
# FIGURE E: Moran's I before/after
# =========================================================================
cat("  Figure E: Moran's I before/after...\n")

moran_data <- data.frame(
  Species = character(), Model = character(),
  Moran_I = numeric(), stringsAsFactors = FALSE
)
for (sp in species_codes) {
  comp <- comp_df[[sp]]
  moran_data <- rbind(moran_data,
    data.frame(Species = sp_short[sp], Model = "tPGOcc (non-spatial)",
               Moran_I = comp$Moran_I[comp$Model == "tPGOcc"],
               stringsAsFactors = FALSE),
    data.frame(Species = sp_short[sp], Model = "stPGOcc (spatial)",
               Moran_I = comp$Moran_I[grepl("stPGOcc", comp$Model)],
               stringsAsFactors = FALSE)
  )
}

fig_e <- ggplot(moran_data, aes(x = Model, y = Moran_I,
                                 group = Species, colour = Species)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 3) +
  scale_colour_manual(values = sp_colors) +
  labs(title = "Spatial autocorrelation reduction",
       subtitle = "Moran's I on model residuals (0 = no autocorrelation)",
       x = NULL, y = "Moran's I") +
  theme_pub +
  theme(axis.text.x = element_text(size = 9))

ggsave(file.path(figs_dir, "pub_fig_spatial_moran.png"),
       fig_e, width = 8, height = 5, dpi = 300)
cat("    Saved: pub_fig_spatial_moran.png\n")


# =========================================================================
# FIGURE F: colext vs spOccupancy coefficient comparison
# =========================================================================
cat("  Figure F: colext vs spOccupancy comparison...\n")

# Shared covariates present in all species psi formulas
shared_vars <- c("bio2", "tree_cover", "grass_cover")

# colext estimates (frequentist: Estimate +/- 1.96*SE)
colext_shared <- colext_psi %>%
  filter(Parameter %in% shared_vars) %>%
  mutate(Q2.5 = Estimate - 1.96 * SE,
         Q97.5 = Estimate + 1.96 * SE,
         Framework = "colext (frequentist)",
         Species_short = case_when(
           Species == "Otis tarda" ~ "O. tarda",
           Species == "Pterocles alchata" ~ "P. alchata",
           Species == "Pterocles orientalis" ~ "P. orientalis",
           Species == "Tetrax tetrax" ~ "T. tetrax"
         )) %>%
  select(Species = Species_short, Parameter, Mean = Estimate,
         Q2.5, Q97.5, Framework)

# stPGOcc estimates (Bayesian)
spocc_shared <- list()
for (sp in species_codes) {
  post <- extract_posterior(fits_st[[sp]], "beta")
  if (!is.null(post)) {
    post_filt <- post[post$Parameter %in% shared_vars, ]
    if (nrow(post_filt) > 0) {
      post_filt$Species <- sp_short[sp]
      post_filt$Framework <- "stPGOcc (Bayesian)"
      spocc_shared[[sp]] <- post_filt[, c("Species", "Parameter", "Mean",
                                           "Q2.5", "Q97.5", "Framework")]
    }
  }
}
spocc_shared_df <- do.call(rbind, spocc_shared)
rownames(spocc_shared_df) <- NULL

comp_betas <- rbind(colext_shared, spocc_shared_df)

fig_f <- ggplot(comp_betas, aes(x = Mean, y = Species,
                                 colour = Framework, shape = Framework)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_pointrange(aes(xmin = Q2.5, xmax = Q97.5),
                  position = position_dodge(width = 0.5),
                  size = 0.5, linewidth = 0.6) +
  facet_wrap(~ Parameter, scales = "free_x", ncol = 3) +
  scale_colour_manual(values = c("colext (frequentist)" = "#756BB1",
                                  "stPGOcc (Bayesian)" = "#E6550D")) +
  scale_shape_manual(values = c("colext (frequentist)" = 2,
                                 "stPGOcc (Bayesian)" = 16)) +
  labs(title = "Coefficient comparison: colext vs stPGOcc",
       subtitle = paste("Shared occupancy covariates.",
                        "colext = initial occupancy (\u03C8\u2081);",
                        "stPGOcc = multi-year occupancy"),
       x = "Standardised coefficient", y = NULL) +
  theme_pub +
  theme(axis.text.y = element_text(face = "italic"),
        legend.box = "vertical",
        strip.text = element_text(face = "bold"))

ggsave(file.path(figs_dir, "pub_fig_spatial_colext_comparison.png"),
       fig_f, width = 10, height = 7, dpi = 300)
cat("    Saved: pub_fig_spatial_colext_comparison.png\n")


# =========================================================================
# FIGURE G: Occupancy trends from z.samples
# =========================================================================
cat("  Figure G: Occupancy trends...\n")

trend_all <- list()
for (sp in species_codes) {
  fit <- fits_st[[sp]]
  z_samps <- fit$z.samples
  if (inherits(z_samps, "mcmc.list")) z_samps <- do.call(rbind, z_samps)

  # z_samps: [n_samples x (J*T)] or [n_samples x J x T]
  # spOccupancy stores z as [n_samples, J, T]
  if (length(dim(z_samps)) == 3) {
    # Mean occupancy per year per MCMC sample
    psi_year <- apply(z_samps, c(1, 3), mean)  # [n_samples x T]
  } else {
    # Reshape: assume J sites * T years flattened
    J_sp <- sp_data[[sp]]$J
    n_samps <- nrow(z_samps)
    z_arr <- array(z_samps, dim = c(n_samps, J_sp, n_years))
    psi_year <- apply(z_arr, c(1, 3), mean)
  }

  for (t in seq_along(years)) {
    trend_all[[paste(sp, t)]] <- data.frame(
      Species = sp_short[sp],
      Year = years[t],
      Mean = mean(psi_year[, t]),
      Q2.5 = quantile(psi_year[, t], 0.025, names = FALSE),
      Q97.5 = quantile(psi_year[, t], 0.975, names = FALSE),
      stringsAsFactors = FALSE
    )
  }
}
trend_df <- do.call(rbind, trend_all)
rownames(trend_df) <- NULL

fig_g <- ggplot(trend_df, aes(x = Year, y = Mean, colour = Species,
                               fill = Species)) +
  geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), alpha = 0.15,
              colour = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  facet_wrap(~ Species, scales = "free_y", ncol = 2) +
  scale_colour_manual(values = sp_colors) +
  scale_fill_manual(values = sp_colors) +
  scale_x_continuous(breaks = years) +
  labs(title = "Occupancy trends (stPGOcc)",
       subtitle = "Mean proportion of occupied sites with 95% credible interval",
       x = NULL, y = "Proportion of sites occupied") +
  theme_pub +
  theme(legend.position = "none",
        strip.text = element_text(face = "italic"),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(figs_dir, "pub_fig_spatial_occupancy_trends.png"),
       fig_g, width = 10, height = 8, dpi = 300)
cat("    Saved: pub_fig_spatial_occupancy_trends.png\n")


# =========================================================================
# FIGURE H: Derived colonization/extinction rates
# =========================================================================
cat("  Figure H: Derived colonization/extinction rates...\n")

colext_rates_all <- list()
for (sp in species_codes) {
  fit <- fits_st[[sp]]
  z_samps <- fit$z.samples
  if (inherits(z_samps, "mcmc.list")) z_samps <- do.call(rbind, z_samps)

  J_sp <- sp_data[[sp]]$J
  n_samps <- nrow(z_samps)

  if (length(dim(z_samps)) == 3) {
    z_arr <- z_samps
  } else {
    z_arr <- array(z_samps, dim = c(n_samps, J_sp, n_years))
  }

  for (t in 2:n_years) {
    # Colonization: z[t]=1 & z[t-1]=0
    col_events <- z_arr[, , t] * (1 - z_arr[, , t - 1])
    ext_events <- (1 - z_arr[, , t]) * z_arr[, , t - 1]
    empty_prev <- 1 - z_arr[, , t - 1]
    occ_prev   <- z_arr[, , t - 1]

    # Per-sample rates (avoid division by zero)
    gamma_s <- rowSums(col_events) / pmax(rowSums(empty_prev), 1)
    eps_s   <- rowSums(ext_events) / pmax(rowSums(occ_prev), 1)

    colext_rates_all[[paste(sp, "col", t)]] <- data.frame(
      Species = sp_short[sp], Year = years[t], Rate = "Colonization",
      Mean = mean(gamma_s),
      Q2.5 = quantile(gamma_s, 0.025, names = FALSE),
      Q97.5 = quantile(gamma_s, 0.975, names = FALSE),
      stringsAsFactors = FALSE
    )
    colext_rates_all[[paste(sp, "ext", t)]] <- data.frame(
      Species = sp_short[sp], Year = years[t], Rate = "Extinction",
      Mean = mean(eps_s),
      Q2.5 = quantile(eps_s, 0.025, names = FALSE),
      Q97.5 = quantile(eps_s, 0.975, names = FALSE),
      stringsAsFactors = FALSE
    )
  }
}
rates_df <- do.call(rbind, colext_rates_all)
rownames(rates_df) <- NULL

fig_h <- ggplot(rates_df, aes(x = Year, y = Mean, colour = Rate,
                               fill = Rate)) +
  geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  facet_wrap(~ Species, scales = "free_y", ncol = 2) +
  scale_colour_manual(values = c("Colonization" = "#2CA02C",
                                  "Extinction" = "#D62728")) +
  scale_fill_manual(values = c("Colonization" = "#2CA02C",
                                "Extinction" = "#D62728")) +
  scale_x_continuous(breaks = years[-1]) +
  labs(title = "Derived turnover rates (stPGOcc)",
       subtitle = "Colonization and extinction rates derived from latent occupancy states",
       x = NULL, y = "Rate") +
  theme_pub +
  theme(strip.text = element_text(face = "italic"),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(figs_dir, "pub_fig_spatial_col_ext_rates.png"),
       fig_h, width = 10, height = 8, dpi = 300)
cat("    Saved: pub_fig_spatial_col_ext_rates.png\n")


# =========================================================================
# FIGURE I: Colonization/extinction hotspot maps
# =========================================================================
cat("  Figure I: Colonization/extinction hotspot maps...\n")

hotspot_plots <- list()
for (sp in species_codes) {
  fit <- fits_st[[sp]]
  dat <- sp_data[[sp]]
  z_samps <- fit$z.samples
  if (inherits(z_samps, "mcmc.list")) z_samps <- do.call(rbind, z_samps)

  J_sp <- dat$J
  n_samps <- nrow(z_samps)

  if (length(dim(z_samps)) == 3) {
    z_arr <- z_samps
  } else {
    z_arr <- array(z_samps, dim = c(n_samps, J_sp, n_years))
  }

  # Cumulative colonization/extinction probability per site
  col_site <- matrix(0, nrow = n_samps, ncol = J_sp)
  ext_site <- matrix(0, nrow = n_samps, ncol = J_sp)
  for (t in 2:n_years) {
    col_site <- col_site + z_arr[, , t] * (1 - z_arr[, , t - 1])
    ext_site <- ext_site + (1 - z_arr[, , t]) * z_arr[, , t - 1]
  }
  col_mean <- colMeans(col_site) / (n_years - 1)  # average annual col prob
  ext_mean <- colMeans(ext_site) / (n_years - 1)

  coords_use <- dat$coords_lonlat
  if (length(col_mean) != nrow(coords_use)) {
    n_use <- min(length(col_mean), nrow(coords_use))
    col_mean <- col_mean[1:n_use]
    ext_mean <- ext_mean[1:n_use]
    coords_use <- coords_use[1:n_use, ]
  }

  map_df <- data.frame(
    lon = coords_use[, 1], lat = coords_use[, 2],
    col = col_mean, ext = ext_mean
  )

  # Colonization map
  hotspot_plots[[paste0(sp, "_col")]] <- ggplot() +
    geom_sf(data = spain, fill = "grey95", colour = "grey60", linewidth = 0.3) +
    geom_point(data = map_df, aes(x = lon, y = lat, colour = col),
               size = 0.2, alpha = 0.7) +
    scale_colour_viridis_c(option = "C", name = "Col.\nrate",
                           limits = c(0, max(map_df$col, na.rm = TRUE))) +
    coord_sf(xlim = c(-10, 5), ylim = c(35, 44)) +
    labs(title = paste0(sp_short[sp], " \u2014 Colonization")) +
    theme_pub +
    theme(axis.title = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(face = "italic", size = 10))

  # Extinction map
  hotspot_plots[[paste0(sp, "_ext")]] <- ggplot() +
    geom_sf(data = spain, fill = "grey95", colour = "grey60", linewidth = 0.3) +
    geom_point(data = map_df, aes(x = lon, y = lat, colour = ext),
               size = 0.2, alpha = 0.7) +
    scale_colour_viridis_c(option = "B", name = "Ext.\nrate",
                           limits = c(0, max(map_df$ext, na.rm = TRUE))) +
    coord_sf(xlim = c(-10, 5), ylim = c(35, 44)) +
    labs(title = paste0(sp_short[sp], " \u2014 Extinction")) +
    theme_pub +
    theme(axis.title = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(face = "italic", size = 10))
}

# Arrange: 4 rows (species) x 2 columns (col, ext)
fig_i_grobs <- list()
for (sp in species_codes) {
  fig_i_grobs <- c(fig_i_grobs,
                    list(hotspot_plots[[paste0(sp, "_col")]],
                         hotspot_plots[[paste0(sp, "_ext")]]))
}

fig_i <- gridExtra::grid.arrange(
  grobs = fig_i_grobs, ncol = 2,
  top = grid::textGrob("Colonization & extinction hotspots (stPGOcc, 2017\u20132023)",
                        gp = grid::gpar(fontsize = 14, fontface = "bold"))
)
ggsave(file.path(figs_dir, "pub_fig_spatial_col_ext_maps.png"),
       fig_i, width = 12, height = 16, dpi = 300)
cat("    Saved: pub_fig_spatial_col_ext_maps.png\n")


# =========================================================================
# FIGURE J: Occupancy change maps (delta psi)
# =========================================================================
cat("  Figure J: Occupancy change maps (delta psi)...\n")

dpsi_plots <- list()
for (sp in species_codes) {
  fit <- fits_st[[sp]]
  dat <- sp_data[[sp]]
  z_samps <- fit$z.samples
  if (inherits(z_samps, "mcmc.list")) z_samps <- do.call(rbind, z_samps)

  J_sp <- dat$J
  n_samps <- nrow(z_samps)

  if (length(dim(z_samps)) == 3) {
    z_arr <- z_samps
  } else {
    z_arr <- array(z_samps, dim = c(n_samps, J_sp, n_years))
  }

  # Delta psi = mean(z_2023) - mean(z_2017)
  delta_psi <- colMeans(z_arr[, , n_years]) - colMeans(z_arr[, , 1])

  coords_use <- dat$coords_lonlat
  if (length(delta_psi) != nrow(coords_use)) {
    n_use <- min(length(delta_psi), nrow(coords_use))
    delta_psi <- delta_psi[1:n_use]
    coords_use <- coords_use[1:n_use, ]
  }

  dpsi_lim <- max(abs(quantile(delta_psi, c(0.01, 0.99))))

  dpsi_df <- data.frame(
    lon = coords_use[, 1], lat = coords_use[, 2], dpsi = delta_psi
  )

  dpsi_plots[[sp]] <- ggplot() +
    geom_sf(data = spain, fill = "grey95", colour = "grey60", linewidth = 0.3) +
    geom_point(data = dpsi_df, aes(x = lon, y = lat, colour = dpsi),
               size = 0.3, alpha = 0.7) +
    scale_colour_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                           midpoint = 0, limits = c(-dpsi_lim, dpsi_lim),
                           name = "\u0394\u03C8") +
    coord_sf(xlim = c(-10, 5), ylim = c(35, 44)) +
    labs(title = paste0(panel_labels[sp], ") ", sp_short[sp]),
         subtitle = paste0("Mean \u0394\u03C8 = ",
                           sprintf("%+.3f", mean(delta_psi)))) +
    theme_pub +
    theme(axis.title = element_blank(),
          plot.title = element_text(face = "italic"))
}

fig_j <- gridExtra::grid.arrange(
  grobs = dpsi_plots, ncol = 2,
  top = grid::textGrob(
    "Occupancy change 2017\u20132023 (stPGOcc)",
    gp = grid::gpar(fontsize = 14, fontface = "bold"))
)
ggsave(file.path(figs_dir, "pub_fig_spatial_delta_psi_maps.png"),
       fig_j, width = 12, height = 11, dpi = 300)
cat("    Saved: pub_fig_spatial_delta_psi_maps.png\n")


###############################################################################
# 6. SUMMARY
###############################################################################

cat("\n")
cat(strrep("=", 70), "\n")
cat("  SPATIAL RESULTS — SUMMARY\n")
cat(strrep("=", 70), "\n\n")

cat("--- Model comparison ---\n")
print(table1, row.names = FALSE)

cat("\n--- Output files ---\n")
cat("Tables:\n")
cat("  ", table1_path, "\n")
cat("  ", table2_path, "\n")
cat("\nFigures:\n")
for (f in list.files(figs_dir, pattern = "pub_fig_spatial_", full.names = TRUE)) {
  cat("  ", f, "\n")
}

cat("\n")
cat(strrep("=", 70), "\n")
cat("  DONE\n")
cat(strrep("=", 70), "\n")
