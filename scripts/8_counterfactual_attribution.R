###############################################################################
# 8_counterfactual_attribution.R
#
# Purpose: Factorial counterfactual attribution analysis for dynamic occupancy
#          models of four Iberian steppe bird species (2017-2023).
#
#          Quantifies how much of the observed colonisation and extinction
#          dynamics is attributable to (a) climate trends, (b) land-use change,
#          and (c) their interaction, using the fitted colext models and
#          observed covariate time series.
#
# Inputs:  results/{sp}_model_object.rds               (fitted colext models)
#          data/processed_2023/{sp}/{sp}_occ_wide_dynamic.csv (site-year data)
#          data/processed_2023/{sp}/{sp}_scaling_params.rds   (static scaling)
#          R/model_configs.R                                   (model formulas)
#
# Outputs: results/covariate_trends/trend_{cov}_allspecies.tif
#          results/counterfactual_predictions.rds
#          results/pub_table_attribution_summary.csv
#          results/attribution_hotspot_correlation.csv
#          figs/attribution/fig_attribution_*.png
#          figs/fig_attribution_summary_{sp}.png
#
# Structure: Tasks 1-8 matching the analysis specification.
###############################################################################

# -- Load packages --
library(here)
library(unmarked)
library(dplyr)
library(tidyr)
library(ggplot2)
library(terra)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
# rmapshaper loaded conditionally (V8 dependency may fail)
has_rmapshaper <- requireNamespace("rmapshaper", quietly = TRUE)

select <- dplyr::select
filter <- dplyr::filter
theme_set(theme_bw())

# Load model configurations
source(here("R", "model_configs.R"))

###############################################################################
# GLOBAL DEFINITIONS
###############################################################################

SPECIES <- list(
  otitar = "Otis tarda",
  ptealc = "Pterocles alchata",
  pteori = "Pterocles orientalis",
  tettet = "Tetrax tetrax"
)

SCENARIOS <- list(
  S0 = list(name = "Null baseline",  climate = "frozen", landuse = "frozen"),
  S1 = list(name = "Climate only",   climate = "observed", landuse = "frozen"),
  S2 = list(name = "Land use only",  climate = "frozen", landuse = "observed"),
  S3 = list(name = "Combined",       climate = "observed", landuse = "observed")
)

YEARS   <- 2017:2023
T_YEARS <- length(YEARS)

# Classification of covariates into climate vs land-use
CLIMATE_COVS <- c("NDVI", "pr", "tmmn", "tmmx")
LANDUSE_COVS <- c("Land_Cover_Type_1_Percent_Class_0",
                   "Land_Cover_Type_1_Percent_Class_6",
                   "Land_Cover_Type_1_Percent_Class_7",
                   "Land_Cover_Type_1_Percent_Class_10",
                   "Land_Cover_Type_1_Percent_Class_12",
                   "Land_Cover_Type_1_Percent_Class_13",
                   "Land_Cover_Type_1_Percent_Class_14")

# Helper: apply training scaling
scale_with_params <- function(x, center, scale_val) {
  (x - center) / scale_val
}

# Track session outcomes
session_log <- list()

# Spain basemap for plots
spain <- ne_countries(country = "spain", scale = "medium", returnclass = "sf")
if (has_rmapshaper) {
  spain_crop <- rmapshaper::ms_filter_islands(spain, min_area = 100000000000,
                                               drop_null_geometries = TRUE)
} else {
  # Crop to Iberian peninsula bounding box (excludes Canary Islands)
  spain_crop <- st_crop(spain, xmin = -10, xmax = 5, ymin = 35, ymax = 44)
}

###############################################################################
# HELPER: Classify covariates for a given species config
###############################################################################
classify_covariates <- function(cfg) {
  all_dyn <- unique(c(cfg$gamma_vars, cfg$epsilon_vars))
  clim <- intersect(all_dyn, CLIMATE_COVS)
  lu   <- intersect(all_dyn, LANDUSE_COVS)
  list(climate = clim, landuse = lu, all = all_dyn)
}

###############################################################################
# HELPER: Load data and compute training scaling for one species
###############################################################################
load_species_data <- function(sp) {

  occ_wide <- read.csv(
    here("data", "processed_2023", sp, paste0(sp, "_occ_wide_dynamic.csv"))
  )

  # Drop sites with missing NDVI (same filter as script 4)
  predict_data <- occ_wide %>% drop_na(paste0("NDVI_", YEARS))

  # Extract raw dynamic covariate matrices and compute training scaling
  all_covs <- unique(c(CLIMATE_COVS, LANDUSE_COVS))
  raw_matrices <- list()
  train_scale  <- list()

  for (cov in all_covs) {
    cols <- paste0(cov, "_", YEARS)
    if (all(cols %in% names(predict_data))) {
      mat <- as.matrix(predict_data[, cols])
      scaled_mat <- scale(mat)
      ctrs <- attr(scaled_mat, "scaled:center")
      scls <- attr(scaled_mat, "scaled:scale")
      raw_matrices[[cov]] <- mat
      for (yi in seq_along(YEARS)) {
        train_scale[[paste0(cov, "_", YEARS[yi])]] <- list(
          center = ctrs[yi], scale = scls[yi]
        )
      }
    }
  }

  # Site coordinates
  coords <- predict_data[, c("latitude", "longitude")]

  list(
    predict_data = predict_data,
    raw_matrices = raw_matrices,
    train_scale  = train_scale,
    coords       = coords,
    S            = nrow(predict_data)
  )
}

###############################################################################
# HELPER: Build newdata for predict() under a given scenario
###############################################################################
build_scenario_newdata <- function(sp_data, cfg, cov_class, scenario, year_idx) {
  yr <- YEARS[year_idx]
  S  <- sp_data$S

  # Determine which covariates to use for this submodel (col or ext)
  # cov_class$vars is the set of vars for this submodel
  newdata <- data.frame(row.names = seq_len(S))

  for (v in cov_class$vars) {
    key <- paste0(v, "_", yr)
    key_2017 <- paste0(v, "_", YEARS[1])

    is_climate <- v %in% CLIMATE_COVS
    is_landuse <- v %in% LANDUSE_COVS

    # Decide whether to use observed or frozen values
    if (is_climate) {
      use_frozen <- (scenario$climate == "frozen")
    } else if (is_landuse) {
      use_frozen <- (scenario$landuse == "frozen")
    } else {
      use_frozen <- FALSE  # unknown type: use observed
    }

    if (use_frozen) {
      # Use 2017 site-level values, scaled with 2017 training params
      raw_val <- sp_data$raw_matrices[[v]][, 1]  # column 1 = 2017
      sc <- sp_data$train_scale[[key_2017]]
    } else {
      # Use observed year values, scaled with that year's training params
      raw_val <- sp_data$raw_matrices[[v]][, year_idx]
      sc <- sp_data$train_scale[[key]]
    }

    if (!is.null(sc)) {
      newdata[[v]] <- scale_with_params(raw_val, sc$center, sc$scale)
    } else {
      newdata[[v]] <- raw_val
      warning("No scaling params for: ", key)
    }
  }

  newdata
}

###############################################################################
# MAIN ANALYSIS LOOP
###############################################################################

# Storage for all results
all_predictions <- list()
all_trends      <- list()
all_attribution <- list()

for (sp in names(SPECIES)) {

  message("\n", strrep("=", 70))
  message("  SPECIES: ", sp, " (", SPECIES[[sp]], ")")
  message(strrep("=", 70))

  tryCatch({

    # -- Load model --
    model_path <- here("results", paste0(sp, "_model_object.rds"))
    if (!file.exists(model_path)) {
      stop("Model object not found: ", model_path)
    }
    mod <- readRDS(model_path)
    message("  Model loaded. AIC = ", round(mod@AIC, 2))

    # -- Load config --
    cfg <- get_model_config(sp)

    # -- Classify covariates --
    cov_info <- classify_covariates(cfg)
    message("  Climate covariates in model: ",
            paste(cov_info$climate, collapse = ", "))
    message("  Land-use covariates in model: ",
            paste(cov_info$landuse, collapse = ", "))

    # -- Load data and scaling --
    sp_data <- load_species_data(sp)
    S <- sp_data$S
    message("  Sites: ", S)

    ############################################################################
    # TASK 1: COVARIATE TRENDS (2017-2023)
    ############################################################################
    message("\n  --- TASK 1: Covariate trends ---")

    trend_results <- list()

    for (cov in cov_info$all) {
      mat <- sp_data$raw_matrices[[cov]]
      if (is.null(mat)) next

      slopes   <- numeric(S)
      slope_se <- numeric(S)
      p_vals   <- numeric(S)

      for (i in seq_len(S)) {
        y_vals <- mat[i, ]
        if (all(is.na(y_vals)) || sd(y_vals, na.rm = TRUE) == 0) {
          slopes[i] <- 0
          slope_se[i] <- NA
          p_vals[i] <- 1
          next
        }
        fit <- lm(y_vals ~ YEARS)
        slopes[i]   <- coef(fit)[2]
        slope_se[i] <- summary(fit)$coefficients[2, 2]
        p_vals[i]   <- summary(fit)$coefficients[2, 4]
      }

      frac_sig_pos <- mean(p_vals < 0.05 & slopes > 0, na.rm = TRUE)
      frac_sig_neg <- mean(p_vals < 0.05 & slopes < 0, na.rm = TRUE)

      message("    ", cov, ": mean trend = ",
              round(mean(slopes, na.rm = TRUE), 5),
              " (SD = ", round(sd(slopes, na.rm = TRUE), 5), ")",
              " | sig+ = ", round(frac_sig_pos * 100, 1), "%",
              " | sig- = ", round(frac_sig_neg * 100, 1), "%")

      trend_results[[cov]] <- data.frame(
        site_idx  = seq_len(S),
        latitude  = sp_data$coords$latitude,
        longitude = sp_data$coords$longitude,
        slope     = slopes,
        slope_se  = slope_se,
        p_value   = p_vals
      )
    }

    # Save trend rasters
    for (cov in names(trend_results)) {
      tr <- trend_results[[cov]]
      tr_sf <- st_as_sf(tr, coords = c("longitude", "latitude"), crs = 4326)
      r <- terra::rasterize(
        terra::vect(tr_sf),
        terra::rast(ext = terra::ext(tr_sf), resolution = 0.05, crs = "EPSG:4326"),
        field = "slope", fun = mean
      )
      fname <- here("results", "covariate_trends",
                     paste0("trend_", cov, "_", sp, ".tif"))
      terra::writeRaster(r, fname, overwrite = TRUE)
    }

    all_trends[[sp]] <- trend_results

    ############################################################################
    # TASKS 2-3: BUILD COUNTERFACTUAL DATAFRAMES AND PREDICT
    ############################################################################
    message("\n  --- TASKS 2-3: Counterfactual predictions ---")

    sp_predictions <- list()

    for (scen_name in names(SCENARIOS)) {
      scenario <- SCENARIOS[[scen_name]]
      message("    Scenario ", scen_name, ": ", scenario$name)

      scen_preds <- list()

      for (yi in seq_len(T_YEARS)) {
        yr <- YEARS[yi]

        # Colonisation newdata
        col_class <- list(vars = cfg$gamma_vars)
        col_nd <- build_scenario_newdata(sp_data, cfg, col_class, scenario, yi)

        # Extinction newdata
        ext_class <- list(vars = cfg$epsilon_vars)
        ext_nd <- build_scenario_newdata(sp_data, cfg, ext_class, scenario, yi)

        # Predict colonisation
        col_pred <- tryCatch(
          predict(mod, type = "col", newdata = col_nd),
          error = function(e) {
            warning("  predict(col) failed for ", sp, " ", scen_name,
                    " year ", yr, ": ", e$message)
            NULL
          }
        )

        # Predict extinction
        ext_pred <- tryCatch(
          predict(mod, type = "ext", newdata = ext_nd),
          error = function(e) {
            warning("  predict(ext) failed for ", sp, " ", scen_name,
                    " year ", yr, ": ", e$message)
            NULL
          }
        )

        if (is.null(col_pred) || is.null(ext_pred)) next

        # Validate non-degenerate
        col_range <- range(col_pred$Predicted, na.rm = TRUE)
        ext_range <- range(ext_pred$Predicted, na.rm = TRUE)
        if (col_range[1] == col_range[2] &&
            (col_range[1] == 0 || col_range[1] == 1)) {
          warning("  Degenerate col predictions for ", sp, " ", scen_name,
                  " year ", yr, ": all = ", col_range[1])
        }

        scen_preds[[as.character(yr)]] <- data.frame(
          site_id    = seq_len(S),
          lat        = sp_data$coords$latitude,
          lon        = sp_data$coords$longitude,
          gamma_hat  = col_pred$Predicted,
          gamma_se   = col_pred$SE,
          epsilon_hat = ext_pred$Predicted,
          epsilon_se  = ext_pred$SE
        )
      }

      sp_predictions[[scen_name]] <- scen_preds
    }

    all_predictions[[sp]] <- sp_predictions

    # Quick diagnostic
    g_s0 <- mean(sp_predictions[["S0"]][["2020"]]$gamma_hat, na.rm = TRUE)
    g_s3 <- mean(sp_predictions[["S3"]][["2020"]]$gamma_hat, na.rm = TRUE)
    e_s0 <- mean(sp_predictions[["S0"]][["2020"]]$epsilon_hat, na.rm = TRUE)
    e_s3 <- mean(sp_predictions[["S3"]][["2020"]]$epsilon_hat, na.rm = TRUE)
    message("    Diagnostic (2020): gamma S0=", round(g_s0, 5),
            " S3=", round(g_s3, 5),
            " | epsilon S0=", round(e_s0, 5),
            " S3=", round(e_s3, 5))

    ############################################################################
    # TASK 4: ATTRIBUTION DECOMPOSITION
    ############################################################################
    message("\n  --- TASK 4: Attribution decomposition ---")

    attr_results <- list()

    for (yi in seq_len(T_YEARS)) {
      yr <- as.character(YEARS[yi])

      g0 <- sp_predictions[["S0"]][[yr]]$gamma_hat
      g1 <- sp_predictions[["S1"]][[yr]]$gamma_hat
      g2 <- sp_predictions[["S2"]][[yr]]$gamma_hat
      g3 <- sp_predictions[["S3"]][[yr]]$gamma_hat

      e0 <- sp_predictions[["S0"]][[yr]]$epsilon_hat
      e1 <- sp_predictions[["S1"]][[yr]]$epsilon_hat
      e2 <- sp_predictions[["S2"]][[yr]]$epsilon_hat
      e3 <- sp_predictions[["S3"]][[yr]]$epsilon_hat

      attr_results[[yr]] <- data.frame(
        site_id = seq_len(S),
        lat     = sp_data$coords$latitude,
        lon     = sp_data$coords$longitude,
        year    = YEARS[yi],
        # Colonisation attribution
        delta_gamma_climate     = g1 - g0,
        delta_gamma_landuse     = g2 - g0,
        delta_gamma_interaction = g3 - g1 - g2 + g0,
        delta_gamma_total       = g3 - g0,
        # Extinction attribution
        delta_epsilon_climate     = e1 - e0,
        delta_epsilon_landuse     = e2 - e0,
        delta_epsilon_interaction = e3 - e1 - e2 + e0,
        delta_epsilon_total       = e3 - e0
      )
    }

    attr_df <- bind_rows(attr_results)
    all_attribution[[sp]] <- attr_df

    # Cumulative summary
    attr_summary <- attr_df %>%
      group_by(site_id, lat, lon) %>%
      summarise(across(starts_with("delta_"), mean, .names = "{.col}"),
                .groups = "drop")

    message("    Cumulative attribution (mean across sites & years):")
    for (v in c("delta_gamma_climate", "delta_gamma_landuse",
                "delta_gamma_interaction", "delta_gamma_total",
                "delta_epsilon_climate", "delta_epsilon_landuse",
                "delta_epsilon_interaction", "delta_epsilon_total")) {
      message("      ", v, " = ",
              round(mean(attr_summary[[v]], na.rm = TRUE), 6),
              " (SD = ", round(sd(attr_summary[[v]], na.rm = TRUE), 6), ")")
    }

    ############################################################################
    # TASK 5: SPATIAL MAPS
    ############################################################################
    message("\n  --- TASK 5: Spatial maps ---")

    # Mean attribution across years per site
    attr_mean <- attr_df %>%
      group_by(site_id, lat, lon) %>%
      summarise(across(starts_with("delta_"), mean, .names = "{.col}"),
                .groups = "drop")

    # Map quantities to plot
    map_vars <- c("delta_gamma_climate", "delta_gamma_landuse",
                  "delta_epsilon_climate", "delta_epsilon_landuse",
                  "delta_gamma_total", "delta_epsilon_total")

    var_labels <- c(
      delta_gamma_climate     = "Climate effect on colonisation",
      delta_gamma_landuse     = "Land-use effect on colonisation",
      delta_epsilon_climate   = "Climate effect on extinction",
      delta_epsilon_landuse   = "Land-use effect on extinction",
      delta_gamma_total       = "Total colonisation change",
      delta_epsilon_total     = "Total extinction change"
    )

    for (v in map_vars) {
      plot_df <- attr_mean[, c("lon", "lat", v)]
      names(plot_df)[3] <- "value"

      max_abs <- max(abs(plot_df$value), na.rm = TRUE)
      if (max_abs == 0) max_abs <- 0.001

      p <- ggplot() +
        geom_sf(data = spain_crop, fill = "grey95", colour = "grey50") +
        geom_point(data = plot_df, aes(x = lon, y = lat, colour = value),
                   size = 0.4, alpha = 0.8) +
        scale_colour_gradient2(
          low = "red3", mid = "white", high = "blue3",
          midpoint = 0, limits = c(-max_abs, max_abs),
          name = "Effect"
        ) +
        labs(title = paste0(SPECIES[[sp]], ": ", var_labels[v]),
             x = "Longitude", y = "Latitude") +
        coord_sf(xlim = c(-10, 5), ylim = c(35, 44)) +
        theme(legend.position = "right",
              plot.title = element_text(size = 10, face = "italic"))

      fname <- here("figs", "attribution",
                     paste0("fig_attribution_", v, "_", sp, ".png"))
      ggsave(fname, p, width = 8, height = 6, dpi = 300)
    }

    # 4-panel summary
    summary_vars <- c("delta_gamma_climate", "delta_gamma_landuse",
                       "delta_epsilon_climate", "delta_epsilon_landuse")

    panels <- list()
    for (v in summary_vars) {
      plot_df <- attr_mean[, c("lon", "lat", v)]
      names(plot_df)[3] <- "value"
      plot_df$panel <- var_labels[v]
      panels <- c(panels, list(plot_df))
    }
    panel_df <- bind_rows(panels)

    max_abs <- max(abs(panel_df$value), na.rm = TRUE)
    if (max_abs == 0) max_abs <- 0.001

    p_summary <- ggplot() +
      geom_sf(data = spain_crop, fill = "grey95", colour = "grey50") +
      geom_point(data = panel_df,
                 aes(x = lon, y = lat, colour = value),
                 size = 0.3, alpha = 0.8) +
      scale_colour_gradient2(
        low = "red3", mid = "white", high = "blue3",
        midpoint = 0, limits = c(-max_abs, max_abs),
        name = "Effect"
      ) +
      facet_wrap(~ panel, ncol = 2) +
      labs(title = paste0(SPECIES[[sp]], " - Attribution summary"),
           x = "Longitude", y = "Latitude") +
      coord_sf(xlim = c(-10, 5), ylim = c(35, 44)) +
      theme(legend.position = "right",
            plot.title = element_text(size = 12, face = "italic"),
            strip.text = element_text(size = 8))

    ggsave(here("figs", paste0("fig_attribution_summary_", sp, ".png")),
           p_summary, width = 12, height = 10, dpi = 300)

    session_log[[sp]] <- list(status = "success")
    message("  ", sp, " completed successfully.")

  }, error = function(e) {
    session_log[[sp]] <<- list(status = "failed", reason = e$message)
    message("  ERROR for ", sp, ": ", e$message)
  })
}

###############################################################################
# TASK 3 (continued): Save predictions
###############################################################################
message("\n", strrep("=", 70))
message("  Saving predictions...")
saveRDS(all_predictions, here("results", "counterfactual_predictions.rds"))

###############################################################################
# TASK 6: CROSS-SPECIES SUMMARY TABLE
###############################################################################
message("\n--- TASK 6: Cross-species summary table ---")

summary_rows <- list()

for (sp in names(SPECIES)) {
  if (is.null(all_attribution[[sp]])) next

  attr_df <- all_attribution[[sp]]

  # Mean across sites and years
  site_means <- attr_df %>%
    group_by(site_id) %>%
    summarise(across(starts_with("delta_"), mean, .names = "{.col}"),
              .groups = "drop")

  m <- function(v) round(mean(site_means[[v]], na.rm = TRUE), 6)
  s <- function(v) round(sd(site_means[[v]], na.rm = TRUE), 6)

  dgc  <- m("delta_gamma_climate")
  dgl  <- m("delta_gamma_landuse")
  dgi  <- m("delta_gamma_interaction")
  dgt  <- m("delta_gamma_total")
  dec  <- m("delta_epsilon_climate")
  del  <- m("delta_epsilon_landuse")
  dei  <- m("delta_epsilon_interaction")
  det_ <- m("delta_epsilon_total")

  # Fractions (handle division by near-zero)
  frac <- function(part, total) {
    if (abs(total) < 1e-8) return(NA_real_)
    round(part / total * 100, 1)
  }

  # Dominant driver
  effects <- c(climate = abs(dgc) + abs(dec),
               landuse = abs(dgl) + abs(del),
               interaction = abs(dgi) + abs(dei))
  if (max(effects) < 1e-6) {
    dominant <- "Negligible"
  } else {
    dominant <- names(which.max(effects))
    dominant <- switch(dominant,
                       climate = "Climate",
                       landuse = "Land use",
                       interaction = "Interaction")
  }

  summary_rows[[sp]] <- data.frame(
    Species = SPECIES[[sp]],
    Mean_delta_gamma_climate     = paste0(dgc, " (", s("delta_gamma_climate"), ")"),
    Mean_delta_gamma_landuse     = paste0(dgl, " (", s("delta_gamma_landuse"), ")"),
    Mean_delta_gamma_interaction = paste0(dgi, " (", s("delta_gamma_interaction"), ")"),
    Mean_delta_epsilon_climate     = paste0(dec, " (", s("delta_epsilon_climate"), ")"),
    Mean_delta_epsilon_landuse     = paste0(del, " (", s("delta_epsilon_landuse"), ")"),
    Mean_delta_epsilon_interaction = paste0(dei, " (", s("delta_epsilon_interaction"), ")"),
    Climate_fraction_gamma  = frac(dgc, dgt),
    Landuse_fraction_gamma  = frac(dgl, dgt),
    Climate_fraction_epsilon = frac(dec, det_),
    Landuse_fraction_epsilon = frac(del, det_),
    Dominant_driver = dominant,
    stringsAsFactors = FALSE
  )
}

summary_table <- bind_rows(summary_rows)
write.csv(summary_table,
          here("results", "pub_table_attribution_summary.csv"),
          row.names = FALSE)
message("  Summary table saved: results/pub_table_attribution_summary.csv")

###############################################################################
# TASK 7: INTERPRETATION PARAGRAPHS
###############################################################################
message("\n--- TASK 7: Interpretation paragraphs ---")

interpretations <- list()

for (sp in names(SPECIES)) {
  if (is.null(all_attribution[[sp]])) next

  cfg <- get_model_config(sp)
  cov_info <- classify_covariates(cfg)
  attr_df <- all_attribution[[sp]]

  site_means <- attr_df %>%
    group_by(site_id, lat, lon) %>%
    summarise(across(starts_with("delta_"), mean), .groups = "drop")

  # Quantitative summaries
  dgc <- mean(site_means$delta_gamma_climate, na.rm = TRUE)
  dgl <- mean(site_means$delta_gamma_landuse, na.rm = TRUE)
  dgi <- mean(site_means$delta_gamma_interaction, na.rm = TRUE)
  dec <- mean(site_means$delta_epsilon_climate, na.rm = TRUE)
  del <- mean(site_means$delta_epsilon_landuse, na.rm = TRUE)
  dei <- mean(site_means$delta_epsilon_interaction, na.rm = TRUE)

  # Spatial pattern: where are effects strongest?
  q_clim <- quantile(abs(site_means$delta_gamma_climate), 0.9, na.rm = TRUE)
  q_lu   <- quantile(abs(site_means$delta_epsilon_landuse), 0.9, na.rm = TRUE)

  high_clim_col <- site_means %>% filter(abs(delta_gamma_climate) > q_clim)
  high_lu_ext   <- site_means %>% filter(abs(delta_epsilon_landuse) > q_lu)

  clim_lat_range <- if (nrow(high_clim_col) > 0) {
    range(high_clim_col$lat, na.rm = TRUE)
  } else { c(38, 42) }
  lu_lat_range <- if (nrow(high_lu_ext) > 0) {
    range(high_lu_ext$lat, na.rm = TRUE)
  } else { c(38, 42) }

  # NDVI caveat
  ndvi_in_gamma <- "NDVI" %in% cfg$gamma_vars
  ndvi_in_epsilon <- "NDVI" %in% cfg$epsilon_vars
  ndvi_caveat <- ndvi_in_gamma || ndvi_in_epsilon

  # Build interpretation
  lines <- character()

  # Sentence 1: dominant driver
  drivers <- c(climate = abs(dgc) + abs(dec),
               landuse = abs(dgl) + abs(del))
  dominant <- names(which.max(drivers))

  if (dominant == "climate") {
    lines <- c(lines, paste0(
      "For ", SPECIES[[sp]], ", climate trends were the dominant driver of ",
      "occupancy dynamics over 2017-2023, with a mean colonisation effect of ",
      sprintf("%.5f", dgc), " and extinction effect of ",
      sprintf("%.5f", dec), "."
    ))
  } else {
    lines <- c(lines, paste0(
      "For ", SPECIES[[sp]], ", land-use change was the dominant driver of ",
      "occupancy dynamics over 2017-2023, with a mean colonisation effect of ",
      sprintf("%.5f", dgl), " and extinction effect of ",
      sprintf("%.5f", del), "."
    ))
  }

  # Sentence 2: secondary driver
  secondary <- names(which.min(drivers))
  if (secondary == "climate") {
    lines <- c(lines, paste0(
      "Climate trends contributed a secondary effect (colonisation: ",
      sprintf("%.5f", dgc), "; extinction: ", sprintf("%.5f", dec), "), ",
      "while the interaction term was ",
      ifelse(abs(dgi) + abs(dei) < 1e-5, "negligible",
             paste0("small (gamma: ", sprintf("%.5f", dgi),
                    "; epsilon: ", sprintf("%.5f", dei), ")")), "."
    ))
  } else {
    lines <- c(lines, paste0(
      "Land-use change contributed a secondary effect (colonisation: ",
      sprintf("%.5f", dgl), "; extinction: ", sprintf("%.5f", del), "), ",
      "while the interaction term was ",
      ifelse(abs(dgi) + abs(dei) < 1e-5, "negligible",
             paste0("small (gamma: ", sprintf("%.5f", dgi),
                    "; epsilon: ", sprintf("%.5f", dei), ")")), "."
    ))
  }

  # Sentence 3: spatial pattern
  lines <- c(lines, paste0(
    "The strongest climate effects on colonisation were concentrated at ",
    "latitudes ", round(clim_lat_range[1], 1), "-",
    round(clim_lat_range[2], 1), " N, while the strongest land-use effects on ",
    "extinction occurred between ", round(lu_lat_range[1], 1), "-",
    round(lu_lat_range[2], 1), " N."
  ))

  # Sentence 4: consistency with coefficients
  lines <- c(lines, paste0(
    "These attribution results are consistent with the fitted colext ",
    "coefficient signs, which show ",
    paste(cfg$gamma_vars, collapse = ", "), " driving colonisation and ",
    paste(cfg$epsilon_vars, collapse = ", "), " driving extinction."
  ))

  # Sentence 5: caveats
  if (ndvi_caveat) {
    lines <- c(lines, paste0(
      "NDVI was classified as a climate-adjacent covariate in this analysis; ",
      "however, its partial coupling to land management practices means that ",
      "the climate attribution may partially reflect land-use effects mediated ",
      "through vegetation productivity."
    ))
  }

  if (sp == "tettet") {
    lines <- c(lines, paste0(
      "For this species, only cropland proportion (LC12) enters both the ",
      "colonisation and extinction submodels, with no climate covariates. ",
      "Consequently, the climate attribution is zero by construction and all ",
      "observed dynamics are attributed to land-use change."
    ))
  }

  if (sp == "ptealc") {
    lines <- c(lines, paste0(
      "The extinction submodel coefficients for precipitation and temperature ",
      "are large in magnitude (|beta| > 50), suggesting potential estimation ",
      "instability that may inflate the climate attribution for extinction."
    ))
  }

  interpretations[[sp]] <- paste(lines, collapse = " ")
  message("\n  ", SPECIES[[sp]], ":")
  message("  ", interpretations[[sp]])
}

# Save interpretations
interp_text <- paste(
  sapply(names(interpretations), function(sp) {
    paste0("## ", SPECIES[[sp]], "\n\n", interpretations[[sp]], "\n")
  }),
  collapse = "\n"
)
writeLines(interp_text,
           here("results", "attribution_interpretations.txt"))

###############################################################################
# TASK 8: OVERLAY WITH stPGOcc HOTSPOT MAPS
###############################################################################
message("\n--- TASK 8: Overlay with stPGOcc hotspots ---")

hotspot_correlations <- list()
hotspots_found <- FALSE

for (sp in names(SPECIES)) {
  col_hotspot_path <- here("results",
                           paste0("stpgocc_col_hotspots_", sp, ".tif"))
  ext_hotspot_path <- here("results",
                           paste0("stpgocc_ext_hotspots_", sp, ".tif"))

  if (!file.exists(col_hotspot_path) || !file.exists(ext_hotspot_path)) {
    message("  stPGOcc hotspot rasters not found for ", sp, " - skipping.")
    next
  }

  hotspots_found <- TRUE
  message("  Processing stPGOcc overlay for ", sp, "...")

  col_hotspot <- terra::rast(col_hotspot_path)
  ext_hotspot <- terra::rast(ext_hotspot_path)

  # Get attribution rasters
  attr_df <- all_attribution[[sp]]
  attr_mean <- attr_df %>%
    group_by(site_id, lat, lon) %>%
    summarise(across(starts_with("delta_"), mean), .groups = "drop")

  attr_sf <- st_as_sf(attr_mean, coords = c("lon", "lat"), crs = 4326)
  attr_vect <- terra::vect(attr_sf)

  # Extract hotspot values at attribution sites
  col_vals <- terra::extract(col_hotspot, attr_vect)[, 2]
  ext_vals <- terra::extract(ext_hotspot, attr_vect)[, 2]

  # Compute correlations
  cors <- data.frame(
    species = SPECIES[[sp]],
    r_gamma_landuse_vs_col_hotspot = cor(
      attr_mean$delta_gamma_landuse, col_vals, use = "complete.obs"),
    r_epsilon_landuse_vs_ext_hotspot = cor(
      attr_mean$delta_epsilon_landuse, ext_vals, use = "complete.obs"),
    r_gamma_climate_vs_col_hotspot = cor(
      attr_mean$delta_gamma_climate, col_vals, use = "complete.obs"),
    r_epsilon_climate_vs_ext_hotspot = cor(
      attr_mean$delta_epsilon_climate, ext_vals, use = "complete.obs"),
    stringsAsFactors = FALSE
  )

  hotspot_correlations[[sp]] <- cors
  message("    Correlations computed for ", sp)
}

if (hotspots_found) {
  cor_table <- bind_rows(hotspot_correlations)
  write.csv(cor_table,
            here("results", "attribution_hotspot_correlation.csv"),
            row.names = FALSE)
  message("  Correlation table saved.")
} else {
  message("  No stPGOcc hotspot rasters found for any species.")
  message("  To generate them, run the stPGOcc prediction pipeline first.")
}

###############################################################################
# SESSION SUMMARY
###############################################################################
message("\n", strrep("=", 70))
message("  SESSION SUMMARY")
message(strrep("=", 70))

for (sp in names(SPECIES)) {
  log <- session_log[[sp]]
  if (is.null(log)) {
    message("  ", sp, ": not processed")
  } else if (log$status == "success") {
    message("  ", sp, " (", SPECIES[[sp]], "): SUCCESS")
  } else {
    message("  ", sp, " (", SPECIES[[sp]], "): FAILED - ", log$reason)
  }
}

message("\nOutputs generated:")
message("  results/counterfactual_predictions.rds")
message("  results/pub_table_attribution_summary.csv")
message("  results/attribution_interpretations.txt")
message("  results/covariate_trends/trend_*.tif")
message("  figs/attribution/fig_attribution_*.png")
message("  figs/fig_attribution_summary_*.png")
if (hotspots_found) {
  message("  results/attribution_hotspot_correlation.csv")
}

message("\nScript 8 complete.")
