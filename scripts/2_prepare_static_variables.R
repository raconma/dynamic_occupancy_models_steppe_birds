###############################################################################
# 2_prepare_static_variables.R
#
# Purpose: For each species, filter repeat visits, extract bioclimatic and
#          topographic covariates, standardise, convert to wide format.
#          Produces the dataset needed before GEE dynamic variable extraction.
#
# Inputs:  data/processed/{sp}/ebd_{sp}_breeding_spain_zf.csv  (from step 1)
#          data/raw/environmental_data/environmental_data_occ/variables_spain.grd
#          data/raw/topology_data/topo_aspect.asc
#          data/raw/topology_data/topo_elev.asc
#          data/raw/topology_data/topo_slope.asc
#
# Outputs: data/processed/{sp}/{sp}_occ_species_observed.csv
#          data/processed/{sp}/{sp}_occ_wide_static.csv
#          data/processed/{sp}/{sp}_occ_wide_latlong.csv
#          data/processed/{sp}/{sp}_scaling_params.rds  <- NEW: for reproducible prediction
#
# Replaces: 2_prepare_static_variables_{otitar,ptealc,pteori,tettet}.R
###############################################################################

# -- Load packages --
library(here)
library(unmarked)
library(auk)
library(sf)
library(dplyr)
library(tidyr)
library(raster)
library(terra)

# Resolve namespace conflicts
select <- dplyr::select
filter <- dplyr::filter

set.seed(1)

# -- Species to process (run all, or subset) --
species_codes <- c("otitar", "ptealc", "pteori", "tettet")

# -- Variables to standardise --
vars_to_scale <- c("bio1", "tree_cover", "bio2", "grass_cover",
                    "topo_aspect", "topo_elev")

# -- Observation covariates for wide format --
obs_covs <- c("time_observations_started", "duration_minutes",
              "effort_distance_km", "number_observers", "protocol_type")

# -- Site covariates for wide format --
site_covs_wide <- c("locality_id", "n_observations", "cells",
                     "latitude", "longitude", vars_to_scale)

# -- Years covered --
years <- 2017:2022

# -- Load environmental rasters (shared across species) --
message("Loading environmental rasters...")
variables_environmental <- raster::stack(
  here("data", "raw", "environmental_data", "environmental_data_occ", "variables_spain.grd"))
variables_aspect <- raster::stack(
  here("data", "raw", "topology_data", "topo_aspect.asc"))
variables_elev <- raster::stack(
  here("data", "raw", "topology_data", "topo_elev.asc"))
variables_slope <- raster::stack(
  here("data", "raw", "topology_data", "topo_slope.asc"))


###############################################################################
# MAIN LOOP: process each species
###############################################################################

for (sp in species_codes) {
  message("\n=== Processing: ", sp, " ===")

  # -- Read zero-filled data from step 1 --
  occ_raw <- read.csv(here("data", "processed", sp,
                            paste0("ebd_", sp, "_breeding_spain_zf.csv")))
  occ_raw <- occ_raw %>% filter(year >= 2017, year <= 2022)
  message("  Raw records: ", nrow(occ_raw))

  # -- Filter repeat visits (auk) --
  occ <- filter_repeat_visits(
    occ_raw,
    min_obs = 3, max_obs = 10,
    annual_closure = TRUE,
    date_var = "observation_date",
    site_vars = c("locality_id", "observer_id")
  )
  message("  Sites after repeat-visit filter: ", n_distinct(occ$site))

  # -- Save observed-only records --
  occ_species_observed <- occ %>% filter(species_observed == TRUE)
  write.csv(occ_species_observed,
            here("data", "processed", sp, paste0(sp, "_occ_species_observed.csv")),
            row.names = FALSE)

  # -- Extract environmental variables at site locations --
  occ_var <- occ %>%
    cbind(as.data.frame(terra::extract(variables_environmental,
                                       occ[, c("longitude", "latitude")],
                                       cellnumbers = TRUE))) %>%
    cbind(as.data.frame(terra::extract(variables_aspect,
                                       occ[, c("longitude", "latitude")]))) %>%
    cbind(as.data.frame(terra::extract(variables_elev,
                                       occ[, c("longitude", "latitude")]))) %>%
    cbind(as.data.frame(terra::extract(variables_slope,
                                       occ[, c("longitude", "latitude")])))

  # Drop NA for key environmental variables
  occ_var <- occ_var %>% drop_na(bio12, tree_cover)
  message("  Sites after NA drop: ", n_distinct(occ_var$site))

  # -- Standardise covariates and SAVE scaling parameters --
  scaling_params <- list()
  occ_var_std <- occ_var
  for (v in vars_to_scale) {
    vals <- occ_var_std[[v]]
    scaling_params[[v]] <- list(center = mean(vals, na.rm = TRUE),
                                 scale  = sd(vals, na.rm = TRUE))
    occ_var_std[[v]] <- (vals - scaling_params[[v]]$center) / scaling_params[[v]]$scale
  }

  # Save scaling parameters for reuse in prediction (CRITICAL for consistency)
  saveRDS(scaling_params,
          here("data", "processed", sp, paste0(sp, "_scaling_params.rds")))
  message("  Scaling parameters saved to: ", sp, "_scaling_params.rds")

  # -- Convert from long to wide format, year by year --
  occ_wide_list <- list()
  for (yr in years) {
    occ_yr <- occ_var_std %>% filter(year == yr)
    if (nrow(occ_yr) == 0) {
      message("  WARNING: No data for year ", yr)
      next
    }
    occ_wide_yr <- format_unmarked_occu(
      occ_yr,
      site_id = "site",
      response = "species_observed",
      site_covs = site_covs_wide,
      obs_covs = obs_covs
    ) %>%
      rename_with(~ paste0(., ".", yr), -locality_id)
    occ_wide_list[[as.character(yr)]] <- occ_wide_yr
  }

  # Merge all years by locality_id
  occ_wide <- occ_wide_list[[1]]
  for (k in 2:length(occ_wide_list)) {
    occ_wide <- occ_wide %>%
      full_join(occ_wide_list[[k]], by = "locality_id", multiple = "all")
  }

  # Remove duplicated rows by locality_id
  occ_wide <- occ_wide %>% distinct(locality_id, .keep_all = TRUE)

  # Convert logical detection columns to numeric (TRUE/FALSE -> 1/0)
  cols_logical <- sapply(occ_wide, is.logical)
  occ_wide[, cols_logical] <- lapply(occ_wide[, cols_logical], as.numeric)

  # Consolidate site covariates: take the first non-NA value across years
  first_year <- years[1]
  for (v in c(vars_to_scale, "cells", "latitude", "longitude")) {
    base_col <- paste0(v, ".", first_year)
    if (base_col %in% names(occ_wide)) {
      occ_wide <- occ_wide %>% rename(!!v := !!base_col)
      # Coalesce with other years
      year_cols <- paste0(v, ".", years[-1])
      existing_cols <- year_cols[year_cols %in% names(occ_wide)]
      if (length(existing_cols) > 0) {
        occ_wide <- occ_wide %>%
          mutate(!!v := coalesce(!!sym(v), !!!syms(existing_cols)))
      }
    }
  }

  # Remove duplicated cells and drop NA
  occ_wide_clean <- occ_wide
  if ("cells" %in% names(occ_wide_clean)) {
    occ_wide_clean <- occ_wide_clean[!duplicated(occ_wide_clean$cells), ]
    occ_wide_clean <- occ_wide_clean %>% drop_na(cells)
  }

  pct_removed <- round((1 - nrow(occ_wide_clean) / nrow(occ_wide)) * 100, 1)
  message("  Rows after dedup: ", nrow(occ_wide_clean),
          " (removed ", pct_removed, "%)")

  # -- Save outputs --
  # Sites for GEE upload
  occ_wide_latlong <- occ_wide_clean %>% select(cells, latitude, longitude)
  write.csv(occ_wide_latlong,
            here("data", "processed", sp, paste0(sp, "_occ_wide_latlong.csv")),
            row.names = FALSE)

  # Full wide dataset with static variables
  write.csv(occ_wide_clean,
            here("data", "processed", sp, paste0(sp, "_occ_wide_static.csv")),
            row.names = FALSE)

  message("  Step 2 complete for ", sp)
}

message("\nStep 2 complete: static variables prepared for all species.")
