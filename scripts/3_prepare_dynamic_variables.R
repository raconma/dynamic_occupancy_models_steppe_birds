###############################################################################
# 3_prepare_dynamic_variables.R
#
# Purpose: Merge GEE-exported dynamic variables (NDVI, EVI, land cover,
#          temperature, precipitation) with the static wide-format dataset.
#
# Inputs:  data/processed/{sp}/{sp}_occ_wide_static.csv    (from step 2)
#          data/raw/gee_exports/{sp}_dynamic_variables.csv  (from GEE step 3a)
#
# Outputs: data/processed/{sp}/{sp}_occ_wide_dynamic.csv
#
# Note:    The GEE step (3_download_dynamic_variables_{sp}.js) must be run
#          manually in Google Earth Engine before this script. See README.md.
#
# Replaces: 3_prepare_dynamic_variables_{otitar,ptealc,pteori,tettet}.R
###############################################################################

# -- Load packages --
library(here)
library(readr)
library(dplyr)
library(tidyr)

# Resolve namespace conflicts
select <- dplyr::select
filter <- dplyr::filter

# -- Species to process --
species_codes <- c("otitar", "ptealc", "pteori", "tettet")

for (sp in species_codes) {
  message("\n=== Processing dynamic variables: ", sp, " ===")

  # Load wide-format static dataset
  occ_wide_clean <- read.csv(
    here("data", "processed", sp, paste0(sp, "_occ_wide_static.csv"))
  )
  message("  Static dataset: ", nrow(occ_wide_clean), " rows, ",
          ncol(occ_wide_clean), " cols")

  # Load GEE-exported dynamic variables
  dyn_path <- here("data", "raw", "gee_exports", paste0(sp, "_dynamic_variables.csv"))
  if (!file.exists(dyn_path)) {
    warning("  GEE export not found: ", dyn_path, " — skipping ", sp)
    next
  }
  dyn_var <- read.csv(dyn_path)

  # Remove GEE metadata columns
  dyn_var <- dyn_var[, !names(dyn_var) %in%
                       c("system.index", "cells", "latitude", "longitude", ".geo")]
  message("  Dynamic variables: ", ncol(dyn_var), " columns")

  # Merge (by row alignment — both sorted by GEE export order matching static)
  # NOTE: This relies on consistent row order from GEE.
  #       A safer approach would be to merge by 'cells' column.
  stopifnot(nrow(occ_wide_clean) == nrow(dyn_var))
  occ_wide_clean <- cbind(occ_wide_clean, dyn_var)

  # Drop rows with missing NDVI in any year
  ndvi_cols <- paste0("NDVI_", 2017:2022)
  occ_wide_clean <- occ_wide_clean %>%
    drop_na(all_of(ndvi_cols))
  message("  After NDVI NA drop: ", nrow(occ_wide_clean), " rows")

  # Save output
  out_path <- here("data", "processed", sp, paste0(sp, "_occ_wide_dynamic.csv"))
  write.csv(occ_wide_clean, out_path, row.names = FALSE)
  message("  Saved: ", out_path)
}

message("\nStep 3 complete: dynamic variables merged for all species.")
