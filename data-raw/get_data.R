###############################################################################
# data-raw/get_data.R
#
# Purpose: Document how to obtain all raw data inputs.
#          Some data cannot be downloaded programmatically (eBird requires
#          a signed license agreement). This script downloads what it can
#          and provides clear instructions for the rest.
#
# Run this BEFORE the analysis pipeline.
###############################################################################

library(here)

cat("
==============================================================================
  DATA ACQUISITION GUIDE
  Dynamic Occupancy Models for Steppe Birds in Spain
==============================================================================

This project requires several datasets. Below are instructions to obtain each.

--- 1. eBird Data (MANUAL DOWNLOAD REQUIRED) ---

  Source: https://ebird.org/data/download
  License: eBird Basic Dataset Terms of Use (cannot be redistributed)

  Steps:
    a) Create an account at https://ebird.org
    b) Request access to the eBird Basic Dataset (EBD)
    c) Download the Spain sampling dataset (ebd_ES_smp_relMar-2024.txt)
       and the sampling events file (ebd_sampling_relMar-2024.txt)
    d) Place both files in:
         data/raw/ebird_raw_mar2024/ebd_ES_smp_relMar-2024.txt
         data/raw/ebird_raw_mar2024/ebd_sampling_relMar-2024.txt

  Approximate size: ~38 GB (uncompressed)

--- 2. Bioclimatic Variables ---

  Source: WorldClim v2.1 (https://www.worldclim.org/data/worldclim21.html)
  These should be the BioClim variables for Spain, stacked as a RasterStack.

  Place in: data/raw/environmental_data/environmental_data_occ/variables_spain.grd
            data/raw/environmental_data/environmental_data_occ/variables_spain.gri

--- 3. Topographic Variables ---

  Source: SRTM-derived (slope, aspect, elevation) for Spain

  Place in:
    data/raw/topology_data/topo_aspect.asc
    data/raw/topology_data/topo_elev.asc
    data/raw/topology_data/topo_slope.asc
    data/raw/topology_data/topo_aspect_masked.tif
    data/raw/topology_data/topo_elev_masked.tif

--- 4. Google Earth Engine Dynamic Variables (MANUAL STEP) ---

  After running Step 2 (which produces {species}_occ_wide_latlong.csv),
  you must run the GEE scripts to extract dynamic variables:

  a) Go to https://code.earthengine.google.com/
  b) Upload each {species}_occ_wide_latlong.csv as a GEE asset
  c) Run scripts/3_download_dynamic_variables_{species}.js for each species
     (update the asset path in the JS file to match your upload)
  d) Export results from GEE to Google Drive
  e) Download and place in:
       data/raw/gee_exports/otitar_dynamic_variables.csv
       data/raw/gee_exports/ptealc_dynamic_variables.csv
       data/raw/gee_exports/pteori_dynamic_variables.csv
       data/raw/gee_exports/tettet_dynamic_variables.csv

--- 5. Spanish Biodiversity Atlas (for validation) ---

  Source: MITECO / Inventario Espanol del Patrimonio Natural y la Biodiversidad

  Place in: data/raw/validation/atlas_biodiversidad/aves_spain.shp
            (with companion .dbf, .shx, .prj files)

==============================================================================
  After placing all data, run the pipeline:

    Rscript R/run_all.R

  Or verify your setup first:

    Rscript R/check_repro.R
==============================================================================
\n")
