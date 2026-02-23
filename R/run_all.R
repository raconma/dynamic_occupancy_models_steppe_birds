###############################################################################
# R/run_all.R
#
# Purpose: ONE COMMAND to reproduce the full analysis pipeline.
#
# Usage:   Rscript R/run_all.R
#
# Prerequisites:
#   1. Raw data placed in data/raw/ (see data-raw/get_data.R)
#   2. GEE dynamic variables exported (see README.md, step 3a)
#   3. R packages installed via renv::restore()
#
# Pipeline:
#   Step 1: Filter eBird data per species
#   Step 2: Prepare static variables (bioclimatic, topographic)
#   Step 3: Merge GEE-exported dynamic variables
#   Step 4: Fit colext dynamic occupancy models + maps + simulations
#   Step 5: Validate predictions against Spanish Biodiversity Atlas
#
# All outputs go to: figs/, results/, data/processed/
###############################################################################

library(here)

# -- Record start time --
t_start <- Sys.time()

# -- Ensure output directories exist --
dir.create(here("figs"), showWarnings = FALSE)
dir.create(here("results"), showWarnings = FALSE)

cat("\n")
cat(strrep("=", 70), "\n")
cat("  DYNAMIC OCCUPANCY MODELS FOR STEPPE BIRDS — FULL PIPELINE\n")
cat("  Started:", format(t_start, "%Y-%m-%d %H:%M:%S"), "\n")
cat(strrep("=", 70), "\n\n")

# ---- Step 1: Prepare eBird database ----
cat(">>> STEP 1: Preparing eBird database...\n")
source(here("scripts", "1_prepare_database_ebird.R"), local = new.env())

# ---- Step 2: Prepare static variables ----
cat("\n>>> STEP 2: Preparing static variables...\n")
source(here("scripts", "2_prepare_static_variables.R"), local = new.env())

# ---- Step 3: Prepare dynamic variables ----
cat("\n>>> STEP 3: Merging dynamic variables from GEE...\n")
source(here("scripts", "3_prepare_dynamic_variables.R"), local = new.env())

# ---- Step 4: Fit occupancy models ----
cat("\n>>> STEP 4: Fitting dynamic occupancy models...\n")
source(here("scripts", "4_occupancy_models.R"), local = new.env())

# ---- Step 5: Validation ----
cat("\n>>> STEP 5: Validating model predictions...\n")
source(here("scripts", "5_validation.R"), local = new.env())

# ---- Save session info ----
session_file <- here("results", "session_info.txt")
sink(session_file)
cat("Session Info\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
sessionInfo()
sink()
cat("\nSession info saved to:", session_file, "\n")

# ---- Summary ----
t_end <- Sys.time()
elapsed <- difftime(t_end, t_start, units = "mins")

cat("\n")
cat(strrep("=", 70), "\n")
cat("  PIPELINE COMPLETE\n")
cat("  Elapsed time:", round(as.numeric(elapsed), 1), "minutes\n")
cat("  Outputs:\n")
cat("    figs/      — all figures\n")
cat("    results/   — model summaries, GOF, validation, session info\n")
cat("    data/processed/ — intermediate data files\n")
cat(strrep("=", 70), "\n")
