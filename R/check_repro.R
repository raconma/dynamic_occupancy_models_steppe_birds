###############################################################################
# R/check_repro.R
#
# Purpose: Smoke test — verify that the environment is ready to run the
#          analysis pipeline. Checks packages, data files, paths, and
#          runs a minimal model fit.
#
# Usage:   Rscript R/check_repro.R
#
# Exit codes: 0 = all checks passed, 1 = one or more checks failed
###############################################################################

library(here)

cat("\n=== REPRODUCIBILITY CHECK ===\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Working directory:", getwd(), "\n")
cat("here() resolves to:", here(), "\n\n")

n_errors <- 0
n_warnings <- 0

check <- function(description, condition) {
  if (condition) {
    cat("  [PASS] ", description, "\n")
  } else {
    cat("  [FAIL] ", description, "\n")
    n_errors <<- n_errors + 1
  }
}

warn_check <- function(description, condition) {
  if (condition) {
    cat("  [PASS] ", description, "\n")
  } else {
    cat("  [WARN] ", description, "\n")
    n_warnings <<- n_warnings + 1
  }
}

# ---- 1. Check R version ----
cat("--- R Version ---\n")
check("R >= 4.0", as.numeric(R.version$major) >= 4)
cat("  Actual:", R.version.string, "\n\n")

# ---- 2. Check required packages ----
cat("--- Required Packages ---\n")
required_pkgs <- c(
  "here", "unmarked", "auk", "terra", "raster", "sf",
  "dplyr", "tidyr", "readr", "ggplot2", "purrr",
  "AICcmodavg", "MuMIn", "rnaturalearth", "rnaturalearthdata",
  "rmapshaper", "lubridate", "pROC", "Metrics", "spdep", "blockCV",
  "psych", "HH", "pander", "knitr", "stringr", "ncdf4", "gridExtra"
)

for (pkg in required_pkgs) {
  installed <- requireNamespace(pkg, quietly = TRUE)
  check(paste0("Package '", pkg, "' installed"), installed)
}

# ---- 3. Check system dependencies (GDAL/GEOS/PROJ) ----
cat("\n--- System Dependencies (sf/terra) ---\n")
if (requireNamespace("sf", quietly = TRUE)) {
  sf_ext <- sf::sf_extSoftVersion()
  cat("  GDAL:", sf_ext["GDAL"], "\n")
  cat("  GEOS:", sf_ext["GEOS"], "\n")
  cat("  PROJ:", sf_ext["proj.4"], "\n")
  check("GDAL available", nchar(sf_ext["GDAL"]) > 0)
  check("GEOS available", nchar(sf_ext["GEOS"]) > 0)
  check("PROJ available", nchar(sf_ext["proj.4"]) > 0)
}

# ---- 4. Check directory structure ----
cat("\n--- Directory Structure ---\n")
required_dirs <- c(
  "scripts", "R", "data", "data/raw", "data/processed",
  "figs", "results", "data-raw"
)
for (d in required_dirs) {
  check(paste0("Directory '", d, "/' exists"), dir.exists(here(d)))
}

# ---- 5. Check data files ----
cat("\n--- Data Files ---\n")
data_files <- list(
  "eBird EBD file" = here("data", "raw", "ebird_raw_mar2024",
                           "ebd_ES_smp_relMar-2024.txt"),
  "eBird sampling file" = here("data", "raw", "ebird_raw_mar2024",
                                "ebd_sampling_relMar-2024.txt"),
  "Bioclimatic raster" = here("data", "raw", "environmental_data",
                               "environmental_data_occ", "variables_spain.grd"),
  "Topography (aspect)" = here("data", "raw", "topology_data", "topo_aspect.asc"),
  "Topography (elevation)" = here("data", "raw", "topology_data", "topo_elev.asc"),
  "Topography (slope)" = here("data", "raw", "topology_data", "topo_slope.asc"),
  "Validation atlas" = here("data", "raw", "validation", "atlas_biodiversidad",
                             "aves_spain.shp")
)

for (name in names(data_files)) {
  warn_check(paste0(name, ": ", basename(data_files[[name]])),
             file.exists(data_files[[name]]))
}

# Check GEE exports
cat("\n--- GEE Dynamic Variable Exports ---\n")
for (sp in c("otitar", "ptealc", "pteori", "tettet")) {
  gee_file <- here("data", "raw", "gee_exports",
                    paste0(sp, "_dynamic_variables.csv"))
  warn_check(paste0("GEE export: ", sp), file.exists(gee_file))
}

# ---- 6. Check key scripts exist ----
cat("\n--- Pipeline Scripts ---\n")
scripts <- c(
  "scripts/1_prepare_database_ebird.R",
  "scripts/2_prepare_static_variables.R",
  "scripts/3_prepare_dynamic_variables.R",
  "scripts/4_occupancy_models.R",
  "scripts/5_validation.R",
  "R/model_configs.R",
  "R/run_all.R"
)
for (s in scripts) {
  check(paste0("Script: ", s), file.exists(here(s)))
}

# ---- 7. Quick functional test (if unmarked available) ----
cat("\n--- Quick Functional Test (unmarked) ---\n")
if (requireNamespace("unmarked", quietly = TRUE)) {
  tryCatch({
    library(unmarked)
    # Create a tiny test dataset
    set.seed(1)
    y_test <- matrix(rbinom(30, 1, 0.5), nrow = 5, ncol = 6)
    umf_test <- unmarkedMultFrame(y = y_test, numPrimary = 3)
    m_test <- colext(~ 1, ~ 1, ~ 1, ~ 1, data = umf_test)
    check("unmarked::colext() runs successfully", !is.null(m_test))
  }, error = function(e) {
    cat("  [FAIL] colext() test failed:", e$message, "\n")
    n_errors <<- n_errors + 1
  })
} else {
  cat("  [SKIP] unmarked not installed\n")
}

# ---- Summary ----
cat("\n", strrep("=", 50), "\n")
cat("  RESULTS:  ", n_errors, " errors  |  ", n_warnings, " warnings\n")
if (n_errors == 0 && n_warnings == 0) {
  cat("  STATUS:   ALL CHECKS PASSED\n")
} else if (n_errors == 0) {
  cat("  STATUS:   READY (with warnings — likely missing data files)\n")
} else {
  cat("  STATUS:   NOT READY — fix errors above before running pipeline\n")
}
cat(strrep("=", 50), "\n")

# Exit with appropriate code
if (n_errors > 0) quit(status = 1, save = "no")
