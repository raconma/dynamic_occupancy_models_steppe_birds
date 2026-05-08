###############################################################################
# 3b_add_stepRep.R
#
# Purpose: Add yearly steppe-representativeness columns to each species'
#          {sp}_occ_wide_dynamic.csv so that the stepRep covariate enters
#          the colext / stPGOcc fits in step 4 with no data wrangling
#          beyond a 4-line edit to scripts/4_occupancy_models.R.
#
# Inputs:  data/processed_2023/{sp}/{sp}_occ_wide_dynamic.csv  (from step 3)
#          data/derived/stepRep_cellyear_{sp}.csv               (build_stepRep.R)
#
# Outputs: data/processed_2023/{sp}/{sp}_occ_wide_dynamic.csv (overwritten in
#          place; idempotent: re-running with stepRep_*_<year> columns already
#          present is a no-op).
#
# Columns added: 4 variants x 7 years = 28
#   stepRep_strict_500m_2017..2023   (PRIMARY recommended for fits)
#   stepRep_strict_1km_2017..2023    (sensitivity: larger buffer)
#   stepRep_broad_500m_2017..2023    (sensitivity: includes dehesa)
#   stepRep_broad_1km_2017..2023     (sensitivity: dehesa + larger buffer)
#
# NA handling: cell-years missing from stepRep_cellyear (i.e. cells that
#   have no checklists in some year but are still in the model table because
#   they have visits in other years) are imputed first by row mean across
#   the years the cell does have data, then by global median if a row is
#   entirely NA (very rare; should not occur given how the pipeline
#   constructs sites).
###############################################################################

suppressPackageStartupMessages({
  library(here)
  library(data.table)
})

YEARS         <- 2017:2023
species_codes <- c("otitar", "ptealc", "pteori", "tettet")
stepRep_vars  <- c("stepRep_strict_500m", "stepRep_strict_1km",
                   "stepRep_broad_500m",  "stepRep_broad_1km")

for (sp in species_codes) {
  message("\n=== ", sp, " ===")

  in_csv   <- here("data", "processed_2023", sp,
                   paste0(sp, "_occ_wide_dynamic.csv"))
  step_csv <- here("data", "derived",
                   paste0("stepRep_cellyear_", sp, ".csv"))

  if (!file.exists(in_csv) || !file.exists(step_csv)) {
    warning("missing input for ", sp, " — skipping"); next
  }

  dt   <- fread(in_csv)
  step <- fread(step_csv)

  expected <- as.vector(outer(stepRep_vars, YEARS,
                              function(v, y) paste0(v, "_", y)))
  if (all(expected %in% names(dt))) {
    message("  stepRep columns already present (", length(expected),
            " cols) — no-op."); next
  }

  # Drop any partial set of stepRep columns so we always rebuild cleanly
  drop_cols <- intersect(expected, names(dt))
  if (length(drop_cols) > 0) {
    message("  dropping ", length(drop_cols),
            " stale stepRep columns before rebuild")
    dt[, (drop_cols) := NULL]
  }

  for (v in stepRep_vars) {
    w <- dcast(step, cells ~ year, value.var = v)
    setnames(w, as.character(YEARS), paste0(v, "_", YEARS))
    dt <- merge(dt, w, by = "cells", all.x = TRUE, sort = FALSE)

    cy  <- paste0(v, "_", YEARS)
    mat <- as.matrix(dt[, ..cy])
    n_pre <- sum(is.na(mat))
    rm    <- rowMeans(mat, na.rm = TRUE)
    for (j in seq_len(ncol(mat))) {
      m <- is.na(mat[, j]); mat[m, j] <- rm[m]
    }
    n_mid <- sum(is.na(mat))
    if (n_mid > 0) {
      mat[is.na(mat)] <- median(step[[v]], na.rm = TRUE)
    }
    n_post <- sum(is.na(mat))
    dt[, (cy) := as.data.frame(mat)]
    message(sprintf("  %-20s : %d NA -> %d after row-mean -> %d after global median",
                    v, n_pre, n_mid, n_post))
  }

  fwrite(dt, in_csv)
  message("  saved: ", in_csv,
          "  (", nrow(dt), " sites x ", ncol(dt), " cols)")
}

message("\nStep 3b complete: stepRep yearly site covariates merged into the model tables.")
