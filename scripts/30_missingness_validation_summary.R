###############################################################################
# 30_missingness_validation_summary.R
#
# Two checks requested in the internal review of the Biological Conservation
# draft:
#
# 1. Informative missingness. Anchored transitions require surveys in both
#    years. If observers stop visiting a cell after the species disappears,
#    anchored transitions under-represent extinctions and bias net change
#    towards zero. In cells with >= 1 detection, compare the probability that
#    a cell surveyed in year t is surveyed again in t+1 when the species was
#    recorded in t versus when it was not (recent transitions, 2018 -> 2023).
#    Proportions with 95% Wilson-type intervals from prop.test().
#
# 2. Validation summary. Means over the five spatial folds of the
#    atlas validation written by scripts/5_validation.R
#    (results/{sp}_validation_cv.csv; II Atlas, aves_spain.shp; 270 km
#    blocks, k = 5).
#
# Outputs: results/manuscript_tables/tableS_missingness.csv
#          results/manuscript_tables/tableS_validation_summary.csv
###############################################################################

library(here)
library(unmarked)

source(here("R", "colext_utils.R"))
set.seed(42)  # no stochastic steps

species <- c("otitar", "ptealc", "pteori", "tettet")
tab_dir <- here("results", "manuscript_tables")
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

miss_rows <- list()
for (sp in species) {
  D <- colext_design(readRDS(colext_model_path(sp)))
  n_t <- D$n_t
  naive <- sapply(1:n_t, function(t) {
    yt <- D$y[, ((t - 1) * D$J + 1):(t * D$J), drop = FALSE]
    ifelse(rowSums(!is.na(yt)) > 0, as.integer(rowSums(yt, na.rm = TRUE) > 0),
           NA_integer_)
  })
  det <- rowSums(naive, na.rm = TRUE) > 0
  recent <- 2:(n_t - 1)
  z0 <- naive[det, recent, drop = FALSE]
  s1 <- !is.na(naive[det, recent + 1, drop = FALSE])
  surveyed0 <- !is.na(z0)
  for (st in c(1, 0)) {
    w <- surveyed0 & z0 == st
    k <- sum(s1[w]); n <- sum(w)
    ci <- suppressWarnings(prop.test(k, n)$conf.int)
    miss_rows[[paste(sp, st)]] <- data.frame(
      species = sp, state_t = if (st == 1) "recorded" else "not recorded",
      n_cell_years = n, n_surveyed_next = k, p_surveyed_next = k / n,
      lo = ci[1], hi = ci[2])
  }
}
miss <- do.call(rbind, miss_rows); rownames(miss) <- NULL
write.csv(miss, file.path(tab_dir, "tableS_missingness.csv"), row.names = FALSE)

val <- do.call(rbind, lapply(species, function(sp) {
  v <- read.csv(here("results", paste0(sp, "_validation_cv.csv")))
  data.frame(species = sp, n_folds = nrow(v),
             AUC = mean(v$AUC), uAUC = mean(v$uAUC), TSS = mean(v$TSS),
             AUC_min = min(v$AUC), AUC_max = max(v$AUC))
}))
write.csv(val, file.path(tab_dir, "tableS_validation_summary.csv"), row.names = FALSE)

options(width = 200)
print(miss, digits = 3)
print(val, digits = 3)
sessionInfo()
