###############################################################################
# 27_spa_losses_table.R
#
# Which Special Protection Areas hold the recent losses? For each species and
# SPA: SPA cells as defined in script 23 (cells with >= 1 eBird detection and
# >= 50 % of their area inside SPAs), each assigned to the SPA covering the
# largest part of the cell, with their
# posterior net change 2018-2023 on data-anchored transitions (net_both,
# script 21), all recent transitions (net_recent), observed (naive)
# extinctions and colonisations, and mean CORINE pseudo-steppe cover
# (script 23). Intended as input for the ecological interpretation of the
# SPA results; cell counts per SPA are small, so it describes where the
# signal sits rather than testing SPA-level effects.
#
# Inputs:  results/net_change_maps/{sp}_cells.csv (21)
#          results/core_spa/cell_habitat.csv      (23)
#          SPAs via R/spatial_layers.R
# Outputs: results/core_spa/spa_losses_by_site.csv
###############################################################################

library(here)
library(sf)

source(here("R", "spatial_layers.R"))
set.seed(42)  # no stochastic steps

species <- c("otitar", "ptealc", "pteori", "tettet")
spa <- load_spa_sites()
hab <- read.csv(here("results", "core_spa", "cell_habitat.csv"))

rows <- list()
for (sp in species) {
  cells <- read.csv(here("results", "net_change_maps", paste0(sp, "_cells.csv")))
  cells <- cells[cells$ever_detected, ]
  cells$steppe_strict <- hab$steppe_strict[match(cells$cells, hab$cells)]
  cells$spa_frac <- hab$spa_frac[match(cells$cells, hab$cells)]
  cells <- cells[cells$spa_frac >= 0.5, ]   # same SPA cells as class_rates.csv
  polys <- wc_cell_polygons(cells$cells)
  ov <- suppressWarnings(st_intersection(polys, spa[, "localId"]))
  ov$area <- as.numeric(st_area(ov))
  ov <- ov[order(ov$cells, -ov$area), ]
  best <- ov[!duplicated(ov$cells), ]
  cells$spa_id <- best$localId[match(cells$cells, best$cells)]
  in_spa <- cells[!is.na(cells$spa_id), ]
  if (!nrow(in_spa)) next
  agg <- do.call(rbind, lapply(split(in_spa, in_spa$spa_id), function(d) {
    data.frame(species = sp, spa_code = d$spa_id[1],
               spa_name = spa$SOName[match(d$spa_id[1], spa$localId)],
               n_detected_cells = nrow(d),
               net_both = sum(d$net_both), net_recent = sum(d$net_recent),
               obs_ext = sum(d$obs_ext), obs_col = sum(d$obs_col),
               mean_steppe_strict = mean(d$steppe_strict, na.rm = TRUE))
  }))
  rows[[sp]] <- agg[order(agg$net_both), ]
}

out <- do.call(rbind, rows); rownames(out) <- NULL
write.csv(out, here("results", "core_spa", "spa_losses_by_site.csv"), row.names = FALSE)
options(width = 200)
for (sp in species) {
  cat("\n==", sp, "== (10 SPAs with the largest data-anchored net loss)\n")
  print(head(out[out$species == sp, -1], 10), digits = 2, row.names = FALSE)
}

sessionInfo()
