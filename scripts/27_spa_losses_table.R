###############################################################################
# 27_spa_losses_table.R
#
# Which Special Protection Areas hold the recent losses? For each species and
# SPA: cells with >= 1 eBird detection whose centroid lies in the SPA, their
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
  pts <- st_as_sf(cells, coords = c("x", "y"), crs = UTM30)
  hit <- st_intersects(pts, spa)
  cells$spa_id <- vapply(hit, function(i) if (length(i)) spa$localId[i[1]] else NA_character_, "")
  cells$steppe_strict <- hab$steppe_strict[match(cells$cells, hab$cells)]
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
