###############################################################################
# 37_checklist_dates_by_year.R
#
# Closure check: did the seasonal timing of checklists change over the study
# period? Checklist numbers rose by 200-250 % between 2017 and 2023; if the
# added checklists fell disproportionately late (after 1 June, when male great
# bustards start leaving breeding areas and sandgrouse broods disperse),
# within-season movements would be confounded with year and could mimic
# local extinction. The detection model has no day-of-year term, so this is
# checked descriptively.
#
# For each species group, complete checklists within its breeding window,
# with the effort filters of the paper (<= 5 h, <= 5 km, <= 10 observers),
# are summarised by year: median and quartiles of day of year and the share
# after 1 June (day 152). Also reported for checklists in cells where the
# species was recorded, the cells that carry the anchored transitions.
#
# Inputs:  data/processed_2023/{sp}/ebd_{sp}_breeding_spain_zf.csv
#          data/processed_2023/{sp}/{sp}_occ_wide_dynamic.csv (recorded cells)
# Output:  results/manuscript_tables/tableS_checklist_dates.csv
###############################################################################

library(here)
library(terra)

source(here("R", "spatial_layers.R"))
set.seed(42)  # no stochastic steps

species <- c("otitar", "tettet", "pteori", "ptealc")
grid <- wc_grid()
rows <- list()
for (sp in species) {
  message("=== ", sp)
  d <- read.csv(here("data", "processed_2023", sp, paste0("ebd_", sp, "_breeding_spain_zf.csv")))
  d <- d[d$all_species_reported & d$year >= 2017 & d$year <= 2023 &
           !is.na(d$duration_minutes) & d$duration_minutes <= 300 &
           (is.na(d$effort_distance_km) | d$effort_distance_km <= 5) &
           d$number_observers <= 10, ]
  d$cells <- terra::cellFromXY(grid, cbind(d$longitude, d$latitude))
  occ <- read.csv(here("data", "processed_2023", sp, paste0(sp, "_occ_wide_dynamic.csv")))
  y <- as.matrix(occ[, grep("^y\\.[0-9]+\\.[0-9]{4}$", names(occ))])
  rec <- occ$cells[rowSums(y == TRUE | y == "TRUE", na.rm = TRUE) > 0]
  for (set in c("all", "recorded_cells")) {
    x <- if (set == "all") d[d$cells %in% occ$cells, ] else d[d$cells %in% rec, ]
    for (yr in sort(unique(x$year))) {
      doy <- x$day_of_year[x$year == yr]
      rows[[paste(sp, set, yr)]] <- data.frame(
        species = sp, cells = set, year = yr, n_checklists = length(doy),
        doy_q25 = quantile(doy, .25, names = FALSE), doy_median = median(doy),
        doy_q75 = quantile(doy, .75, names = FALSE), share_after_1june = mean(doy >= 152))
    }
  }
}
out <- do.call(rbind, rows)
write.csv(out, here("results", "manuscript_tables", "tableS_checklist_dates.csv"), row.names = FALSE)
print(out, digits = 3)

sessionInfo()
