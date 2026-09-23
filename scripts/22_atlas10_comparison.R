###############################################################################
# 22_atlas10_comparison.R
#
# Aggregate the 5x5 km posterior quantities of script 21 to the 10x10 km UTM
# atlas grid and compare 2017-2023 status with the II Breeding Bird Atlas
# (fieldwork 1998-2002; aves_spain.shp, SEO98 grid). The file matches the
# MITECO "Inventario Nacional de Biodiversidad 2007, Aves Reproductoras"
# dataset in GBIF exactly for O. tarda (528 squares); it is NOT the III Atlas
# (2014-2018), whose square-level data are not public.
#
# Per 10x10 square (only squares containing >= 1 modelled 5x5 cell):
#   p_occ_2023   1 - prod(1 - P(z = 1 in 2023 | data)) over modelled subcells
#                (a lower bound: unmodelled subcells contribute nothing)
#   net_recent   posterior colonisations - extinctions 2018-2023, cells with
#                >= 1 detection only (see script 21 on never-detected cells)
#   effort       surveyed cell-years summed over modelled subcells
# Status vs Atlas II, for squares with effort >= MIN_EFFORT:
#   retained   Atlas II present, p_occ_2023 >= 0.5
#   uncertain  Atlas II present, 0.1 <= p_occ_2023 < 0.5
#   lost       Atlas II present, p_occ_2023 < 0.1
#   new        Atlas II absent,  p_occ_2023 >= 0.5
# Sensitivity: MIN_EFFORT = 3 vs 10 surveyed cell-years.
#
# Inputs:
#   results/net_change_maps/{sp}_cells.csv   (script 21)
#   data-raw/data/validation/atlas_biodiversidad/atlas_biodiversidad/aves_spain.shp
#
# Outputs (new):
#   results/atlas10/{sp}_squares.csv
#   results/atlas10/status_vs_atlas2.csv
#   results/atlas10/fig_status_vs_atlas2.png
#   results/atlas10/fig_net_change_10km.png
###############################################################################

library(here)
library(sf)
library(ggplot2)

set.seed(42)  # no stochastic steps; kept for workflow consistency

species <- c("otitar", "ptealc", "pteori", "tettet")
species_names <- c(otitar = "O. tarda", ptealc = "P. alchata",
                   pteori = "P. orientalis", tettet = "T. tetrax")
atlas_col <- c(otitar = "OTITAR", ptealc = "PTEALC", pteori = "PTEORI",
               tettet = "TETTET")
UTM30 <- 25830
MIN_EFFORT <- c(3, 10)
out_dir <- here("results", "atlas10")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

atlas2 <- st_transform(st_read(here("data-raw", "data", "validation",
                                    "atlas_biodiversidad", "atlas_biodiversidad",
                                    "aves_spain.shp"), quiet = TRUE), UTM30)
atlas2$sq_id <- seq_len(nrow(atlas2))

spain <- NULL
if (requireNamespace("rnaturalearth", quietly = TRUE)) {
  spain <- st_transform(rnaturalearth::ne_countries(
    scale = 50, country = "Spain", returnclass = "sf"), UTM30)
}

classify <- function(a2, p, effort, min_effort) {
  s <- ifelse(effort < min_effort, "insufficient effort",
       ifelse(a2 & p >= 0.5, "retained",
       ifelse(a2 & p >= 0.1, "uncertain",
       ifelse(a2, "lost",
       ifelse(p >= 0.5, "new", "absent")))))
  factor(s, levels = c("retained", "uncertain", "lost", "new", "absent",
                       "insufficient effort"))
}

sq_all <- list()
status_rows <- list()

for (sp in species) {
  cells <- read.csv(here("results", "net_change_maps", paste0(sp, "_cells.csv")))
  pts <- st_as_sf(cells, coords = c("x", "y"), crs = UTM30, remove = FALSE)
  cells$sq_id <- st_join(pts, atlas2[, "sq_id"], join = st_within)$sq_id
  cells <- cells[!is.na(cells$sq_id), ]

  agg <- do.call(rbind, lapply(split(cells, cells$sq_id), function(d) {
    data.frame(sq_id = d$sq_id[1], n_cells = nrow(d),
               any_detected = any(d$ever_detected),
               effort = sum(d$n_years_surveyed),
               p_occ_2023 = 1 - prod(1 - d$z1_2023),
               net_recent = sum(d$net_recent[d$ever_detected]),
               ext_full = sum(d$ext_full[d$ever_detected]),
               col_full = sum(d$col_full[d$ever_detected]),
               obs_ext = sum(d$obs_ext), obs_col = sum(d$obs_col))
  }))
  agg$atlas2 <- as.numeric(st_drop_geometry(atlas2)[agg$sq_id, atlas_col[sp]]) > 0
  for (m in MIN_EFFORT) {
    agg[[paste0("status_eff", m)]] <- classify(agg$atlas2, agg$p_occ_2023,
                                               agg$effort, m)
  }
  agg$species <- sp
  write.csv(agg, file.path(out_dir, paste0(sp, "_squares.csv")), row.names = FALSE)
  sq_all[[sp]] <- agg

  n_a2_total <- sum(as.numeric(st_drop_geometry(atlas2)[[atlas_col[sp]]]) > 0)
  for (m in MIN_EFFORT) {
    st <- agg[[paste0("status_eff", m)]]
    a2_cl <- st[agg$atlas2 & st != "insufficient effort"]
    status_rows[[paste(sp, m)]] <- data.frame(
      species = sp, min_effort = m,
      atlas2_squares_total = n_a2_total,
      atlas2_squares_modelled = sum(agg$atlas2),
      atlas2_squares_classified = length(a2_cl),
      retained = sum(a2_cl == "retained"), uncertain = sum(a2_cl == "uncertain"),
      lost = sum(a2_cl == "lost"),
      pct_lost = round(100 * mean(a2_cl == "lost"), 1),
      pct_retained = round(100 * mean(a2_cl == "retained"), 1),
      new = sum(st == "new"),
      # Recent (2018-23) posterior net change split by Atlas II status
      net_recent_in_atlas2 = round(sum(agg$net_recent[agg$atlas2]), 2),
      net_recent_outside_atlas2 = round(sum(agg$net_recent[!agg$atlas2]), 2),
      obs_col_in_atlas2 = sum(agg$obs_col[agg$atlas2]),
      obs_col_outside_atlas2 = sum(agg$obs_col[!agg$atlas2])
    )
  }
}

status_df <- do.call(rbind, status_rows); rownames(status_df) <- NULL
write.csv(status_df, file.path(out_dir, "status_vs_atlas2.csv"), row.names = FALSE)
options(width = 220)
print(status_df)

###############################################################################
# MAPS
###############################################################################

map_df <- do.call(rbind, lapply(names(sq_all), function(sp) {
  d <- merge(atlas2[, "sq_id"], sq_all[[sp]], by = "sq_id")
  d$species <- factor(species_names[sp], levels = species_names)
  d
}))
a2_outline <- do.call(rbind, lapply(species, function(sp) {
  d <- atlas2[as.numeric(st_drop_geometry(atlas2)[[atlas_col[sp]]]) > 0, "sq_id"]
  d$species <- factor(species_names[sp], levels = species_names)
  d
}))

frame <- function(title) {
  list(coord_sf(crs = UTM30, xlim = c(-60000, 1130000),
                ylim = c(3980000, 4880000), expand = FALSE),
       facet_wrap(~ species), labs(title = title, x = NULL, y = NULL),
       theme_bw(base_size = 9),
       theme(strip.text = element_text(face = "italic"),
             axis.text = element_blank(), axis.ticks = element_blank(),
             panel.grid = element_blank(), legend.position = "bottom"))
}
base <- ggplot()
if (!is.null(spain)) base <- base + geom_sf(data = spain, fill = "grey97",
                                            colour = "grey60", linewidth = 0.2)

status_cols <- c(retained = "#1B7837", uncertain = "#F1B6DA", lost = "#C51B7D",
                 new = "#2166AC", absent = "grey93",
                 "insufficient effort" = "grey80")
p_status <- base +
  geom_sf(data = map_df, aes(fill = status_eff3), colour = NA) +
  scale_fill_manual(values = status_cols, name = NULL, drop = FALSE) +
  frame(paste("Status 2023 vs II Atlas (1998-2002), 10x10 km;",
                    "squares with >= 1 modelled cell and >= 3 surveyed cell-years"))
ggsave(file.path(out_dir, "fig_status_vs_atlas2.png"), p_status,
       width = 8, height = 7.5, dpi = 300)

det_sq <- map_df[map_df$any_detected, ]
lim <- max(abs(det_sq$net_recent), na.rm = TRUE)
p_net <- base +
  geom_sf(data = a2_outline, fill = "grey85", colour = NA) +
  geom_sf(data = det_sq, aes(fill = net_recent), colour = "grey40",
          linewidth = 0.1) +
  scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC",
                       midpoint = 0, limits = c(-lim, lim),
                       name = "Posterior colonisations - extinctions, 2018-2023") +
  frame(paste("Recent net change in 10x10 km squares with eBird detections",
                    "(grey: II Atlas range)"))
ggsave(file.path(out_dir, "fig_net_change_10km.png"), p_net,
       width = 8, height = 7.5, dpi = 300)

sessionInfo()
