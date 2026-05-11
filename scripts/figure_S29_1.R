###############################################################################
# Figure S29.1 — two-panel figure for Appendix S29
#
# Panel A: Annual checklist-weighted mean of stepRep_strict_500m, by species,
#          for "all peninsular cells" and "focal cells (>=1 detection)".
# Panel B: Annual difference (focal - peninsular) per species — makes
#          explicit the divergent temporal trends.
###############################################################################

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})

WT  <- "/Users/gfandos/Documents/GitHub/dynamic_occupancy_models_steppe_birds/.claude/worktrees/loving-einstein-aa930a"
OUT <- "/Users/gfandos/Documents/GitHub/dynamic_occupancy_models_steppe_birds/results/stepRep_v2_run/figs"

sp_lookup <- c(otitar = "Otis tarda", ptealc = "Pterocles alchata",
               pteori = "Pterocles orientalis", tettet = "Tetrax tetrax")

trend <- rbindlist(lapply(names(sp_lookup), function(sp) {
  chk <- readRDS(file.path(WT, "data", "derived",
                            paste0("checklist_stepRep_", sp, ".rds")))
  setDT(chk)
  # focal cells = cells with >=1 detection of the species
  focal_cells <- chk[species_observed == TRUE, unique(cells)]

  # Per-year weighted mean (each checklist counts equally)
  pen <- chk[, .(weighted_mean = mean(step_strict_500m, na.rm = TRUE),
                  n_checklists  = .N), by = year]
  pen[, scope := "peninsular"]

  fcl <- chk[cells %in% focal_cells,
             .(weighted_mean = mean(step_strict_500m, na.rm = TRUE),
               n_checklists  = .N), by = year]
  fcl[, scope := "focal"]

  out <- rbind(pen, fcl)
  out[, species := sp][]
}))
trend[, species_name := sp_lookup[species]]

# Save the underlying data
fwrite(trend, file.path(OUT, "..", "fig_s29_1_trend_data.csv"))

# Panel A: dual line per species
pA <- ggplot(trend, aes(x = year, y = weighted_mean,
                          colour = scope, linetype = scope,
                          shape = scope, group = scope)) +
  geom_line(linewidth = 0.85) +
  geom_point(size = 2.4) +
  facet_wrap(~ species_name, ncol = 2) +
  scale_colour_manual(values = c(peninsular = "grey45",
                                  focal = "#c0392b"),
                      breaks = c("peninsular", "focal"),
                      labels = c(peninsular = "All peninsular cells",
                                  focal = "Cells with >=1 detection"),
                      name = NULL) +
  scale_linetype_manual(values = c(peninsular = "dashed",
                                    focal = "solid"),
                        breaks = c("peninsular", "focal"),
                        labels = c(peninsular = "All peninsular cells",
                                    focal = "Cells with >=1 detection"),
                        name = NULL) +
  scale_shape_manual(values = c(peninsular = 1, focal = 16),
                     breaks = c("peninsular", "focal"),
                     labels = c(peninsular = "All peninsular cells",
                                 focal = "Cells with >=1 detection"),
                     name = NULL) +
  scale_x_continuous(breaks = 2017:2023) +
  ylim(0.1, 0.7) +
  labs(x = NULL, y = "stepRep_strict_500m\n(annual checklist-weighted mean)",
       tag = "A") +
  theme_bw() +
  theme(legend.position = "top",
        strip.background = element_rect(fill = "grey95"),
        plot.tag = element_text(face = "bold", size = 14))

# Panel B: delta focal - peninsular per year, faceted by species
delta_wide <- dcast(trend[, .(year, scope, species_name, weighted_mean)],
                    species_name + year ~ scope, value.var = "weighted_mean")
delta_wide[, delta := focal - peninsular]

pB <- ggplot(delta_wide, aes(x = year, y = delta, fill = species_name)) +
  geom_col(width = 0.7, colour = "black", linewidth = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  facet_wrap(~ species_name, ncol = 2) +
  scale_x_continuous(breaks = 2017:2023) +
  scale_fill_manual(values = c("Otis tarda" = "#e08214",
                                 "Pterocles alchata" = "#5e3c99",
                                 "Pterocles orientalis" = "#1b7837",
                                 "Tetrax tetrax" = "#4393c3"),
                    guide = "none") +
  labs(x = "Year",
       y = "Delta stepRep (focal cells - peninsular)",
       tag = "B") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "grey95"),
        plot.tag = element_text(face = "bold", size = 14))

combined <- pA / pB + plot_layout(heights = c(1, 0.8))
ggsave(file.path(OUT, "figure_S29_1.png"), combined,
       width = 11, height = 11, dpi = 300)
ggsave(file.path(OUT, "figure_S29_1.pdf"), combined,
       width = 11, height = 11, device = cairo_pdf)

cat("Figure S29.1 written. Delta range per species:\n")
print(delta_wide[, .(min_delta = min(delta), max_delta = max(delta),
                      delta_2017 = delta[year == 2017],
                      delta_2023 = delta[year == 2023]),
                  by = species_name])
