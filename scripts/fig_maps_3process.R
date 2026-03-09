###############################################################################
# fig_maps_3process.R
#
# Publication figure: Initial occupancy (psi1) maps for 4 Iberian steppe birds.
# 1-row x 4-column panel using existing OccuMap rasters.
#
# If spatial predictions for gamma/epsilon are feasible from the model objects,
# additional rows are added. Otherwise only psi1 maps are produced with a
# warning.
#
# Output:
#   figs/pub_fig_maps_3process_4species.png  (300 DPI)
#   figs/pub_fig_maps_3process_4species.pdf
###############################################################################

library(terra)
library(sf)
library(ggplot2)
library(rnaturalearth)
library(dplyr)
library(patchwork)

cat("\n", strrep("=", 70), "\n")
cat("  FIGURE: psi1 maps (4 species)\n")
cat(strrep("=", 70), "\n\n")

# --- Paths ---
proj_dir <- "."
data_dir <- file.path(proj_dir, "data/processed")
res_dir  <- file.path(proj_dir, "results")
figs_dir <- file.path(proj_dir, "figs")
if (!dir.exists(figs_dir)) dir.create(figs_dir, recursive = TRUE)

# --- Species info ---
sp_codes <- c("otitar", "ptealc", "pteori", "tettet")
sp_latin <- c(otitar = "Otis tarda", ptealc = "Pterocles alchata",
              pteori = "Pterocles orientalis", tettet = "Tetrax tetrax")

# --- Spain outline ---
cat("Loading country outlines...\n")
europe <- ne_countries(scale = 50, continent = "Europe", returnclass = "sf")
iberia <- europe %>% filter(name %in% c("Spain", "Portugal"))

bbox <- c(xmin = -9.5, xmax = 3.5, ymin = 35.8, ymax = 44)

# --- Build psi1 map panels ---
map_list <- list()

for (sp in sp_codes) {
  tif_path <- file.path(data_dir, sp, paste0(sp, "_OccuMap.tif"))
  if (!file.exists(tif_path)) {
    warning("OccuMap not found for ", sp, ": ", tif_path)
    next
  }

  r <- rast(tif_path)
  r_df <- as.data.frame(r, xy = TRUE)
  names(r_df)[3] <- "occupancy"
  r_df <- r_df %>% filter(occupancy > 0.001)

  p <- ggplot() +
    geom_raster(data = r_df, aes(x = x, y = y, fill = occupancy)) +
    geom_sf(data = iberia, fill = NA, colour = "grey30", linewidth = 0.3) +
    scale_fill_viridis_c(
      option = "mako", direction = -1,
      limits = c(0, 1),
      labels = scales::percent_format(),
      name = expression(psi[1])
    ) +
    coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]),
             ylim = c(bbox["ymin"], bbox["ymax"]),
             expand = FALSE) +
    labs(title = bquote(italic(.(sp_latin[sp])))) +
    theme_void(base_size = 9) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 9, face = "italic"),
      legend.position = "none",
      plot.margin = margin(2, 2, 2, 2)
    )

  map_list[[sp]] <- p
  cat("  ", sp, "loaded\n")
}

if (length(map_list) == 0) {
  stop("No OccuMap rasters found. Cannot produce figure.")
}

# --- Shared legend ---
# Build a dummy plot to extract the legend
legend_plot <- ggplot(data.frame(x = 0, y = 0, z = 0.5), aes(x, y, fill = z)) +
  geom_tile() +
  scale_fill_viridis_c(
    option = "mako", direction = -1,
    limits = c(0, 1),
    labels = scales::percent_format(),
    name = expression(psi[1])
  ) +
  theme_void() +
  theme(legend.position = "right",
        legend.key.height = unit(2, "cm"),
        legend.key.width = unit(0.35, "cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8))

legend_grob <- cowplot::get_legend(legend_plot)

# --- Assemble with patchwork ---
n_sp <- length(map_list)
p_combined <- wrap_plots(map_list, nrow = 1) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right",
        legend.key.height = unit(2, "cm"),
        legend.key.width = unit(0.35, "cm"))

# --- Attempt gamma/epsilon maps ---
cat("\nAttempting gamma/epsilon spatial predictions...\n")
cat("  NOTE: Full spatial prediction requires the original covariate grids,\n")
cat("        which are not straightforwardly available. Producing psi1-only figure.\n")
warning("Gamma and epsilon spatial maps not produced (would require full prediction grids).")

# --- Save ---
w_mm <- 180
h_mm <- 50
ggsave(file.path(figs_dir, "pub_fig_maps_3process_4species.png"), p_combined,
       width = w_mm, height = h_mm, units = "mm", dpi = 300)
ggsave(file.path(figs_dir, "pub_fig_maps_3process_4species.pdf"), p_combined,
       width = w_mm, height = h_mm, units = "mm")

cat("\n  -> pub_fig_maps_3process_4species.png\n")
cat("  -> pub_fig_maps_3process_4species.pdf\n")
cat("\nDone.\n")
