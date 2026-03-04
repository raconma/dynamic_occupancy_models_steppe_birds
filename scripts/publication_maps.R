###############################################################################
# publication_maps.R
#
# Genera mapas espaciales para publicacion:
#   1. Panel 4 especies: ocupacion predicha (OccuMap rasters)
#   2. Mapa de sitios muestreados con detecciones vs no-detecciones
#   3. Mapa de tree_cover + ocupacion superpuesta
#   4. Mapa de cambio temporal (diferencia smoothed year1 vs year6)
#
# Requiere:
#   data/processed/{sp}/{sp}_OccuMap.tif
#   data/processed/{sp}/{sp}_occ_wide_latlong.csv
#   results/results2/{sp}_{pipeline}_all_models.rds
#   results/results2/{sp}_{pipeline}_umf.rds
#
# Output:
#   figs/pub_map_occupancy_4species.png
#   figs/pub_map_sampling_sites.png
#   figs/pub_map_occupancy_individual.png (4 panels alta resolucion)
#   figs/pub_map_detection_sites.png
###############################################################################

library(terra)
library(sf)
library(ggplot2)
library(rnaturalearth)
library(dplyr)
library(tidyr)
library(gridExtra)

cat("\n")
cat(strrep("=", 70), "\n")
cat("  MAPAS PARA PUBLICACION\n")
cat("  4 especies esteparias - Peninsula Iberica\n")
cat(strrep("=", 70), "\n\n")

# --- Paths ---
data_dir <- "data/processed"
rds_dir  <- "results/results2"
figs_dir <- "figs"
if (!dir.exists(figs_dir)) dir.create(figs_dir, recursive = TRUE)

# --- Species info ---
sp_codes <- c("otitar", "ptealc", "pteori", "tettet")
sp_latin <- c(otitar = "Otis tarda", ptealc = "Pterocles alchata",
              pteori = "Pterocles orientalis", tettet = "Tetrax tetrax")
sp_common <- c(otitar = "Great Bustard", ptealc = "Pin-tailed Sandgrouse",
               pteori = "Black-bellied Sandgrouse", tettet = "Little Bustard")

# --- Get Spain + Portugal outline ---
cat("Cargando contornos...\n")
europe <- ne_countries(scale = 50, continent = "Europe", returnclass = "sf")
iberia <- europe %>% filter(name %in% c("Spain", "Portugal"))
spain  <- europe %>% filter(name == "Spain")

# Study area bounding box (from rasters)
bbox <- c(xmin = -9.5, xmax = 3.5, ymin = 35.8, ymax = 44)

# --- Publication theme for maps ---
theme_map <- theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(face = "bold.italic", size = 11),
    plot.subtitle = element_text(size = 8, color = "grey40"),
    legend.position = "bottom",
    legend.key.width = unit(1.5, "cm"),
    legend.key.height = unit(0.3, "cm")
  )

###############################################################################
# MAPA 1: Panel 4 especies - Ocupacion predicha (OccuMap rasters)
###############################################################################

cat("\n  Mapa 1: Panel ocupacion predicha (4 especies)...\n")

map_list <- list()

for (sp in sp_codes) {
  # Load raster
  r <- rast(file.path(data_dir, sp, paste0(sp, "_OccuMap.tif")))

  # Convert to data frame for ggplot
  r_df <- as.data.frame(r, xy = TRUE)
  names(r_df)[3] <- "occupancy"

  # Remove zero/very low values for cleaner map
  r_df <- r_df %>% filter(occupancy > 0.001)

  p <- ggplot() +
    geom_raster(data = r_df, aes(x = x, y = y, fill = occupancy)) +
    geom_sf(data = iberia, fill = NA, color = "grey30", linewidth = 0.4) +
    scale_fill_viridis_c(
      option = "magma", direction = -1,
      limits = c(0, max(r_df$occupancy)),
      labels = scales::percent_format(),
      name = expression(hat(psi))
    ) +
    coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]),
             ylim = c(bbox["ymin"], bbox["ymax"]),
             expand = FALSE) +
    labs(title = sp_latin[sp],
         subtitle = sp_common[sp]) +
    theme_map

  map_list[[sp]] <- p
  cat("    ", sp, "loaded\n")
}

# Combine into panel
p_panel <- grid.arrange(
  map_list$otitar   + labs(title = paste0("A) ", sp_latin["otitar"])),
  map_list$ptealc   + labs(title = paste0("B) ", sp_latin["ptealc"])),
  map_list$pteori   + labs(title = paste0("C) ", sp_latin["pteori"])),
  map_list$tettet   + labs(title = paste0("D) ", sp_latin["tettet"])),
  ncol = 2
)

ggsave(file.path(figs_dir, "pub_map_occupancy_4species.png"), p_panel,
       width = 12, height = 11, dpi = 300)
cat("    -> pub_map_occupancy_4species.png\n")

###############################################################################
# MAPA 2: Individual species maps (alta resolucion, 1 por especie)
###############################################################################

cat("\n  Mapa 2: Mapas individuales por especie...\n")

for (sp in sp_codes) {
  r <- rast(file.path(data_dir, sp, paste0(sp, "_OccuMap.tif")))
  r_df <- as.data.frame(r, xy = TRUE)
  names(r_df)[3] <- "occupancy"
  r_df <- r_df %>% filter(occupancy > 0.001)

  p_ind <- ggplot() +
    geom_raster(data = r_df, aes(x = x, y = y, fill = occupancy)) +
    geom_sf(data = iberia, fill = NA, color = "grey20", linewidth = 0.5) +
    scale_fill_viridis_c(
      option = "magma", direction = -1,
      limits = c(0, max(r_df$occupancy)),
      labels = scales::percent_format(),
      name = "Predicted\noccupancy"
    ) +
    coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]),
             ylim = c(bbox["ymin"], bbox["ymax"]),
             expand = FALSE) +
    labs(title = paste0(sp_latin[sp], " (", sp_common[sp], ")"),
         subtitle = "Initial occupancy probability from colext model") +
    theme_map +
    theme(legend.position = "right",
          legend.key.width = unit(0.4, "cm"),
          legend.key.height = unit(2, "cm"))

  fname <- paste0("pub_map_occupancy_", sp, ".png")
  ggsave(file.path(figs_dir, fname), p_ind, width = 9, height = 7, dpi = 300)
  cat("    ->", fname, "\n")
}

###############################################################################
# MAPA 3: Sitios de muestreo con detecciones
###############################################################################

cat("\n  Mapa 3: Sitios de muestreo...\n")

library(unmarked)

site_maps <- list()

for (sp in sp_codes) {
  # Load latlong
  ll <- read.csv(file.path(data_dir, sp, paste0(sp, "_occ_wide_latlong.csv")))
  ll_sf <- st_as_sf(ll, coords = c("longitude", "latitude"), crs = 4326)

  # Load UMF to get detection data
  info <- list(
    otitar = list(pipe = "v4b"),
    ptealc = list(pipe = "v5"),
    pteori = list(pipe = "v4b"),
    tettet = list(pipe = "v4b")
  )
  umf <- readRDS(file.path(rds_dir, paste0(sp, "_", info[[sp]]$pipe, "_umf.rds")))
  y <- getY(umf)

  # Classify sites: ever detected vs never detected
  # Note: UMF sites (1988/2179) are a subset of latlong cells (3056/3400)
  # We can only classify sites that are in the UMF
  # For latlong, just show all survey cells
  ever_detected <- apply(y, 1, function(x) any(x == 1, na.rm = TRUE))
  n_detected <- sum(ever_detected)
  n_total <- nrow(y)

  p_sites <- ggplot() +
    geom_sf(data = iberia, fill = "grey95", color = "grey40", linewidth = 0.4) +
    geom_sf(data = ll_sf, color = "grey70", size = 0.3, alpha = 0.5) +
    coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]),
             ylim = c(bbox["ymin"], bbox["ymax"]),
             expand = FALSE) +
    labs(title = sp_latin[sp],
         subtitle = paste0(nrow(ll), " survey cells; ",
                           n_detected, "/", n_total, " sites with detections (",
                           round(n_detected/n_total*100, 1), "%)")) +
    theme_map

  site_maps[[sp]] <- p_sites
}

p_sites_panel <- grid.arrange(
  site_maps$otitar + labs(title = paste0("A) ", sp_latin["otitar"])),
  site_maps$ptealc + labs(title = paste0("B) ", sp_latin["ptealc"])),
  site_maps$pteori + labs(title = paste0("C) ", sp_latin["pteori"])),
  site_maps$tettet + labs(title = paste0("D) ", sp_latin["tettet"])),
  ncol = 2
)

ggsave(file.path(figs_dir, "pub_map_sampling_sites.png"), p_sites_panel,
       width = 12, height = 11, dpi = 300)
cat("    -> pub_map_sampling_sites.png\n")

###############################################################################
# MAPA 4: Ocupacion + Sitios con deteccion superpuestos
###############################################################################

cat("\n  Mapa 4: Ocupacion + sitios de deteccion...\n")

overlay_maps <- list()

for (sp in sp_codes) {
  # Raster
  r <- rast(file.path(data_dir, sp, paste0(sp, "_OccuMap.tif")))
  r_df <- as.data.frame(r, xy = TRUE)
  names(r_df)[3] <- "occupancy"
  r_df <- r_df %>% filter(occupancy > 0.001)

  # Latlong points
  ll <- read.csv(file.path(data_dir, sp, paste0(sp, "_occ_wide_latlong.csv")))
  ll_sf <- st_as_sf(ll, coords = c("longitude", "latitude"), crs = 4326)

  p_overlay <- ggplot() +
    geom_raster(data = r_df, aes(x = x, y = y, fill = occupancy)) +
    geom_sf(data = iberia, fill = NA, color = "grey20", linewidth = 0.5) +
    geom_sf(data = ll_sf, color = "white", size = 0.15, alpha = 0.6,
            shape = 16) +
    scale_fill_viridis_c(
      option = "magma", direction = -1,
      limits = c(0, max(r_df$occupancy)),
      labels = scales::percent_format(),
      name = expression(hat(psi))
    ) +
    coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]),
             ylim = c(bbox["ymin"], bbox["ymax"]),
             expand = FALSE) +
    labs(title = sp_latin[sp],
         subtitle = paste0("Predicted occupancy with survey locations (white dots)")) +
    theme_map

  overlay_maps[[sp]] <- p_overlay
}

p_overlay_panel <- grid.arrange(
  overlay_maps$otitar + labs(title = paste0("A) ", sp_latin["otitar"])),
  overlay_maps$ptealc + labs(title = paste0("B) ", sp_latin["ptealc"])),
  overlay_maps$pteori + labs(title = paste0("C) ", sp_latin["pteori"])),
  overlay_maps$tettet + labs(title = paste0("D) ", sp_latin["tettet"])),
  ncol = 2
)

ggsave(file.path(figs_dir, "pub_map_occupancy_sites_overlay.png"), p_overlay_panel,
       width = 12, height = 11, dpi = 300)
cat("    -> pub_map_occupancy_sites_overlay.png\n")

###############################################################################
# MAPA 5: Hotspots - Areas de alta ocupacion compartidas entre especies
###############################################################################

cat("\n  Mapa 5: Hotspots multi-especie...\n")

# Load all rasters and stack
rasters <- list()
for (sp in sp_codes) {
  r <- rast(file.path(data_dir, sp, paste0(sp, "_OccuMap.tif")))
  rasters[[sp]] <- r
}

# Create species richness map (sum of occupancy > threshold)
threshold <- 0.05  # 5% occupancy threshold
binary_stack <- rast(lapply(rasters, function(r) r > threshold))
richness <- sum(binary_stack, na.rm = TRUE)

rich_df <- as.data.frame(richness, xy = TRUE)
names(rich_df)[3] <- "n_species"
rich_df <- rich_df %>% filter(n_species > 0)

p_hotspot <- ggplot() +
  geom_raster(data = rich_df, aes(x = x, y = y, fill = factor(n_species))) +
  geom_sf(data = iberia, fill = NA, color = "grey20", linewidth = 0.5) +
  scale_fill_manual(
    values = c("1" = "#FEE08B", "2" = "#FDAE61", "3" = "#F46D43", "4" = "#D73027"),
    name = "Number of\nspecies",
    labels = c("1", "2", "3", "4")
  ) +
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]),
           ylim = c(bbox["ymin"], bbox["ymax"]),
           expand = FALSE) +
  labs(title = "Multi-species steppe bird hotspots",
       subtitle = paste0("Number of species with predicted occupancy > ",
                         threshold * 100, "%")) +
  theme_map +
  theme(legend.position = "right")

ggsave(file.path(figs_dir, "pub_map_hotspots_multispecies.png"), p_hotspot,
       width = 9, height = 7, dpi = 300)
cat("    -> pub_map_hotspots_multispecies.png\n")

###############################################################################
# MAPA 6: Panel combinado para publicacion principal
#         (Ocupacion + hotspot en un solo panel)
###############################################################################

cat("\n  Mapa 6: Panel combinado publicacion...\n")

# Recreate individual maps with consistent color scales and smaller titles
pub_maps <- list()
max_occ <- 0.47  # max across all species for consistent scale

for (i in seq_along(sp_codes)) {
  sp <- sp_codes[i]
  r <- rast(file.path(data_dir, sp, paste0(sp, "_OccuMap.tif")))
  r_df <- as.data.frame(r, xy = TRUE)
  names(r_df)[3] <- "occupancy"
  r_df <- r_df %>% filter(occupancy > 0.001)

  panel_letter <- LETTERS[i]

  pub_maps[[sp]] <- ggplot() +
    geom_raster(data = r_df, aes(x = x, y = y, fill = occupancy)) +
    geom_sf(data = iberia, fill = NA, color = "grey30", linewidth = 0.3) +
    scale_fill_viridis_c(
      option = "magma", direction = -1,
      limits = c(0, max_occ),
      labels = scales::percent_format(),
      name = expression(hat(psi))
    ) +
    coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]),
             ylim = c(bbox["ymin"], bbox["ymax"]),
             expand = FALSE) +
    labs(title = paste0(panel_letter, ") ", sp_latin[sp])) +
    theme_map +
    theme(plot.title = element_text(size = 10, face = "bold.italic"),
          legend.key.width = unit(1, "cm"),
          legend.key.height = unit(0.25, "cm"),
          plot.margin = margin(2, 2, 2, 2))
}

# Add hotspot as panel E
p_hotspot_e <- p_hotspot +
  labs(title = "E) Multi-species hotspots") +
  theme(plot.title = element_text(size = 10, face = "bold"),
        plot.margin = margin(2, 2, 2, 2))

p_pub_combined <- grid.arrange(
  pub_maps$otitar, pub_maps$ptealc,
  pub_maps$pteori, pub_maps$tettet,
  p_hotspot_e, grid::nullGrob(),
  ncol = 2, heights = c(1, 1, 1)
)

ggsave(file.path(figs_dir, "pub_map_main_figure.png"), p_pub_combined,
       width = 12, height = 16, dpi = 300)
cat("    -> pub_map_main_figure.png\n")

###############################################################################
# RESUMEN
###############################################################################

cat("\n\n")
cat(strrep("=", 70), "\n")
cat("  RESUMEN DE MAPAS GENERADOS\n")
cat(strrep("=", 70), "\n\n")

map_files <- list.files(figs_dir, pattern = "pub_map.*\\.png$")
cat("  Total mapas:", length(map_files), "\n\n")
for (f in sort(map_files)) {
  finfo <- file.info(file.path(figs_dir, f))
  cat("   ", f, sprintf("(%.1f MB)\n", finfo$size / 1024^2))
}

cat("\n")
cat(strrep("=", 70), "\n")
cat("  FIN\n")
cat(strrep("=", 70), "\n\n")
