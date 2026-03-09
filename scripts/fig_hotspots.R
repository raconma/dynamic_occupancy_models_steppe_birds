###############################################################################
# fig_hotspots.R
#
# Publication Figure 1: Hotspot maps for 4 Iberian steppe birds.
# 4-column (species) x 3-row (psi1, gamma, epsilon) panel.
#
# Binary hotspot maps: top 25% of values = hotspot (coloured), rest = grey.
# Row colours: psi1 = navy, gamma = orange, epsilon = red.
#
# Data sources:
#   psi1:    data/processed/{sp}/{sp}_OccuMap.tif
#   gamma:   results/stpgocc_col_hotspots_{sp}.tif
#   epsilon: results/stpgocc_ext_hotspots_{sp}.tif
#
# Output:
#   figs/pub_fig1_hotspots_4species.png  (300 DPI, ~180 mm wide)
#   figs/pub_fig1_hotspots_4species.pdf
###############################################################################

library(terra)
library(sf)
library(ggplot2)
library(rnaturalearth)
library(dplyr)
library(patchwork)

cat("\n", strrep("=", 70), "\n")
cat("  FIGURE 1: Hotspot maps (4 species x 3 processes)\n")
cat(strrep("=", 70), "\n\n")

# --- Paths ---
proj_dir <- "."
data_dir <- file.path(proj_dir, "data/processed")
res_dir  <- file.path(proj_dir, "results")
figs_dir <- file.path(proj_dir, "figs")
if (!dir.exists(figs_dir)) dir.create(figs_dir, recursive = TRUE)

# --- Species info ---
sp_codes <- c("otitar", "ptealc", "pteori", "tettet")
sp_italic <- c(
  otitar = "O. tarda",
  ptealc = "P. alchata",
  pteori = "P. orientalis",
  tettet = "T. tetrax"
)

# --- Row definitions ---
row_processes <- c("psi1", "gamma", "epsilon")
row_colours <- c(psi1 = "#1B4F8A", gamma = "#E07B39", epsilon = "#B03A2E")
row_labels  <- c(
  psi1    = expression(psi[1]~"hotspots"),
  gamma   = expression(gamma~"hotspots"),
  epsilon = expression(epsilon~"hotspots")
)

# --- Spain outline ---
cat("Loading country outlines...\n")
europe <- ne_countries(scale = 50, continent = "Europe", returnclass = "sf")
iberia <- europe %>% filter(name %in% c("Spain", "Portugal"))

bbox <- c(xmin = -9.5, xmax = 3.5, ymin = 35.8, ymax = 44)

# --- Helper: build one hotspot panel ---
make_hotspot_panel <- function(raster_path, hotspot_colour, sp_code, process,
                               show_col_title = FALSE, show_row_label = FALSE) {

  if (!file.exists(raster_path)) {
    warning("Raster not found: ", raster_path)
    return(ggplot() + theme_void())
  }

  r <- rast(raster_path)
  r_df <- as.data.frame(r, xy = TRUE)
  names(r_df)[3] <- "value"

  # Remove zero / near-zero cells (no data)
  r_df <- r_df %>% filter(value > 1e-10)

  if (nrow(r_df) == 0) {
    warning("No data in raster: ", raster_path)
    return(ggplot() + theme_void())
  }

  # Threshold: top 25% = hotspot
  q75 <- quantile(r_df$value, 0.75, na.rm = TRUE)
  r_df$hotspot <- ifelse(r_df$value >= q75, "hot", "not")

  p <- ggplot() +
    # Non-hotspot cells (grey)
    geom_raster(data = r_df %>% filter(hotspot == "not"),
                aes(x = x, y = y), fill = "#EEEEEE") +
    # Hotspot cells (coloured)
    geom_raster(data = r_df %>% filter(hotspot == "hot"),
                aes(x = x, y = y), fill = hotspot_colour) +
    # Country outline
    geom_sf(data = iberia, fill = NA, colour = "grey40", linewidth = 0.25) +
    coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]),
             ylim = c(bbox["ymin"], bbox["ymax"]),
             expand = FALSE) +
    theme_void(base_size = 8) +
    theme(
      plot.margin = margin(1, 1, 1, 1),
      plot.title = element_text(hjust = 0.5, size = 8, face = "italic",
                                margin = margin(b = 1))
    )

  # Column title (species name) only on top row
  if (show_col_title) {
    p <- p + labs(title = bquote(italic(.(sp_italic[sp_code]))))
  } else {
    p <- p + labs(title = NULL)
  }

  return(p)
}

# --- Build all 12 panels ---
cat("Building panels...\n")
panel_list <- list()

for (i_proc in seq_along(c("psi1", "gamma", "epsilon"))) {
  process <- c("psi1", "gamma", "epsilon")[i_proc]
  colour  <- row_colours[process]

  for (j_sp in seq_along(sp_codes)) {
    sp <- sp_codes[j_sp]

    # Determine raster path
    if (process == "psi1") {
      rpath <- file.path(data_dir, sp, paste0(sp, "_OccuMap.tif"))
    } else if (process == "gamma") {
      rpath <- file.path(res_dir, paste0("stpgocc_col_hotspots_", sp, ".tif"))
    } else {
      rpath <- file.path(res_dir, paste0("stpgocc_ext_hotspots_", sp, ".tif"))
    }

    show_title <- (i_proc == 1)  # column headers only on first row

    panel <- make_hotspot_panel(rpath, colour, sp, process,
                                show_col_title = show_title)

    panel_list[[paste0(process, "_", sp)]] <- panel
    cat("  ", process, "-", sp, "\n")
  }
}

# --- Row labels as text panels ---
make_row_label <- function(process) {
  label_expr <- switch(process,
    psi1    = 'psi[1]~"hotspots"',
    gamma   = 'gamma~"hotspots"',
    epsilon = 'epsilon~"hotspots"'
  )
  colour <- row_colours[process]

  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = label_expr,
             colour = colour, size = 3, fontface = "bold", angle = 90,
             parse = TRUE) +
    xlim(0, 1) + ylim(0, 1) +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
}

# --- Assemble with patchwork: 3 rows x (label + 4 maps) ---
# We'll use a design layout

row1 <- make_row_label("psi1") +
  panel_list[["psi1_otitar"]] + panel_list[["psi1_ptealc"]] +
  panel_list[["psi1_pteori"]] + panel_list[["psi1_tettet"]] +
  plot_layout(widths = c(0.06, 1, 1, 1, 1), nrow = 1)

row2 <- make_row_label("gamma") +
  panel_list[["gamma_otitar"]] + panel_list[["gamma_ptealc"]] +
  panel_list[["gamma_pteori"]] + panel_list[["gamma_tettet"]] +
  plot_layout(widths = c(0.06, 1, 1, 1, 1), nrow = 1)

row3 <- make_row_label("epsilon") +
  panel_list[["epsilon_otitar"]] + panel_list[["epsilon_ptealc"]] +
  panel_list[["epsilon_pteori"]] + panel_list[["epsilon_tettet"]] +
  plot_layout(widths = c(0.06, 1, 1, 1, 1), nrow = 1)

combined <- row1 / row2 / row3

# --- Save ---
cat("\nSaving figure...\n")

w_mm <- 180
h_mm <- 135  # 3 rows

ggsave(file.path(figs_dir, "pub_fig1_hotspots_4species.png"), combined,
       width = w_mm, height = h_mm, units = "mm", dpi = 300, bg = "white")
cat("  Saved: figs/pub_fig1_hotspots_4species.png\n")

ggsave(file.path(figs_dir, "pub_fig1_hotspots_4species.pdf"), combined,
       width = w_mm, height = h_mm, units = "mm", bg = "white")
cat("  Saved: figs/pub_fig1_hotspots_4species.pdf\n")

cat("\nDone.\n")
