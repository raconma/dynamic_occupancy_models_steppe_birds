###############################################################################
# fig_conceptual_framework.R
#
# Conceptual Figure: Demographic regimes of range dynamics
# Three-panel figure (A-C) with coherent patch narrative across panels.
# Bird silhouettes excluded — to be added manually in PowerPoint/Illustrator.
#
# Output: figs/pub_fig_conceptual_framework.pdf / .png
###############################################################################

library(ggplot2)
library(grid)
library(gridExtra)

# ---- Colour palette ----
col_occupied    <- "#A8B07A"
col_empty       <- "#E8E0D0"
col_extinct_bg  <- "#D2D5C2"
col_border      <- "#4A4A4A"
col_border_dash <- "#9A9A8A"
col_arrow       <- "#5B6B3C"
col_extinct_x   <- "#9E5A3C"
col_ground      <- "#F7F3EB"
col_text        <- "#333333"

# ---- Patch layout: 3x3, same positions in all panels ----
# Patches are numbered 1-9 (top-left to bottom-right)
#   1  2  3
#   4  5  6
#   7  8  9
patch_x <- c(1.8, 5.0, 8.2, 1.8, 5.0, 8.2, 1.8, 5.0, 8.2)
patch_y <- c(5.2, 5.2, 5.2, 3.1, 3.1, 3.1, 1.0, 1.0, 1.0)

set.seed(42)
patch_x <- patch_x + runif(9, -0.12, 0.12)
patch_y <- patch_y + runif(9, -0.08, 0.08)

patch_rx <- c(0.70, 0.75, 0.65, 0.72, 0.68, 0.74, 0.70, 0.66, 0.71)
patch_ry <- patch_rx * 0.60

patches <- data.frame(id = 1:9, x = patch_x, y = patch_y,
                      rx = patch_rx, ry = patch_ry)

# ---- Coherent patch narrative A -> B -> C ----
# Panel A: 7 occupied, 1 recently extinct (patch 6), 1 being recolonised (patch 2)
# Panel B: patches 6 and 3 lost (extinct), patch 2 recolonised successfully,
#          patch 9 now empty (long-empty), one arrow tries to reach patch 3
# Panel C: only patches 4,5,8 remain occupied. 1,6,3,9 empty. 2,7 recently extinct.
#          No arrows at all.

# The story:
# Patch 1: occupied -> occupied -> empty (long gone)
# Patch 2: empty(recolonising) -> occupied (success!) -> recently extinct
# Patch 3: occupied -> recently extinct -> empty (long gone)
# Patch 4: occupied -> occupied -> occupied (core)
# Patch 5: occupied -> occupied -> occupied (core)
# Patch 6: recently extinct -> empty -> empty
# Patch 7: occupied -> occupied -> recently extinct
# Patch 8: occupied -> occupied -> occupied (core, last refuge)
# Patch 9: occupied -> empty -> empty

states_A <- c("occupied",      "recolonising", "occupied",
              "occupied",      "occupied",     "extinct",
              "occupied",      "occupied",     "occupied")

states_B <- c("occupied",      "occupied",     "extinct",
              "occupied",      "occupied",     "empty",
              "occupied",      "occupied",     "empty")

states_C <- c("empty",         "extinct",      "empty",
              "occupied",      "occupied",     "empty",
              "extinct",       "occupied",     "empty")

# ---- Colonisation arrows ----
# Panel A: 3 arrows (active recolonisation system)
#   patch 1 -> patch 2 (recolonising)
#   patch 5 -> patch 6 (attempting, towards extinct patch — shows effort)
#   patch 8 -> patch 9 (maintaining connectivity)
arrows_A <- data.frame(
  x    = c(patches$x[1], patches$x[7], patches$x[4]),
  y    = c(patches$y[1], patches$y[7], patches$y[4]),
  xend = c(patches$x[2], patches$x[4], patches$x[1]),
  yend = c(patches$y[2], patches$y[4], patches$y[1])
)

# Panel B: 1 arrow (weakened recolonisation)
#   patch 5 -> patch 3 (trying to recolonise recently extinct patch — will fail)
arrows_B <- data.frame(
  x    = patches$x[5],
  y    = patches$y[5],
  xend = patches$x[3],
  yend = patches$y[3]
)

# Panel C: NO arrows (recolonisation has ceased)
arrows_C <- data.frame(x = numeric(0), y = numeric(0),
                        xend = numeric(0), yend = numeric(0))

# ---- Extinct positions ----
get_extinct_pos <- function(states) {
  idx <- which(states == "extinct")
  if (length(idx) == 0) return(data.frame(x = numeric(0), y = numeric(0)))
  data.frame(x = patches$x[idx], y = patches$y[idx])
}

# ---- Helpers ----
ellipse_poly <- function(cx, cy, rx, ry, n = 50) {
  t <- seq(0, 2 * pi, length.out = n)
  data.frame(x = cx + rx * cos(t), y = cy + ry * sin(t))
}

build_patches <- function(states) {
  dfs <- lapply(1:9, function(i) {
    ep <- ellipse_poly(patches$x[i], patches$y[i],
                       patches$rx[i], patches$ry[i])
    ep$id    <- i
    ep$state <- states[i]
    ep
  })
  do.call(rbind, dfs)
}

# ---- Panel builder ----
make_panel <- function(states, arrows_df, label, subtitle, annotation,
                       sub_annotation = NULL) {

  patch_data <- build_patches(states)
  extinct_df <- get_extinct_pos(states)

  fill_map   <- c(occupied = col_occupied, empty = col_empty,
                  extinct = col_extinct_bg, recolonising = col_empty)
  border_map <- c(occupied = col_border, empty = col_border_dash,
                  extinct = col_border_dash, recolonising = col_border_dash)
  lty_map    <- c(occupied = "solid", empty = "dashed",
                  extinct = "dashed", recolonising = "dashed")

  p <- ggplot() +
    annotate("rect", xmin = -0.3, xmax = 10.3, ymin = -0.3, ymax = 6.6,
             fill = col_ground, colour = NA) +
    annotate("segment", x = -0.3, xend = 10.3, y = -0.3, yend = -0.3,
             colour = "#D0C8B8", linewidth = 0.3)

  # Patches
  for (s in unique(states)) {
    sub_d <- patch_data[patch_data$state == s, ]
    if (nrow(sub_d) == 0) next
    p <- p + geom_polygon(
      data = sub_d,
      aes(x = x, y = y, group = id),
      fill = fill_map[s],
      colour = border_map[s],
      linewidth = 0.35,
      linetype = lty_map[s]
    )
  }

  # Arrows
  if (nrow(arrows_df) > 0) {
    p <- p + geom_curve(
      data = arrows_df,
      aes(x = x, y = y, xend = xend, yend = yend),
      colour = col_arrow,
      linewidth = 0.5,
      curvature = 0.25,
      arrow = arrow(length = unit(1.8, "mm"), type = "closed"),
      alpha = 0.80
    )
  }

  # X marks
  if (nrow(extinct_df) > 0) {
    xs <- 0.16
    for (i in 1:nrow(extinct_df)) {
      ex <- extinct_df$x[i]
      ey <- extinct_df$y[i]
      p <- p +
        annotate("segment",
                 x = ex - xs, xend = ex + xs,
                 y = ey - xs * 0.65, yend = ey + xs * 0.65,
                 colour = col_extinct_x, linewidth = 0.65) +
        annotate("segment",
                 x = ex + xs, xend = ex - xs,
                 y = ey - xs * 0.65, yend = ey + xs * 0.65,
                 colour = col_extinct_x, linewidth = 0.65)
    }
  }

  # Labels
  p <- p +
    annotate("text", x = 0.2, y = 7.8, label = label,
             fontface = "bold", size = 4.5, hjust = 0, colour = col_text) +
    annotate("text", x = 1.3, y = 7.8, label = subtitle,
             fontface = "plain", size = 3, hjust = 0, colour = col_text)

  # Greek annotation
  p <- p +
    annotate("text", x = 5, y = -0.9, label = annotation,
             size = 3.2, colour = col_text, parse = TRUE)

  # Sub-annotation
  if (!is.null(sub_annotation)) {
    p <- p +
      annotate("text", x = 5, y = -1.5, label = sub_annotation,
               fontface = "italic", size = 2.3, colour = "#777777")
  }

  p + coord_fixed(ratio = 1, xlim = c(-0.5, 10.5),
                  ylim = c(-1.8, 8.3), clip = "off") +
    theme_void() +
    theme(plot.margin = margin(1, 3, 3, 3, "mm"),
          plot.background = element_rect(fill = "white", colour = NA))
}

# ---- Three panels ----
panel_A <- make_panel(states_A, arrows_A,
                      "A", "Metapopulation equilibrium",
                      "gamma %~~% epsilon")

panel_B <- make_panel(states_B, arrows_B,
                      "B", "Declining but recoverable",
                      "epsilon > gamma")

panel_C <- make_panel(states_C, arrows_C,
                      "C", "Demographic trap",
                      "epsilon ~'>>>'~ gamma",
                      "Suitable habitat persists, but recolonisation has ceased")

# ---- Legend (with proper spacing) ----
legend_panel <- ggplot() +
  # 1. Occupied
  geom_polygon(data = ellipse_poly(1.0, 0.5, 0.25, 0.15),
               aes(x = x, y = y), fill = col_occupied,
               colour = col_border, linewidth = 0.3) +
  annotate("text", x = 1.5, y = 0.5, label = "Occupied patch",
           size = 2.3, hjust = 0, colour = col_text) +

  # 2. Empty suitable
  geom_polygon(data = ellipse_poly(3.5, 0.5, 0.25, 0.15),
               aes(x = x, y = y), fill = col_empty,
               colour = col_border_dash, linewidth = 0.3,
               linetype = "dashed") +
  annotate("text", x = 4.0, y = 0.5, label = "Empty suitable patch",
           size = 2.3, hjust = 0, colour = col_text) +

  # 3. Local extinction
  geom_polygon(data = ellipse_poly(6.6, 0.5, 0.25, 0.15),
               aes(x = x, y = y), fill = col_extinct_bg,
               colour = col_border_dash, linewidth = 0.3,
               linetype = "dashed") +
  annotate("segment", x = 6.48, xend = 6.72, y = 0.4, yend = 0.6,
           colour = col_extinct_x, linewidth = 0.5) +
  annotate("segment", x = 6.72, xend = 6.48, y = 0.4, yend = 0.6,
           colour = col_extinct_x, linewidth = 0.5) +
  annotate("text", x = 7.1, y = 0.5, label = "Recent extinction",
           size = 2.3, hjust = 0, colour = col_text) +

  # 4. Colonisation arrow
  annotate("segment", x = 9.0, xend = 9.5, y = 0.5, yend = 0.5,
           colour = col_arrow, linewidth = 0.5,
           arrow = arrow(length = unit(1.5, "mm"), type = "closed")) +
  annotate("text", x = 9.7, y = 0.5, label = "Colonisation",
           size = 2.3, hjust = 0, colour = col_text) +

  coord_cartesian(xlim = c(0.5, 11), ylim = c(0, 1)) +
  theme_void() +
  theme(plot.margin = margin(0, 5, 2, 5, "mm"))

# ---- Compose ----
fig <- grid.arrange(
  arrangeGrob(panel_A, panel_B, panel_C, ncol = 3),
  legend_panel,
  nrow = 2, heights = c(10, 1.2)
)

# ---- Save ----
out_dir <- file.path(getwd(), "figs")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ggsave(file.path(out_dir, "pub_fig_conceptual_framework.pdf"),
       fig, width = 190, height = 90, units = "mm", dpi = 300)

ggsave(file.path(out_dir, "pub_fig_conceptual_framework.png"),
       fig, width = 190, height = 90, units = "mm", dpi = 600)

cat("Saved: figs/pub_fig_conceptual_framework.pdf\n")
cat("Saved: figs/pub_fig_conceptual_framework.png\n")
