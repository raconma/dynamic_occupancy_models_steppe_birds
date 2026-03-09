## ------------------------------------------------------------------
## fig_bivariate_risk_map.R
## 4-panel bivariate risk map: log10(epsilon / gamma) across Iberia
## + 4-panel extinction-only map
## ------------------------------------------------------------------

library(terra)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(patchwork)
library(dplyr)

# --- Species codes and labels ---
spp   <- c("otitar", "ptealc", "pteori", "tettet")
labels <- c(
  otitar = "O. tarda",
  ptealc = "P. alchata",
  pteori = "P. orientalis",
  tettet = "T. tetrax"
)

sp_colours <- c(otitar = "#D55E00", ptealc = "#0072B2",
                pteori = "#009E73", tettet = "#CC79A7")

# --- Country outline (WGS84) ---
iberia <- ne_countries(scale = 50,
                       country = c("Spain", "Portugal"),
                       returnclass = "sf")

# --- Build combined data.frame ---
all_df <- list()
ext_df <- list()

for (sp in spp) {
  col_r <- rast(sprintf("results/stpgocc_col_hotspots_%s.tif", sp))
  ext_r <- rast(sprintf("results/stpgocc_ext_hotspots_%s.tif", sp))

  # --- Extinction-only data ---
  edf <- as.data.frame(ext_r, xy = TRUE)
  names(edf)[3] <- "epsilon"
  edf <- edf[!is.na(edf$epsilon), ]
  edf$species <- labels[sp]
  ext_df[[sp]] <- edf

  # --- Bivariate ratio ---
  ratio_r <- ext_r / col_r
  df <- as.data.frame(ratio_r, xy = TRUE)
  names(df)[3] <- "ratio"
  df <- df[!is.na(df$ratio), ]

  df$log_ratio <- log10(df$ratio)

  # Cap where gamma was near-zero
  col_df <- as.data.frame(col_r, xy = TRUE)
  names(col_df)[3] <- "gamma"
  df <- left_join(df, col_df, by = c("x", "y"))

  finite_vals <- df$log_ratio[is.finite(df$log_ratio) & df$gamma >= 1e-6]
  cap_val <- quantile(finite_vals, 0.99, na.rm = TRUE)
  low_cap <- quantile(finite_vals, 0.01, na.rm = TRUE)

  df$log_ratio[df$gamma < 1e-6 | !is.finite(df$log_ratio)] <- cap_val
  df$log_ratio[df$log_ratio > cap_val] <- cap_val
  df$log_ratio[df$log_ratio < low_cap] <- low_cap

  df$species <- labels[sp]
  all_df[[sp]] <- df
}

plot_df <- bind_rows(all_df)
plot_df$species <- factor(plot_df$species, levels = labels)

ext_plot_df <- bind_rows(ext_df)
ext_plot_df$species <- factor(ext_plot_df$species, levels = labels)

# =====================================================================
# FIGURE A: Bivariate risk map log10(epsilon/gamma) -- viridis palette
# =====================================================================

# --- Build plot with viridis ---
p_biv <- ggplot(plot_df) +
  geom_raster(aes(x = x, y = y, fill = log_ratio)) +
  geom_sf(data = iberia, fill = NA, colour = "grey20", linewidth = 0.35,
          inherit.aes = FALSE) +
  scale_fill_viridis_c(
    option   = "inferno",
    direction = 1,
    name     = expression(log[10](epsilon / gamma)),
    guide    = guide_colorbar(
      title.position = "top",
      barheight = unit(50, "mm"),
      barwidth  = unit(4, "mm")
    )
  ) +
  facet_wrap(~ species, ncol = 2) +
  coord_sf(xlim = c(-10, 5), ylim = c(35, 44), crs = 4326, expand = FALSE) +
  labs(caption = "Low values = growth (\u03b5 \u2248 \u03b3)      High values = decline (\u03b5 \u226b \u03b3)") +
  theme_void(base_size = 10) +
  theme(
    strip.text   = element_text(face = "italic", size = 10, margin = margin(b = 2)),
    legend.position = "right",
    legend.title = element_text(size = 9),
    legend.text  = element_text(size = 8),
    plot.caption = element_text(hjust = 0.5, size = 8, colour = "grey40",
                                margin = margin(t = 5)),
    plot.margin  = margin(5, 5, 5, 5)
  )

ggsave("figs/pub_fig3_bivariate_risk_map.png", p_biv,
       width = 180, height = 160, units = "mm", dpi = 300)
ggsave("figs/pub_fig3_bivariate_risk_map.pdf", p_biv,
       width = 180, height = 160, units = "mm")
message("Done: bivariate risk map")

# =====================================================================
# FIGURE B: Extinction-only map
# =====================================================================

p_ext <- ggplot(ext_plot_df) +
  geom_raster(aes(x = x, y = y, fill = epsilon)) +
  geom_sf(data = iberia, fill = NA, colour = "grey20", linewidth = 0.35,
          inherit.aes = FALSE) +
  scale_fill_viridis_c(
    option    = "rocket",
    direction = -1,
    name      = expression(bar(epsilon)),
    guide     = guide_colorbar(
      title.position = "top",
      barheight = unit(50, "mm"),
      barwidth  = unit(4, "mm")
    )
  ) +
  facet_wrap(~ species, ncol = 2) +
  coord_sf(xlim = c(-10, 5), ylim = c(35, 44), crs = 4326, expand = FALSE) +
  labs(caption = "Mean predicted local extinction rate (2017-2023)") +
  theme_void(base_size = 10) +
  theme(
    strip.text   = element_text(face = "italic", size = 10, margin = margin(b = 2)),
    legend.position = "right",
    legend.title = element_text(size = 9),
    legend.text  = element_text(size = 8),
    plot.caption = element_text(hjust = 0.5, size = 8, colour = "grey40",
                                margin = margin(t = 5)),
    plot.margin  = margin(5, 5, 5, 5)
  )

ggsave("figs/pub_fig3b_extinction_only_map.png", p_ext,
       width = 180, height = 160, units = "mm", dpi = 300)
ggsave("figs/pub_fig3b_extinction_only_map.pdf", p_ext,
       width = 180, height = 160, units = "mm")
message("Done: extinction-only map")
