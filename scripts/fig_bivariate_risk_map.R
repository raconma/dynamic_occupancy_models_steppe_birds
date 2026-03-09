## ------------------------------------------------------------------
## fig_bivariate_risk_map.R
## 4-panel bivariate risk map: log10(epsilon / gamma) across Iberia
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

# --- Country outline (WGS84) ---
iberia <- ne_countries(scale = 50,
                       country = c("Spain", "Portugal"),
                       returnclass = "sf")

# --- Build combined data.frame ---
all_df <- list()

for (sp in spp) {
  col_r <- rast(sprintf("results/stpgocc_col_hotspots_%s.tif", sp))
  ext_r <- rast(sprintf("results/stpgocc_ext_hotspots_%s.tif", sp))

  # ratio = epsilon / gamma

  ratio_r <- ext_r / col_r

  # Convert to data.frame

  df <- as.data.frame(ratio_r, xy = TRUE)
  names(df)[3] <- "ratio"
  df <- df[!is.na(df$ratio), ]

  # log10 transform

  df$log_ratio <- log10(df$ratio)

  # Cap where gamma was near-zero (ratio explodes)
  # Identify cells with gamma < 1e-6
  col_df <- as.data.frame(col_r, xy = TRUE)
  names(col_df)[3] <- "gamma"
  df <- left_join(df, col_df, by = c("x", "y"))

  # Compute 99th percentile of non-capped distribution
  finite_vals <- df$log_ratio[is.finite(df$log_ratio) & df$gamma >= 1e-6]
  cap_val <- quantile(finite_vals, 0.99, na.rm = TRUE)

  # Cap: where gamma < 1e-6 or log_ratio > cap_val or is infinite

  df$log_ratio[df$gamma < 1e-6 | !is.finite(df$log_ratio)] <- cap_val
  df$log_ratio[df$log_ratio > cap_val] <- cap_val

  # Also cap the lower tail symmetrically

  low_cap <- quantile(finite_vals, 0.01, na.rm = TRUE)
  df$log_ratio[df$log_ratio < low_cap] <- low_cap

  df$species <- labels[sp]
  all_df[[sp]] <- df
}

plot_df <- bind_rows(all_df)
plot_df$species <- factor(plot_df$species, levels = labels)

# --- Symmetric colour limits ---
max_abs <- max(abs(plot_df$log_ratio), na.rm = TRUE)

# --- Build plot ---
p <- ggplot(plot_df) +
  geom_raster(aes(x = x, y = y, fill = log_ratio)) +
  geom_sf(data = iberia, fill = NA, colour = "grey30", linewidth = 0.3,
          inherit.aes = FALSE) +
  scale_fill_gradient2(
    low      = "#2166AC",
    mid      = "grey95",
    high     = "#B2182B",
    midpoint = 0,
    limits   = c(-max_abs, max_abs),
    name     = expression(log[10](epsilon / gamma))
  ) +
  facet_wrap(~ species, ncol = 2) +
  coord_sf(
    xlim = c(-10, 5),
    ylim = c(35, 44),
    crs  = 4326,
    expand = FALSE
  ) +
  theme_void(base_size = 10) +
  theme(
    strip.text = element_text(face = "italic", size = 10, margin = margin(b = 2)),
    legend.position = "right",
    legend.title = element_text(size = 9),
    legend.text  = element_text(size = 8),
    plot.margin  = margin(5, 5, 5, 5)
  )

# Add Growth / Decline annotations to the legend via guide
p <- p +
  guides(fill = guide_colorbar(
    title.position = "top",
    barheight = unit(50, "mm"),
    barwidth  = unit(4, "mm"),
    label.position = "right",
    ticks = TRUE
  )) +
  annotate("text", x = -Inf, y = -Inf, label = "") # placeholder; we add via patchwork below

# --- Attempt contour at log10(ratio) = 2 ---
# We try adding contours; wrap in tryCatch to skip on error
p_final <- tryCatch({
  p + geom_contour(aes(x = x, y = y, z = log_ratio),
                   breaks = 2, colour = "black", linewidth = 0.3,
                   linetype = "dashed")
}, error = function(e) {
  message("Contour at log10(ratio)=2 skipped: ", e$message)
  p
})

# Add Growth / Decline text annotations at legend ends using patchwork caption
# Instead, add them as custom labels in the plot margins
p_final <- p_final +
  labs(caption = "Blue = Growth (\u03b5 < \u03b3)          Red = Decline (\u03b5 > \u03b3)") +
  theme(
    plot.caption = element_text(hjust = 0.5, size = 8, colour = "grey40",
                                margin = margin(t = 5))
  )

# --- Save ---
ggsave("figs/pub_fig3_bivariate_risk_map.png", p_final,
       width = 180, height = 160, units = "mm", dpi = 300)

ggsave("figs/pub_fig3_bivariate_risk_map.pdf", p_final,
       width = 180, height = 160, units = "mm")

message("Done: figs/pub_fig3_bivariate_risk_map.{png,pdf}")
