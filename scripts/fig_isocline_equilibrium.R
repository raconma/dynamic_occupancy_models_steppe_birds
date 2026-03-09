###############################################################################
# fig_isocline_equilibrium.R
#
# Key summary figure for GCB: colonisation–extinction isocline plot.
#
# Panel (a): γ vs ε scatter on log–log axes. The γ = ε diagonal separates
#            population growth (below) from decline (above). Dashed contours
#            show equilibrium occupancy ψ* = γ/(γ+ε). Species are labelled
#            directly. All four species fall deep in the decline zone.
#
# Panel (b): Extinction debt. Stacked horizontal bars showing current
#            occupancy vs equilibrium ψ*. Debt = (current − ψ*) / current.
#
# Output: figs/pub_fig_isocline_equilibrium.png  (300 DPI, 12 × 6.5 in)
#         figs/pub_fig_isocline_equilibrium.pdf
###############################################################################

library(here)
library(ggplot2)
library(dplyr)
library(scales)
library(patchwork)

theme_set(theme_bw(base_size = 12))

###############################################################################
# 1. Load data
###############################################################################

iso_data  <- read.csv(here("results", "isocline_plot_data.csv"))
debt_data <- read.csv(here("results", "extinction_debt_table.csv"))
eq_data   <- read.csv(here("results", "equilibrium_occupancy_table.csv"))

# Merge
plot_data <- iso_data %>%
  left_join(debt_data, by = "species") %>%
  left_join(eq_data %>% select(species, psi_star_med, psi_star_lo, psi_star_hi,
                                recol_yr_med),
            by = "species")

# Species labels (italic for plot)
sp_labels <- c(
  "otitar" = "O. tarda",
  "ptealc" = "P. alchata",
  "pteori" = "P. orientalis",
  "tettet" = "T. tetrax"
)
plot_data$label <- sp_labels[plot_data$species]

###############################################################################
# 2. Colour palette (colourblind-safe, Wong 2011)
###############################################################################

sp_colours <- c(
  "otitar" = "#D55E00",   # vermillion
  "ptealc" = "#0072B2",   # blue
  "pteori" = "#009E73",   # bluish green
  "tettet" = "#CC79A7"    # reddish purple
)

sp_shapes <- c(
  "otitar" = 16,   # circle
  "ptealc" = 17,   # triangle
  "pteori" = 15,   # square
  "tettet" = 18    # diamond
)

###############################################################################
# 3. Build Panel (a): γ vs ε isocline plot
###############################################################################
message("Building isocline plot...")

# Axis ranges (log10 scale)
# x-axis starts at 1e-8 to include pteori (γ_baseline ≈ 3.8 × 10⁻⁷)
x_lo <- 1e-8
x_hi <- 1
y_lo <- 0.01    # extinction rates are all > 0.01
y_hi <- 1

# --- Isocline contours: ψ* = γ/(γ+ε)  →  ε = γ(1−ψ*)/ψ* ---
psi_levels <- c(0.0001, 0.001, 0.01, 0.1, 0.5)
psi_labels <- c("0.01%", "0.1%", "1%", "10%", "50%")
# Plotmath-safe labels for isocline annotation
psi_plotmath <- c(
  'psi*"*"=="0.01%"', 'psi*"*"=="0.1%"', 'psi*"*"=="1%"',
  'psi*"*"=="10%"', 'psi*"*"=="50%"'
)

gamma_seq <- 10^seq(log10(x_lo), log10(x_hi), length.out = 400)

iso_list <- lapply(seq_along(psi_levels), function(i) {
  eps <- gamma_seq * (1 - psi_levels[i]) / psi_levels[i]
  df <- data.frame(gamma = gamma_seq, epsilon = eps,
                   psi_label = psi_labels[i],
                   psi_expr  = psi_plotmath[i],
                   stringsAsFactors = FALSE)
  df[df$epsilon >= y_lo & df$epsilon <= y_hi, ]
})
iso_df <- do.call(rbind, iso_list)

# Label positions: place at the right end of each visible contour
iso_label_pos <- iso_df %>%
  group_by(psi_label) %>%
  # pick a point near the right edge of each contour within the plot area
  filter(epsilon <= 0.9 * y_hi, epsilon >= 1.2 * y_lo) %>%
  slice_max(gamma, n = 1) %>%
  ungroup()

# --- Zone shading polygons (data coordinates, log-transformed by ggplot) ---
# Decline zone: upper-left triangle where ε > γ
decline_poly <- data.frame(
  x = c(x_lo, y_lo, y_hi, x_lo),
  y = c(y_lo, y_lo, y_hi, y_hi)
)
# Growth zone: lower-right triangle where γ > ε
growth_poly <- data.frame(
  x = c(y_lo, x_hi, x_hi),
  y = c(y_lo, y_lo, y_hi)
)

# --- Build plot ---
p_a <- ggplot() +

  # Zone shading (drawn first, behind everything)
  geom_polygon(data = decline_poly, aes(x = x, y = y),
               fill = "#FFEBEE", alpha = 0.5) +
  geom_polygon(data = growth_poly, aes(x = x, y = y),
               fill = "#E8F5E9", alpha = 0.5) +

  # Zone labels (clean text, no arrows, no boxes)
  # "decline" label: upper-left, away from P. orientalis (γ≈3.8e-7, ε≈0.44)
  annotate("text", x = 1e-7, y = 0.8,
           label = "Decline",
           colour = "#666666", size = 4.5, fontface = "bold.italic",
           alpha = 0.7) +
  # "growth" label: bottom-right, below and right of legend
  annotate("text", x = 0.5, y = 0.013,
           label = "Growth",
           colour = "#666666", size = 4.5, fontface = "bold.italic",
           alpha = 0.7) +

  # ψ* isocline contours
  geom_line(data = iso_df,
            aes(x = gamma, y = epsilon, group = psi_label),
            colour = "grey55", linewidth = 0.35, linetype = "dashed") +

  # Isocline labels (ψ* values) — use parse = TRUE for Greek letters
  geom_text(data = iso_label_pos,
            aes(x = gamma, y = epsilon, label = psi_expr),
            colour = "grey45", size = 2.5, hjust = 1, vjust = -0.5,
            parse = TRUE) +

  # γ = ε diagonal (net zero line)
  geom_abline(slope = 1, intercept = 0,
              colour = "grey30", linewidth = 0.7) +
  annotate("text", x = 0.055, y = 0.12,
           label = expression(gamma == epsilon),
           colour = "grey30", size = 3.5, angle = 42) +

  # 95% CI crosshairs
  geom_errorbar(data = plot_data,
                aes(x = gamma_baseline,
                    ymin = epsilon_lo, ymax = epsilon_hi,
                    colour = species),
                width = 0, linewidth = 0.45, alpha = 0.45) +
  geom_errorbar(data = plot_data,
                aes(y = epsilon_baseline,
                    xmin = gamma_lo, xmax = gamma_hi,
                    colour = species),
                width = 0, linewidth = 0.45, alpha = 0.45,
                orientation = "y") +

  # Species points
  geom_point(data = plot_data,
             aes(x = gamma_baseline, y = epsilon_baseline,
                 colour = species, shape = species),
             size = 4, stroke = 0.9) +

  # Scales
  scale_x_log10(
    name = expression("Colonisation rate,"~gamma~"(yr"^{-1}*")"),
    limits = c(x_lo, x_hi),
    breaks = 10^seq(-8, 0, by = 2),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  scale_y_log10(
    name = expression("Extinction rate,"~epsilon~"(yr"^{-1}*")"),
    limits = c(y_lo, y_hi),
    breaks = 10^seq(-2, 0),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  scale_colour_manual(
    values = sp_colours,
    labels = c("otitar" = expression(italic("O. tarda")),
               "ptealc" = expression(italic("P. alchata")),
               "pteori" = expression(italic("P. orientalis")),
               "tettet" = expression(italic("T. tetrax")))
  ) +
  scale_shape_manual(
    values = sp_shapes,
    labels = c("otitar" = expression(italic("O. tarda")),
               "ptealc" = expression(italic("P. alchata")),
               "pteori" = expression(italic("P. orientalis")),
               "tettet" = expression(italic("T. tetrax")))
  ) +

  # Theme — clean, legend in lower-right (empty growth zone)
  theme(
    legend.position = c(0.18, 0.82),
    legend.justification = c(0.5, 0.5),
    legend.background = element_rect(fill = "white", color = "grey40", linewidth = 0.4),
    legend.margin = margin(4, 6, 4, 6),
    legend.key.size = unit(0.45, "cm"),
    legend.text = element_text(size = 9.5),
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
    plot.margin = margin(8, 12, 8, 8)
  )


###############################################################################
# 4. Build Panel (b): Extinction debt
###############################################################################
message("Creating extinction debt panel...")

# Prepare data
debt_plot <- debt_data %>%
  mutate(
    label = sp_labels[species],
    label = factor(label, levels = rev(sp_labels)),
    # Debt region = current − ψ*
    debt_pct = current_occupancy_pct - psi_star_pct
  )

p_b <- ggplot(debt_plot) +
  # Full bar = current occupancy (light grey = "debt" portion)
  geom_col(aes(x = current_occupancy_pct, y = label),
           fill = "#FFCDD2", width = 0.55) +
  # Equilibrium portion (dark green = "secure")
  geom_col(aes(x = psi_star_pct, y = label),
           fill = "#66BB6A", width = 0.55) +
  # Debt fraction text
  geom_text(aes(x = current_occupancy_pct + 0.015, y = label,
                label = ifelse(debt_fraction >= 100,
                               paste0(round(debt_fraction), "%"),
                               paste0(round(debt_fraction), "%"))),
            hjust = 0, size = 3.3, colour = "#B71C1C", fontface = "bold") +
  # Vertical line at zero
  geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3) +
  # Scales
  scale_x_continuous(
    name = "Site occupancy (%)",
    limits = c(0, 1.6),
    expand = c(0, 0),
    breaks = seq(0, 1.5, 0.5)
  ) +
  labs(y = NULL,
       subtitle = "Extinction debt") +
  # Legend annotation inside plot
  annotate("rect", xmin = 0.85, xmax = 0.97, ymin = 0.55, ymax = 0.72,
           fill = "#66BB6A") +
  annotate("text", x = 1.0, y = 0.64, label = expression(psi*"*"),
           hjust = 0, size = 2.8) +
  annotate("rect", xmin = 0.85, xmax = 0.97, ymin = 0.35, ymax = 0.52,
           fill = "#FFCDD2") +
  annotate("text", x = 1.0, y = 0.44, label = "Transient",
           hjust = 0, size = 2.8) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
    plot.subtitle = element_text(size = 9.5, colour = "grey30"),
    axis.text.y = element_text(face = "italic", size = 10),
    plot.margin = margin(8, 12, 8, 8)
  )


###############################################################################
# 5. Combine panels
###############################################################################
message("Combining panels...")

combined <- p_a + p_b +
  plot_layout(widths = c(3, 2)) +
  plot_annotation(
    tag_levels = "a",
    tag_prefix = "(",
    tag_suffix = ")",
    theme = theme(
      plot.tag = element_text(size = 13, face = "bold")
    )
  )


###############################################################################
# 6. Save
###############################################################################
message("Saving figure...")

ggsave(here("figs", "pub_fig_isocline_equilibrium.png"),
       combined, width = 12, height = 5.5, dpi = 300, bg = "white")
message("  Saved: figs/pub_fig_isocline_equilibrium.png")

ggsave(here("figs", "pub_fig_isocline_equilibrium.pdf"),
       combined, width = 12, height = 5.5, bg = "white")
message("  Saved: figs/pub_fig_isocline_equilibrium.pdf")

message("\nFigure complete.")
