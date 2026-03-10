#!/usr/bin/env Rscript
# =============================================================================
# fig_naive_vs_corrected.R
#
# Horizontal arrow plot: naive vs detection-corrected transition rates
# Species on y-axis, arrows from naive → corrected on log x-axis
# Two panels: Colonisation (gamma) | Extinction (epsilon) with free x-scales
# Fold-reduction annotated on each arrow
# =============================================================================

library(ggplot2)
library(dplyr)
library(patchwork)

# ── Read data ──
d <- read.csv("results/naive_vs_corrected_full.csv", stringsAsFactors = FALSE)

# ── Species labels & colours ──
sp_order <- c("otitar", "ptealc", "pteori", "tettet")
sp_labels <- c(otitar = "O. tarda", ptealc = "P. alchata",
               pteori = "P. orientalis", tettet = "T. tetrax")
sp_colours <- c(otitar = "#D55E00", ptealc = "#0072B2",
                pteori = "#009E73", tettet = "#CC79A7")

# Convert proportions to percentages
d <- d %>%
  mutate(across(c(naive_gamma, corrected_gamma_median,
                  corrected_gamma_lo, corrected_gamma_hi,
                  naive_epsilon, corrected_epsilon_median,
                  corrected_epsilon_lo, corrected_epsilon_hi),
                ~ . * 100))

# ── Detection floor for gamma ──
detection_floor <- 1e-3  # 0.001 %

d <- d %>%
  mutate(
    gamma_below_floor = corrected_gamma_median < detection_floor,
    corrected_gamma_plot    = pmax(corrected_gamma_median, detection_floor),
    corrected_gamma_lo_plot = pmax(corrected_gamma_lo, detection_floor),
    corrected_gamma_hi_plot = pmax(corrected_gamma_hi, detection_floor),
    species_label = factor(sp_labels[species], levels = rev(sp_labels))
  )

# ── Format fold-change labels ──
format_ratio <- function(r) {
  ifelse(r >= 1000, paste0(format(round(r), big.mark = ","), "\u00d7"),
  ifelse(r >= 10, paste0(round(r), "\u00d7"),
         paste0(round(r, 1), "\u00d7")))
}

# ── Panel A: Colonisation (gamma) ──
p_gamma <- ggplot(d) +
  # Arrow from naive to corrected
  geom_segment(aes(x = naive_gamma, xend = corrected_gamma_plot,
                   y = species_label, yend = species_label,
                   colour = species),
               arrow = arrow(length = unit(2, "mm"), type = "closed"),
               linewidth = 0.7) +
  # Error bar on corrected
  geom_errorbar(orientation = "y",
aes(xmin = corrected_gamma_lo_plot,
                     xmax = corrected_gamma_hi_plot,
                     y = species_label),
                 width = 0.15, colour = "grey40", linewidth = 0.35) +
  # Naive point (open circle)
  geom_point(aes(x = naive_gamma, y = species_label),
             shape = 21, size = 3.5, fill = "white", colour = "grey50",
             stroke = 0.8) +
  # Corrected point (filled, species colour)
  geom_point(aes(x = corrected_gamma_plot, y = species_label,
                 colour = species),
             shape = 16, size = 3) +
  # Fold-change annotation
  geom_text(aes(x = sqrt(naive_gamma * corrected_gamma_plot),
                y = species_label,
                label = format_ratio(ratio_gamma)),
            vjust = -1.0, size = 3, colour = "grey30") +
  # "< detection floor" annotation
  geom_text(data = d %>% filter(gamma_below_floor),
            aes(x = corrected_gamma_plot, y = species_label,
                label = "< floor"),
            hjust = -0.2, vjust = 1.8,
            fontface = "italic", colour = "grey50", size = 2.3) +
  scale_colour_manual(values = sp_colours, guide = "none") +
  scale_x_log10(labels = function(x) {
    ifelse(x >= 0.01, format(x, drop0trailing = TRUE),
           format(x, scientific = TRUE, digits = 1))
  }) +
  labs(title = expression(Colonisation~(gamma)),
       x = "Transition rate (%, log scale)", y = NULL) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.y = element_text(face = "italic", size = 9),
    plot.title = element_text(size = 10, hjust = 0.5),
    plot.margin = margin(5, 10, 5, 5)
  )

# ── Panel B: Extinction (epsilon) ──
p_epsilon <- ggplot(d) +
  # Arrow from naive to corrected
  geom_segment(aes(x = naive_epsilon, xend = corrected_epsilon_median,
                   y = species_label, yend = species_label,
                   colour = species),
               arrow = arrow(length = unit(2, "mm"), type = "closed"),
               linewidth = 0.7) +
  # Error bar on corrected
  geom_errorbar(orientation = "y",
aes(xmin = corrected_epsilon_lo,
                     xmax = corrected_epsilon_hi,
                     y = species_label),
                 width = 0.15, colour = "grey40", linewidth = 0.35) +
  # Naive point (open circle)
  geom_point(aes(x = naive_epsilon, y = species_label),
             shape = 21, size = 3.5, fill = "white", colour = "grey50",
             stroke = 0.8) +
  # Corrected point (filled, species colour)
  geom_point(aes(x = corrected_epsilon_median, y = species_label,
                 colour = species),
             shape = 16, size = 3) +
  # Fold-change annotation
  geom_text(aes(x = sqrt(naive_epsilon * corrected_epsilon_median),
                y = species_label,
                label = format_ratio(ratio_epsilon)),
            vjust = -1.0, size = 3, colour = "grey30") +
  scale_colour_manual(values = sp_colours, guide = "none") +
  scale_x_log10(labels = function(x) format(x, drop0trailing = TRUE)) +
  labs(title = expression(Extinction~(epsilon)),
       x = "Transition rate (%, log scale)", y = NULL) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.y = element_text(face = "italic", size = 9),
    plot.title = element_text(size = 10, hjust = 0.5),
    plot.margin = margin(5, 10, 5, 5)
  )

# ── Combine ──
p_final <- p_gamma + p_epsilon +
  plot_annotation(
    caption = "Open circles = naive rates; filled points = detection-corrected (bootstrap median + 95% CI)"
  ) &
  theme(plot.caption = element_text(size = 7.5, colour = "grey50", hjust = 0.5))

# ── Save ──
dir.create("figs", showWarnings = FALSE)

ggsave("figs/pub_fig1_naive_vs_corrected.png", p_final,
       width = 200, height = 100, units = "mm", dpi = 300)
ggsave("figs/pub_fig1_naive_vs_corrected.pdf", p_final,
       width = 200, height = 100, units = "mm")

message("Done: saved pub_fig1_naive_vs_corrected.png and .pdf")
