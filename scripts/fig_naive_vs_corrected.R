#!/usr/bin/env Rscript
# fig_naive_vs_corrected.R
# 2-row x 4-column dot plot: naive vs detection-corrected transition rates

library(ggplot2)
library(dplyr)
library(tidyr)

# ── Read data ──
d <- read.csv("results/naive_vs_corrected_full.csv", stringsAsFactors = FALSE)

# ── Species labels & colours ──
sp_labels <- c(otitar = "italic('O. tarda')", ptealc = "italic('P. alchata')",
               pteori = "italic('P. orientalis')", tettet = "italic('T. tetrax')")
sp_colours <- c(otitar = "#D55E00", ptealc = "#0072B2",
                pteori = "#009E73", tettet = "#CC79A7")

# Convert proportions to percentages
d <- d %>%
  mutate(across(c(naive_gamma, corrected_gamma_median,
                  corrected_gamma_lo, corrected_gamma_hi,
                  naive_epsilon, corrected_epsilon_median,
                  corrected_epsilon_lo, corrected_epsilon_hi),
                ~ . * 100))

# ── Detection floor ──
detection_floor <- 1e-3  # 0.001 %

d <- d %>%
  mutate(gamma_below_floor = corrected_gamma_median < detection_floor,
         corrected_gamma_plot    = pmax(corrected_gamma_median, detection_floor),
         corrected_gamma_lo_plot = pmax(corrected_gamma_lo, detection_floor),
         corrected_gamma_hi_plot = pmax(corrected_gamma_hi, detection_floor))

# ── Reshape to long format for faceting ──
gamma_df <- d %>%
  select(species, naive_gamma, corrected_gamma_plot,
         corrected_gamma_lo_plot, corrected_gamma_hi_plot,
         ratio_gamma, gamma_below_floor) %>%
  rename(naive = naive_gamma, corrected = corrected_gamma_plot,
         lo = corrected_gamma_lo_plot, hi = corrected_gamma_hi_plot,
         ratio = ratio_gamma, below_floor = gamma_below_floor) %>%
  mutate(process = "Colonisation~(gamma)")

epsilon_df <- d %>%
  select(species, naive_epsilon, corrected_epsilon_median,
         corrected_epsilon_lo, corrected_epsilon_hi, ratio_epsilon) %>%
  rename(naive = naive_epsilon, corrected = corrected_epsilon_median,
         lo = corrected_epsilon_lo, hi = corrected_epsilon_hi,
         ratio = ratio_epsilon) %>%
  mutate(below_floor = FALSE, process = "Extinction~(epsilon)")

plot_df <- bind_rows(gamma_df, epsilon_df) %>%
  mutate(species_label = factor(sp_labels[species], levels = sp_labels),
         process = factor(process, levels = c("Colonisation~(gamma)",
                                              "Extinction~(epsilon)")),
         ratio_label = ifelse(ratio >= 10,
                              paste0(round(ratio), "x"),
                              paste0(round(ratio, 1), "x")))

# ── Use numeric y: 1 = Corrected, 2 = Naive ──

# ── Plot ──
p <- ggplot(plot_df) +
  # Dashed connector between naive (y=2) and corrected (y=1)
  geom_segment(aes(x = naive, xend = corrected, y = 2, yend = 1),
               linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  # Error bars on corrected (y=1)
  geom_errorbar(aes(xmin = lo, xmax = hi, y = 1),
                width = 0.15, colour = "grey40", linewidth = 0.35,
                orientation = "y") +
  # Naive point (open circle, grey50)
  geom_point(aes(x = naive, y = 2),
             shape = 21, size = 3, fill = "white", colour = "grey50") +
  # Corrected point (filled, species colour)
  geom_point(aes(x = corrected, y = 1, colour = species),
             shape = 16, size = 3) +
  scale_colour_manual(values = sp_colours, guide = "none") +
  # "< detection floor" annotation
  geom_text(data = plot_df %>% filter(below_floor),
            aes(x = corrected, y = 1, label = "< detection floor"),
            hjust = -0.1, vjust = -1.2,
            fontface = "italic", colour = "grey50", size = 2.5) +
  # Ratio annotation (upper-right)
  geom_text(aes(x = Inf, y = Inf, label = ratio_label),
            hjust = 1.15, vjust = 1.6,
            colour = "grey40", size = 3.5) +
  facet_grid(process ~ species_label,
             labeller = label_parsed) +
  scale_x_log10() +
  scale_y_continuous(breaks = c(1, 2), labels = c("Corrected", "Naive"),
                     limits = c(0.5, 2.5)) +
  labs(x = "Transition rate (%, log scale)", y = NULL) +
  theme_classic(base_size = 10) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 9),
    panel.spacing = unit(0.6, "lines"),
    axis.text.y = element_text(size = 8)
  )

# ── Save ──
dir.create("figs", showWarnings = FALSE)

ggsave("figs/pub_fig1_naive_vs_corrected.png", p,
       width = 180, height = 120, units = "mm", dpi = 300)
ggsave("figs/pub_fig1_naive_vs_corrected.pdf", p,
       width = 180, height = 120, units = "mm")

message("Done: saved pub_fig1_naive_vs_corrected.png and .pdf")
