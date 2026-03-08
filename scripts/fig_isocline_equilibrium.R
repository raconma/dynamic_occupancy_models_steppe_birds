###############################################################################
# fig_isocline_equilibrium.R
#
# Creates the key summary figure for GCB: an isocline plot showing
# colonisation vs extinction rates for all 4 species, with the γ=ε
# diagonal separating growth from decline. Uses PhyloPic bird
# silhouettes and shows equilibrium occupancy (ψ*) + extinction debt.
#
# Output: figs/pub_fig_isocline_equilibrium.png
#         figs/pub_fig_isocline_equilibrium.pdf
###############################################################################

library(here)
library(ggplot2)
library(dplyr)
library(rphylopic)
library(cowplot)
library(grid)
library(scales)

theme_set(theme_bw(base_size = 13))

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

# Species labels
plot_data$label <- c(
  "otitar" = "O. tarda",
  "ptealc" = "P. alchata",
  "pteori" = "P. orientalis",
  "tettet" = "T. tetrax"
)[plot_data$species]

# Formatted labels with psi* and debt
plot_data$annotation <- paste0(
  plot_data$label, "\n",
  "psi* = ", sprintf("%.2f", plot_data$psi_star_med), "%\n",
  "debt = ", sprintf("%.0f", plot_data$debt_fraction), "%"
)

###############################################################################
# 2. Fetch PhyloPic silhouettes
###############################################################################
message("Fetching PhyloPic silhouettes...")

# UUIDs found via search
phylo_uuids <- c(
  otitar = "4291fa38-2fbd-4b3e-970c-ee8a937e4f4d",  # Otididae (family)
  ptealc = "ca5a230c-93a6-4e49-a4d6-82af5d140f3f",  # Pterocles alchata
  pteori = "7d6a29b5-ff04-45da-b625-4b6534f8725b",  # Pterocles orientalis
  tettet = "52ac389c-7936-48f0-86cd-cc03b651a81a"    # Tetrax tetrax
)

# Pre-fetch images
phylo_imgs <- list()
for (sp in names(phylo_uuids)) {
  tryCatch({
    phylo_imgs[[sp]] <- get_phylopic(uuid = phylo_uuids[[sp]])
    message("  Got silhouette for ", sp)
  }, error = function(e) {
    message("  Failed for ", sp, ": ", e$message)
    phylo_imgs[[sp]] <<- NULL
  })
}

###############################################################################
# 3. Species-specific colours (colourblind-safe palette)
###############################################################################

sp_colours <- c(
  "otitar" = "#D55E00",   # vermillion
  "ptealc" = "#0072B2",   # blue
  "pteori" = "#009E73",   # bluish green
  "tettet" = "#CC79A7"    # reddish purple
)

sp_shapes <- c(
  "otitar" = 16,
  "ptealc" = 17,
  "pteori" = 15,
  "tettet" = 18
)

###############################################################################
# 4. Build main isocline plot
###############################################################################
message("Building isocline plot...")

# Use mean (site-averaged) rates for position, baseline for CIs
# For the isocline, we plot baseline (intercept) rates with bootstrap CIs

# Define axis range (log10 scale)
# x-axis extends to 1e-8 to include pteori (gamma_baseline ~ 3.8e-7)
axis_min_x <- 1e-8
axis_min_y <- 1e-3
axis_max   <- 1

# Generate isocline contours (ψ* = γ/(γ+ε) = constant)
psi_star_levels <- c(0.0001, 0.001, 0.01, 0.1, 0.5)
iso_lines <- expand.grid(
  gamma = 10^seq(log10(axis_min_x), log10(axis_max), length.out = 300)
)
iso_contours <- list()
for (psi in psi_star_levels) {
  # ψ* = γ/(γ+ε) → ε = γ(1-ψ*)/ψ*
  df <- data.frame(
    gamma = iso_lines$gamma,
    epsilon = iso_lines$gamma * (1 - psi) / psi,
    psi_star_label = paste0(psi * 100, "%")
  )
  df <- df[df$epsilon >= axis_min_y & df$epsilon <= axis_max, ]
  iso_contours <- c(iso_contours, list(df))
}
iso_contour_df <- do.call(rbind, iso_contours)

# Main plot
p_main <- ggplot() +
  # Shade region above γ=ε diagonal (decline zone: ε > γ)
  annotate("polygon",
           x = c(axis_min_x, axis_max, axis_max, axis_min_x),
           y = c(axis_min_x, axis_max, axis_max, axis_max),
           fill = "#FFE0E0", alpha = 0.4) +
  # Shade region below γ=ε diagonal (growth zone: γ > ε)
  annotate("polygon",
           x = c(axis_min_x, axis_max, axis_min_x, axis_min_x),
           y = c(axis_min_x, axis_max, axis_min_x, axis_min_x),
           fill = "#E0FFE0", alpha = 0.4) +

  # ψ* isoclines
  geom_line(data = iso_contour_df,
            aes(x = gamma, y = epsilon, group = psi_star_label),
            colour = "grey60", linewidth = 0.3, linetype = "dashed") +

  # γ = ε diagonal
  geom_abline(slope = 1, intercept = 0,
              colour = "grey30", linewidth = 0.7, linetype = "solid") +

  # Error bars (bootstrap CIs)
  geom_errorbar(data = plot_data,
                aes(x = gamma_baseline,
                    ymin = epsilon_lo, ymax = epsilon_hi,
                    colour = species),
                width = 0, linewidth = 0.5, alpha = 0.6) +
  geom_errorbar(data = plot_data,
                aes(y = epsilon_baseline,
                    xmin = gamma_lo, xmax = gamma_hi,
                    colour = species),
                width = 0, linewidth = 0.5, alpha = 0.6,
                orientation = "y") +

  # Species points
  geom_point(data = plot_data,
             aes(x = gamma_baseline, y = epsilon_baseline,
                 colour = species, shape = species),
             size = 4, stroke = 1.2) +

  # Scales
  scale_x_log10(
    name = expression(paste("Colonisation rate (", gamma, ")")),
    limits = c(1e-8, 1),
    breaks = 10^(seq(-8, 0, by = 2)),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  scale_y_log10(
    name = expression(paste("Extinction rate (", epsilon, ")")),
    limits = c(1e-3, 1),
    breaks = 10^(-3:0),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  scale_colour_manual(values = sp_colours,
                      labels = c("otitar" = expression(italic("O. tarda")),
                                 "ptealc" = expression(italic("P. alchata")),
                                 "pteori" = expression(italic("P. orientalis")),
                                 "tettet" = expression(italic("T. tetrax")))) +
  scale_shape_manual(values = sp_shapes,
                     labels = c("otitar" = expression(italic("O. tarda")),
                                "ptealc" = expression(italic("P. alchata")),
                                "pteori" = expression(italic("P. orientalis")),
                                "tettet" = expression(italic("T. tetrax")))) +

  # Zone labels
  annotate("text", x = 0.1, y = 0.003,
           label = expression(paste("Growth (", gamma, " > ", epsilon, ")")),
           colour = "#2E7D32",
           size = 3.5, fontface = "italic", alpha = 0.7) +
  annotate("text", x = 1e-6, y = 0.7,
           label = expression(paste("Decline (", epsilon, " > ", gamma, ")")),
           colour = "#C62828",
           size = 3.5, fontface = "italic", alpha = 0.7) +

  # Isocline labels (positioned along the contour lines)
  annotate("text", x = 2e-6, y = 0.6,
           label = expression(paste(psi, "* = 0.01%")), colour = "grey50",
           size = 2.5, angle = 30) +
  annotate("text", x = 2e-4, y = 0.6,
           label = expression(paste(psi, "* = 1%")), colour = "grey50",
           size = 2.5, angle = 30) +
  annotate("text", x = 0.03, y = 0.6,
           label = expression(paste(psi, "* = 10%")), colour = "grey50",
           size = 2.5, angle = 30) +

  # Theme
  labs(colour = "Species", shape = "Species") +
  theme(
    legend.position = c(0.82, 0.22),
    legend.background = element_rect(fill = alpha("white", 0.9), colour = "grey80"),
    legend.key.size = unit(0.8, "lines"),
    legend.text.align = 0,
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey92"),
    plot.margin = margin(10, 15, 10, 10)
  )


###############################################################################
# 5. Add PhyloPic silhouettes
###############################################################################
message("Adding PhyloPic silhouettes...")

# Position silhouettes near each species point (adjusted for 1e-8 to 1 x-axis)
silhouette_positions <- data.frame(
  species = c("otitar", "ptealc", "pteori", "tettet"),
  x = c(5e-4, 0.01, 3e-8, 3e-3),
  y = c(0.5, 0.8, 0.2, 0.05),
  stringsAsFactors = FALSE
)

for (i in 1:nrow(silhouette_positions)) {
  sp <- silhouette_positions$species[i]
  if (!is.null(phylo_imgs[[sp]])) {
    p_main <- p_main +
      add_phylopic(
        img = phylo_imgs[[sp]],
        x = silhouette_positions$x[i],
        y = silhouette_positions$y[i],
        height = 0.08,
        color = sp_colours[sp],
        alpha = 0.35
      )
  }
}


###############################################################################
# 6. Create extinction debt bar chart (inset or panel B)
###############################################################################
message("Creating extinction debt panel...")

debt_plot_data <- debt_data %>%
  mutate(
    label = c("otitar" = "O. tarda", "ptealc" = "P. alchata",
              "pteori" = "P. orientalis", "tettet" = "T. tetrax")[species],
    label = factor(label, levels = rev(c("O. tarda", "P. alchata",
                                          "P. orientalis", "T. tetrax")))
  )

p_debt <- ggplot(debt_plot_data) +
  # Current occupancy (full bar)
  geom_col(aes(x = current_occupancy_pct, y = label),
           fill = "grey80", width = 0.6) +
  # Equilibrium (the part that will persist)
  geom_col(aes(x = psi_star_pct, y = label),
           fill = "#4CAF50", width = 0.6, alpha = 0.8) +
  # Debt annotation
  geom_text(aes(x = current_occupancy_pct + 0.02, y = label,
                label = paste0(round(debt_fraction), "% debt")),
            hjust = 0, size = 3.2, colour = "grey30") +
  # Vertical line at zero
  geom_vline(xintercept = 0, colour = "grey40") +
  # Scale
  scale_x_continuous(
    name = "Site occupancy (%)",
    limits = c(0, 1.7),
    expand = c(0, 0)
  ) +
  labs(y = NULL,
       title = "Extinction debt") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 11, face = "bold", hjust = 0),
    axis.text.y = element_text(face = "italic"),
    plot.margin = margin(5, 15, 5, 5)
  )


###############################################################################
# 7. Create recolonisation timescale panel
###############################################################################
message("Creating recolonisation timescale panel...")

recol_data <- eq_data %>%
  filter(vcov_ok == TRUE) %>%
  mutate(
    label = c("otitar" = "O. tarda", "ptealc" = "P. alchata",
              "pteori" = "P. orientalis", "tettet" = "T. tetrax")[species],
    label = factor(label, levels = rev(c("O. tarda", "P. alchata",
                                          "P. orientalis", "T. tetrax")))
  )

p_recol <- ggplot(recol_data) +
  geom_point(aes(x = recol_yr_med, y = label, colour = species),
             size = 3, show.legend = FALSE) +
  geom_errorbar(aes(xmin = recol_yr_lo, xmax = pmin(recol_yr_hi, 1e7),
                     y = label, colour = species),
                width = 0.2, show.legend = FALSE,
                orientation = "y") +
  scale_x_log10(
    name = "Years to recolonise (1/gamma)",
    breaks = c(100, 1000, 10000, 100000, 1000000),
    labels = c("100", "1K", "10K", "100K", "1M")
  ) +
  scale_colour_manual(values = sp_colours) +
  # Reference lines
  geom_vline(xintercept = 100, linetype = "dotted", colour = "grey60") +
  annotate("text", x = 80, y = 0.5, label = "100 yr", colour = "grey50",
           size = 2.5, angle = 90) +
  labs(y = NULL,
       title = "Recolonisation timescale") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 11, face = "bold", hjust = 0),
    axis.text.y = element_text(face = "italic"),
    plot.margin = margin(5, 15, 5, 5)
  )


###############################################################################
# 8. Combine panels
###############################################################################
message("Combining panels...")

library(patchwork)

combined <- p_main +
  (p_debt / p_recol) +
  plot_layout(widths = c(3, 2)) +
  plot_annotation(
    tag_levels = list(c("a", "b", "c")),
    title = "Colonisation–extinction dynamics reveal extinction debt in Iberian steppe birds",
    subtitle = "Baseline transition rates (intercepts ± 95% CI) with equilibrium occupancy contours",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10, colour = "grey40")
    )
  )


###############################################################################
# 9. Save
###############################################################################
message("Saving figure...")

ggsave(here("figs", "pub_fig_isocline_equilibrium.png"),
       combined, width = 14, height = 7, dpi = 300, bg = "white")
message("  Saved: figs/pub_fig_isocline_equilibrium.png")

ggsave(here("figs", "pub_fig_isocline_equilibrium.pdf"),
       combined, width = 14, height = 7, bg = "white")
message("  Saved: figs/pub_fig_isocline_equilibrium.pdf")

# Print PhyloPic attributions
cat("\n=== PhyloPic credits ===\n")
for (sp in names(phylo_uuids)) {
  tryCatch({
    attr_info <- get_attribution(uuid = phylo_uuids[[sp]])
    cat(sp, ":", attr_info$contributor, "(", attr_info$license, ")\n")
  }, error = function(e) {
    cat(sp, ": attribution unavailable\n")
  })
}

message("\nFigure complete.")
