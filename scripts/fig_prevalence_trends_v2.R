#!/usr/bin/env Rscript
# ============================================================================
# fig_prevalence_trends_v2.R
# Option A: Model-simulated prevalence trends with parametric bootstrap CIs
# Option B (conditional): Compact bar chart if all CIs are flat
# ============================================================================

library(unmarked)
library(ggplot2)
library(dplyr)
library(patchwork)

# ---- Configuration --------------------------------------------------------
sp_codes  <- c("otitar", "ptealc", "pteori", "tettet")
sp_latin  <- c(otitar = "Otis tarda", ptealc = "Pterocles alchata",
               pteori = "Pterocles orientalis", tettet = "Tetrax tetrax")
sp_colours <- c(otitar = "#D55E00", ptealc = "#0072B2",
                pteori = "#009E73", tettet = "#CC79A7")
years <- 2017:2023
n_years <- length(years)

# ---- Load data -------------------------------------------------------------
boot <- readRDS("results/bootstrap_draws_5000.rds")

# Extract psi1 (intercept-only, logistic) from each model
psi1_vals <- setNames(numeric(length(sp_codes)), sp_codes)
gam_vals  <- setNames(numeric(length(sp_codes)), sp_codes)
eps_vals  <- setNames(numeric(length(sp_codes)), sp_codes)

for (sp in sp_codes) {
  m <- readRDS(paste0("results/", sp, "_model_object.rds"))
  psi1_vals[sp] <- plogis(coef(m, type = "psi")[1])
  gam_vals[sp]  <- plogis(coef(m, type = "col")[1])
  eps_vals[sp]  <- plogis(coef(m, type = "ext")[1])
}

# ---- Forward simulation function ------------------------------------------
forward_sim <- function(psi1, gamma_vec, epsilon_vec, n_years) {
  # gamma_vec and epsilon_vec are vectors of length n_draws
  n_draws <- length(gamma_vec)
  psi_mat <- matrix(NA, nrow = n_draws, ncol = n_years)
  psi_mat[, 1] <- psi1
  for (t in 2:n_years) {
    psi_mat[, t] <- psi_mat[, t - 1] * (1 - epsilon_vec) +
                    (1 - psi_mat[, t - 1]) * gamma_vec
  }
  psi_mat
}

# ---- Simulate and build data for plotting ----------------------------------
plot_data <- list()

for (sp in sp_codes) {
  gam_draws <- boot[[sp]]$gamma_draws
  eps_draws <- boot[[sp]]$epsilon_draws
  psi1 <- psi1_vals[sp]

  psi_mat <- forward_sim(psi1, gam_draws, eps_draws, n_years)

  # Summaries per year
  mean_psi <- colMeans(psi_mat)
  ci_lo    <- apply(psi_mat, 2, quantile, probs = 0.025)
  ci_hi    <- apply(psi_mat, 2, quantile, probs = 0.975)

  plot_data[[sp]] <- data.frame(
    species = sp,
    year    = years,
    mean    = mean_psi * 100,
    lo      = ci_lo * 100,
    hi      = ci_hi * 100
  )
}

df <- bind_rows(plot_data)

# ---- Equilibrium occupancy (psi*) -----------------------------------------
psi_star <- gam_vals / (gam_vals + eps_vals) * 100

# ---- Build Option A panels -------------------------------------------------
panels <- list()
flat_flags <- logical(length(sp_codes))
names(flat_flags) <- sp_codes

for (i in seq_along(sp_codes)) {
  sp <- sp_codes[i]
  d  <- filter(df, species == sp)

  # Check flatness: CI width > 2x absolute trend (year 1 to year 7)
  mean_ci_width <- mean(d$hi - d$lo)
  abs_trend     <- abs(d$mean[n_years] - d$mean[1])
  flat_flags[sp] <- mean_ci_width > 2 * abs_trend

  y_max <- max(d$hi) * 1.10

  # Label for psi*
  psi_star_label <- sprintf("psi* = %.2f%%", psi_star[sp])

  panels[[sp]] <- ggplot(d, aes(x = year)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = sp_colours[sp], alpha = 0.25) +
    geom_line(aes(y = mean), colour = sp_colours[sp], linewidth = 0.8) +
    geom_point(aes(y = mean), colour = sp_colours[sp], size = 2) +
    geom_hline(yintercept = psi_star[sp], linetype = "dashed", colour = "grey40") +
    annotate("text", x = 2020, y = psi_star[sp],
             label = psi_star_label, vjust = -0.6, size = 3, colour = "grey30") +
    scale_x_continuous(breaks = years) +
    coord_cartesian(ylim = c(0, y_max)) +
    labs(title = bquote(italic(.(sp_latin[sp]))),
         x = NULL, y = "% cells occupied") +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(size = 10))
}

fig_a <- (panels[[1]] + panels[[2]]) / (panels[[3]] + panels[[4]])

# Save Option A
dir.create("figs", showWarnings = FALSE)
ggsave("figs/pub_fig7_prevalence_trends_A.png", fig_a,
       width = 8, height = 6, dpi = 300)
ggsave("figs/pub_fig7_prevalence_trends_A.pdf", fig_a,
       width = 8, height = 6)
cat("Option A saved.\n")

# ---- Check if Option B needed ---------------------------------------------
all_flat <- all(flat_flags)
cat("Flat CI flags:", paste(names(flat_flags), flat_flags, sep = "=", collapse = ", "), "\n")

if (all_flat) {
  cat("WARNING: All four species have flat CI bands. Generating Option B.\n")

  # Build summary data for Option B
  bar_data <- data.frame()
  for (sp in sp_codes) {
    gam_draws <- boot[[sp]]$gamma_draws
    eps_draws <- boot[[sp]]$epsilon_draws

    bar_data <- bind_rows(bar_data, data.frame(
      species  = sp,
      latin    = sp_latin[sp],
      param    = c("psi1", "gamma_bar", "epsilon_bar"),
      value    = c(psi1_vals[sp], mean(gam_draws), mean(eps_draws)) * 100,
      lo       = c(NA, quantile(gam_draws, 0.025), quantile(eps_draws, 0.025)) * 100,
      hi       = c(NA, quantile(gam_draws, 0.975), quantile(eps_draws, 0.975)) * 100
    ))
  }

  bar_data$param <- factor(bar_data$param,
    levels = c("psi1", "gamma_bar", "epsilon_bar"),
    labels = c(expression(psi[1]), expression(bar(gamma)), expression(bar(epsilon))))

  # Use species latin name as factor
  bar_data$latin <- factor(bar_data$latin, levels = sp_latin)

  # Param as character labels for dodging
  bar_data$param_label <- rep(c("psi1", "gamma_bar", "epsilon_bar"),
                              times = length(sp_codes))
  bar_data$param_label <- factor(bar_data$param_label,
    levels = c("psi1", "gamma_bar", "epsilon_bar"),
    labels = c("\u03C8\u2081", "\u03B3\u0305", "\u03B5\u0305"))

  fig_b <- ggplot(bar_data, aes(x = latin, y = value, fill = param_label)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.8) +
    geom_errorbar(aes(ymin = lo, ymax = hi),
                  position = position_dodge(width = 0.7), width = 0.2,
                  na.rm = TRUE) +
    scale_fill_manual(
      values = c("\u03C8\u2081" = "grey50",
                 "\u03B3\u0305" = "#4DAF4A",
                 "\u03B5\u0305" = "#E41A1C"),
      name = "Parameter"
    ) +
    labs(x = NULL, y = "Probability (%)") +
    theme_classic(base_size = 10) +
    theme(axis.text.x = element_text(face = "italic", angle = 20, hjust = 1))

  ggsave("figs/pub_fig7_prevalence_trends_B.png", fig_b,
         width = 7, height = 4.5, dpi = 300)
  ggsave("figs/pub_fig7_prevalence_trends_B.pdf", fig_b,
         width = 7, height = 4.5)
  cat("Option B saved.\n")
} else {
  cat("Option B not needed (not all CIs are flat).\n")
}

cat("Done.\n")
