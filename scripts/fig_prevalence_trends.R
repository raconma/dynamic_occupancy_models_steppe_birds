###############################################################################
# fig_prevalence_trends.R
#
# Publication figure: Occupancy prevalence trends 2017-2023 for 4 Iberian
# steppe bird species.
#
# Data sources (tried in order):
#   1. results/{sp}_simulation_prevalence.csv  (pre-computed)
#   2. Forward simulation from model intercepts + bootstrap CI
#
# Output:
#   figs/pub_fig_prevalence_trends.png  (300 DPI)
#   figs/pub_fig_prevalence_trends.pdf
###############################################################################

library(unmarked)
library(ggplot2)
library(dplyr)
library(patchwork)
library(MASS)  # mvrnorm

cat("\n", strrep("=", 70), "\n")
cat("  FIGURE: Occupancy prevalence trends 2017-2023\n")
cat(strrep("=", 70), "\n\n")

# --- Paths ---
res_dir  <- "results"
figs_dir <- "figs"
if (!dir.exists(figs_dir)) dir.create(figs_dir, recursive = TRUE)

# --- Species info ---
sp_codes <- c("otitar", "ptealc", "pteori", "tettet")
sp_latin <- c(otitar = "Otis tarda", ptealc = "Pterocles alchata",
              pteori = "Pterocles orientalis", tettet = "Tetrax tetrax")
sp_colours <- c(otitar = "#D55E00", ptealc = "#0072B2",
                pteori = "#009E73", tettet = "#CC79A7")

years <- 2017:2023

# --- Load bootstrap draws ---
boot_path <- file.path(res_dir, "bootstrap_draws_5000.rds")
boot_draws <- NULL
if (file.exists(boot_path)) {
  boot_draws <- readRDS(boot_path)
  cat("  Bootstrap draws loaded (5000 per species)\n")
} else {
  warning("bootstrap_draws_5000.rds not found; will use mvrnorm for CIs")
}

# --- Forward simulation function ---
simulate_prevalence <- function(psi1, gam, eps, years) {
  n <- length(years)
  psi <- numeric(n)
  psi[1] <- psi1
  for (t in 2:n) {
    psi[t] <- psi[t - 1] * (1 - eps) + (1 - psi[t - 1]) * gam
  }
  data.frame(year = years, prevalence = psi * 100)
}

# --- Process each species ---
all_trends <- list()

for (sp in sp_codes) {
  cat("\n  Processing", sp, "...\n")

  # Try pre-computed CSV
  csv_path <- file.path(res_dir, paste0(sp, "_simulation_prevalence.csv"))

  if (file.exists(csv_path)) {
    cat("    Using pre-computed CSV\n")
    sim <- read.csv(csv_path)

    # The CSV has: year, mean_occ, sd_occ, mean_col, sd_col, mean_ext, sd_ext
    # mean_occ is mean prevalence across sites (proportion)
    trend_df <- data.frame(
      year     = sim$year,
      mean_prev = sim$mean_occ * 100,
      stringsAsFactors = FALSE
    )

    # Compute CI from bootstrap draws if available
    if (!is.null(boot_draws) && sp %in% names(boot_draws)) {
      bd <- boot_draws[[sp]]
      gam_draws <- bd$gamma_draws
      eps_draws <- bd$epsilon_draws

      # Load model to get psi1 intercept
      mod_path <- file.path(res_dir, paste0(sp, "_model_object.rds"))
      if (file.exists(mod_path)) {
        mod <- readRDS(mod_path)
        cc <- coef(mod)
        psi1_base <- plogis(cc["psi(Int)"])

        # Simulate forward for each bootstrap draw
        n_boot <- length(gam_draws)
        boot_mat <- matrix(NA, nrow = n_boot, ncol = length(years))
        for (b in seq_len(n_boot)) {
          sim_b <- simulate_prevalence(psi1_base, gam_draws[b], eps_draws[b], years)
          boot_mat[b, ] <- sim_b$prevalence
        }
        trend_df$lo_prev <- apply(boot_mat, 2, quantile, probs = 0.025, na.rm = TRUE)
        trend_df$hi_prev <- apply(boot_mat, 2, quantile, probs = 0.975, na.rm = TRUE)
      } else {
        # Fallback: use sd_occ for CI
        trend_df$lo_prev <- (sim$mean_occ - 1.96 * sim$sd_occ) * 100
        trend_df$hi_prev <- (sim$mean_occ + 1.96 * sim$sd_occ) * 100
      }
    } else {
      # Use sd from CSV
      trend_df$lo_prev <- (sim$mean_occ - 1.96 * sim$sd_occ) * 100
      trend_df$hi_prev <- (sim$mean_occ + 1.96 * sim$sd_occ) * 100
    }
  } else {
    # Generate from model intercepts
    cat("    CSV not found; generating from model intercepts\n")
    mod_path <- file.path(res_dir, paste0(sp, "_model_object.rds"))
    if (!file.exists(mod_path)) {
      warning("Model not found for ", sp)
      next
    }

    mod <- readRDS(mod_path)
    cc <- coef(mod)
    V  <- vcov(mod)

    psi1 <- plogis(cc["psi(Int)"])
    gam  <- plogis(cc["col(Int)"])
    eps  <- plogis(cc["ext(Int)"])

    sim_base <- simulate_prevalence(psi1, gam, eps, years)
    trend_df <- data.frame(year = years, mean_prev = sim_base$prevalence)

    # Bootstrap CI with mvrnorm
    n_boot <- 500
    draws <- mvrnorm(n_boot, cc, V)
    boot_mat <- matrix(NA, nrow = n_boot, ncol = length(years))
    for (b in seq_len(n_boot)) {
      psi1_b <- plogis(draws[b, "psi(Int)"])
      gam_b  <- plogis(draws[b, "col(Int)"])
      eps_b  <- plogis(draws[b, "ext(Int)"])
      sim_b  <- simulate_prevalence(psi1_b, gam_b, eps_b, years)
      boot_mat[b, ] <- sim_b$prevalence
    }
    trend_df$lo_prev <- apply(boot_mat, 2, quantile, probs = 0.025, na.rm = TRUE)
    trend_df$hi_prev <- apply(boot_mat, 2, quantile, probs = 0.975, na.rm = TRUE)

    # Save generated CSV
    out_csv <- file.path(res_dir, paste0(sp, "_simulation_prevalence_generated.csv"))
    write.csv(trend_df, out_csv, row.names = FALSE)
    cat("    Saved:", out_csv, "\n")
  }

  # Clamp lower CI at 0
  trend_df$lo_prev <- pmax(trend_df$lo_prev, 0)

  trend_df$species <- sp
  all_trends[[sp]] <- trend_df
}

trends <- bind_rows(all_trends)

# --- Compute equilibrium occupancy per species ---
eq_occ <- data.frame()
for (sp in sp_codes) {
  mod_path <- file.path(res_dir, paste0(sp, "_model_object.rds"))
  if (!file.exists(mod_path)) next
  mod <- readRDS(mod_path)
  cc <- coef(mod)
  gam <- plogis(cc["col(Int)"])
  eps <- plogis(cc["ext(Int)"])
  psi_star <- gam / (gam + eps) * 100
  eq_occ <- rbind(eq_occ, data.frame(species = sp, psi_star = psi_star))
}

# --- Build 2x2 panel figure ---
panel_list <- list()

for (sp in sp_codes) {
  df_sp <- trends %>% filter(species == sp)
  if (nrow(df_sp) == 0) next

  psi_star_val <- eq_occ$psi_star[eq_occ$species == sp]
  psi_star_label <- sprintf("psi* = %.2f%%", psi_star_val)

  p <- ggplot(df_sp, aes(x = year, y = mean_prev)) +
    geom_ribbon(aes(ymin = lo_prev, ymax = hi_prev),
                fill = sp_colours[sp], alpha = 0.2) +
    geom_line(colour = sp_colours[sp], linewidth = 1) +
    geom_point(colour = sp_colours[sp], size = 2) +
    geom_hline(yintercept = psi_star_val, linetype = "dashed",
               colour = "grey40", linewidth = 0.5) +
    annotate("text", x = max(years), y = psi_star_val,
             label = psi_star_label,
             hjust = 1, vjust = -0.5, size = 3, colour = "grey30") +
    scale_x_continuous(breaks = years) +
    labs(title = bquote(italic(.(sp_latin[sp]))),
         x = "Year", y = "% cells occupied") +
    theme_classic(base_size = 10) +
    theme(
      plot.title = element_text(face = "italic", size = 10),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 9),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  panel_list[[sp]] <- p
}

p_final <- wrap_plots(panel_list, ncol = 2)

# --- Save ---
ggsave(file.path(figs_dir, "pub_fig_prevalence_trends.png"), p_final,
       width = 180, height = 160, units = "mm", dpi = 300)
ggsave(file.path(figs_dir, "pub_fig_prevalence_trends.pdf"), p_final,
       width = 180, height = 160, units = "mm")

cat("\n  -> pub_fig_prevalence_trends.png\n")
cat("  -> pub_fig_prevalence_trends.pdf\n")
cat("\nDone.\n")
