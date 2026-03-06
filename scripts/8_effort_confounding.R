###############################################################################
# 8_effort_confounding.R
#
# Purpose: Diagnose whether apparent occupancy trends (2017-2023) are
#          confounded with increasing eBird survey effort over time.
#
# Tasks:
#   1. Effort trend diagnostics (raw metrics + Spearman correlations)
#   2. Model-estimated detection trends (p_hat, psi_hat, naive/model ratio)
#   3. Residual effort confounding test (binomial GLM)
#   4. Cross-species summary table + interpretation
#
# Inputs:
#   data/processed_2023/{sp}/{sp}_occ_wide_dynamic.csv
#   results/results2/{sp}_{pipeline}_all_models.rds
#   results/results2/{sp}_{pipeline}_umf.rds
#
# Outputs:
#   figs/pub_fig_effort_diagnostics_{sp}.png     (per-species, 4 panels)
#   figs/pub_fig_effort_cross_species.png         (cross-species summary)
#   results/pub_effort_confounding_summary.csv    (summary table)
#   results/pub_effort_interpretation.txt         (plain-text paragraphs)
#
###############################################################################

# -- Load packages --
library(here)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(unmarked)

select  <- dplyr::select
filter  <- dplyr::filter

# Load model configurations (for detection formula)
source(here("R", "model_configs.R"))

# -- Species info (consistent with publication_figures_models.R) --
sp_codes <- c("otitar", "ptealc", "pteori", "tettet")
sp_latin <- c(otitar = "Otis tarda", ptealc = "Pterocles alchata",
              pteori = "Pterocles orientalis", tettet = "Tetrax tetrax")
sp_short <- c(otitar = "O. tarda", ptealc = "P. alchata",
              pteori = "P. orientalis", tettet = "T. tetrax")
sp_colors <- c("O. tarda" = "#1B9E77", "P. alchata" = "#D95F02",
               "P. orientalis" = "#7570B3", "T. tetrax" = "#E7298A")

# Best model info (from publication_figures_models.R)
best_info <- list(
  otitar = list(pipe = "v4b", name = "m9: ambos lc_crop"),
  ptealc = list(pipe = "v5",  name = "m6: gam~pr_lag"),
  pteori = list(pipe = "v4b", name = "m0: Baseline estatico"),
  tettet = list(pipe = "v4b", name = "m7: gam~grass+crop")
)

# -- Constants --
YEARS_ALL   <- 2017:2023   # full data range
YEARS_MODEL <- 2017:2022   # range of colext models
J_reps      <- 10           # max secondary periods

# -- Publication theme --
theme_pub <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(face = "italic"),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    axis.title = element_text(size = 10)
  )

# -- Helper: extract observation-level columns in year-major, rep-minor order --
get_obs_cols <- function(df, prefix, years, reps = 1:10) {
  cols <- c()
  for (yr in years) {
    for (k in reps) {
      cn <- paste0(prefix, ".", k, ".", yr)
      if (cn %in% names(df)) cols <- c(cols, cn)
    }
  }
  cols
}

# -- Directories --
figs_dir    <- here("figs")
results_dir <- here("results")
rds_dir     <- here("results", "results2")

###############################################################################
# LOAD COLEXT MODELS
###############################################################################
cat("\n", strrep("=", 70), "\n")
cat("  EFFORT CONFOUNDING DIAGNOSTICS\n")
cat(strrep("=", 70), "\n\n")

cat("Loading colext models...\n")
best_mods <- list()
best_umfs <- list()
for (sp in sp_codes) {
  info <- best_info[[sp]]
  mods <- readRDS(file.path(rds_dir, paste0(sp, "_", info$pipe, "_all_models.rds")))
  best_mods[[sp]] <- mods[[info$name]]
  best_umfs[[sp]] <- readRDS(file.path(rds_dir, paste0(sp, "_", info$pipe, "_umf.rds")))
  cat("  ", sp, " (", info$pipe, "): ", info$name,
      " — AIC: ", round(best_mods[[sp]]@AIC, 1), "\n")
}

###############################################################################
# TASK 1: EFFORT TREND DIAGNOSTICS
###############################################################################
cat("\n--- Task 1: Effort trend diagnostics ---\n")

effort_metrics <- data.frame()

for (sp in sp_codes) {
  cat("  Processing", sp, "...\n")

  # Load wide-format data
  occ_wide <- read.csv(
    here("data", "processed_2023", sp, paste0(sp, "_occ_wide_dynamic.csv"))
  )

  for (yr in YEARS_ALL) {

    # Detection history columns for this year
    y_cols <- paste0("y.", 1:J_reps, ".", yr)
    y_cols <- y_cols[y_cols %in% names(occ_wide)]
    if (length(y_cols) == 0) next

    y_mat <- as.matrix(occ_wide[, y_cols])

    # Duration columns for this year
    dur_cols <- paste0("duration_minutes.", 1:J_reps, ".", yr)
    dur_cols <- dur_cols[dur_cols %in% names(occ_wide)]
    dur_mat <- if (length(dur_cols) > 0) as.matrix(occ_wide[, dur_cols]) else NULL

    # 1. Total checklists (non-NA observations)
    n_checklists <- sum(!is.na(y_mat))

    # 2. Number of sites with at least one visit
    sites_visited <- sum(rowSums(!is.na(y_mat)) > 0)

    # 3. Mean visits per visited site
    visits_per_site <- rowSums(!is.na(y_mat))
    mean_visits <- mean(visits_per_site[visits_per_site > 0])

    # 4. Mean checklist duration (across all non-NA observations)
    mean_duration <- if (!is.null(dur_mat)) mean(dur_mat, na.rm = TRUE) else NA

    # 5. Number of sites with at least 1 detection (naive occupied)
    any_det <- rowSums(y_mat, na.rm = TRUE) > 0 & rowSums(!is.na(y_mat)) > 0
    n_naive_occ <- sum(any_det)

    # 6. Naive occupancy (among visited sites)
    naive_occ <- n_naive_occ / sites_visited

    # 7. Proportion unsampled
    n_sites <- nrow(occ_wide)
    prop_unsampled <- sum(visits_per_site == 0) / n_sites

    effort_metrics <- rbind(effort_metrics, data.frame(
      species       = sp,
      species_short = sp_short[sp],
      year          = yr,
      n_checklists  = n_checklists,
      sites_visited = sites_visited,
      mean_visits   = round(mean_visits, 2),
      mean_duration = round(mean_duration, 1),
      n_naive_occ   = n_naive_occ,
      naive_occ     = round(naive_occ, 4),
      prop_unsampled = round(prop_unsampled, 4),
      n_sites       = n_sites,
      stringsAsFactors = FALSE
    ))
  }
}

effort_metrics$species_short <- factor(effort_metrics$species_short,
  levels = c("O. tarda", "P. alchata", "P. orientalis", "T. tetrax"))

cat("  Effort metrics computed for", nrow(effort_metrics), "species-years.\n")

# -- Spearman correlations --
cat("\n  Spearman correlations:\n")
spearman_results <- data.frame()

for (sp in sp_codes) {
  d <- effort_metrics %>% filter(species == sp)

  metrics <- c("n_checklists", "sites_visited", "mean_visits",
               "mean_duration", "naive_occ", "prop_unsampled")
  labels  <- c("Total checklists", "Sites visited", "Mean visits/site",
               "Mean duration (min)", "Naive occupancy", "Prop. unsampled")

  for (i in seq_along(metrics)) {
    ct <- cor.test(d$year, d[[metrics[i]]], method = "spearman", exact = FALSE)
    spearman_results <- rbind(spearman_results, data.frame(
      species = sp,
      species_short = sp_short[sp],
      metric  = labels[i],
      rho     = round(ct$estimate, 3),
      p_value = round(ct$p.value, 4),
      stringsAsFactors = FALSE
    ))
    cat("    ", sp, " — ", labels[i], ": rho=", round(ct$estimate, 3),
        ", p=", round(ct$p.value, 4), "\n")
  }
}

###############################################################################
# TASK 2: MODEL-ESTIMATED DETECTION TRENDS
###############################################################################
cat("\n--- Task 2: Model-estimated detection & occupancy trends ---\n")

model_trends <- data.frame()

for (sp in sp_codes) {
  cat("  Processing", sp, "...\n")

  mod  <- best_mods[[sp]]
  umf  <- best_umfs[[sp]]
  nY   <- umf@numPrimary

  # --- Detection probability p_hat per year ---
  # Compute manually from coefficients to avoid UMF version issues.
  # Strategy: extract detection beta coefficients, compute mean covariate
  # values per year from the raw data, then compute mean p per year.

  occ_wide <- read.csv(
    here("data", "processed_2023", sp, paste0(sp, "_occ_wide_dynamic.csv"))
  )
  cfg <- get_model_config(sp)
  det_vars <- all.vars(cfg$p_formula)
  det_coefs <- coef(mod, type = "det")

  # Map formula variable names to raw data prefixes
  var_prefix_map <- c(
    effort    = "effort_distance_km",
    duration  = "duration_minutes",
    observers = "number_observers",
    time      = "time_observations_started"
  )

  p_annual <- numeric(nY)
  for (t in seq_len(nY)) {
    yr <- YEARS_MODEL[t]
    logit_p <- det_coefs["p(Int)"]  # intercept

    for (v in det_vars) {
      # Determine source of covariate
      if (v %in% names(var_prefix_map)) {
        # Observation-level: get raw values for this year
        prefix <- var_prefix_map[v]
        v_cols <- paste0(prefix, ".", 1:J_reps, ".", yr)
        v_cols <- v_cols[v_cols %in% names(occ_wide)]
        if (length(v_cols) > 0) {
          raw_vals <- as.matrix(occ_wide[, v_cols])
          # Compute global mean and sd for standardisation (as in model fitting)
          # For obs-level covariates, scale() was applied across ALL years
          all_cols <- get_obs_cols(occ_wide, prefix, YEARS_MODEL)
          all_vals <- as.matrix(occ_wide[, all_cols])
          global_mean <- mean(all_vals, na.rm = TRUE)
          global_sd   <- sd(as.vector(all_vals), na.rm = TRUE)
          v_scaled_mean <- mean((raw_vals - global_mean) / global_sd, na.rm = TRUE)
        } else {
          v_scaled_mean <- 0
        }
      } else if (grepl("_obs$", v)) {
        # Expanded yearly/static covariate (e.g., NDVI_obs, pr_obs, topo_aspect_obs)
        base_var <- sub("_obs$", "", v)
        if (base_var %in% c("NDVI", "pr", "tmmn", "tmmx", "EVI")) {
          # Dynamic yearly: get this year's values
          col_name <- paste0(base_var, "_", yr)
          if (col_name %in% names(occ_wide)) {
            raw_vals <- occ_wide[[col_name]]
            # Scale using all years
            all_yearly <- occ_wide[, paste0(base_var, "_", YEARS_MODEL)]
            all_yearly_mat <- as.matrix(all_yearly)
            v_scaled_mean <- mean(
              (raw_vals - mean(all_yearly_mat[, t])) /
                sd(all_yearly_mat[, t]),
              na.rm = TRUE)
            # NOTE: scale() on a matrix gives per-column center/scale
            # Mean of scaled values for a column is ~0 by definition
            v_scaled_mean <- 0  # by construction, mean of scaled values ≈ 0
          } else {
            v_scaled_mean <- 0
          }
        } else {
          # Static expanded (topo_aspect_obs, topo_elev_obs)
          v_scaled_mean <- 0  # mean of scaled static covariate ≈ 0
        }
      } else {
        v_scaled_mean <- 0
      }

      coef_name <- paste0("p(", v, ")")
      if (coef_name %in% names(det_coefs)) {
        logit_p <- logit_p + det_coefs[coef_name] * v_scaled_mean
      }
    }

    p_annual[t] <- plogis(logit_p)
  }

  # --- Model-estimated occupancy psi_hat per year (smoothed) ---
  sm <- smoothed(mod)
  psi_hat <- as.numeric(sm["occupied", ])

  # --- Naive occupancy from effort_metrics (matching model years) ---
  d_eff <- effort_metrics %>%
    filter(species == sp, year %in% YEARS_MODEL[1:nY])

  model_trends <- rbind(model_trends, data.frame(
    species       = sp,
    species_short = sp_short[sp],
    year          = YEARS_MODEL[1:nY],
    p_hat         = round(p_annual, 4),
    psi_hat       = round(psi_hat, 4),
    naive_occ     = d_eff$naive_occ,
    ratio_naive_model = round(d_eff$naive_occ / psi_hat, 4),
    stringsAsFactors = FALSE
  ))
}

model_trends$species_short <- factor(model_trends$species_short,
  levels = c("O. tarda", "P. alchata", "P. orientalis", "T. tetrax"))

cat("  Model trends computed.\n")

###############################################################################
# TASK 3: RESIDUAL EFFORT CONFOUNDING TEST (binomial GLM)
###############################################################################
cat("\n--- Task 3: Binomial GLM — naive_occ ~ year + log(n_checklists) ---\n")

glm_results <- data.frame()

for (sp in sp_codes) {
  d <- effort_metrics %>% filter(species == sp)

  # Binomial GLM: prop sites detected ~ year + log(checklists)
  glm_fit <- glm(
    cbind(n_naive_occ, sites_visited - n_naive_occ) ~ year + log(n_checklists),
    data = d, family = binomial(link = "logit")
  )

  s <- summary(glm_fit)
  coefs <- s$coefficients

  # Year coefficient
  year_est <- coefs["year", "Estimate"]
  year_se  <- coefs["year", "Std. Error"]
  year_p   <- coefs["year", "Pr(>|z|)"]
  year_lo  <- year_est - 1.96 * year_se
  year_hi  <- year_est + 1.96 * year_se

  # log(checklists) coefficient
  effort_est <- coefs["log(n_checklists)", "Estimate"]
  effort_se  <- coefs["log(n_checklists)", "Std. Error"]
  effort_p   <- coefs["log(n_checklists)", "Pr(>|z|)"]

  glm_results <- rbind(glm_results, data.frame(
    species       = sp,
    species_short = sp_short[sp],
    year_coef     = round(year_est, 4),
    year_ci_lo    = round(year_lo, 4),
    year_ci_hi    = round(year_hi, 4),
    year_p        = round(year_p, 4),
    effort_coef   = round(effort_est, 4),
    effort_p      = round(effort_p, 4),
    stringsAsFactors = FALSE
  ))

  cat("  ", sp, ": year=", round(year_est, 4),
      " [", round(year_lo, 4), ",", round(year_hi, 4), "] p=", round(year_p, 4),
      " | log(checklists)=", round(effort_est, 4), " p=", round(effort_p, 4), "\n")
}

###############################################################################
# TASK 4: SUMMARY TABLE
###############################################################################
cat("\n--- Task 4: Summary table ---\n")

summary_table <- data.frame()

for (sp in sp_codes) {
  # Spearman: year ~ naive_occ
  sp_naive <- spearman_results %>%
    filter(species == sp, metric == "Naive occupancy")
  # Spearman: year ~ total checklists
  sp_effort <- spearman_results %>%
    filter(species == sp, metric == "Total checklists")

  # GLM results
  glm_sp <- glm_results %>% filter(species == sp)

  # Verdict logic
  # If year coef is significant AND effort coef is not → likely real

  # If year coef is NOT significant AND effort coef is significant → likely confounded
  # Otherwise → ambiguous
  if (glm_sp$year_p < 0.05 & glm_sp$effort_p >= 0.05) {
    verdict <- "likely real trend"
  } else if (glm_sp$year_p >= 0.05 & glm_sp$effort_p < 0.05) {
    verdict <- "likely confounded"
  } else if (glm_sp$year_p >= 0.05 & glm_sp$effort_p >= 0.05) {
    # Neither significant — check Spearman
    if (abs(sp_naive$rho) > 0.7 & abs(sp_effort$rho) > 0.7) {
      verdict <- "ambiguous (both correlated)"
    } else {
      verdict <- "ambiguous"
    }
  } else {
    # Both significant
    if (abs(glm_sp$year_coef) > abs(glm_sp$effort_coef)) {
      verdict <- "likely real trend (effort also significant)"
    } else {
      verdict <- "ambiguous (both significant)"
    }
  }

  summary_table <- rbind(summary_table, data.frame(
    Species              = sp_short[sp],
    rho_naive_occ        = sp_naive$rho,
    p_naive_occ          = sp_naive$p_value,
    rho_checklists       = sp_effort$rho,
    p_checklists         = sp_effort$p_value,
    GLM_year_coef        = glm_sp$year_coef,
    GLM_year_CI          = paste0("[", glm_sp$year_ci_lo, ", ", glm_sp$year_ci_hi, "]"),
    GLM_year_p           = glm_sp$year_p,
    GLM_log_effort_coef  = glm_sp$effort_coef,
    GLM_log_effort_p     = glm_sp$effort_p,
    Verdict              = verdict,
    stringsAsFactors = FALSE
  ))
}

write.csv(summary_table,
          here("results", "pub_effort_confounding_summary.csv"),
          row.names = FALSE)
cat("  Summary table saved: results/pub_effort_confounding_summary.csv\n")

###############################################################################
# FIGURES: Per-species diagnostics (4 panels each)
###############################################################################
cat("\n--- Generating per-species diagnostic figures ---\n")

for (sp in sp_codes) {
  cat("  Plotting", sp, "...\n")

  d_eff <- effort_metrics %>% filter(species == sp)
  d_mod <- model_trends %>% filter(species == sp)
  col <- sp_colors[sp_short[sp]]

  # Panel A: Total checklists over time
  pA <- ggplot(d_eff, aes(x = year, y = n_checklists)) +
    geom_line(color = "grey40", linewidth = 0.8) +
    geom_point(color = "grey40", size = 2.5) +
    labs(x = NULL, y = "Total checklists",
         title = paste0(sp_latin[sp], " (", sp, ")"),
         subtitle = "A) eBird survey effort") +
    scale_x_continuous(breaks = YEARS_ALL) +
    theme_pub

  # Panel B: Naive occupancy + checklists (dual axis)
  scale_factor <- max(d_eff$naive_occ) / max(d_eff$n_checklists)

  pB <- ggplot(d_eff, aes(x = year)) +
    geom_line(aes(y = naive_occ), color = col, linewidth = 1) +
    geom_point(aes(y = naive_occ), color = col, size = 2.5) +
    geom_line(aes(y = n_checklists * scale_factor),
              color = "grey50", linewidth = 0.8, linetype = "dashed") +
    scale_y_continuous(
      name = "Naive occupancy",
      sec.axis = sec_axis(~ . / scale_factor, name = "Total checklists")
    ) +
    scale_x_continuous(breaks = YEARS_ALL) +
    labs(x = NULL,
         subtitle = "B) Naive occupancy (solid) vs. effort (dashed)") +
    theme_pub +
    theme(axis.title.y.right = element_text(color = "grey50"),
          axis.text.y.right = element_text(color = "grey50"))

  # Panel C: Detection probability p_hat over time (model years only)
  pC <- ggplot(d_mod, aes(x = year, y = p_hat)) +
    geom_line(color = "steelblue", linewidth = 0.8) +
    geom_point(color = "steelblue", size = 2.5) +
    labs(x = NULL, y = expression(hat(p)),
         subtitle = "C) Mean detection probability (colext model)") +
    scale_x_continuous(breaks = YEARS_MODEL) +
    coord_cartesian(ylim = c(0, max(d_mod$p_hat) * 1.3)) +
    theme_pub

  # Panel D: Ratio naive / model occupancy
  pD <- ggplot(d_mod, aes(x = year, y = ratio_naive_model)) +
    geom_line(color = "darkorange", linewidth = 0.8) +
    geom_point(color = "darkorange", size = 2.5) +
    geom_hline(yintercept = 1, linetype = "dotted", color = "grey50") +
    labs(x = "Year", y = expression(psi[naive] / hat(psi)),
         subtitle = "D) Naive / model-estimated occupancy ratio") +
    scale_x_continuous(breaks = YEARS_MODEL) +
    theme_pub

  # Combine
  p_combined <- pA / pB / pC / pD +
    plot_annotation(
      title = paste0("Effort confounding diagnostics: ", sp_latin[sp]),
      theme = theme(plot.title = element_text(face = "bold", size = 14))
    )

  ggsave(here("figs", paste0("pub_fig_effort_diagnostics_", sp, ".png")),
         p_combined, width = 8, height = 12, dpi = 300)
  cat("    -> pub_fig_effort_diagnostics_", sp, ".png\n")
}

###############################################################################
# FIGURE: Cross-species summary
###############################################################################
cat("\n--- Generating cross-species summary figure ---\n")

# Panel 1: Checklists trend all species
p1 <- ggplot(effort_metrics, aes(x = year, y = n_checklists,
                                  color = species_short)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_color_manual(values = sp_colors) +
  scale_x_continuous(breaks = YEARS_ALL) +
  labs(x = NULL, y = "Total checklists",
       subtitle = "A) eBird survey effort by species",
       color = "Species") +
  theme_pub

# Panel 2: Naive occupancy all species
p2 <- ggplot(effort_metrics, aes(x = year, y = naive_occ,
                                  color = species_short)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_color_manual(values = sp_colors) +
  scale_x_continuous(breaks = YEARS_ALL) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = NULL, y = "Naive occupancy",
       subtitle = "B) Proportion of visited sites with detection",
       color = "Species") +
  theme_pub

# Panel 3: Detection probability all species
p3 <- ggplot(model_trends, aes(x = year, y = p_hat,
                                color = species_short)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_color_manual(values = sp_colors) +
  scale_x_continuous(breaks = YEARS_MODEL) +
  labs(x = NULL, y = expression("Mean " * hat(p)),
       subtitle = "C) Model-estimated detection probability",
       color = "Species") +
  theme_pub

# Panel 4: GLM coefficients (year effect controlling for effort)
glm_plot_data <- glm_results
glm_plot_data$species_short <- factor(glm_plot_data$species_short,
  levels = c("O. tarda", "P. alchata", "P. orientalis", "T. tetrax"))

p4 <- ggplot(glm_plot_data, aes(x = species_short, y = year_coef,
                                 color = species_short)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_pointrange(aes(ymin = year_ci_lo, ymax = year_ci_hi), size = 0.8) +
  scale_color_manual(values = sp_colors) +
  labs(x = NULL, y = "Year coefficient (logit scale)",
       subtitle = "D) GLM year effect controlling for log(checklists)") +
  theme_pub +
  theme(legend.position = "none",
        axis.text.x = element_text(face = "italic"))

# Combine
p_cross <- (p1 + p2) / (p3 + p4) +
  plot_annotation(
    title = "Cross-species effort confounding assessment",
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  ) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(here("figs", "pub_fig_effort_cross_species.png"),
       p_cross, width = 12, height = 9, dpi = 300)
cat("  -> pub_fig_effort_cross_species.png\n")

###############################################################################
# INTERPRETATION PARAGRAPHS
###############################################################################
cat("\n--- Generating interpretation ---\n")

sink(here("results", "pub_effort_interpretation.txt"))

cat("EFFORT CONFOUNDING ASSESSMENT\n")
cat("Generated:", as.character(Sys.time()), "\n")
cat(strrep("=", 70), "\n\n")

for (sp in sp_codes) {
  st <- summary_table %>% filter(Species == sp_short[sp])
  d_eff <- effort_metrics %>% filter(species == sp)
  d_mod <- model_trends %>% filter(species == sp)

  cat(sp_latin[sp], " (", sp_short[sp], ")\n", sep = "")
  cat(strrep("-", 40), "\n")

  # Compute key diagnostics for interpretation
  chk_change <- round((tail(d_eff$n_checklists, 1) / head(d_eff$n_checklists, 1)
                        - 1) * 100, 0)
  naive_change <- round((tail(d_eff$naive_occ, 1) - head(d_eff$naive_occ, 1)) * 100, 1)
  p_trend <- cor.test(d_mod$year, d_mod$p_hat, method = "spearman",
                       exact = FALSE)

  cat("eBird checklists changed by ", ifelse(chk_change > 0, "+", ""),
      chk_change, "% between ", min(d_eff$year), " and ",
      max(d_eff$year), ". ", sep = "")

  cat("Naive occupancy ",
      ifelse(naive_change > 0, "increased", "decreased"),
      " by ", abs(naive_change), " percentage points ",
      "(Spearman rho = ", st$rho_naive_occ, ", p = ", st$p_naive_occ, "). ", sep = "")

  cat("The Spearman correlation between year and total checklists was ",
      st$rho_checklists, " (p = ", st$p_checklists, "). ", sep = "")

  cat("Mean model-estimated detection probability ",
      ifelse(p_trend$estimate > 0.3, "showed an increasing trend",
             ifelse(p_trend$estimate < -0.3, "showed a decreasing trend",
                    "remained relatively stable")),
      " over the study period (rho = ", round(p_trend$estimate, 3), "). ", sep = "")

  cat("After controlling for log(checklists) in a binomial GLM, the year effect was ",
      round(st$GLM_year_coef, 4), " (95% CI: ", st$GLM_year_CI,
      ", p = ", st$GLM_year_p, "). ", sep = "")

  cat("VERDICT: ", st$Verdict, ".\n\n", sep = "")
}

sink()
cat("  -> results/pub_effort_interpretation.txt\n")

###############################################################################
# DONE
###############################################################################
cat("\n", strrep("=", 70), "\n")
cat("  EFFORT CONFOUNDING DIAGNOSTICS COMPLETE\n")
cat("  Figures: figs/pub_fig_effort_diagnostics_{sp}.png (x4)\n")
cat("           figs/pub_fig_effort_cross_species.png\n")
cat("  Table:   results/pub_effort_confounding_summary.csv\n")
cat("  Text:    results/pub_effort_interpretation.txt\n")
cat(strrep("=", 70), "\n")
