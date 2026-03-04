###############################################################################
# publication_figures_models.R
#
# Genera figuras avanzadas para publicacion usando los objetos modelo (.rds)
# de los pipelines v4b y v5 para las 4 especies esteparias.
#
# Requiere:
#   results/results2/{sp}_{pipeline}_all_models.rds
#   results/results2/{sp}_{pipeline}_umf.rds
#
# Output (en figs/):
#   pub_fig_occupancy_trends.png        - Tendencias de ocupacion 2017-2022
#   pub_fig_col_ext_rates.png           - Tasas colonizacion/extincion
#   pub_fig_response_psi_tree.png       - Curva respuesta psi ~ tree_cover
#   pub_fig_response_eps_tree.png       - Curva respuesta epsilon ~ tree_cover
#   pub_fig_response_gam_dynamic.png    - Curvas respuesta colonizacion dinamicas
#   pub_fig_coef_all_submodels.png      - Forest plot todos los coeficientes
#   pub_fig_model_selection_heatmap.png - Heatmap seleccion de modelos
#   pub_fig_detection_covariates.png    - Efecto covariables de deteccion
###############################################################################

library(unmarked)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(tidyr)

cat("\n")
cat(strrep("=", 70), "\n")
cat("  FIGURAS AVANZADAS PARA PUBLICACION\n")
cat("  Usando objetos modelo (.rds)\n")
cat(strrep("=", 70), "\n\n")

# --- Paths ---
rds_dir <- "results/results2"
results_dir <- "results"
figs_dir <- "figs"
if (!dir.exists(figs_dir)) dir.create(figs_dir, recursive = TRUE)

# --- Species info ---
sp_codes <- c("otitar", "ptealc", "pteori", "tettet")
sp_latin <- c(otitar = "Otis tarda", ptealc = "Pterocles alchata",
              pteori = "Pterocles orientalis", tettet = "Tetrax tetrax")
sp_common <- c(otitar = "Great Bustard", ptealc = "Pin-tailed Sandgrouse",
               pteori = "Black-bellied Sandgrouse", tettet = "Little Bustard")
sp_short <- c(otitar = "O. tarda", ptealc = "P. alchata",
              pteori = "P. orientalis", tettet = "T. tetrax")
sp_colors <- c("O. tarda" = "#1B9E77", "P. alchata" = "#D95F02",
               "P. orientalis" = "#7570B3", "T. tetrax" = "#E7298A")

years <- 2017:2022

# --- Best model per species ---
best_info <- list(
  otitar = list(pipe = "v4b", name = "m9: ambos lc_crop"),
  ptealc = list(pipe = "v5",  name = "m6: gam~pr_lag"),
  pteori = list(pipe = "v4b", name = "m0: Baseline estatico"),
  tettet = list(pipe = "v4b", name = "m7: gam~grass+crop")
)

# --- Load all best models ---
cat("Cargando modelos...\n")
best_mods <- list()
best_umfs <- list()
for (sp in sp_codes) {
  info <- best_info[[sp]]
  mods <- readRDS(file.path(rds_dir, paste0(sp, "_", info$pipe, "_all_models.rds")))
  best_mods[[sp]] <- mods[[info$name]]
  best_umfs[[sp]] <- readRDS(file.path(rds_dir, paste0(sp, "_", info$pipe, "_umf.rds")))
  cat("  ", sp, "(", info$pipe, "):", info$name, "- AIC:", round(best_mods[[sp]]@AIC, 1), "\n")
}

# --- Publication theme ---
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

###############################################################################
# FIGURA 1: Tendencias de ocupacion (projected + smoothed)
###############################################################################

cat("\n  Figura 1: Tendencias de ocupacion...\n")

trend_data <- data.frame()
for (sp in sp_codes) {
  mod <- best_mods[[sp]]
  proj <- projected(mod)
  sm <- smoothed(mod)

  trend_data <- rbind(trend_data,
    data.frame(
      Species = sp_short[sp],
      Year = years,
      Projected = as.numeric(proj["occupied", ]),
      Smoothed = as.numeric(sm["occupied", ]),
      stringsAsFactors = FALSE
    )
  )
}

trend_long <- trend_data %>%
  pivot_longer(cols = c(Projected, Smoothed),
               names_to = "Method", values_to = "Occupancy")

trend_long$Species <- factor(trend_long$Species,
                              levels = c("O. tarda", "P. alchata",
                                         "P. orientalis", "T. tetrax"))

p_trend <- ggplot(trend_long, aes(x = Year, y = Occupancy, color = Species,
                                   linetype = Method)) +
  geom_line(linewidth = 1) +
  geom_point(data = trend_long %>% filter(Method == "Smoothed"),
             size = 2, show.legend = FALSE) +
  scale_color_manual(values = sp_colors) +
  scale_linetype_manual(values = c("Projected" = "dashed", "Smoothed" = "solid")) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  scale_x_continuous(breaks = years) +
  labs(title = "Occupancy trends for four steppe bird species (2017-2022)",
       subtitle = "Solid: smoothed estimates; Dashed: projected from transition probabilities",
       x = "Year", y = "Occupancy probability",
       color = "Species", linetype = "Method") +
  theme_pub +
  theme(legend.box = "vertical")

ggsave(file.path(figs_dir, "pub_fig_occupancy_trends.png"), p_trend,
       width = 9, height = 5.5, dpi = 300)
cat("    -> pub_fig_occupancy_trends.png\n")

# --- Per-species panel version ---
p_trend_panel <- ggplot(trend_long %>% filter(Method == "Smoothed"),
                         aes(x = Year, y = Occupancy, color = Species)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.5) +
  scale_color_manual(values = sp_colors) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  scale_x_continuous(breaks = years) +
  facet_wrap(~Species, scales = "free_y", ncol = 2) +
  labs(title = "Smoothed occupancy trends by species (2017-2022)",
       x = "Year", y = "Occupancy probability") +
  theme_pub +
  theme(legend.position = "none",
        strip.text = element_text(face = "italic", size = 11))

ggsave(file.path(figs_dir, "pub_fig_occupancy_trends_panel.png"), p_trend_panel,
       width = 9, height = 7, dpi = 300)
cat("    -> pub_fig_occupancy_trends_panel.png\n")

###############################################################################
# FIGURA 2: Tasas de colonizacion y extincion por ano
###############################################################################

cat("\n  Figura 2: Tasas de colonizacion/extincion...\n")

rate_data <- data.frame()
for (sp in sp_codes) {
  mod <- best_mods[[sp]]
  umf <- best_umfs[[sp]]
  nYears <- mod@data@numPrimary

  # Get mean col/ext rates per transition (nYears-1 transitions)
  col_coefs <- coef(mod, type = "col")
  ext_coefs <- coef(mod, type = "ext")

  # For models with dynamic covariates, we need site-level predictions per year
  # Use mean covariate values per year from yearlySiteCovs
  yscovs <- yearlySiteCovs(umf)

  # Get gamma/epsilon formulas
  gam_form <- mod@gamformula
  eps_form <- mod@epsformula

  for (t in 1:(nYears - 1)) {
    # For each transition, compute mean prediction across sites
    # Simple approach: use overall coefs and mean covariate values per year

    # Gamma (colonization)
    gam_int <- col_coefs[1]
    gam_logit <- gam_int
    if (length(col_coefs) > 1) {
      # Dynamic covariates in gamma - get mean value for this transition
      gam_vars <- names(col_coefs)[-1]
      for (v in gam_vars) {
        v_clean <- sub("^col\\((.+)\\)$", "\\1", v)
        # Check yearlySiteCovs (dynamic) and siteCovs (static)
        nSites <- numSites(umf)
        if (v_clean %in% names(yscovs)) {
          start_row <- (t - 1) * nSites + 1
          end_row <- t * nSites
          mean_val <- mean(yscovs[[v_clean]][start_row:end_row], na.rm = TRUE)
        } else if (v_clean %in% names(siteCovs(umf))) {
          mean_val <- mean(siteCovs(umf)[[v_clean]], na.rm = TRUE)
        } else {
          mean_val <- 0
        }
        gam_logit <- gam_logit + col_coefs[v] * mean_val
      }
    }

    # Epsilon (extinction)
    ext_int <- ext_coefs[1]
    ext_logit <- ext_int
    if (length(ext_coefs) > 1) {
      ext_vars <- names(ext_coefs)[-1]
      for (v in ext_vars) {
        v_clean <- sub("^ext\\((.+)\\)$", "\\1", v)
        nSites <- numSites(umf)
        # Check siteCovs (static) first, then yearlySiteCovs
        if (v_clean %in% names(siteCovs(umf))) {
          mean_val <- mean(siteCovs(umf)[[v_clean]], na.rm = TRUE)
        } else if (v_clean %in% names(yscovs)) {
          start_row <- (t - 1) * nSites + 1
          end_row <- t * nSites
          mean_val <- mean(yscovs[[v_clean]][start_row:end_row], na.rm = TRUE)
        } else {
          mean_val <- 0
        }
        ext_logit <- ext_logit + ext_coefs[v] * mean_val
      }
    }

    rate_data <- rbind(rate_data, data.frame(
      Species = sp_short[sp],
      Transition = paste0(years[t], "-", years[t + 1]),
      Transition_num = t,
      Colonization = plogis(gam_logit),
      Extinction = plogis(ext_logit),
      stringsAsFactors = FALSE
    ))
  }
}

rate_long <- rate_data %>%
  pivot_longer(cols = c(Colonization, Extinction),
               names_to = "Rate", values_to = "Probability")

rate_long$Species <- factor(rate_long$Species,
                             levels = c("O. tarda", "P. alchata",
                                        "P. orientalis", "T. tetrax"))

p_rates <- ggplot(rate_long, aes(x = Transition, y = Probability,
                                  color = Species, group = interaction(Species, Rate))) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  scale_color_manual(values = sp_colors) +
  facet_wrap(~Rate, scales = "free_y") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "Estimated colonization and extinction rates (2017-2022)",
       subtitle = "Predictions at mean covariate values per transition period",
       x = "Transition period", y = "Probability",
       color = "Species") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(figs_dir, "pub_fig_col_ext_rates.png"), p_rates,
       width = 10, height = 5.5, dpi = 300)
cat("    -> pub_fig_col_ext_rates.png\n")

###############################################################################
# FIGURA 3: Curva de respuesta psi ~ tree_cover (4 especies)
###############################################################################

cat("\n  Figura 3: Curva respuesta psi ~ tree_cover...\n")

tree_seq <- seq(-2, 3, length.out = 200)
psi_curves <- data.frame()

for (sp in sp_codes) {
  mod <- best_mods[[sp]]
  psi_coefs <- coef(mod, type = "psi")

  # Build prediction: fix all covariates at 0 (mean), vary tree_cover
  psi_logit <- psi_coefs[1]  # intercept
  # Add tree_cover effect (coef name is psi(tree_cover))
  tc_name <- "psi(tree_cover)"
  if (tc_name %in% names(psi_coefs)) {
    psi_pred <- psi_logit + psi_coefs[tc_name] * tree_seq
  } else {
    psi_pred <- rep(psi_logit, length(tree_seq))
  }

  psi_curves <- rbind(psi_curves, data.frame(
    tree_cover = tree_seq,
    psi = plogis(psi_pred),
    Species = sp_short[sp],
    stringsAsFactors = FALSE
  ))
}

psi_curves$Species <- factor(psi_curves$Species,
                              levels = c("O. tarda", "P. alchata",
                                         "P. orientalis", "T. tetrax"))

p_psi_tree <- ggplot(psi_curves, aes(x = tree_cover, y = psi, color = Species)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = sp_colors) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  labs(title = "Initial occupancy probability as a function of tree cover",
       subtitle = "Other covariates held at mean values (standardized = 0)",
       x = "Tree cover (standardized)", y = expression(psi ~ "(initial occupancy)"),
       color = "Species") +
  theme_pub

ggsave(file.path(figs_dir, "pub_fig_response_psi_tree.png"), p_psi_tree,
       width = 8, height = 5, dpi = 300)
cat("    -> pub_fig_response_psi_tree.png\n")

###############################################################################
# FIGURA 4: Curva de respuesta epsilon ~ tree_cover (4 especies)
###############################################################################

cat("\n  Figura 4: Curva respuesta epsilon ~ tree_cover...\n")

eps_curves <- data.frame()

for (sp in sp_codes) {
  mod <- best_mods[[sp]]
  ext_coefs <- coef(mod, type = "ext")

  ext_logit <- ext_coefs[1]  # intercept
  tc_name <- "ext(tree_cover)"
  if (tc_name %in% names(ext_coefs)) {
    ext_pred <- ext_logit + ext_coefs[tc_name] * tree_seq
  } else {
    ext_pred <- rep(ext_logit, length(tree_seq))
  }

  eps_curves <- rbind(eps_curves, data.frame(
    tree_cover = tree_seq,
    epsilon = plogis(ext_pred),
    Species = sp_short[sp],
    stringsAsFactors = FALSE
  ))
}

eps_curves$Species <- factor(eps_curves$Species,
                              levels = c("O. tarda", "P. alchata",
                                         "P. orientalis", "T. tetrax"))

p_eps_tree <- ggplot(eps_curves, aes(x = tree_cover, y = epsilon, color = Species)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = sp_colors) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  geom_hline(yintercept = 0.5, linetype = "dotted", alpha = 0.4) +
  labs(title = "Local extinction probability as a function of tree cover",
       subtitle = "Other covariates held at mean values (standardized = 0)",
       x = "Tree cover (standardized)", y = expression(epsilon ~ "(extinction probability)"),
       color = "Species") +
  theme_pub

ggsave(file.path(figs_dir, "pub_fig_response_eps_tree.png"), p_eps_tree,
       width = 8, height = 5, dpi = 300)
cat("    -> pub_fig_response_eps_tree.png\n")

###############################################################################
# FIGURA 5: Curvas respuesta colonizacion ~ covariable dinamica clave
###############################################################################

cat("\n  Figura 5: Curvas respuesta colonizacion...\n")

# Each species has a different key dynamic predictor for colonization
# otitar: lc_crop (v4b), ptealc: pr_lag (v5), pteori: intercept-only, tettet: lc_grass+lc_crop (v4b)
cov_seq <- seq(-2, 3, length.out = 200)

gam_curves <- data.frame()

# otitar: gam ~ lc_crop
mod <- best_mods[["otitar"]]
gam_coefs <- coef(mod, type = "col")
gam_pred <- gam_coefs[1] + gam_coefs["col(lc_crop)"] * cov_seq
gam_curves <- rbind(gam_curves, data.frame(
  Covariate_value = cov_seq, Colonization = plogis(gam_pred),
  Species = "O. tarda", Covariate = "Land cover: cropland",
  stringsAsFactors = FALSE
))

# ptealc: gam ~ pr_lag
mod <- best_mods[["ptealc"]]
gam_coefs <- coef(mod, type = "col")
gam_pred <- gam_coefs[1] + gam_coefs["col(pr_lag)"] * cov_seq
gam_curves <- rbind(gam_curves, data.frame(
  Covariate_value = cov_seq, Colonization = plogis(gam_pred),
  Species = "P. alchata", Covariate = "Precipitation (lag-1)",
  stringsAsFactors = FALSE
))

# tettet: gam ~ lc_grass + lc_crop (vary lc_grass, hold lc_crop at 0)
mod <- best_mods[["tettet"]]
gam_coefs <- coef(mod, type = "col")
# Curve for lc_grass
gam_pred_grass <- gam_coefs[1] + gam_coefs["col(lc_grass)"] * cov_seq
gam_curves <- rbind(gam_curves, data.frame(
  Covariate_value = cov_seq, Colonization = plogis(gam_pred_grass),
  Species = "T. tetrax", Covariate = "Land cover: grassland",
  stringsAsFactors = FALSE
))
# Curve for lc_crop
gam_pred_crop <- gam_coefs[1] + gam_coefs["col(lc_crop)"] * cov_seq
gam_curves <- rbind(gam_curves, data.frame(
  Covariate_value = cov_seq, Colonization = plogis(gam_pred_crop),
  Species = "T. tetrax", Covariate = "Land cover: cropland",
  stringsAsFactors = FALSE
))

gam_curves$Species <- factor(gam_curves$Species,
                               levels = c("O. tarda", "P. alchata", "T. tetrax"))

p_gam <- ggplot(gam_curves, aes(x = Covariate_value, y = Colonization,
                                  color = Species, linetype = Covariate)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = sp_colors[c("O. tarda", "P. alchata", "T. tetrax")]) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, NA)) +
  geom_hline(yintercept = 0.5, linetype = "dotted", alpha = 0.3) +
  labs(title = "Colonization probability as a function of dynamic covariates",
       subtitle = "Key predictors per species; other covariates at mean values. P. orientalis excluded (intercept-only)",
       x = "Covariate value (standardized)", y = expression(gamma ~ "(colonization probability)"),
       color = "Species", linetype = "Dynamic covariate") +
  theme_pub +
  theme(legend.box = "vertical")

ggsave(file.path(figs_dir, "pub_fig_response_gam_dynamic.png"), p_gam,
       width = 9, height = 5.5, dpi = 300)
cat("    -> pub_fig_response_gam_dynamic.png\n")

###############################################################################
# FIGURA 6: Forest plot todos los coeficientes (4 submodelos x 4 especies)
###############################################################################

cat("\n  Figura 6: Forest plot completo...\n")

all_coef_data <- data.frame()

for (sp in sp_codes) {
  mod <- best_mods[[sp]]
  info <- best_info[[sp]]

  # Read from CSV for consistency (has SE and p-values)
  f <- file.path(results_dir, paste0(sp, "_", info$pipe, "_best_model_coefficients.csv"))
  if (!file.exists(f)) next
  coefs_df <- read.csv(f)

  for (i in 1:nrow(coefs_df)) {
    all_coef_data <- rbind(all_coef_data, data.frame(
      Species = sp_short[sp],
      Submodel = coefs_df$Submodel[i],
      Parameter = coefs_df$Parameter[i],
      Estimate = coefs_df$Estimate[i],
      SE = coefs_df$SE[i],
      p_value = coefs_df$P...z..[i],
      stringsAsFactors = FALSE
    ))
  }
}

# Clean submodel names
all_coef_data$Submodel <- recode(all_coef_data$Submodel,
                                  "psi" = "psi (initial occupancy)",
                                  "col" = "gamma (colonization)",
                                  "ext" = "epsilon (extinction)",
                                  "det" = "p (detection)")

# Remove intercepts for cleaner plot
plot_coefs <- all_coef_data %>%
  filter(Parameter != "(Intercept)") %>%
  mutate(
    lower = Estimate - 1.96 * SE,
    upper = Estimate + 1.96 * SE,
    Significant = ifelse(p_value < 0.05, "p < 0.05", "p >= 0.05"),
    Species = factor(Species, levels = rev(c("O. tarda", "P. alchata",
                                              "P. orientalis", "T. tetrax"))),
    Submodel = factor(Submodel, levels = c("psi (initial occupancy)",
                                            "gamma (colonization)",
                                            "epsilon (extinction)",
                                            "p (detection)"))
  )

p_forest <- ggplot(plot_coefs, aes(x = Estimate, y = Species, color = Significant)) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.25, linewidth = 0.6) +
  geom_point(size = 2) +
  scale_color_manual(values = c("p < 0.05" = "#D32F2F", "p >= 0.05" = "#9E9E9E")) +
  facet_grid(Parameter ~ Submodel, scales = "free_x", space = "free_y",
             switch = "y") +
  labs(title = "Coefficient estimates across species and submodels",
       subtitle = "95% confidence intervals from best model per species",
       x = "Estimate (logit scale)", y = "", color = "") +
  theme_pub +
  theme(
    strip.text.y.left = element_text(angle = 0, hjust = 1, size = 8),
    strip.text.x = element_text(size = 9),
    panel.spacing.x = unit(0.8, "lines"),
    strip.placement = "outside"
  )

ggsave(file.path(figs_dir, "pub_fig_coef_all_submodels.png"), p_forest,
       width = 14, height = 10, dpi = 300)
cat("    -> pub_fig_coef_all_submodels.png\n")

###############################################################################
# FIGURA 7: Heatmap de seleccion de modelos (dAIC)
###############################################################################

cat("\n  Figura 7: Heatmap seleccion de modelos...\n")

heatmap_data <- data.frame()

for (sp in sp_codes) {
  for (pipe in c("v4b", "v5")) {
    f <- file.path(results_dir, paste0(sp, "_", pipe, "_model_selection.csv"))
    if (!file.exists(f)) next
    ms <- read.csv(f)
    baseline_aic <- ms$AIC[grep("Baseline", ms$Model)]
    if (length(baseline_aic) == 0) baseline_aic <- max(ms$AIC)

    for (i in 1:nrow(ms)) {
      heatmap_data <- rbind(heatmap_data, data.frame(
        Species = sp_short[sp],
        Pipeline = pipe,
        Model = ms$Model[i],
        AIC = ms$AIC[i],
        dAIC = round(ms$AIC[i] - baseline_aic, 1),
        Boundary = grepl("Boundary", ms$Status),
        stringsAsFactors = FALSE
      ))
    }
  }
}

# Keep only v4b for heatmap (cleaner)
heat_v4b <- heatmap_data %>%
  filter(Pipeline == "v4b") %>%
  mutate(
    Species = factor(Species, levels = c("O. tarda", "P. alchata",
                                          "P. orientalis", "T. tetrax")),
    # Simplify model names
    Model_short = gsub("^m[0-9]+: ", "", Model)
  )

p_heat <- ggplot(heat_v4b, aes(x = Model_short, y = Species, fill = dAIC)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(dAIC, 1)), size = 3) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, name = expression(Delta*"AIC")) +
  labs(title = "Model comparison: dAIC vs static baseline (v4b pipeline)",
       subtitle = "Negative values indicate improvement over baseline",
       x = "", y = "") +
  theme_pub +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(face = "italic"),
    legend.position = "right"
  )

ggsave(file.path(figs_dir, "pub_fig_model_selection_heatmap.png"), p_heat,
       width = 12, height = 5, dpi = 300)
cat("    -> pub_fig_model_selection_heatmap.png\n")

###############################################################################
# FIGURA 8: Efecto de covariables de deteccion
###############################################################################

cat("\n  Figura 8: Covariables de deteccion...\n")

det_cov_seq <- seq(-2, 3, length.out = 200)
det_curves <- data.frame()
det_vars <- c("effort", "duration", "observers", "time")
det_labels <- c(effort = "Search effort (km)",
                duration = "Duration (min)",
                observers = "Number of observers",
                time = "Time of day")

for (sp in sp_codes) {
  mod <- best_mods[[sp]]
  det_coefs <- coef(mod, type = "det")

  for (dvar in det_vars) {
    # coef names have p() wrapper: p(effort), p(duration), etc.
    dvar_wrapped <- paste0("p(", dvar, ")")
    if (!(dvar_wrapped %in% names(det_coefs))) next
    det_logit <- det_coefs[1] + det_coefs[dvar_wrapped] * det_cov_seq
    det_curves <- rbind(det_curves, data.frame(
      Covariate_value = det_cov_seq,
      Detection = plogis(det_logit),
      Species = sp_short[sp],
      Covariate = det_labels[dvar],
      stringsAsFactors = FALSE
    ))
  }
}

det_curves$Species <- factor(det_curves$Species,
                              levels = c("O. tarda", "P. alchata",
                                         "P. orientalis", "T. tetrax"))

p_det <- ggplot(det_curves, aes(x = Covariate_value, y = Detection, color = Species)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = sp_colors) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 0.8)) +
  facet_wrap(~Covariate, scales = "free_x") +
  labs(title = "Detection probability as a function of survey covariates",
       subtitle = "Other detection covariates held at mean values",
       x = "Covariate value (standardized)", y = "Detection probability",
       color = "Species") +
  theme_pub

ggsave(file.path(figs_dir, "pub_fig_detection_covariates.png"), p_det,
       width = 10, height = 7, dpi = 300)
cat("    -> pub_fig_detection_covariates.png\n")

###############################################################################
# FIGURA 9: Panel resumen combinado (ocupacion + tree_cover effect)
###############################################################################

cat("\n  Figura 9: Panel resumen combinado...\n")

# Combine occupancy trend + tree cover effects in one panel
p_combined <- grid.arrange(
  p_trend + labs(title = "A) Occupancy trends (2017-2022)") +
    theme(legend.position = "right", legend.box = "vertical",
          plot.title = element_text(size = 11)),
  p_eps_tree + labs(title = "B) Tree cover effect on extinction") +
    theme(legend.position = "right",
          plot.title = element_text(size = 11)),
  p_gam + labs(title = "C) Dynamic covariates effect on colonization") +
    theme(legend.position = "right", legend.box = "vertical",
          plot.title = element_text(size = 11)),
  ncol = 1, heights = c(1, 1, 1.1)
)

ggsave(file.path(figs_dir, "pub_fig_summary_panel.png"), p_combined,
       width = 10, height = 14, dpi = 300)
cat("    -> pub_fig_summary_panel.png\n")

###############################################################################
# TABLA EXTRA: Tendencias numericas para el texto
###############################################################################

cat("\n  Tabla extra: Tendencias numericas...\n")

trend_summary <- trend_data %>%
  group_by(Species) %>%
  summarise(
    Psi_2017 = Smoothed[Year == 2017],
    Psi_2022 = Smoothed[Year == 2022],
    Change_abs = Psi_2022 - Psi_2017,
    Change_pct = round((Psi_2022 - Psi_2017) / Psi_2017 * 100, 1),
    .groups = "drop"
  )

write.csv(trend_summary, file.path(results_dir, "pub_occupancy_trend_summary.csv"),
          row.names = FALSE)
cat("    -> pub_occupancy_trend_summary.csv\n")
print(trend_summary)

###############################################################################
# RESUMEN
###############################################################################

cat("\n\n")
cat(strrep("=", 70), "\n")
cat("  RESUMEN DE FIGURAS GENERADAS\n")
cat(strrep("=", 70), "\n\n")

figs_generated <- list.files(figs_dir, pattern = "pub_fig.*\\.png$")
cat("  Total figuras:", length(figs_generated), "\n\n")
for (f in sort(figs_generated)) {
  cat("   ", f, "\n")
}

cat("\n")
cat(strrep("=", 70), "\n")
cat("  FIN\n")
cat(strrep("=", 70), "\n\n")
