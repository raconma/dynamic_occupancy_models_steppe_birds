###############################################################################
# publication_response_curves.R
#
# Curvas de respuesta con intervalo de confianza al 95% para TODAS las
# variables significativas de cada submodelo (psi, gamma, epsilon, det)
# y las 4 especies esteparias.
#
# Metodo: Delta method sobre escala logit
#   logit(y) = beta0 + beta_j * x_j   (resto covariables fijadas en 0 = media)
#   var(logit) = V[0,0] + x^2*V[j,j] + 2*x*V[0,j]
#   CI = plogis(logit +/- 1.96 * sqrt(var))
#
# Requiere:
#   results/results2/{sp}_{pipeline}_all_models.rds
#
# Output (en figs/):
#   pub_fig_response_psi_all.png       - Curvas psi (todas las variables)
#   pub_fig_response_eps_all.png       - Curvas epsilon
#   pub_fig_response_gam_all.png       - Curvas gamma
#   pub_fig_response_det_all.png       - Curvas deteccion
#   pub_fig_response_panel_main.png    - Panel resumen para publicacion
###############################################################################

library(unmarked)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

cat("\n")
cat(strrep("=", 70), "\n")
cat("  CURVAS DE RESPUESTA CON 95% CI\n")
cat("  4 especies - Todos los submodelos\n")
cat(strrep("=", 70), "\n\n")

# --- Paths ---
rds_dir  <- "results/results2"
figs_dir <- "figs"
if (!dir.exists(figs_dir)) dir.create(figs_dir, recursive = TRUE)

# --- Species info ---
sp_codes <- c("otitar", "ptealc", "pteori", "tettet")
sp_latin <- c(otitar = "Otis tarda", ptealc = "Pterocles alchata",
              pteori = "Pterocles orientalis", tettet = "Tetrax tetrax")
sp_short <- c(otitar = "O. tarda", ptealc = "P. alchata",
              pteori = "P. orientalis", tettet = "T. tetrax")
sp_colors <- c("O. tarda" = "#1B9E77", "P. alchata" = "#D95F02",
               "P. orientalis" = "#7570B3", "T. tetrax" = "#E7298A")

# Pretty variable names for plots
var_labels <- c(
  tree_cover = "Tree cover", grass_cover = "Grass cover",
  bio1 = "Mean temperature (BIO1)", bio2 = "Diurnal range (BIO2)",
  topo_elev = "Elevation", topo_aspect = "Aspect",
  lc_grass = "Grassland cover (dynamic)", lc_crop = "Cropland cover (dynamic)",
  lc_shrub = "Shrubland cover (dynamic)",
  pr_lag = "Precipitation (lag-1)", NDVI_lag = "NDVI (lag-1)",
  tmmx_lag = "Max temperature (lag-1)", lc_grass_lag = "Grassland (lag-1)",
  effort = "Search effort", duration = "Duration",
  observers = "Number of observers", time = "Time of day"
)

# --- Best model per species ---
best_info <- list(
  otitar = list(pipe = "v4b", name = "m9: ambos lc_crop"),
  ptealc = list(pipe = "v5",  name = "m6: gam~pr_lag"),
  pteori = list(pipe = "v4b", name = "m0: Baseline estatico"),
  tettet = list(pipe = "v4b", name = "m7: gam~grass+crop")
)

# --- Load models ---
cat("Cargando modelos...\n")
best_mods <- list()
for (sp in sp_codes) {
  info <- best_info[[sp]]
  mods <- readRDS(file.path(rds_dir, paste0(sp, "_", info$pipe, "_all_models.rds")))
  best_mods[[sp]] <- mods[[info$name]]
  cat("  ", sp, ":", info$name, "\n")
}

# --- Publication theme ---
theme_pub <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(size = 9),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    axis.title = element_text(size = 10)
  )

###############################################################################
# FUNCION: Generar curva de respuesta con CI para una variable
###############################################################################

make_response_curve <- function(mod, submodel, var_name, sp_label,
                                 x_seq = seq(-2.5, 3, length.out = 200)) {
  # submodel: "psi", "col", "ext", "det"
  # var_name: nombre limpio (sin wrapper), e.g. "tree_cover"

  # Prefix mapping
  prefix <- switch(submodel,
                   psi = "psi", col = "col", ext = "ext", det = "p")

  # Get coefficients and vcov
  all_coefs <- coef(mod)
  V <- vcov(mod)

  int_name <- paste0(prefix, "(Int)")
  var_full <- paste0(prefix, "(", var_name, ")")

  # Check variable exists
  if (!(var_full %in% names(all_coefs))) return(NULL)

  beta0 <- all_coefs[int_name]
  beta_j <- all_coefs[var_full]

  # Variance components from full vcov
  v00 <- V[int_name, int_name]
  vjj <- V[var_full, var_full]
  v0j <- V[int_name, var_full]

  # Predictions on logit scale (other covariates = 0 = mean)
  logit_pred <- beta0 + beta_j * x_seq
  var_logit  <- v00 + x_seq^2 * vjj + 2 * x_seq * v0j
  se_logit   <- sqrt(pmax(var_logit, 0))  # pmax to avoid numerical negatives

  # Back-transform to probability scale
  pred  <- plogis(logit_pred)
  lower <- plogis(logit_pred - 1.96 * se_logit)
  upper <- plogis(logit_pred + 1.96 * se_logit)

  data.frame(
    x = x_seq,
    pred = pred,
    lower = lower,
    upper = upper,
    Species = sp_label,
    Variable = ifelse(var_name %in% names(var_labels),
                      var_labels[var_name], var_name),
    Var_raw = var_name,
    Submodel = submodel,
    stringsAsFactors = FALSE
  )
}

###############################################################################
# GENERAR TODAS LAS CURVAS
###############################################################################

cat("\nGenerando curvas de respuesta...\n")

all_curves <- data.frame()
x_seq <- seq(-2.5, 3, length.out = 200)

for (sp in sp_codes) {
  mod <- best_mods[[sp]]
  all_coefs <- coef(mod)

  # For each submodel, extract variable names
  for (submodel in c("psi", "col", "ext", "det")) {
    prefix <- switch(submodel, psi = "psi", col = "col", ext = "ext", det = "p")
    # Get all coefs for this submodel (excluding intercept)
    pattern <- paste0("^", prefix, "\\((?!Int)")
    sub_vars <- names(all_coefs)[grepl(pattern, names(all_coefs), perl = TRUE)]

    if (length(sub_vars) == 0) next

    for (v_full in sub_vars) {
      # Extract clean variable name
      v_clean <- sub(paste0("^", prefix, "\\((.+)\\)$"), "\\1", v_full)
      curve <- make_response_curve(mod, submodel, v_clean, sp_short[sp], x_seq)
      if (!is.null(curve)) {
        all_curves <- rbind(all_curves, curve)
      }
    }
  }
  cat("  ", sp, "- done\n")
}

# Factor species
all_curves$Species <- factor(all_curves$Species,
                              levels = c("O. tarda", "P. alchata",
                                         "P. orientalis", "T. tetrax"))

cat("  Total curvas:", nrow(all_curves) / 200, "\n")

###############################################################################
# FIGURA 1: Curvas psi (ocupacion inicial) - todas las variables
###############################################################################

cat("\n  Figura 1: Curvas psi...\n")

psi_curves <- all_curves %>% filter(Submodel == "psi")

p_psi <- ggplot(psi_curves, aes(x = x, y = pred, color = Species, fill = Species)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = sp_colors) +
  scale_fill_manual(values = sp_colors) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  facet_wrap(~Variable, scales = "free_x", ncol = 3) +
  labs(title = "Initial occupancy (psi) response curves",
       subtitle = "95% CI via delta method; other covariates held at mean (standardized = 0)",
       x = "Covariate value (standardized)", y = expression(psi),
       color = "Species", fill = "Species") +
  theme_pub

ggsave(file.path(figs_dir, "pub_fig_response_psi_all.png"), p_psi,
       width = 12, height = 8, dpi = 300)
cat("    -> pub_fig_response_psi_all.png\n")

###############################################################################
# FIGURA 2: Curvas epsilon (extincion) - todas las variables
###############################################################################

cat("\n  Figura 2: Curvas epsilon...\n")

eps_curves <- all_curves %>% filter(Submodel == "ext")

p_eps <- ggplot(eps_curves, aes(x = x, y = pred, color = Species, fill = Species)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 0.5, linetype = "dotted", alpha = 0.3) +
  scale_color_manual(values = sp_colors) +
  scale_fill_manual(values = sp_colors) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  facet_wrap(~Variable, scales = "free_x", ncol = 3) +
  labs(title = "Local extinction (epsilon) response curves",
       subtitle = "95% CI via delta method; other covariates held at mean",
       x = "Covariate value (standardized)", y = expression(epsilon),
       color = "Species", fill = "Species") +
  theme_pub

ggsave(file.path(figs_dir, "pub_fig_response_eps_all.png"), p_eps,
       width = 12, height = 6, dpi = 300)
cat("    -> pub_fig_response_eps_all.png\n")

###############################################################################
# FIGURA 3: Curvas gamma (colonizacion) - todas las variables
###############################################################################

cat("\n  Figura 3: Curvas gamma...\n")

gam_curves <- all_curves %>% filter(Submodel == "col")

if (nrow(gam_curves) > 0) {
  p_gam <- ggplot(gam_curves, aes(x = x, y = pred, color = Species, fill = Species)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, color = NA) +
    geom_line(linewidth = 0.9) +
    scale_color_manual(values = sp_colors) +
    scale_fill_manual(values = sp_colors) +
    scale_y_continuous(labels = scales::percent_format(), limits = c(0, NA)) +
    facet_wrap(~Variable, scales = "free", ncol = 3) +
    labs(title = "Colonization (gamma) response curves",
         subtitle = "95% CI via delta method; P. orientalis excluded (intercept-only colonization)",
         x = "Covariate value (standardized)", y = expression(gamma),
         color = "Species", fill = "Species") +
    theme_pub

  ggsave(file.path(figs_dir, "pub_fig_response_gam_all.png"), p_gam,
         width = 12, height = 6, dpi = 300)
  cat("    -> pub_fig_response_gam_all.png\n")
} else {
  cat("    (sin variables de colonizacion)\n")
}

###############################################################################
# FIGURA 4: Curvas deteccion - todas las variables
###############################################################################

cat("\n  Figura 4: Curvas deteccion...\n")

det_curves <- all_curves %>% filter(Submodel == "det")

p_det <- ggplot(det_curves, aes(x = x, y = pred, color = Species, fill = Species)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = sp_colors) +
  scale_fill_manual(values = sp_colors) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 0.8)) +
  facet_wrap(~Variable, scales = "free_x", ncol = 2) +
  labs(title = "Detection probability response curves",
       subtitle = "95% CI via delta method; other detection covariates held at mean",
       x = "Covariate value (standardized)", y = "p (detection)",
       color = "Species", fill = "Species") +
  theme_pub

ggsave(file.path(figs_dir, "pub_fig_response_det_all.png"), p_det,
       width = 10, height = 8, dpi = 300)
cat("    -> pub_fig_response_det_all.png\n")

###############################################################################
# FIGURA 5: Panel principal para publicacion
#   A) psi ~ tree_cover  B) psi ~ grass_cover
#   C) epsilon ~ tree_cover  D) gamma ~ dynamic covariates
###############################################################################

cat("\n  Figura 5: Panel principal publicacion...\n")

# Helper to make a focused panel
make_panel <- function(data, title_text, ylab_expr, ylimits = c(0, 1)) {
  ggplot(data, aes(x = x, y = pred, color = Species, fill = Species)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.12, color = NA) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = sp_colors) +
    scale_fill_manual(values = sp_colors) +
    scale_y_continuous(labels = scales::percent_format(), limits = ylimits) +
    labs(title = title_text, x = "Covariate value (standardized)",
         y = ylab_expr, color = "Species", fill = "Species") +
    theme_pub +
    theme(legend.position = "none",
          plot.title = element_text(size = 10))
}

# A) psi ~ tree_cover
pA <- make_panel(
  psi_curves %>% filter(Var_raw == "tree_cover"),
  "A) Initial occupancy ~ Tree cover",
  expression(psi)
)

# B) psi ~ grass_cover
pB <- make_panel(
  psi_curves %>% filter(Var_raw == "grass_cover"),
  "B) Initial occupancy ~ Grass cover",
  expression(psi)
)

# C) epsilon ~ tree_cover
pC <- make_panel(
  eps_curves %>% filter(Var_raw == "tree_cover"),
  "C) Extinction ~ Tree cover",
  expression(epsilon)
) + geom_hline(yintercept = 0.5, linetype = "dotted", alpha = 0.3)

# D) epsilon ~ grass_cover
pD <- make_panel(
  eps_curves %>% filter(Var_raw == "grass_cover"),
  "D) Extinction ~ Grass cover",
  expression(epsilon)
) + geom_hline(yintercept = 0.5, linetype = "dotted", alpha = 0.3)

# E) gamma ~ key dynamic covariates (faceted)
gam_key <- gam_curves %>%
  filter(Var_raw %in% c("lc_crop", "lc_grass", "pr_lag"))

if (nrow(gam_key) > 0) {
  pE <- ggplot(gam_key, aes(x = x, y = pred, color = Species, fill = Species)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.12, color = NA) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = sp_colors) +
    scale_fill_manual(values = sp_colors) +
    scale_y_continuous(labels = scales::percent_format(), limits = c(0, NA)) +
    facet_wrap(~Variable, ncol = 3) +
    labs(title = "E) Colonization ~ Dynamic covariates",
         x = "Covariate value (standardized)",
         y = expression(gamma),
         color = "Species", fill = "Species") +
    theme_pub +
    theme(plot.title = element_text(size = 10),
          legend.position = "bottom")
} else {
  pE <- grid::nullGrob()
}

# F) detection ~ effort
pF <- make_panel(
  det_curves %>% filter(Var_raw == "effort"),
  "F) Detection ~ Search effort",
  "p (detection)",
  ylimits = c(0, 0.7)
)

# Add legend back to individual panels for assembly
pA <- pA + theme(legend.position = "right")
pF <- pF + theme(legend.position = "right")

# Simple assembly with gridExtra
p_main <- grid.arrange(
  pA, pB,
  pC, pD,
  pE,
  pF, grid::nullGrob(),
  layout_matrix = rbind(
    c(1, 2),
    c(3, 4),
    c(5, 5),
    c(6, 7)
  ),
  heights = c(1, 1, 1.1, 0.8)
)

ggsave(file.path(figs_dir, "pub_fig_response_panel_main.png"), p_main,
       width = 11, height = 16, dpi = 300)
cat("    -> pub_fig_response_panel_main.png\n")

###############################################################################
# FIGURA 6: Curva individual por especie (una pagina por especie)
###############################################################################

cat("\n  Figura 6: Curvas por especie individual...\n")

for (sp in sp_codes) {
  sp_data <- all_curves %>% filter(Species == sp_short[sp])
  if (nrow(sp_data) == 0) next

  sp_data$Submodel_label <- recode(sp_data$Submodel,
                                    psi = "psi (initial occupancy)",
                                    col = "gamma (colonization)",
                                    ext = "epsilon (extinction)",
                                    det = "p (detection)")

  p_sp <- ggplot(sp_data, aes(x = x, y = pred)) +
    geom_ribbon(aes(ymin = lower, ymax = upper),
                fill = sp_colors[sp_short[sp]], alpha = 0.2) +
    geom_line(color = sp_colors[sp_short[sp]], linewidth = 1) +
    geom_hline(yintercept = 0.5, linetype = "dotted", alpha = 0.2) +
    scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
    facet_grid(Submodel_label ~ Variable, scales = "free") +
    labs(title = paste0(sp_latin[sp], " (", best_info[[sp]]$pipe,
                        ": ", best_info[[sp]]$name, ")"),
         subtitle = "Response curves with 95% CI; other covariates at mean values",
         x = "Covariate value (standardized)", y = "Probability") +
    theme_pub +
    theme(strip.text = element_text(size = 8))

  fname <- paste0("pub_fig_response_", sp, ".png")
  ggsave(file.path(figs_dir, fname), p_sp,
         width = 14, height = 10, dpi = 300)
  cat("    ->", fname, "\n")
}

###############################################################################
# RESUMEN
###############################################################################

cat("\n\n")
cat(strrep("=", 70), "\n")
cat("  RESUMEN\n")
cat(strrep("=", 70), "\n\n")

fig_files <- list.files(figs_dir, pattern = "pub_fig_response.*\\.png$")
cat("  Figuras de curvas de respuesta:", length(fig_files), "\n\n")
for (f in sort(fig_files)) {
  cat("   ", f, "\n")
}

cat("\n  Todas las curvas incluyen 95% CI (delta method)\n")
cat("  Covariables no mostradas fijadas en su media (0 estandarizado)\n")
cat("\n")
cat(strrep("=", 70), "\n")
cat("  FIN\n")
cat(strrep("=", 70), "\n\n")
