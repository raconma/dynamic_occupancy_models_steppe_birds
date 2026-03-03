###############################################################################
# publication_tables_figures.R
#
# Genera tablas y figuras para publicacion a partir de los resultados
# de los pipelines v4b y v5 para las 4 especies esteparias.
#
# PARTE 1: Tablas comparativas (directamente desde CSVs de resultados)
# PARTE 2: Figuras de tendencia de ocupacion (requiere re-ejecutar modelos)
# PARTE 3: Mapas espaciales (requiere re-ejecutar modelos)
#
# Output:
#   Tablas (CSV listos para publicacion):
#     results/pub_table1_model_comparison.csv
#     results/pub_table2_coefficients_psi.csv
#     results/pub_table3_coefficients_gamma_epsilon.csv
#     results/pub_table4_detection.csv
#     results/pub_table5_dynamic_effects.csv
#
#   Figuras (PNG):
#     figs/pub_fig1_aic_comparison.png
#     figs/pub_fig2_coefficient_forest_plot.png
#     figs/pub_fig3_tree_cover_effect.png
#     figs/pub_fig4_occupancy_trends.png (requiere objetos modelo)
#     figs/pub_fig5_spatial_maps.png (requiere objetos modelo)
###############################################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)

cat("\n")
cat(strrep("=", 70), "\n")
cat("  TABLAS Y FIGURAS PARA PUBLICACION\n")
cat("  4 especies esteparias - Pipelines v4b y v5\n")
cat(strrep("=", 70), "\n\n")

results_dir <- "/Users/gfandos/Documents/GitHub/dynamic_occupancy_models_steppe_birds/results"
figs_dir <- "/Users/gfandos/Documents/GitHub/dynamic_occupancy_models_steppe_birds/figs"
if (!dir.exists(figs_dir)) dir.create(figs_dir, recursive = TRUE)

# Nombres comunes
sp_names <- c(
  otitar = "Otis tarda",
  ptealc = "Pterocles alchata",
  pteori = "Pterocles orientalis",
  tettet = "Tetrax tetrax"
)
sp_common <- c(
  otitar = "Great Bustard",
  ptealc = "Pin-tailed Sandgrouse",
  pteori = "Black-bellied Sandgrouse",
  tettet = "Little Bustard"
)
species <- names(sp_names)

###############################################################################
# PARTE 1: TABLAS
###############################################################################

cat("--- PARTE 1: Tablas para publicacion ---\n\n")

# =========================================================================
# TABLA 1: Comparacion de modelos por especie y pipeline
# =========================================================================

cat("  Tabla 1: Comparacion de modelos...\n")

tab1_rows <- list()

for (sp in species) {
  for (pipeline in c("v4b", "v5")) {
    f <- file.path(results_dir, paste0(sp, "_", pipeline, "_model_selection.csv"))
    if (!file.exists(f)) next
    ms <- read.csv(f)

    # Mejor modelo (menor AIC)
    best <- ms[which.min(ms$AIC), ]
    baseline <- ms[grep("Baseline", ms$Model), ]

    # dAIC vs baseline
    if (nrow(baseline) > 0 && !is.na(baseline$AIC[1])) {
      daic_vs_bl <- round(best$AIC - baseline$AIC[1], 1)
    } else {
      daic_vs_bl <- NA
    }

    # Boundary?
    has_boundary <- grepl("Boundary", best$Status)

    tab1_rows[[length(tab1_rows) + 1]] <- data.frame(
      Species = sp_names[sp],
      Pipeline = pipeline,
      Best_Model = best$Model,
      AIC = round(best$AIC, 1),
      nPars = best$nPars,
      dAIC_vs_Baseline = daic_vs_bl,
      Weight = ifelse("Weight" %in% names(best) && !is.na(best$Weight), best$Weight, NA),
      Boundary = has_boundary,
      n_Models = nrow(ms),
      stringsAsFactors = FALSE
    )
  }
}

tab1 <- do.call(rbind, tab1_rows)
write.csv(tab1, file.path(results_dir, "pub_table1_model_comparison.csv"), row.names = FALSE)
cat("    -> pub_table1_model_comparison.csv\n")
print(tab1)

# =========================================================================
# TABLA 2: Coeficientes de psi (ocupacion inicial) - 4 especies
# =========================================================================

cat("\n  Tabla 2: Coeficientes de psi...\n")

tab2_rows <- list()

for (sp in species) {
  # Usar v4b como referencia (o v5 si v4b no tiene)
  for (pipeline in c("v4b", "v5")) {
    f <- file.path(results_dir, paste0(sp, "_", pipeline, "_best_model_coefficients.csv"))
    if (!file.exists(f)) next
    coefs <- read.csv(f)
    psi_coefs <- coefs[coefs$Submodel == "psi", ]
    if (nrow(psi_coefs) == 0) next

    for (i in 1:nrow(psi_coefs)) {
      tab2_rows[[length(tab2_rows) + 1]] <- data.frame(
        Species = sp_names[sp],
        Pipeline = pipeline,
        Parameter = psi_coefs$Parameter[i],
        Estimate = round(psi_coefs$Estimate[i], 3),
        SE = round(psi_coefs$SE[i], 3),
        z = round(psi_coefs$z[i], 2),
        p_value = signif(psi_coefs$P...z..[i], 3),
        Significance = ifelse(psi_coefs$P...z..[i] < 0.001, "***",
                       ifelse(psi_coefs$P...z..[i] < 0.01, "**",
                       ifelse(psi_coefs$P...z..[i] < 0.05, "*",
                       ifelse(psi_coefs$P...z..[i] < 0.1, ".", "")))),
        stringsAsFactors = FALSE
      )
    }
    break  # Solo usar el primer pipeline disponible
  }
}

tab2 <- do.call(rbind, tab2_rows)
write.csv(tab2, file.path(results_dir, "pub_table2_coefficients_psi.csv"), row.names = FALSE)
cat("    -> pub_table2_coefficients_psi.csv\n")

# =========================================================================
# TABLA 3: Coeficientes de gamma y epsilon - 4 especies (mejor modelo)
# =========================================================================

cat("\n  Tabla 3: Coeficientes de gamma/epsilon...\n")

tab3_rows <- list()

for (sp in species) {
  # Determinar mejor pipeline para esta especie
  best_pipeline <- NULL
  best_aic <- Inf

  for (pipeline in c("v4b", "v5")) {
    f_ms <- file.path(results_dir, paste0(sp, "_", pipeline, "_model_selection.csv"))
    if (!file.exists(f_ms)) next
    ms <- read.csv(f_ms)
    # Mejor modelo valido (OK status o el de menor AIC)
    ok_rows <- grep("^OK", ms$Status)
    if (length(ok_rows) > 0) {
      min_aic_ok <- min(ms$AIC[ok_rows])
      if (min_aic_ok < best_aic) {
        best_aic <- min_aic_ok
        best_pipeline <- pipeline
      }
    }
  }

  if (is.null(best_pipeline)) {
    # Si ningun modelo es OK, usar v4b
    best_pipeline <- "v4b"
  }

  f <- file.path(results_dir, paste0(sp, "_", best_pipeline, "_best_model_coefficients.csv"))
  if (!file.exists(f)) next
  coefs <- read.csv(f)

  for (subm in c("col", "ext")) {
    sm_coefs <- coefs[coefs$Submodel == subm, ]
    if (nrow(sm_coefs) == 0) next

    for (i in 1:nrow(sm_coefs)) {
      tab3_rows[[length(tab3_rows) + 1]] <- data.frame(
        Species = sp_names[sp],
        Best_Pipeline = best_pipeline,
        Submodel = ifelse(subm == "col", "Colonization (gamma)", "Extinction (epsilon)"),
        Parameter = sm_coefs$Parameter[i],
        Estimate = round(sm_coefs$Estimate[i], 3),
        SE = round(sm_coefs$SE[i], 3),
        z = round(sm_coefs$z[i], 2),
        p_value = signif(sm_coefs$P...z..[i], 3),
        Prob = round(sm_coefs$Prob[i], 4),
        Significance = ifelse(sm_coefs$P...z..[i] < 0.001, "***",
                       ifelse(sm_coefs$P...z..[i] < 0.01, "**",
                       ifelse(sm_coefs$P...z..[i] < 0.05, "*",
                       ifelse(sm_coefs$P...z..[i] < 0.1, ".", "")))),
        stringsAsFactors = FALSE
      )
    }
  }
}

tab3 <- do.call(rbind, tab3_rows)
write.csv(tab3, file.path(results_dir, "pub_table3_coefficients_gamma_epsilon.csv"), row.names = FALSE)
cat("    -> pub_table3_coefficients_gamma_epsilon.csv\n")

# =========================================================================
# TABLA 4: Parametros de deteccion - 4 especies
# =========================================================================

cat("\n  Tabla 4: Parametros de deteccion...\n")

tab4_rows <- list()

for (sp in species) {
  for (pipeline in c("v4b", "v5")) {
    f <- file.path(results_dir, paste0(sp, "_", pipeline, "_best_model_coefficients.csv"))
    if (!file.exists(f)) next
    coefs <- read.csv(f)
    det_coefs <- coefs[coefs$Submodel == "det", ]
    if (nrow(det_coefs) == 0) next

    for (i in 1:nrow(det_coefs)) {
      tab4_rows[[length(tab4_rows) + 1]] <- data.frame(
        Species = sp_names[sp],
        Parameter = det_coefs$Parameter[i],
        Estimate = round(det_coefs$Estimate[i], 3),
        SE = round(det_coefs$SE[i], 3),
        p_value = signif(det_coefs$P...z..[i], 3),
        Prob = round(det_coefs$Prob[i], 3),
        stringsAsFactors = FALSE
      )
    }
    break  # Solo primer pipeline
  }
}

tab4 <- do.call(rbind, tab4_rows)
write.csv(tab4, file.path(results_dir, "pub_table4_detection.csv"), row.names = FALSE)
cat("    -> pub_table4_detection.csv\n")

# =========================================================================
# TABLA 5: Resumen de efectos dinamicos por especie
# =========================================================================

cat("\n  Tabla 5: Resumen efectos dinamicos...\n")

tab5 <- data.frame(
  Species = c("Otis tarda", "Pterocles alchata", "Pterocles orientalis", "Tetrax tetrax"),
  Common_Name = c("Great Bustard", "Pin-tailed Sandgrouse", "Black-bellied Sandgrouse", "Little Bustard"),
  v4b_Best = c("gam~lc_crop, eps~tree+lc_crop", "Baseline (boundary)",
               "Baseline (boundary)", "gam~lc_grass+lc_crop"),
  v4b_dAIC = c(-1.3, NA, NA, -48.8),
  v4b_Boundary_Resolved = c("N/A (no boundary)", "No", "No", "Yes"),
  v5_Best = c("eps~tree+pr_lag", "gam~pr_lag", "Baseline (boundary)",
              "m11 (no converge)"),
  v5_dAIC = c(-1.5, -19.4, NA, NA),
  v5_Boundary_Resolved = c("N/A", "Yes", "No", "No"),
  Overall_Best_Pipeline = c("v3 static", "v5 (pr_lag)", "v4b static", "v4b (lc dynamic)"),
  Key_Dynamic_Predictor = c("None (tree_cover static)", "pr_lag (colonization)",
                            "None", "lc_grass + lc_crop (colonization)"),
  tree_cover_eps_p = c(0.0025, 0.073, 0.479, 0.0002),
  Detection_p = c(0.253, 0.250, 0.133, 0.253),
  stringsAsFactors = FALSE
)

write.csv(tab5, file.path(results_dir, "pub_table5_dynamic_effects.csv"), row.names = FALSE)
cat("    -> pub_table5_dynamic_effects.csv\n")
print(tab5[, c("Species", "v4b_dAIC", "v5_dAIC", "Key_Dynamic_Predictor", "tree_cover_eps_p")])

###############################################################################
# PARTE 2: FIGURAS
###############################################################################

cat("\n\n--- PARTE 2: Figuras ---\n\n")

# Tema comun para publicacion
theme_pub <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95"),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 12),
    axis.title = element_text(size = 10)
  )

# =========================================================================
# FIGURA 1: Comparacion AIC por especie (barplot dAIC vs baseline)
# =========================================================================

cat("  Figura 1: Comparacion AIC...\n")

fig1_data <- data.frame(
  Species = rep(c("O. tarda", "P. alchata", "P. orientalis", "T. tetrax"), each = 2),
  Pipeline = rep(c("v4b\n(land cover)", "v5\n(lag-1)"), 4),
  dAIC = c(-1.3, -1.5,   # otitar
           NA, -19.4,     # ptealc
           NA, NA,         # pteori
           -48.8, NA),     # tettet
  stringsAsFactors = FALSE
)
fig1_data$Species <- factor(fig1_data$Species,
                             levels = c("O. tarda", "P. alchata", "P. orientalis", "T. tetrax"))

p1 <- ggplot(fig1_data, aes(x = Species, y = dAIC, fill = Pipeline)) +
  geom_col(position = position_dodge(0.8), width = 0.7, na.rm = TRUE) +
  geom_hline(yintercept = -2, linetype = "dashed", color = "red", alpha = 0.7) +
  annotate("text", x = 4.5, y = -2, label = "dAIC = -2\n(substantial)", size = 3,
           hjust = 1, vjust = -0.3, color = "red") +
  scale_fill_manual(values = c("v4b\n(land cover)" = "#2166AC", "v5\n(lag-1)" = "#B2182B")) +
  labs(title = "Dynamic covariate improvement over static baseline",
       subtitle = "dAIC < -2 indicates substantial improvement (Burnham & Anderson 2002)",
       x = "", y = expression(Delta*"AIC vs static baseline"),
       fill = "Pipeline") +
  coord_flip() +
  theme_pub +
  theme(legend.position = "right")

ggsave(file.path(figs_dir, "pub_fig1_aic_comparison.png"), p1,
       width = 8, height = 5, dpi = 300)
cat("    -> pub_fig1_aic_comparison.png\n")

# =========================================================================
# FIGURA 2: Forest plot de coeficientes clave (tree_cover en psi y eps)
# =========================================================================

cat("  Figura 2: Forest plot de tree_cover...\n")

fig2_data <- data.frame(
  Species = c("O. tarda", "P. alchata", "P. orientalis", "T. tetrax",
              "O. tarda", "P. alchata", "P. orientalis", "T. tetrax"),
  Submodel = c(rep("psi (initial occupancy)", 4), rep("epsilon (extinction)", 4)),
  Estimate = c(-11.44, -3.32, -4.82, -3.29,   # psi
               2.65, 1.50, 0.72, 3.06),          # epsilon
  SE = c(1.85, 0.59, 0.64, 0.69,               # psi
         0.88, 0.84, 1.01, 0.82),                # epsilon
  p_value = c(1e-10, 1e-8, 6e-14, 2e-6,        # psi
              0.0025, 0.073, 0.479, 0.0002),      # epsilon
  stringsAsFactors = FALSE
)

fig2_data$lower <- fig2_data$Estimate - 1.96 * fig2_data$SE
fig2_data$upper <- fig2_data$Estimate + 1.96 * fig2_data$SE
fig2_data$Significant <- ifelse(fig2_data$p_value < 0.05, "p < 0.05", "p >= 0.05")
fig2_data$Species <- factor(fig2_data$Species,
                             levels = rev(c("O. tarda", "P. alchata", "P. orientalis", "T. tetrax")))

p2 <- ggplot(fig2_data, aes(x = Estimate, y = Species, color = Significant)) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, linewidth = 0.8) +
  geom_point(size = 3) +
  scale_color_manual(values = c("p < 0.05" = "#D32F2F", "p >= 0.05" = "#9E9E9E")) +
  facet_wrap(~Submodel, scales = "free_x") +
  labs(title = "Tree cover effect on occupancy dynamics",
       subtitle = "Logit-scale coefficients with 95% CI",
       x = "Coefficient estimate (logit scale)", y = "",
       color = "") +
  theme_pub +
  theme(legend.position = "bottom")

ggsave(file.path(figs_dir, "pub_fig2_coefficient_forest_plot.png"), p2,
       width = 9, height = 5, dpi = 300)
cat("    -> pub_fig2_coefficient_forest_plot.png\n")

# =========================================================================
# FIGURA 3: Efecto de tree_cover en extincion (curva de respuesta)
# =========================================================================

cat("  Figura 3: Curva de respuesta tree_cover en extincion...\n")

# Datos de la curva de respuesta (usando coeficientes del mejor modelo)
tree_seq <- seq(-2, 3, length.out = 100)  # Escala estandarizada

curves <- data.frame()
sp_coefs <- list(
  "O. tarda" = c(int = 0.098, slope = 2.723),  # s4 (sensibilidad)
  "P. alchata" = c(int = -1.735, slope = 1.159),
  "P. orientalis" = c(int = -1.738, slope = 0.716),
  "T. tetrax" = c(int = -0.071, slope = 3.055)
)

for (sp_name in names(sp_coefs)) {
  coefs <- sp_coefs[[sp_name]]
  logit_eps <- coefs["int"] + coefs["slope"] * tree_seq
  prob_eps <- plogis(logit_eps)
  curves <- rbind(curves, data.frame(
    tree_cover_scaled = tree_seq,
    extinction_prob = prob_eps,
    Species = sp_name,
    stringsAsFactors = FALSE
  ))
}

curves$Species <- factor(curves$Species,
                          levels = c("O. tarda", "P. alchata", "P. orientalis", "T. tetrax"))

p3 <- ggplot(curves, aes(x = tree_cover_scaled, y = extinction_prob, color = Species)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("O. tarda" = "#1B9E77", "P. alchata" = "#D95F02",
                                 "P. orientalis" = "#7570B3", "T. tetrax" = "#E7298A")) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  labs(title = "Local extinction probability as a function of tree cover",
       subtitle = "Response curves from best models (logistic link function)",
       x = "Tree cover (standardized)", y = "Extinction probability",
       color = "Species") +
  theme_pub

ggsave(file.path(figs_dir, "pub_fig3_tree_cover_effect.png"), p3,
       width = 8, height = 5, dpi = 300)
cat("    -> pub_fig3_tree_cover_effect.png\n")

# =========================================================================
# FIGURA 4: Comparacion de detectabilidad entre especies
# =========================================================================

cat("  Figura 4: Detectabilidad entre especies...\n")

det_data <- data.frame(
  Species = factor(c("O. tarda", "P. alchata", "P. orientalis", "T. tetrax"),
                    levels = c("O. tarda", "P. alchata", "P. orientalis", "T. tetrax")),
  Detection = c(0.253, 0.250, 0.133, 0.253),
  stringsAsFactors = FALSE
)

p4 <- ggplot(det_data, aes(x = Species, y = Detection, fill = Species)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%", Detection * 100)), vjust = -0.5, size = 4) +
  scale_fill_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A")) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 0.35)) +
  labs(title = "Baseline detection probability by species",
       subtitle = "P. orientalis has notably lower detectability",
       x = "", y = "Detection probability") +
  theme_pub +
  theme(legend.position = "none")

ggsave(file.path(figs_dir, "pub_fig4_detection_comparison.png"), p4,
       width = 7, height = 5, dpi = 300)
cat("    -> pub_fig4_detection_comparison.png\n")

# =========================================================================
# FIGURA 5: Panel de coeficientes psi (4 especies)
# =========================================================================

cat("  Figura 5: Panel coeficientes psi...\n")

psi_data <- tab2[tab2$Parameter != "(Intercept)", ]
psi_data$lower <- psi_data$Estimate - 1.96 * psi_data$SE
psi_data$upper <- psi_data$Estimate + 1.96 * psi_data$SE
psi_data$Significant <- ifelse(psi_data$p_value < 0.05, "p < 0.05", "p >= 0.05")

p5 <- ggplot(psi_data, aes(x = Estimate, y = Parameter, color = Significant)) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
  geom_point(size = 2.5) +
  scale_color_manual(values = c("p < 0.05" = "#D32F2F", "p >= 0.05" = "#9E9E9E")) +
  facet_wrap(~Species, scales = "free_x", ncol = 2) +
  labs(title = "Initial occupancy (psi) coefficients across species",
       x = "Estimate (logit scale)", y = "", color = "") +
  theme_pub +
  theme(legend.position = "bottom")

ggsave(file.path(figs_dir, "pub_fig5_psi_coefficients.png"), p5,
       width = 10, height = 7, dpi = 300)
cat("    -> pub_fig5_psi_coefficients.png\n")

###############################################################################
# PARTE 3: FIGURAS QUE REQUIEREN OBJETOS MODELO
###############################################################################

cat("\n\n--- PARTE 3: Figuras que requieren re-ejecutar modelos ---\n")
cat("  NOTA: Las siguientes figuras requieren guardar los objetos modelo\n")
cat("  con saveRDS() en los scripts de cada pipeline.\n")
cat("  Anade al final de cada pipeline:\n")
cat("    saveRDS(best_mod, 'results/{sp}_{pipeline}_best_model.rds')\n")
cat("    saveRDS(occ_umf, 'results/{sp}_{pipeline}_umf.rds')\n\n")

# Verificar si existen objetos modelo guardados
model_files <- list.files(results_dir, pattern = "_best_model\\.rds$")
if (length(model_files) > 0) {
  cat("  Modelos encontrados: ", paste(model_files, collapse = ", "), "\n")
  cat("  TODO: Implementar figuras de tendencia y mapas espaciales\n")
} else {
  cat("  No se encontraron objetos modelo (.rds)\n")
  cat("  Para generar figuras de tendencia y mapas, ejecutar:\n\n")
  for (sp in species) {
    for (pipeline in c("v4b", "v5")) {
      cat(sprintf("    # %s %s:\n", sp, pipeline))
      cat(sprintf("    # Anade saveRDS(best_mod, 'results/%s_%s_best_model.rds') al final del pipeline\n", sp, pipeline))
    }
  }
}

cat("\n  Figuras pendientes (requieren objetos modelo):\n")
cat("    - pub_fig_occupancy_trends.png: Tendencia de ocupacion 2017-2022 (projected())\n")
cat("    - pub_fig_colonization_extinction_rates.png: Tasas gamma/epsilon por ano\n")
cat("    - pub_fig_spatial_maps.png: Mapas de psi/epsilon para las 4 especies\n")
cat("    - pub_fig_response_curves.png: Curvas de respuesta completas\n")

###############################################################################
# RESUMEN
###############################################################################

cat("\n\n")
cat(strrep("=", 70), "\n")
cat("  RESUMEN DE OUTPUTS GENERADOS\n")
cat(strrep("=", 70), "\n\n")

cat("  TABLAS:\n")
cat("    1. pub_table1_model_comparison.csv - Comparacion de modelos (4 sp x 2 pipelines)\n")
cat("    2. pub_table2_coefficients_psi.csv - Coeficientes psi (4 sp)\n")
cat("    3. pub_table3_coefficients_gamma_epsilon.csv - Coeficientes gamma/epsilon (4 sp)\n")
cat("    4. pub_table4_detection.csv - Parametros de deteccion (4 sp)\n")
cat("    5. pub_table5_dynamic_effects.csv - Resumen efectos dinamicos (4 sp)\n")

cat("\n  FIGURAS:\n")
cat("    1. pub_fig1_aic_comparison.png - Barplot dAIC vs baseline\n")
cat("    2. pub_fig2_coefficient_forest_plot.png - Forest plot tree_cover\n")
cat("    3. pub_fig3_tree_cover_effect.png - Curva respuesta extincion\n")
cat("    4. pub_fig4_detection_comparison.png - Detectabilidad por especie\n")
cat("    5. pub_fig5_psi_coefficients.png - Panel coeficientes psi\n")

cat("\n")
cat(strrep("=", 70), "\n")
cat("  FIN\n")
cat(strrep("=", 70), "\n\n")
