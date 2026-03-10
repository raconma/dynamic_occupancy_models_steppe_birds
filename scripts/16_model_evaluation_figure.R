#!/usr/bin/env Rscript
# =============================================================================
# 16_model_evaluation_figure.R
#
# Model evaluation figure for GCB publication.
#
# Left column:  Spatial validation — predicted psi_1 maps with
#               III Atlas of Breeding Birds in Spain (SEO/BirdLife 2022)
#               presence polygons overlaid (independent validation).
# Right column: Temporal comparison — model-predicted occupancy rate vs
#               naive eBird detection rate per year (descriptive, not
#               independent validation).
#
# Spatial AUC values from spatially-blocked 5-fold CV (scripts/5_validation.R)
# using the Atlas as independent ground truth.
#
# Inputs:
#   data/processed_2023/{sp}/occ_{sp}_prediction.csv   -- spatial predictions
#   data/processed_2023/{sp}/{sp}_occ_wide_dynamic.csv -- detection histories
#   data-raw/.../aves_spain.shp                        -- III Atlas shapefile
#   results/{sp}_model_object.rds                      -- fitted colext models
#   results/{sp}_validation_cv.csv                     -- atlas CV results
#
# Outputs:
#   figs/pub_fig_model_evaluation.png  -- composite figure (300 DPI)
#   figs/pub_fig_model_evaluation.pdf
#   results/model_evaluation_trends.csv
# =============================================================================

library(unmarked)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(sf)
library(rnaturalearth)

cat("\n", strrep("=", 70), "\n")
cat("  MODEL EVALUATION FIGURE (Atlas + Temporal Trends)\n")
cat(strrep("=", 70), "\n\n")

# --- Configuration ---
sp_codes  <- c("otitar", "ptealc", "pteori", "tettet")
sp_latin  <- c(otitar = "O. tarda", ptealc = "P. alchata",
               pteori = "P. orientalis", tettet = "T. tetrax")
# Atlas column names (uppercase species codes)
sp_atlas_col <- c(otitar = "OTITAR", ptealc = "PTEALC",
                  pteori = "PTEORI", tettet = "TETTET")
years <- 2017:2023
n_years <- length(years)

dir.create("figs", showWarnings = FALSE)

# =============================================================================
# STEP 1: Load geographic data
# =============================================================================
cat("Step 1: Loading geographic data...\n")

# Iberian Peninsula boundary
europe <- ne_countries(scale = 50, continent = "Europe", returnclass = "sf")
iberia <- europe[europe$name %in% c("Spain", "Portugal"), ]

# III Atlas of Breeding Birds in Spain (SEO/BirdLife 2022)
atlas_path <- "data-raw/data/validation/atlas_biodiversidad/atlas_biodiversidad/aves_spain.shp"
atlas <- st_read(atlas_path, quiet = TRUE)
atlas_ll <- st_transform(atlas, 4326)  # to WGS84 for plotting

cat("  Atlas loaded:", nrow(atlas), "polygons\n")
for (sp in sp_codes) {
  n_pres <- sum(atlas[[sp_atlas_col[sp]]] == 1, na.rm = TRUE)
  cat("    ", sp, ":", n_pres, "presence polygons\n")
}
cat("\n")

# =============================================================================
# STEP 2: Load atlas-validated AUC (from 5_validation.R spatially-blocked CV)
# =============================================================================
cat("Step 2: Loading atlas-validated spatial AUC...\n")

atlas_auc <- data.frame()
for (sp in sp_codes) {
  cv_path <- paste0("results/", sp, "_validation_cv.csv")
  if (file.exists(cv_path)) {
    cv <- read.csv(cv_path)
    mean_auc <- mean(cv$AUC)
    sd_auc   <- sd(cv$AUC)
    mean_tss <- mean(cv$TSS)
    cat(sprintf("  %-8s: AUC = %.3f +/- %.3f  (TSS = %.3f)\n",
                sp, mean_auc, sd_auc, mean_tss))
    atlas_auc <- rbind(atlas_auc,
      data.frame(species = sp, atlas_auc = mean_auc,
                 atlas_auc_sd = sd_auc, atlas_tss = mean_tss))
  } else {
    cat("  WARNING:", cv_path, "not found\n")
  }
}
cat("\n")

# =============================================================================
# STEP 3: Compute naive occupancy + model predictions per year
# =============================================================================
cat("Step 3: Computing naive occupancy + model predictions per year...\n")

naive_summary <- data.frame()
predicted_summary <- data.frame()

for (sp in sp_codes) {
  cat("  ", sp, "...")

  # --- Naive occupancy from detection histories ---
  dh <- read.csv(paste0("data/processed_2023/", sp, "/", sp,
                         "_occ_wide_dynamic.csv"),
                 stringsAsFactors = FALSE)

  for (yr in years) {
    site_col <- paste0("site.", yr)
    y_cols <- paste0("y.", 1:10, ".", yr)
    surveyed_idx <- which(!is.na(dh[[site_col]]))

    if (length(surveyed_idx) == 0) {
      naive_summary <- rbind(naive_summary,
        data.frame(species = sp, year = yr, n_surveyed = 0,
                   n_detected = 0, naive_occ = NA))
      next
    }

    y_mat <- as.matrix(dh[surveyed_idx, y_cols])
    naive_det <- apply(y_mat, 1, function(x) {
      as.integer(any(x == 1, na.rm = TRUE))
    })

    naive_summary <- rbind(naive_summary,
      data.frame(species = sp, year = yr,
                 n_surveyed = length(surveyed_idx),
                 n_detected = sum(naive_det),
                 naive_occ  = mean(naive_det)))
  }

  # --- Model-predicted occupancy ---
  mod <- readRDS(paste0("results/", sp, "_model_object.rds"))
  sm <- smoothed(mod)

  predicted_summary <- rbind(predicted_summary,
    data.frame(species = sp, year = years,
               psi_smoothed = as.numeric(sm["occupied", ])))

  cat(" done\n")
}

# Save trends
eval_trends <- merge(naive_summary, predicted_summary,
                     by = c("species", "year"))
write.csv(eval_trends, "results/model_evaluation_trends.csv",
          row.names = FALSE)
cat("  Trends saved to results/model_evaluation_trends.csv\n\n")

# =============================================================================
# STEP 4: Create composite figure (GCB publication quality)
# =============================================================================
cat("Step 4: Creating composite figure...\n")

# --- Colour definitions ---
col_naive <- "#636363"   # dark grey for naive/eBird
col_model <- "#B2182B"   # dark red for model predictions

# --- Left column: Spatial maps with Atlas presence polygons ---
map_panels <- list()

for (sp in sp_codes) {
  # Load prediction raster
  pred_path <- paste0("data/processed_2023/", sp, "/occ_", sp,
                       "_prediction.csv")
  if (!file.exists(pred_path)) {
    cat("  WARNING: prediction file not found for", sp, "\n")
    next
  }
  pred_df <- read.csv(pred_path)
  names(pred_df) <- tolower(names(pred_df))

  # Atlas presence polygons for this species (in WGS84)
  atlas_pres <- atlas_ll[atlas_ll[[sp_atlas_col[sp]]] == 1, ]

  # Atlas-validated AUC for annotation
  sp_auc_row <- atlas_auc[atlas_auc$species == sp, ]
  auc_label <- if (nrow(sp_auc_row) > 0) {
    sprintf("AUC = %.2f", sp_auc_row$atlas_auc)
  } else {
    ""
  }

  # Show colorbar + atlas legend only on last species panel
  is_last <- (sp == sp_codes[length(sp_codes)])

  p_map <- ggplot() +
    geom_raster(data = pred_df,
                aes(x = longitude, y = latitude, fill = occ_prob)) +
    scale_fill_distiller(palette = "BuPu", direction = 1,
                         name = expression(paste("Predicted ", psi[1])),
                         limits = c(0, NA),
                         guide = if (is_last) {
                           guide_colorbar(barwidth = 8, barheight = 0.5,
                                          title.position = "top",
                                          title.hjust = 0.5)
                         } else {
                           "none"
                         }) +
    geom_sf(data = iberia, fill = NA, colour = "grey30", linewidth = 0.3) +
    # Atlas presence polygons as semi-transparent fill
    geom_sf(data = atlas_pres, fill = "#E66101", colour = NA,
            alpha = 0.35) +
    coord_sf(xlim = c(-9.5, 3.5), ylim = c(35.8, 44), expand = FALSE) +
    # AUC annotation
    {if (auc_label != "")
      annotate("label", x = -9, y = 36.3, label = auc_label,
               hjust = 0, size = 2.8, fontface = "bold",
               fill = alpha("white", 0.85),
               label.padding = unit(0.15, "lines"))
    } +
    labs(title = bquote(bold(.(LETTERS[which(sp_codes == sp)])) ~
                          italic(.(sp_latin[sp]))),
         x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(size = 10),
      axis.text = element_text(size = 7),
      panel.grid = element_blank(),
      legend.position = if (is_last) "bottom" else "none",
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      legend.box = "horizontal",
      legend.margin = margin(0, 0, 0, 0),
      plot.margin = margin(2, 4, 2, 2)
    )

  map_panels[[sp]] <- p_map
}

# --- Right column: Temporal trend (naive vs model) ---
trend_panels <- list()

for (sp in sp_codes) {
  ns <- naive_summary[naive_summary$species == sp, ]
  ps <- predicted_summary[predicted_summary$species == sp, ]

  trend_df <- data.frame(
    year = rep(years, 2),
    occ = c(ns$naive_occ * 100, ps$psi_smoothed * 100),
    type = factor(rep(c("Observed (naive eBird)", "Predicted (colext)"),
                      each = n_years),
                  levels = c("Predicted (colext)", "Observed (naive eBird)"))
  )

  # Show legend only for the first species
  show_legend <- (sp == sp_codes[1])

  p_trend <- ggplot(trend_df, aes(x = year, y = occ, colour = type,
                                   shape = type, linetype = type)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 3) +
    scale_colour_manual(
      values = c("Predicted (colext)" = col_model,
                 "Observed (naive eBird)" = col_naive),
      name = NULL) +
    scale_shape_manual(
      values = c("Predicted (colext)" = 16,
                 "Observed (naive eBird)" = 1),
      name = NULL) +
    scale_linetype_manual(
      values = c("Predicted (colext)" = "solid",
                 "Observed (naive eBird)" = "dashed"),
      name = NULL) +
    scale_x_continuous(breaks = years) +
    labs(x = NULL, y = "Occupancy (% cells)") +
    theme_classic(base_size = 10) +
    theme(
      legend.position = if (show_legend) c(0.72, 0.97) else "none",
      legend.text = element_text(size = 8),
      legend.key.size = unit(0.5, "cm"),
      legend.key.width = unit(1.2, "cm"),
      legend.background = element_rect(fill = alpha("white", 0.9),
                                        colour = "grey80",
                                        linewidth = 0.3),
      axis.text = element_text(size = 8),
      axis.title.y = element_text(size = 9),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.margin = margin(2, 8, 2, 2)
    )

  trend_panels[[sp]] <- p_trend
}

# --- Assemble composite figure ---
row_panels <- list()
for (sp in sp_codes) {
  row_panels[[sp]] <- map_panels[[sp]] + trend_panels[[sp]] +
    plot_layout(widths = c(1, 1.2))
}

fig_final <- row_panels[[1]] / row_panels[[2]] /
  row_panels[[3]] / row_panels[[4]] +
  plot_layout(heights = c(1, 1, 1, 1.15))

ggsave("figs/pub_fig_model_evaluation.png", fig_final,
       width = 210, height = 290, units = "mm", dpi = 300)
ggsave("figs/pub_fig_model_evaluation.pdf", fig_final,
       width = 210, height = 290, units = "mm")

cat("\n  -> figs/pub_fig_model_evaluation.png\n")
cat("  -> figs/pub_fig_model_evaluation.pdf\n")

# --- Print summary ---
cat("\n", strrep("-", 50), "\n")
cat("  SUMMARY\n")
cat(strrep("-", 50), "\n\n")

cat("Spatial validation (Atlas III, spatially-blocked 5-fold CV):\n")
for (sp in sp_codes) {
  a <- atlas_auc[atlas_auc$species == sp, ]
  if (nrow(a) > 0) {
    cat(sprintf("  %-8s: AUC = %.3f +/- %.3f  (TSS = %.3f)\n",
                sp, a$atlas_auc, a$atlas_auc_sd, a$atlas_tss))
  }
}

cat("\nTemporal trends (descriptive, same eBird data):\n")
for (sp in sp_codes) {
  ns <- naive_summary[naive_summary$species == sp, ]
  ps <- predicted_summary[predicted_summary$species == sp, ]
  cat(sprintf("  %-8s: Naive %.1f%%->%.1f%%  |  Model %.1f%%->%.1f%%\n",
              sp,
              ns$naive_occ[1] * 100, ns$naive_occ[n_years] * 100,
              ps$psi_smoothed[1] * 100, ps$psi_smoothed[n_years] * 100))
}

cat("\nDone.\n")
