###############################################################################
# fig_forest_coefs.R
#
# Publication figure: Forest plot of gamma and epsilon standardised
# coefficients for 4 Iberian steppe bird species.
#
# Two-panel layout (left = gamma, right = epsilon).
# Species distinguished by point shape; non-significant in grey.
# Driver categories colour-coded: Climate, Land use, NDVI.
# P. alchata gamma panel excluded (complete separation).
#
# Output:
#   figs/pub_fig_forest_gamma_epsilon.png  (300 DPI)
#   figs/pub_fig_forest_gamma_epsilon.pdf
###############################################################################

library(unmarked)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

cat("\n", strrep("=", 70), "\n")
cat("  FIGURE: Forest plot of gamma/epsilon coefficients\n")
cat(strrep("=", 70), "\n\n")

# --- Paths ---
res_dir  <- "results"
figs_dir <- "figs"
if (!dir.exists(figs_dir)) dir.create(figs_dir, recursive = TRUE)

# --- Species info ---
sp_codes <- c("otitar", "ptealc", "pteori", "tettet")
sp_short <- c(otitar = "O. tarda", ptealc = "P. alchata",
              pteori = "P. orientalis", tettet = "T. tetrax")

# Wong palette for species
sp_colours <- c("O. tarda" = "#D55E00", "P. alchata" = "#0072B2",
                "P. orientalis" = "#009E73", "T. tetrax" = "#CC79A7")
sp_shapes  <- c("O. tarda" = 16, "P. alchata" = 17,
                "P. orientalis" = 15, "T. tetrax" = 18)

# Short labels for covariates
short_labels <- c(
  "Land_Cover_Type_1_Percent_Class_6"  = "LC6",
  "Land_Cover_Type_1_Percent_Class_7"  = "LC7",
  "Land_Cover_Type_1_Percent_Class_12" = "LC12",
  "Land_Cover_Type_1_Percent_Class_13" = "LC13",
  "NDVI" = "NDVI",
  "pr"   = "Precipitation",
  "tmmn" = "Min. temperature",
  "tmmx" = "Max. temperature"
)

# Driver categories
driver_cat <- c(
  "pr" = "Climate", "tmmn" = "Climate", "tmmx" = "Climate",
  "Land_Cover_Type_1_Percent_Class_6"  = "Land use",
  "Land_Cover_Type_1_Percent_Class_7"  = "Land use",
  "Land_Cover_Type_1_Percent_Class_12" = "Land use",
  "Land_Cover_Type_1_Percent_Class_13" = "Land use",
  "NDVI" = "NDVI (climate-adjacent)"
)
cat_colours <- c("Climate" = "#0072B2", "Land use" = "#D55E00",
                 "NDVI (climate-adjacent)" = "#009E73")

# --- Extract coefficients ---
all_coefs <- data.frame()

for (sp in sp_codes) {
  mod_path <- file.path(res_dir, paste0(sp, "_model_object.rds"))
  if (!file.exists(mod_path)) {
    warning("Model not found: ", mod_path)
    next
  }

  mod <- readRDS(mod_path)
  cc  <- coef(mod)
  se  <- sqrt(diag(vcov(mod)))
  names(se) <- names(cc)

  for (prefix in c("col", "ext")) {
    submodel <- ifelse(prefix == "col", "gamma", "epsilon")
    # Get non-intercept coefficient names
    pattern <- paste0("^", prefix, "\\((?!Int)")
    idx <- grepl(pattern, names(cc), perl = TRUE)
    if (!any(idx)) next

    for (nm in names(cc)[idx]) {
      var_raw <- sub(paste0("^", prefix, "\\((.+)\\)$"), "\\1", nm)
      est <- cc[nm]
      s   <- se[nm]

      short_lab <- if (var_raw %in% names(short_labels)) short_labels[var_raw] else var_raw
      cat_lab   <- if (var_raw %in% names(driver_cat))   driver_cat[var_raw]   else "Other"

      all_coefs <- rbind(all_coefs, data.frame(
        species   = sp_short[sp],
        submodel  = submodel,
        covariate = short_lab,
        cov_raw   = var_raw,
        estimate  = as.numeric(est),
        se        = as.numeric(s),
        lower     = as.numeric(est - 1.96 * s),
        upper     = as.numeric(est + 1.96 * s),
        significant = abs(est / s) > 1.96,
        driver_category = cat_lab,
        stringsAsFactors = FALSE
      ))
    }
  }
  cat("  ", sp, "extracted\n")
}

# --- Remove P. alchata gamma (separation) ---
ptealc_gam <- all_coefs %>% filter(species == "P. alchata" & submodel == "gamma")
all_coefs  <- all_coefs %>% filter(!(species == "P. alchata" & submodel == "gamma"))

cat("\n  P. alchata gamma coefficients excluded (separation):\n")
if (nrow(ptealc_gam) > 0) print(ptealc_gam[, c("covariate", "estimate", "se")])

# Also flag P. alchata epsilon coefficients with huge SE as non-significant
# (separation artefact — |estimate/SE| can be ~1 despite meaningless values)
all_coefs$significant[all_coefs$species == "P. alchata" &
                       all_coefs$submodel == "epsilon" &
                       all_coefs$se > 5] <- FALSE

# Factor levels
all_coefs$species <- factor(all_coefs$species,
                             levels = c("O. tarda", "P. alchata",
                                        "P. orientalis", "T. tetrax"))

# Colour mapping: significant uses driver_category colour; non-sig is grey
all_coefs$point_colour <- ifelse(all_coefs$significant,
                                  cat_colours[all_coefs$driver_category],
                                  "grey60")

# --- Build gamma panel ---
df_gam <- all_coefs %>% filter(submodel == "gamma")
df_eps <- all_coefs %>% filter(submodel == "epsilon")

make_forest_panel <- function(df, panel_title, annotate_text = NULL) {
  p <- ggplot(df, aes(x = estimate, y = covariate, shape = species)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_errorbarh(aes(xmin = lower, xmax = upper),
                   height = 0.2, colour = "grey40", linewidth = 0.4) +
    geom_point(aes(colour = significant, fill = driver_category),
               size = 3, stroke = 0.5) +
    scale_shape_manual(values = c("O. tarda" = 16, "P. alchata" = 17,
                                   "P. orientalis" = 15, "T. tetrax" = 18),
                        name = "Species") +
    scale_colour_manual(values = c("TRUE" = "black", "FALSE" = "grey60"),
                         guide = "none") +
    scale_fill_manual(values = cat_colours, name = "Driver category",
                       guide = guide_legend(override.aes = list(shape = 22, size = 4))) +
    labs(title = panel_title, x = "Standardised coefficient", y = NULL) +
    theme_classic(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      axis.text = element_text(size = 9),
      axis.title = element_text(size = 10),
      legend.position = "none"
    )

  if (!is.null(annotate_text)) {
    p <- p + annotate("text", x = Inf, y = -Inf, label = annotate_text,
                       hjust = 1.1, vjust = -0.5, size = 3, fontface = "italic",
                       colour = "grey40")
  }
  p
}

p_gam <- make_forest_panel(df_gam,
                            expression(bold("Colonisation (" * gamma * ")")),
                            "P. alchata excluded (separation)")
p_eps <- make_forest_panel(df_eps,
                            expression(bold("Extinction (" * epsilon * ")")))

# --- Combined figure with shared legend at bottom ---
p_combined <- (p_gam | p_eps) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom",
        legend.box = "horizontal")

# Override to show both shape and fill legends
p_combined <- p_gam + p_eps +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# Re-add legends properly
p_gam2 <- p_gam + theme(legend.position = "bottom")
p_eps2 <- p_eps + theme(legend.position = "bottom")

p_final <- p_gam2 + p_eps2 + plot_layout(guides = "collect")

# --- Save ---
ggsave(file.path(figs_dir, "pub_fig_forest_gamma_epsilon.png"), p_final,
       width = 220, height = 130, units = "mm", dpi = 300)
ggsave(file.path(figs_dir, "pub_fig_forest_gamma_epsilon.pdf"), p_final,
       width = 220, height = 130, units = "mm")

cat("\n  -> pub_fig_forest_gamma_epsilon.png\n")
cat("  -> pub_fig_forest_gamma_epsilon.pdf\n")
cat("\nDone.\n")
