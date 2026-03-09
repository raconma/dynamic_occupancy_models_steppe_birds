###############################################################################
# fig_response_curves.R
#
# Publication figure: Response curves for 6 key driver x process combinations
# using the delta method for 95% CI.
#
# Panels:
#   1. T. tetrax:  LC12 on gamma (solid) AND epsilon (dashed), same axes
#   2. O. tarda:   LC13 on epsilon
#   3. P. orientalis: NDVI on epsilon
#   4. T. tetrax:  LC12 on gamma with CI ribbon & rug
#   5. O. tarda:   tmmx on gamma
#   6. O. tarda:   tmmx on epsilon
#
# Output:
#   figs/pub_fig_response_curves.png  (300 DPI)
#   figs/pub_fig_response_curves.pdf
###############################################################################

library(unmarked)
library(ggplot2)
library(dplyr)
library(patchwork)

cat("\n", strrep("=", 70), "\n")
cat("  FIGURE: Response curves (6 key driver x process combinations)\n")
cat(strrep("=", 70), "\n\n")

# --- Paths ---
res_dir  <- "results"
figs_dir <- "figs"
if (!dir.exists(figs_dir)) dir.create(figs_dir, recursive = TRUE)

# --- Variable name mapping ---
var_map <- c(
  "LC12" = "Land_Cover_Type_1_Percent_Class_12",
  "LC13" = "Land_Cover_Type_1_Percent_Class_13",
  "LC6"  = "Land_Cover_Type_1_Percent_Class_6",
  "LC7"  = "Land_Cover_Type_1_Percent_Class_7"
)

# --- Load all models ---
sp_codes <- c("otitar", "ptealc", "pteori", "tettet")
sp_latin <- c(otitar = "Otis tarda", ptealc = "Pterocles alchata",
              pteori = "Pterocles orientalis", tettet = "Tetrax tetrax")
models <- list()
for (sp in sp_codes) {
  mod_path <- file.path(res_dir, paste0(sp, "_model_object.rds"))
  if (file.exists(mod_path)) {
    models[[sp]] <- readRDS(mod_path)
    cat("  Loaded model:", sp, "\n")
  } else {
    warning("Model not found: ", mod_path)
  }
}

# --- Delta method response curve function ---
response_curve <- function(mod, prefix, var_name, x_seq = seq(-2.5, 3, length.out = 200)) {
  # prefix: "col" or "ext"
  cc <- coef(mod)
  V  <- vcov(mod)

  int_name <- paste0(prefix, "(Int)")
  var_full <- paste0(prefix, "(", var_name, ")")

  if (!(var_full %in% names(cc))) {
    warning("Coefficient not found: ", var_full)
    return(NULL)
  }

  beta0 <- cc[int_name]
  beta_j <- cc[var_full]

  v00 <- V[int_name, int_name]
  vjj <- V[var_full, var_full]
  v0j <- V[int_name, var_full]

  logit_pred <- beta0 + beta_j * x_seq
  var_logit  <- v00 + x_seq^2 * vjj + 2 * x_seq * v0j
  se_logit   <- sqrt(pmax(var_logit, 0))

  data.frame(
    x     = x_seq,
    pred  = plogis(logit_pred),
    lower = plogis(logit_pred - 1.96 * se_logit),
    upper = plogis(logit_pred + 1.96 * se_logit)
  )
}

x_seq <- seq(-2.5, 3, length.out = 200)

# --- Theme ---
theme_rc <- theme_classic(base_size = 10) +
  theme(
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 10),
    plot.title = element_text(face = "bold", size = 10),
    plot.subtitle = element_text(size = 8, colour = "grey40")
  )

# --- Panel 1: T. tetrax -- LC12 on BOTH gamma and epsilon ---
cat("  Panel 1: T. tetrax LC12 gamma+epsilon\n")
lc12_full <- var_map["LC12"]
rc_gam_tt <- response_curve(models$tettet, "col", lc12_full, x_seq)
rc_ext_tt <- response_curve(models$tettet, "ext", lc12_full, x_seq)

p1 <- NULL
if (!is.null(rc_gam_tt) && !is.null(rc_ext_tt)) {
  rc_gam_tt$process <- "Colonisation (gamma)"
  rc_ext_tt$process <- "Extinction (epsilon)"
  rc_both <- bind_rows(rc_gam_tt, rc_ext_tt)

  p1 <- ggplot(rc_both, aes(x = x, y = pred, colour = process, fill = process)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, colour = NA) +
    geom_line(aes(linetype = process), linewidth = 0.9) +
    scale_colour_manual(values = c("Colonisation (gamma)" = "#0072B2",
                                    "Extinction (epsilon)" = "#D55E00")) +
    scale_fill_manual(values = c("Colonisation (gamma)" = "#0072B2",
                                  "Extinction (epsilon)" = "#D55E00")) +
    scale_linetype_manual(values = c("Colonisation (gamma)" = "solid",
                                      "Extinction (epsilon)" = "dashed")) +
    scale_y_continuous(labels = scales::percent_format(), limits = c(0, NA)) +
    labs(title = expression(italic("T. tetrax") ~ ": LC12 (cropland)"),
         x = "LC12 (standardised)", y = "Probability",
         colour = NULL, fill = NULL, linetype = NULL) +
    theme_rc + theme(legend.position = c(0.65, 0.85),
                      legend.background = element_blank(),
                      legend.text = element_text(size = 8))
}

# --- Panel 2: O. tarda -- LC13 (urban) on epsilon ---
cat("  Panel 2: O. tarda LC13 epsilon\n")
lc13_full <- var_map["LC13"]
rc_ext_ot <- response_curve(models$otitar, "ext", lc13_full, x_seq)

p2 <- NULL
if (!is.null(rc_ext_ot)) {
  p2 <- ggplot(rc_ext_ot, aes(x = x, y = pred)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#D55E00", alpha = 0.2) +
    geom_line(colour = "#D55E00", linewidth = 0.9) +
    scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
    labs(title = expression(italic("O. tarda") ~ ": LC13 (urban) on " * epsilon),
         x = "LC13 (standardised)", y = expression(epsilon)) +
    theme_rc
}

# --- Panel 3: P. orientalis -- NDVI on epsilon ---
cat("  Panel 3: P. orientalis NDVI epsilon\n")
rc_ext_po <- response_curve(models$pteori, "ext", "NDVI", x_seq)

p3_subtitle <- NULL
if (file.exists(file.path(res_dir, "attribution_ndvi_decomposed_table.csv"))) {
  p3_subtitle <- "See NDVI decomposition table for climate vs. land-use attribution"
}

p3 <- NULL
if (!is.null(rc_ext_po)) {
  p3 <- ggplot(rc_ext_po, aes(x = x, y = pred)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#009E73", alpha = 0.2) +
    geom_line(colour = "#009E73", linewidth = 0.9) +
    scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
    labs(title = expression(italic("P. orientalis") ~ ": NDVI on " * epsilon),
         subtitle = p3_subtitle,
         x = "NDVI (standardised)", y = expression(epsilon)) +
    theme_rc
}

# --- Panel 4: T. tetrax -- LC12 on gamma with CI and rug ---
cat("  Panel 4: T. tetrax LC12 gamma with rug\n")

# Try to get observed LC12 values from the model's site covariates
rug_data <- NULL
tryCatch({
  umf <- models$tettet@data
  sc <- siteCovs(umf)
  if (lc12_full %in% names(sc)) {
    rug_data <- data.frame(x_obs = sc[[lc12_full]])
  } else {
    # Try yearly site covariates
    ysc <- yearlySiteCovs(umf)
    if (lc12_full %in% names(ysc)) {
      rug_data <- data.frame(x_obs = ysc[[lc12_full]])
    }
  }
}, error = function(e) {
  cat("    Could not extract rug data:", e$message, "\n")
})

p4 <- NULL
if (!is.null(rc_gam_tt)) {
  p4 <- ggplot(rc_gam_tt, aes(x = x, y = pred)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#CC79A7", alpha = 0.2) +
    geom_line(colour = "#CC79A7", linewidth = 0.9) +
    scale_y_continuous(labels = scales::percent_format(), limits = c(0, NA)) +
    labs(title = expression(italic("T. tetrax") ~ ": LC12 on " * gamma),
         x = "LC12 (standardised)", y = expression(gamma)) +
    theme_rc

  if (!is.null(rug_data)) {
    p4 <- p4 + geom_rug(data = rug_data, aes(x = x_obs), inherit.aes = FALSE,
                          sides = "b", colour = "grey40", alpha = 0.3, length = unit(0.03, "npc"))
  }
}

# --- Panel 5: O. tarda -- tmmx on gamma ---
cat("  Panel 5: O. tarda tmmx gamma\n")
rc_gam_ot <- response_curve(models$otitar, "col", "tmmx", x_seq)

p5 <- NULL
if (!is.null(rc_gam_ot)) {
  p5 <- ggplot(rc_gam_ot, aes(x = x, y = pred)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#D55E00", alpha = 0.2) +
    geom_line(colour = "#D55E00", linewidth = 0.9) +
    scale_y_continuous(labels = scales::percent_format(), limits = c(0, NA)) +
    labs(title = expression(italic("O. tarda") ~ ": Max. temperature on " * gamma),
         x = "tmmx (standardised)", y = expression(gamma)) +
    theme_rc
}

# --- Panel 6: O. tarda -- tmmx on epsilon ---
cat("  Panel 6: O. tarda tmmx epsilon\n")
rc_ext_ot_tmmx <- response_curve(models$otitar, "ext", "tmmx", x_seq)

p6 <- NULL
if (!is.null(rc_ext_ot_tmmx)) {
  p6 <- ggplot(rc_ext_ot_tmmx, aes(x = x, y = pred)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#D55E00", alpha = 0.2) +
    geom_line(colour = "#D55E00", linewidth = 0.9) +
    scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
    labs(title = expression(italic("O. tarda") ~ ": Max. temperature on " * epsilon),
         x = "tmmx (standardised)", y = expression(epsilon)) +
    theme_rc
}

# --- Assemble 3x2 layout ---
panels <- list(p1, p2, p3, p4, p5, p6)
panels <- panels[!sapply(panels, is.null)]

if (length(panels) == 0) {
  stop("No panels could be created. Check model files.")
}

# Pad to 6 if needed
while (length(panels) < 6) {
  panels[[length(panels) + 1]] <- plot_spacer()
}

p_final <- (panels[[1]] + panels[[2]] + panels[[3]]) /
           (panels[[4]] + panels[[5]] + panels[[6]])

# --- Save ---
ggsave(file.path(figs_dir, "pub_fig_response_curves.png"), p_final,
       width = 260, height = 180, units = "mm", dpi = 300)
ggsave(file.path(figs_dir, "pub_fig_response_curves.pdf"), p_final,
       width = 260, height = 180, units = "mm")

cat("\n  -> pub_fig_response_curves.png\n")
cat("  -> pub_fig_response_curves.pdf\n")
cat("\nDone.\n")
