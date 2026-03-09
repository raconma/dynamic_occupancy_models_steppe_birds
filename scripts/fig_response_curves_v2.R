## fig_response_curves_v2.R
## 2-row response-curve figure (Fig 5)
## Row 1: T. tetrax LC12 gamma + epsilon
## Row 2: O. tarda tmmx gamma | P. alchata pr epsilon | P. orientalis NDVI epsilon

library(unmarked)
library(ggplot2)
library(dplyr)
library(patchwork)

# ---------- load models ----------
sp_codes <- c("otitar", "ptealc", "pteori", "tettet")
models <- setNames(
  lapply(sp_codes, function(s) readRDS(paste0("results/", s, "_model_object.rds"))),
  sp_codes
)

# ---------- delta-method response curve ----------
response_curve <- function(mod, prefix, var_name, x_seq) {
  cc <- coef(mod); V <- vcov(mod)
  int_name <- paste0(prefix, "(Int)")
  var_full <- paste0(prefix, "(", var_name, ")")
  if (!(var_full %in% names(cc))) return(NULL)
  beta0 <- cc[int_name]; beta_j <- cc[var_full]
  v00 <- V[int_name, int_name]; vjj <- V[var_full, var_full]; v0j <- V[int_name, var_full]
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

# ---------- common theme ----------
base_theme <- theme_classic(base_size = 10) +
  theme(plot.title = element_text(size = 9, face = "plain"))

# ==========================================================
# ROW 1 – T. tetrax: gamma and epsilon vs LC12 (cropland)
# ==========================================================
x_lc12 <- seq(-2.5, 3, length.out = 300)
lc12_var <- "Land_Cover_Type_1_Percent_Class_12"

df_gam <- response_curve(models$tettet, "col", lc12_var, x_lc12) %>% mutate(process = "gamma")
df_eps <- response_curve(models$tettet, "ext", lc12_var, x_lc12) %>% mutate(process = "epsilon")

# rug data from scaling object (pool all years)
sc <- readRDS("results/tettet_train_dyn_scale.rds")
lc12_keys <- grep("^Land_Cover_Type_1_Percent_Class_12_", names(sc), value = TRUE)
# Each element has $center and $scale; we can approximate observed range
# but actual raw data is not stored. Use center +/- 2*scale as proxy range.
# Better: try to extract siteCovs from model
rug_vals <- tryCatch({
  ysc <- models$tettet@data@yearlySiteCovs
  if (!is.null(ysc) && lc12_var %in% names(ysc)) ysc[[lc12_var]] else NULL
}, error = function(e) NULL)

# crossover point (gamma = epsilon => beta0_col + beta_col*x = beta0_ext + beta_ext*x)
cc_tt <- coef(models$tettet)
b0_col <- cc_tt[paste0("col(Int)")]; b1_col <- cc_tt[paste0("col(", lc12_var, ")")]
b0_ext <- cc_tt[paste0("ext(Int)")]; b1_ext <- cc_tt[paste0("ext(", lc12_var, ")")]
x_cross <- (b0_ext - b0_col) / (b1_col - b1_ext)

p1 <- ggplot() +
  geom_ribbon(data = df_gam, aes(x, ymin = lower, ymax = upper), fill = "#E07B39", alpha = 0.15) +
  geom_ribbon(data = df_eps, aes(x, ymin = lower, ymax = upper), fill = "#B03A2E", alpha = 0.15) +
  geom_line(data = df_gam, aes(x, pred), colour = "#E07B39", linewidth = 0.9) +
  geom_line(data = df_eps, aes(x, pred), colour = "#B03A2E", linewidth = 0.9, linetype = "dashed") +
  labs(
    title = expression(italic("T. tetrax") ~ ": cropland drives both processes"),
    x = "Cropland cover (standardised)",
    y = "Predicted probability"
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  base_theme

# add rug if available
if (!is.null(rug_vals)) {
  rug_df <- data.frame(xr = rug_vals[!is.na(rug_vals)])
  p1 <- p1 + geom_rug(data = rug_df, aes(x = xr), sides = "b", alpha = 0.15, length = unit(0.03, "npc"))
}

# add crossover line if within plot range
if (!is.na(x_cross) && x_cross >= -2.5 && x_cross <= 3) {
  y_cross <- plogis(b0_col + b1_col * x_cross)
  p1 <- p1 +
    geom_vline(xintercept = x_cross, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
    annotate("text", x = x_cross + 0.15, y = y_cross + 0.07, label = expression(gamma == epsilon),
             size = 3, colour = "grey30", hjust = 0)
}

# direct labels on lines
lab_x <- 2.2
p1 <- p1 +
  annotate("text", x = lab_x,
           y = plogis(b0_col + b1_col * lab_x) + 0.04,
           label = expression(gamma ~ "(colonisation)"),
           colour = "#E07B39", size = 3, hjust = 0) +
  annotate("text", x = lab_x,
           y = plogis(b0_ext + b1_ext * lab_x) - 0.04,
           label = expression(epsilon ~ "(extinction)"),
           colour = "#B03A2E", size = 3, hjust = 0)

# ==========================================================
# ROW 2 – three single-driver panels
# ==========================================================
x_std <- seq(-2.5, 3, length.out = 300)

# --- O. tarda: tmmx on gamma ---
df_ot <- response_curve(models$otitar, "col", "tmmx", x_std)
p2 <- ggplot(df_ot, aes(x, pred)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#D55E00", alpha = 0.2) +
  geom_line(colour = "#D55E00", linewidth = 0.8) +
  labs(
    title = expression(italic("O. tarda") ~ ": Max. temp. on " * gamma),
    x = "Max. temperature (std.)", y = "Predicted probability"
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  base_theme

# --- P. alchata: pr on epsilon ---
df_pa <- response_curve(models$ptealc, "ext", "pr", x_std)
p3 <- ggplot(df_pa, aes(x, pred)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#0072B2", alpha = 0.2) +
  geom_line(colour = "#0072B2", linewidth = 0.8) +
  labs(
    title = expression(italic("P. alchata") ~ ": Precipitation on " * epsilon),
    x = "Precipitation (std.)", y = "Predicted probability"
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  base_theme

# --- P. orientalis: NDVI on epsilon ---
df_po <- response_curve(models$pteori, "ext", "NDVI", x_std)
p4 <- ggplot(df_po, aes(x, pred)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#009E73", alpha = 0.2) +
  geom_line(colour = "#009E73", linewidth = 0.8) +
  labs(
    title = expression(italic("P. orientalis") ~ ": NDVI on " * epsilon),
    x = "NDVI (std.)", y = "Predicted probability"
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  base_theme

# ==========================================================
# ASSEMBLY with patchwork
# ==========================================================
layout <- "AAA
BCD"

fig <- p1 + p2 + p3 + p4 +
  plot_layout(design = layout, heights = c(3, 2))   # 60 / 40

dir.create("figs", showWarnings = FALSE)

ggsave("figs/pub_fig5_response_curves.png", fig,
       width = 180, height = 160, units = "mm", dpi = 300)
ggsave("figs/pub_fig5_response_curves.pdf", fig,
       width = 180, height = 160, units = "mm")

cat("Done – saved figs/pub_fig5_response_curves.{png,pdf}\n")
