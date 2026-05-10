###############################################################################
# compare_v1_v2.R
#
# Side-by-side comparison of v1 (without stepRep_obs) and v2 (with
# stepRep_obs in detection sub-model) colext fits for the four steppe
# species. Produces:
#   - comparison_table.csv (all coefficients, both models, plus deltas)
#   - comparison_summary.csv (one row per species: AIC, log-lik,
#     stepRep effect, gamma/epsilon shifts on logit scale)
#   - figs/v1_vs_v2_forest.png   (forest plot of coefficient changes)
#   - figs/v1_vs_v2_aic.png      (delta-AIC bar)
#   - figs/v1_vs_v2_gamma_epsilon.png (annual mean gamma & epsilon)
###############################################################################

suppressPackageStartupMessages({
  library(here); library(unmarked); library(data.table); library(ggplot2)
})

REPO <- "/Users/gfandos/Documents/GitHub/dynamic_occupancy_models_steppe_birds"
RUN  <- file.path(REPO, "results", "stepRep_v2_run")
dir.create(file.path(RUN, "figs"), recursive = TRUE, showWarnings = FALSE)

species <- c(otitar = "Otis tarda", ptealc = "Pterocles alchata",
             pteori = "Pterocles orientalis", tettet = "Tetrax tetrax")

get_coefs <- function(fit, model_label) {
  rbindlist(lapply(c("psi","col","ext","det"), function(type) {
    co <- tryCatch(coef(fit, type = type), error = function(e) NULL)
    if (is.null(co) || length(co) == 0) return(NULL)
    se <- tryCatch(sqrt(diag(vcov(fit, type = type))),
                   error = function(e) rep(NA_real_, length(co)))
    z  <- co / se
    data.table(param = names(co), Estimate = unname(co), SE = unname(se),
               z = unname(z),
               Pval = 2 * pnorm(-abs(unname(z))),
               type = type, model = model_label)
  }))
}

# ---- Per-species AIC / coefs / annual gamma & epsilon ----
all_coefs <- list()
summary_rows <- list()
ge_rows <- list()  # annual gamma, epsilon at the data points

for (sp in names(species)) {
  cat("\n== ", sp, " ==\n", sep = "")

  v1_path <- file.path(REPO, "results", paste0(sp, "_model_object.rds"))
  v2_path <- file.path(RUN,  paste0(sp, "_v2_model_object.rds"))
  if (!file.exists(v1_path) || !file.exists(v2_path)) {
    cat("  missing fit file -- skipping\n"); next
  }
  m1 <- readRDS(v1_path)
  m2 <- readRDS(v2_path)

  c1 <- get_coefs(m1, "v1")
  c2 <- get_coefs(m2, "v2")
  c1[, species := sp]; c2[, species := sp]
  # Strip unmarked's "psi(...)" / "p(...)" / "col(...)" / "ext(...)" wrappers
  for (dt in list(c1, c2))
    dt[, param := gsub("^(psi|p|col|ext)\\(|\\)$", "", param)]
  all_coefs[[sp]] <- rbind(c1, c2)

  # AIC and log-lik
  ll1 <- m1@AIC; ll2 <- m2@AIC
  step_row <- c2[type == "det" & param == "stepRep_obs"]
  step_beta <- if (nrow(step_row) > 0) step_row$Estimate else NA_real_
  step_se   <- if (nrow(step_row) > 0) step_row$SE       else NA_real_

  # Annual mean gamma and epsilon -- unmarked's predict() at each
  # cell-year using model's own data
  pred_gamma_v1 <- tryCatch(predict(m1, type = "col"), error = function(e) NULL)
  pred_gamma_v2 <- tryCatch(predict(m2, type = "col"), error = function(e) NULL)
  pred_ext_v1   <- tryCatch(predict(m1, type = "ext"), error = function(e) NULL)
  pred_ext_v2   <- tryCatch(predict(m2, type = "ext"), error = function(e) NULL)

  ge_rows[[sp]] <- data.table(
    species = sp,
    metric  = c("gamma_mean","gamma_mean","epsilon_mean","epsilon_mean"),
    model   = c("v1","v2","v1","v2"),
    value   = c(
      if (!is.null(pred_gamma_v1)) mean(pred_gamma_v1$Predicted, na.rm = TRUE) else NA,
      if (!is.null(pred_gamma_v2)) mean(pred_gamma_v2$Predicted, na.rm = TRUE) else NA,
      if (!is.null(pred_ext_v1))   mean(pred_ext_v1$Predicted,   na.rm = TRUE) else NA,
      if (!is.null(pred_ext_v2))   mean(pred_ext_v2$Predicted,   na.rm = TRUE) else NA
    )
  )

  summary_rows[[sp]] <- data.table(
    species   = sp,
    species_n = species[[sp]],
    AIC_v1    = round(ll1, 2),
    AIC_v2    = round(ll2, 2),
    delta_AIC = round(ll2 - ll1, 2),
    stepRep_beta = round(step_beta, 3),
    stepRep_SE   = round(step_se,   3),
    stepRep_z    = round(step_beta / step_se, 2),
    gamma_mean_v1 = round(mean(pred_gamma_v1$Predicted, na.rm = TRUE), 5),
    gamma_mean_v2 = round(mean(pred_gamma_v2$Predicted, na.rm = TRUE), 5),
    delta_gamma_pct = round(100 * (mean(pred_gamma_v2$Predicted, na.rm = TRUE) /
                                   mean(pred_gamma_v1$Predicted, na.rm = TRUE) - 1), 1),
    epsilon_mean_v1 = round(mean(pred_ext_v1$Predicted, na.rm = TRUE), 5),
    epsilon_mean_v2 = round(mean(pred_ext_v2$Predicted, na.rm = TRUE), 5),
    delta_epsilon_pct = round(100 * (mean(pred_ext_v2$Predicted, na.rm = TRUE) /
                                     mean(pred_ext_v1$Predicted, na.rm = TRUE) - 1), 1)
  )
  print(summary_rows[[sp]])
}

# ---- Save tables ----
all_coefs_dt   <- rbindlist(all_coefs, fill = TRUE)
all_coefs_dt[, species_name := species[species]]
fwrite(all_coefs_dt, file.path(RUN, "comparison_table.csv"))

summary_dt <- rbindlist(summary_rows, fill = TRUE)
fwrite(summary_dt, file.path(RUN, "comparison_summary.csv"))

ge_dt <- rbindlist(ge_rows, fill = TRUE)
ge_dt[, species_name := species[species]]
fwrite(ge_dt, file.path(RUN, "gamma_epsilon_comparison.csv"))

cat("\n=== Summary table ===\n")
print(summary_dt)

# ---- Figure 1: forest plot of detection coefficients (v1 vs v2) ----
det_coefs <- all_coefs_dt[type == "det"]
det_coefs[, species_name := species[species]]
det_coefs[, model_lab := factor(model, levels = c("v2","v1"))]

p_forest <- ggplot(det_coefs,
                   aes(x = Estimate, y = param, colour = model_lab)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(xmin = Estimate - 1.96*SE, xmax = Estimate + 1.96*SE),
                 height = 0.18, position = position_dodge(0.5), linewidth = 0.6) +
  geom_point(position = position_dodge(0.5), size = 2.4) +
  facet_wrap(~ species_name, scales = "free_y", ncol = 2) +
  scale_colour_manual(values = c(v1 = "grey55", v2 = "#d62728"),
                      breaks = c("v1","v2"),
                      labels = c("v1 (no stepRep)", "v2 (+ stepRep)"),
                      name = NULL) +
  labs(x = "Coefficient (logit-scale)  ± 95% CI", y = NULL,
       title = "Detection sub-model coefficients: v1 vs v2") +
  theme_bw() + theme(legend.position = "top")
ggsave(file.path(RUN, "figs", "v1_vs_v2_detection_forest.png"),
       p_forest, width = 10, height = 7, dpi = 150)

# ---- Figure 2: AIC and gamma/epsilon shifts ----
p_aic <- ggplot(summary_dt, aes(x = reorder(species_n, delta_AIC),
                                 y = delta_AIC, fill = delta_AIC < 0)) +
  geom_col() + geom_hline(yintercept = 0) +
  geom_text(aes(label = round(delta_AIC, 1)), hjust = 1.1, colour = "white") +
  scale_fill_manual(values = c("TRUE" = "#2ca02c", "FALSE" = "#d62728"),
                    labels = c("TRUE" = "v2 better", "FALSE" = "v1 better"),
                    name = NULL) +
  coord_flip() +
  labs(x = NULL, y = "ΔAIC (v2 − v1)",
       title = "Model fit: ΔAIC after adding stepRep_obs") +
  theme_bw() + theme(legend.position = "top")
ggsave(file.path(RUN, "figs", "v1_vs_v2_aic.png"),
       p_aic, width = 8, height = 4, dpi = 150)

# Gamma/epsilon shifts
ge_long <- ge_dt[, .(species_name, metric, model, value)]
ge_long[, metric_lab := factor(metric,
                                levels = c("gamma_mean", "epsilon_mean"),
                                labels = c("Mean colonisation γ",
                                            "Mean extinction ε"))]
p_ge <- ggplot(ge_long, aes(x = species_name, y = value, fill = model)) +
  geom_col(position = position_dodge(0.7), width = 0.6) +
  geom_text(aes(label = sprintf("%.4f", value)),
            position = position_dodge(0.7), vjust = -0.3, size = 3) +
  facet_wrap(~ metric_lab, scales = "free_y") +
  scale_fill_manual(values = c(v1 = "grey55", v2 = "#d62728"),
                    labels = c(v1 = "v1 (no stepRep)",
                                v2 = "v2 (+ stepRep)"),
                    name = NULL) +
  labs(x = NULL, y = "Mean predicted probability") +
  theme_bw() + theme(legend.position = "top",
                     axis.text.x = element_text(angle = 20, hjust = 1))
ggsave(file.path(RUN, "figs", "v1_vs_v2_gamma_epsilon.png"),
       p_ge, width = 10, height = 5, dpi = 150)

cat("\n=== artefacts written to", RUN, "===\n")
