###############################################################################
# 48_figures_dynamics.R
#
# Figures of the dynamic results:
#   Fig. 3  Covariates of initial occupancy, colonisation and extinction in
#           the primary models (script 41): standardised logit-scale
#           coefficients with 95 % Wald intervals; intervals beyond the axis
#           are truncated (arrows), estimates beyond it drawn pale at the edge.
#   Fig. 4  SPAs judged by occupancy and by dynamics: initial occupancy and
#           annual extinction averaged over recorded cells with SPA cover set
#           to 0 and 1 (script 49, PREDICT=1; primary SPA model).
#   Fig. 2  Dynamic model against a map updated yearly, rolling origin
#           (script 54): mean gain in held-out log-likelihood per test year
#           (2021-2023) and the three test years.
#   Fig. S5 The same against all four static references, summed over the
#           test years.
#   Fig. S4 SPA effects on colonisation and extinction in models without SPA
#           cover in initial occupancy (scripts 43-46). Coloured: all
#           controls (script 46); grey: other control sets. Fits with a
#           standard error above 5 are omitted.
#
# Inputs:  results/core_spa/spa_in_dynamics.csv, spa_in_dynamics_checks.csv,
#          spa_neighbourhood.csv, spa_full_model.csv, spa_probabilities.csv;
#          results/robustness/rolling_validation.csv;
#          results/{sp}_model_object.rds
# Outputs: results/manuscript_figures/submission/
#          Fig2_temporal_validation, Fig3_responses, Fig4_spa_static_dynamic,
#          FigS4_spa_controls, FigS5_rolling_all_references (.png, .pdf, .tif)
###############################################################################

library(here)
library(unmarked)
library(ggplot2)

set.seed(42)

species <- c("otitar", "tettet", "pteori", "ptealc")
sp_lab <- c(otitar = "O. tarda", tettet = "T. tetrax", pteori = "P. orientalis", ptealc = "P. alchata")
sp_col <- c(otitar = "#D55E00", tettet = "#CC79A7", pteori = "#009E73", ptealc = "#0072B2")
out_dir <- here("results", "manuscript_figures", "submission")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
theme_fig <- theme_classic(base_size = 8) +
  theme(axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
        strip.background = element_blank(), strip.text = element_text(face = "bold", size = 8, hjust = 0),
        panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
        axis.text = element_text(size = 7), legend.position = "none", plot.margin = margin(4, 6, 4, 4))
minus <- function(x) sub("^-", "\u2212", format(x, trim = TRUE))   # typographic minus on axes
save_fig <- function(p, name, w, h) {
  for (ext in c("png", "pdf")) ggsave(file.path(out_dir, paste0(name, ".", ext)), p, width = w, height = h, units = "mm", dpi = 600)
  ggsave(file.path(out_dir, paste0(name, ".tif")), p, width = w, height = h, units = "mm", dpi = 600, compression = "lzw")
}

# --- Fig. S4: SPA effects without SPA cover in initial occupancy -------------
cd <- here("results", "core_spa")
keep <- c("species", "term", "estimate", "lo", "hi", "max_se", "set")
a <- read.csv(file.path(cd, "spa_in_dynamics.csv")); a <- a[a$model == "spa", ]; a$set <- "other"
b <- read.csv(file.path(cd, "spa_in_dynamics_checks.csv"))
b <- b[b$check %in% c("effort", "det_spa", "coords", "all"), ]; b$set <- "other"
n <- read.csv(file.path(cd, "spa_neighbourhood.csv")); n <- n[n$neighbourhood == "both", ]; n$set <- "other"
f <- read.csv(file.path(cd, "spa_full_model.csv")); f$set <- "main"
d <- rbind(a[, keep], b[, keep], n[, keep], f[, keep])
d <- d[d$term %in% c("col(spa)", "ext(spa)") & is.finite(d$lo) & is.finite(d$max_se) & d$max_se <= 5, ]
d$process <- factor(ifelse(d$term == "col(spa)", "A  Colonisation", "B  Extinction"))
d$sp <- factor(sp_lab[d$species], levels = rev(sp_lab))
LIM <- c(0.05, 100)
d$lo_c <- pmax(exp(d$lo), LIM[1]); d$hi_c <- pmin(exp(d$hi), LIM[2])
d$trunc <- exp(d$hi) > LIM[2] | exp(d$lo) < LIM[1]
oth <- d[d$set == "other", ]; mn <- d[d$set == "main", ]
oth <- oth[order(oth$sp, oth$process), ]
oth$off <- ave(seq_len(nrow(oth)), oth$sp, oth$process, FUN = function(i) seq(-0.32, -0.12, length.out = length(i)))
p2 <- ggplot() +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey55", linewidth = 0.3) +
  geom_segment(data = oth, aes(x = lo_c, xend = hi_c, y = as.numeric(sp) + off, yend = as.numeric(sp) + off),
               colour = "grey70", linewidth = 0.3) +
  geom_point(data = oth, aes(x = exp(estimate), y = as.numeric(sp) + off), colour = "grey60", size = 0.8) +
  geom_segment(data = mn, aes(x = lo_c, xend = hi_c, y = as.numeric(sp) + 0.12, yend = as.numeric(sp) + 0.12,
                              colour = species), linewidth = 0.7) +
  geom_point(data = mn, aes(x = exp(estimate), y = as.numeric(sp) + 0.12, colour = species), size = 2.2) +
  geom_segment(data = oth[exp(oth$hi) > LIM[2], ], aes(x = LIM[2] / 1.6, xend = LIM[2], y = as.numeric(sp) + off,
               yend = as.numeric(sp) + off), colour = "grey60", linewidth = 0.3,
               arrow = arrow(length = unit(1.2, "mm"), type = "closed")) +
  geom_segment(data = mn[exp(mn$hi) > LIM[2], ], aes(x = LIM[2] / 1.6, xend = LIM[2], y = as.numeric(sp) + 0.12,
               yend = as.numeric(sp) + 0.12, colour = species), linewidth = 0.7,
               arrow = arrow(length = unit(1.6, "mm"), type = "closed")) +
  scale_colour_manual(values = sp_col) +
  scale_x_log10(limits = LIM, breaks = c(0.1, 1, 10, 100), labels = c("0.1", "1", "10", "100"), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq_along(levels(d$sp)), labels = levels(d$sp), limits = c(0.5, 4.5), expand = c(0, 0)) +
  facet_wrap(~ process, nrow = 1) +
  labs(x = "Odds ratio, SPA vs unprotected cell", y = NULL) +
  theme_fig + theme(axis.text.y = element_text(face = "italic", size = 8), panel.spacing.x = unit(6, "mm"))
save_fig(p2, "FigS4_spa_controls", 140, 70)

# --- Fig. 3: covariates of initial occupancy, colonisation and extinction ----
lab <- c(pseudo_steppe = "Pseudo-steppe", Land_Cover_Type_1_Percent_Class_12 = "Cropland",
         Land_Cover_Type_1_Percent_Class_7 = "Open shrubland", grass_cover = "Herbaceous cover",
         tree_cover = "Tree cover", bio1 = "Mean temperature", bio2 = "Diurnal temperature range",
         topo_aspect = "Aspect")
rows <- list()
for (sp in species) {
  m <- readRDS(here("results", paste0(sp, "_model_object.rds")))
  cf <- coef(m); se <- sqrt(diag(vcov(m)))
  i <- grepl("^(psi|col|ext)[(]", names(cf)) & !grepl("Int", names(cf))
  rows[[sp]] <- data.frame(species = sp, process = sub("[(].*", "", names(cf)[i]),
                           covariate = sub("^[a-z]+[(](.*)[)]$", "\\1", names(cf)[i]),
                           est = cf[i], lo = cf[i] - 1.96 * se[i], hi = cf[i] + 1.96 * se[i])
}
r <- do.call(rbind, rows)
r$proc <- factor(c(psi = "A  Initial occupancy", col = "B  Colonisation", ext = "C  Extinction")[r$process])
r$cov <- factor(lab[r$covariate], levels = rev(lab))
r$sp <- factor(r$species, levels = species)
XL <- c(-7, 7)
r$lo_c <- pmax(r$lo, XL[1]); r$hi_c <- pmin(r$hi, XL[2]); r$beyond <- abs(r$est) > XL[2]
r$cut_hi <- r$hi > XL[2]; r$cut_lo <- r$lo < XL[1]
# manual vertical offsets: species present in a row are spread around it, O. tarda on top
r <- r[order(r$proc, r$cov, match(r$species, species)), ]
r$y <- as.numeric(r$cov) + ave(match(r$species, species), r$proc, r$cov, FUN = function(k) {
  n <- length(k); if (n == 1) 0 else seq(0.24, -0.24, length.out = n) * min(1, (n - 1) / 2)
})
p3 <- ggplot(r, aes(y = y, colour = sp)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.3) +
  geom_segment(aes(x = lo_c, xend = hi_c, yend = y), linewidth = 0.5) +
  geom_segment(data = r[r$cut_hi, ], aes(x = XL[2] - 0.6, xend = XL[2], yend = y), linewidth = 0.5,
               arrow = arrow(length = unit(1.4, "mm"), type = "closed"), show.legend = FALSE) +
  geom_segment(data = r[r$cut_lo, ], aes(x = XL[1] + 0.6, xend = XL[1], yend = y), linewidth = 0.5,
               arrow = arrow(length = unit(1.4, "mm"), type = "closed"), show.legend = FALSE) +
  geom_point(data = r[!r$beyond, ], aes(x = est), size = 1.8) +
  # estimates beyond the axis: pale marker at the edge, not interpreted (cf. Fig. 4)
  geom_point(data = r[r$beyond, ], aes(x = sign(est) * (XL[2] - 0.25)), size = 1.8, alpha = 0.35,
             show.legend = FALSE) +
  scale_colour_manual(values = sp_col, labels = sp_lab, breaks = species, name = NULL) +
  scale_x_continuous(limits = XL, breaks = seq(-6, 6, 3), labels = minus, expand = c(0, 0)) +
  scale_y_continuous(breaks = seq_along(levels(r$cov)), labels = levels(r$cov),
                     limits = c(0.5, length(levels(r$cov)) + 0.5), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  facet_wrap(~ proc, nrow = 1) +
  labs(x = "Standardised effect (logit scale)", y = NULL) +
  theme_fig + theme(legend.position = "bottom", legend.text = element_text(face = "italic"),
                    legend.key.width = unit(4, "mm"), panel.spacing.x = unit(7, "mm"),
                    legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(1, "mm"),
                    panel.grid.major.y = element_line(colour = "grey96", linewidth = 0.3))
save_fig(p3, "Fig3_responses", 170, 85)
print(r[r$beyond, c("species", "process", "covariate", "est")])

# --- Fig. 4: SPAs judged by occupancy and by dynamics -------------------------
# initial occupancy and annual extinction averaged over recorded cells with SPA
# cover set to 0 (outside) and 1 (inside), primary SPA model (script 49,
# PREDICT=1); P. alchata omitted (SPA model not interpreted)
q <- read.csv(file.path(cd, "spa_probabilities.csv"))
q$panel <- factor(ifelse(q$quantity == "initial_occupancy", "A  Initial occupancy (snapshot)",
                         "B  Annual local extinction (dynamics)"))
q$sp <- factor(sp_lab[q$species], levels = rev(sp_lab[c("otitar", "tettet", "pteori")]))
q$where <- factor(q$where, levels = c("outside", "inside"), labels = c("Outside SPAs", "Inside SPAs"))
q$y <- as.numeric(q$sp) + ifelse(q$where == "Outside SPAs", 0.13, -0.13)
shade <- data.frame(panel = factor("A  Initial occupancy (snapshot)", levels = levels(q$panel)))
pct <- function(x) paste0(round(100 * x), " %")
p4 <- ggplot(q, aes(y = y, colour = species)) +
  geom_rect(data = shade, aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), fill = "grey97", inherit.aes = FALSE) +
  geom_vline(xintercept = c(0.2, 0.4, 0.6), colour = "grey92", linewidth = 0.3) +
  geom_segment(data = reshape(q[, c("species", "panel", "sp", "where", "estimate")], idvar = c("species", "panel", "sp"),
                              timevar = "where", direction = "wide"),
               aes(x = `estimate.Outside SPAs`, xend = `estimate.Inside SPAs`,
                   y = as.numeric(sp) + 0.13, yend = as.numeric(sp) - 0.13, colour = species),
               linewidth = 0.35, alpha = 0.6) +
  geom_segment(aes(x = lo, xend = hi, yend = y), linewidth = 0.5) +
  geom_point(aes(x = estimate, shape = where), size = 2.2, stroke = 0.7, fill = "white") +
  scale_colour_manual(values = sp_col, guide = "none") +
  scale_shape_manual(values = c(`Outside SPAs` = 21, `Inside SPAs` = 16), name = NULL) +
  scale_x_continuous(labels = pct, limits = c(0, 0.62), breaks = seq(0, 0.6, 0.2), expand = c(0, 0)) +
  scale_y_continuous(breaks = 1:3, labels = levels(q$sp), limits = c(0.5, 3.5), expand = c(0, 0)) +
  facet_wrap(~ panel, nrow = 1) +
  labs(x = "Probability for a recorded cell", y = NULL) +
  theme_fig + theme(axis.text.y = element_text(face = "italic", size = 8), panel.spacing.x = unit(6, "mm"),
                    legend.position = "bottom", legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(1, "mm")) +
  guides(shape = guide_legend(override.aes = list(colour = "grey30")))
save_fig(p4, "Fig4_spa_static_dynamic", 150, 62)

# --- Fig. 2: dynamic versus static prediction, rolling origin ----------------
# gain of the dynamic model over four static references, summed over the test
# years 2021-2023 (script 54); filled symbols, maps that use detections
v <- read.csv(here("results", "robustness", "rolling_validation.csv"))
v <- v[v$test_year == "pooled", ]
ref_lab <- c(S1 = "Persistence map", S2 = "Map updated yearly", S3 = "Three-year atlas",
             S4 = "Habitat-only model")
v$ref <- factor(ref_lab[v$reference], levels = ref_lab)
v$panel <- factor(ifelse(v$cells == "all", "A  All surveyed cells", "B  Cells with earlier records"))
v$sp <- factor(sp_lab[v$species], levels = rev(sp_lab))
v$y <- as.numeric(v$sp) + c(0.27, 0.09, -0.09, -0.27)[as.integer(v$ref)]
XL2 <- c(-80, 345)
v$lo_c <- pmax(v$lo, XL2[1]); v$hi_c <- pmin(v$hi, XL2[2])
ann2 <- data.frame(panel = factor("B  Cells with earlier records", levels = levels(v$panel)), x = XL2[2] - 5, y = 0.58,
                   lab = "dynamic model better \u2192")
p5 <- ggplot(v, aes(y = y, colour = species)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.3) +
  geom_segment(aes(x = lo_c, xend = hi_c, yend = y), linewidth = 0.5) +
  geom_point(aes(x = d_loglik, shape = ref), size = 1.9, stroke = 0.6, fill = "white") +
  geom_text(data = ann2, aes(x = x, y = y, label = lab), hjust = 1, size = 2.3, colour = "grey40", inherit.aes = FALSE) +
  scale_colour_manual(values = sp_col, guide = "none") +
  scale_shape_manual(values = c(16, 17, 15, 21), name = NULL) +
  scale_y_continuous(breaks = seq_along(levels(v$sp)), labels = levels(v$sp), limits = c(0.5, 4.5), expand = c(0, 0)) +
  scale_x_continuous(limits = XL2, breaks = seq(0, 300, 100), labels = minus) +
  coord_cartesian(clip = "off") +
  facet_wrap(~ panel, nrow = 1) +
  labs(x = "Log-likelihood of the 2021\u20132023 detections, dynamic minus static model", y = NULL) +
  theme_fig + theme(axis.text.y = element_text(face = "italic", size = 8), panel.spacing.x = unit(6, "mm"),
                    legend.position = "bottom", legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(1, "mm"),
                    panel.grid.major.y = element_line(colour = "grey96", linewidth = 0.3)) +
  guides(shape = guide_legend(nrow = 1, override.aes = list(colour = "grey30", size = 2)))
save_fig(p5, "FigS5_rolling_all_references", 150, 70)

# --- Fig. 2: dynamic model against a map updated yearly ------------------------
# mean gain per test year (pooled sum / 3, with its block-bootstrap interval)
# and the three test years as small points
v2 <- read.csv(here("results", "robustness", "rolling_validation.csv"))
v2 <- v2[v2$reference == "S2", ]
v2$panel <- factor(ifelse(v2$cells == "all", "A  All surveyed cells", "B  Cells with earlier records"))
v2$sp <- factor(sp_lab[v2$species], levels = rev(sp_lab))
pool <- v2[v2$test_year == "pooled", ]; pool[, c("d_loglik", "lo", "hi")] <- pool[, c("d_loglik", "lo", "hi")] / 3
yrs <- v2[v2$test_year != "pooled", ]
yrs$y <- as.numeric(yrs$sp) - 0.22 + 0.11 * (as.integer(yrs$test_year) - 2021)
ann3 <- data.frame(panel = factor("B  Cells with earlier records", levels = levels(v2$panel)), x = 118, y = 0.62,
                   lab = "dynamic model better \u2192")
p6 <- ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.3) +
  geom_point(data = yrs, aes(x = d_loglik, y = y, colour = species), size = 1, alpha = 0.45) +
  geom_segment(data = pool, aes(x = lo, xend = hi, y = as.numeric(sp) + 0.12, yend = as.numeric(sp) + 0.12, colour = species),
               linewidth = 0.6) +
  geom_point(data = pool, aes(x = d_loglik, y = as.numeric(sp) + 0.12, colour = species), size = 2.3) +
  geom_text(data = ann3, aes(x = x, y = y, label = lab), hjust = 1, size = 2.3, colour = "grey40") +
  scale_colour_manual(values = sp_col, guide = "none") +
  scale_y_continuous(breaks = seq_along(levels(v2$sp)), labels = levels(v2$sp), limits = c(0.5, 4.5), expand = c(0, 0)) +
  scale_x_continuous(limits = c(-35, 120), breaks = seq(-25, 100, 25), labels = minus) +
  coord_cartesian(clip = "off") +
  facet_wrap(~ panel, nrow = 1) +
  labs(x = "Gain in log-likelihood per test year, dynamic model over a map updated yearly", y = NULL) +
  theme_fig + theme(axis.text.y = element_text(face = "italic", size = 8), panel.spacing.x = unit(6, "mm"),
                    panel.grid.major.y = element_line(colour = "grey96", linewidth = 0.3))
save_fig(p6, "Fig2_temporal_validation", 140, 58)
print(pool[, c("species", "cells", "d_loglik", "lo", "hi")]); print(range(yrs$d_loglik))

sessionInfo()
