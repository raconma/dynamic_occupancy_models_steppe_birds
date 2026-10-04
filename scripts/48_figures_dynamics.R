###############################################################################
# 48_figures_dynamics.R
#
# Figures of the dynamic results:
#   Fig. 3  Covariates of initial occupancy, colonisation and extinction in
#           the primary models (script 41): standardised logit-scale
#           coefficients with 95 % Wald intervals; intervals beyond the axis
#           are truncated (arrows), estimates beyond it drawn pale at the edge.
#   Fig. 4  SPAs judged by occupancy and by dynamics: odds ratios (cell fully
#           inside vs fully outside SPAs) for initial occupancy, colonisation
#           and extinction in the fixed SPA model (script 49, psi_spa_t).
#   Fig. 2  Dynamic versus static prediction of the 2023 detections (script
#           51): difference in held-out log-likelihood from the static
#           closed model, all surveyed cells and cells recorded in 2017-2022.
#   Fig. S4 SPA effects on colonisation and extinction in models without SPA
#           cover in initial occupancy (scripts 43-46). Coloured: all
#           controls (script 46); grey: other control sets. Fits with a
#           standard error above 5 are omitted.
#
# Inputs:  results/core_spa/spa_in_dynamics.csv, spa_in_dynamics_checks.csv,
#          spa_neighbourhood.csv, spa_full_model.csv, spa_psi_check.csv;
#          results/robustness/temporal_validation_2023.csv;
#          results/{sp}_model_object.rds
# Outputs: results/manuscript_figures/submission/
#          Fig2_temporal_validation, Fig3_responses, Fig4_spa_static_dynamic,
#          FigS4_spa_controls (.png, .pdf, .tif)
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
q <- read.csv(file.path(cd, "spa_psi_check.csv"))
q <- q[q$check == "psi_spa_t" & q$term %in% c("psi(spa)", "col(spa)", "ext(spa)"), ]
q$process <- factor(c("psi(spa)" = "A  Initial occupancy", "col(spa)" = "B  Colonisation",
                      "ext(spa)" = "C  Extinction")[q$term])
q$sp <- factor(sp_lab[q$species], levels = rev(sp_lab))
LIM3 <- c(0.05, 100)
q$lo_c <- pmax(exp(q$lo), LIM3[1]); q$hi_c <- pmin(exp(q$hi), LIM3[2])
q$alpha <- ifelse(q$species == "ptealc", 0.35, 1)          # not interpreted (SE rule)
p4 <- ggplot(q, aes(y = sp, colour = species, alpha = alpha)) +
  # shade the snapshot panel (occupancy) apart from the dynamic ones
  geom_rect(data = data.frame(process = factor("A  Initial occupancy", levels = levels(q$process))),
            aes(xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf), fill = "grey97", inherit.aes = FALSE) +
  geom_vline(xintercept = c(0.1, 10, 100), colour = "grey92", linewidth = 0.3) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey55", linewidth = 0.3) +
  geom_segment(aes(x = lo_c, xend = hi_c, yend = sp), linewidth = 0.7) +
  geom_segment(data = q[exp(q$hi) > LIM3[2], ], aes(x = LIM3[2] / 1.6, xend = LIM3[2], yend = sp),
               linewidth = 0.7, arrow = arrow(length = unit(1.6, "mm"), type = "closed")) +
  geom_point(aes(x = exp(estimate)), size = 2.2) +
  geom_text(data = data.frame(process = factor("C  Extinction", levels = levels(q$process)), x = 0.055, y = 0.62,
                              lab = "\u2190 lower extinction inside SPAs"),
            aes(x = x, y = y, label = lab), hjust = 0, size = 2.3, colour = "grey40", inherit.aes = FALSE) +
  scale_colour_manual(values = sp_col) + scale_alpha_identity() +
  scale_x_log10(limits = LIM3, breaks = c(0.1, 1, 10, 100), labels = c("0.1", "1", "10", "100"), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  facet_wrap(~ process, nrow = 1) +
  labs(x = "Odds ratio, SPA vs unprotected cell", y = NULL) +
  theme_fig + theme(axis.text.y = element_text(face = "italic", size = 8), panel.spacing.x = unit(5, "mm"))
save_fig(p4, "Fig4_spa_static_dynamic", 160, 55)

# --- Fig. 2: dynamic versus static prediction of 2023 -------------------------
# gain of the dynamic model over the two static references: a persistence map
# (fixed-state model, filled) and a habitat-only model (open)
v <- read.csv(here("results", "robustness", "temporal_validation_2023.csv"))
v <- v[v$model == "dynamic" & v$cells %in% c("all", "detected_2017_2022"), ]
v <- rbind(transform(v, ref = "Persistence map", est = d_loglik_vs_closed, lo = d_lo, hi = d_hi),
           transform(v, ref = "Habitat-only model", est = d_loglik_vs_sdm, lo = d_sdm_lo, hi = d_sdm_hi))
v$ref <- factor(v$ref, levels = c("Persistence map", "Habitat-only model"))
ann4 <- data.frame(panel = factor("A  All surveyed cells", levels = c("A  All surveyed cells", "B  Cells with records in 2017\u20132022")),
                   x = 128, y = 0.62, lab = "dynamic model better \u2192")
v$panel <- factor(ifelse(v$cells == "all", "A  All surveyed cells", "B  Cells with records in 2017\u20132022"))
v$sp <- factor(sp_lab[v$species], levels = rev(sp_lab))
v$y <- as.numeric(v$sp) + ifelse(v$ref == "Persistence map", 0.14, -0.14)
p5 <- ggplot(v, aes(y = y, colour = species)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.3) +
  geom_segment(aes(x = lo, xend = hi, yend = y, linetype = ref), linewidth = 0.6) +
  geom_point(aes(x = est, shape = ref, fill = after_scale(ifelse(shape == 21, "white", colour))), size = 2.1, stroke = 0.7) +
  scale_colour_manual(values = sp_col, guide = "none") +
  geom_text(data = ann4, aes(x = x, y = y, label = lab), hjust = 1, size = 2.3, colour = "grey40", inherit.aes = FALSE) +
  scale_shape_manual(values = c(`Persistence map` = 16, `Habitat-only model` = 21), name = "Static reference") +
  scale_linetype_manual(values = c(`Persistence map` = "solid", `Habitat-only model` = "22"), name = "Static reference") +
  scale_y_continuous(breaks = seq_along(levels(v$sp)), labels = levels(v$sp), limits = c(0.5, 4.5), expand = c(0, 0)) +
  scale_x_continuous(limits = c(-20, 130), breaks = seq(0, 120, 40), labels = minus) +
  facet_wrap(~ panel, nrow = 1) +
  labs(x = "Log-likelihood of the 2023 detections, dynamic minus static model", y = NULL) +
  theme_fig + theme(axis.text.y = element_text(face = "italic", size = 8), panel.spacing.x = unit(6, "mm"),
                    legend.position = "bottom", legend.key.width = unit(8, "mm"), legend.margin = margin(0, 0, 0, 0),
                    legend.box.spacing = unit(1, "mm")) +
  guides(shape = guide_legend(override.aes = list(colour = "grey30", fill = c("grey30", "white"))))
save_fig(p5, "Fig2_temporal_validation", 140, 62)

sessionInfo()
