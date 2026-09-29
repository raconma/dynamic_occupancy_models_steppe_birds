###############################################################################
# 48_figures_dynamics.R
#
# Two main-text figures of the dynamic results:
#   Fig. 2  SPA cover and colonisation and extinction: odds ratios (cell fully
#           inside vs fully outside SPAs). Coloured: model with all controls
#           (script 46); grey: the other control sets (scripts 43-45), as a
#           robustness check. Fits with a standard error above 5 are omitted.
#   Fig. 3  Covariates of colonisation and extinction in the primary models
#           (script 41): standardised logit-scale coefficients with 95 % Wald
#           intervals; intervals beyond the axis are truncated (arrows).
#
# Inputs:  results/core_spa/spa_in_dynamics.csv, spa_in_dynamics_checks.csv,
#          spa_neighbourhood.csv, spa_full_model.csv; results/{sp}_model_object.rds
# Outputs: results/manuscript_figures/submission/Fig2_spa_dynamics.{png,pdf,tif}
#          results/manuscript_figures/submission/Fig3_drivers.{png,pdf,tif}
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
        legend.position = "none", plot.margin = margin(4, 6, 4, 4))
save_fig <- function(p, name, w, h) {
  for (ext in c("png", "pdf")) ggsave(file.path(out_dir, paste0(name, ".", ext)), p, width = w, height = h, units = "mm", dpi = 600)
  ggsave(file.path(out_dir, paste0(name, ".tif")), p, width = w, height = h, units = "mm", dpi = 600, compression = "lzw")
}

# --- Fig. 2: SPA cover and colonisation and extinction -----------------------
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
save_fig(p2, "Fig2_spa_dynamics", 140, 70)

# --- Fig. 3: covariates of colonisation and extinction -----------------------
lab <- c(pseudo_steppe = "Pseudo-steppe", Land_Cover_Type_1_Percent_Class_12 = "Cropland",
         Land_Cover_Type_1_Percent_Class_7 = "Open shrubland", grass_cover = "Herbaceous cover",
         tree_cover = "Tree cover", bio1 = "Mean temperature", bio2 = "Temperature range",
         topo_aspect = "Aspect")
rows <- list()
for (sp in species) {
  m <- readRDS(here("results", paste0(sp, "_model_object.rds")))
  cf <- coef(m); se <- sqrt(diag(vcov(m)))
  i <- grepl("^(col|ext)[(]", names(cf)) & !grepl("Int", names(cf))
  rows[[sp]] <- data.frame(species = sp, process = substr(names(cf)[i], 1, 3),
                           covariate = sub("^[a-z]+[(](.*)[)]$", "\\1", names(cf)[i]),
                           est = cf[i], lo = cf[i] - 1.96 * se[i], hi = cf[i] + 1.96 * se[i])
}
r <- do.call(rbind, rows)
r$process <- factor(ifelse(r$process == "col", "A  Colonisation", "B  Extinction"))
r$cov <- factor(lab[r$covariate], levels = rev(lab))
r$sp <- factor(r$species, levels = species)
XL <- c(-6, 6)
r$lo_c <- pmax(r$lo, XL[1]); r$hi_c <- pmin(r$hi, XL[2]); r$est_c <- pmin(pmax(r$est, XL[1]), XL[2])
r$cut_hi <- r$hi > XL[2]; r$cut_lo <- r$lo < XL[1]
# manual vertical offsets: species present in a row are spread around it,
# O. tarda on top (as in the legend)
r <- r[order(r$process, r$cov, match(r$species, species)), ]
r$y <- as.numeric(r$cov) + ave(match(r$species, species), r$process, r$cov, FUN = function(k) {
  n <- length(k); if (n == 1) 0 else seq(0.22, -0.22, length.out = n) * min(1, (n - 1) / 2)
})
p3 <- ggplot(r, aes(y = y, colour = sp)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.3) +
  geom_segment(aes(x = lo_c, xend = hi_c, yend = y), linewidth = 0.5) +
  geom_segment(data = r[r$cut_hi, ], aes(x = XL[2] - 0.6, xend = XL[2], yend = y), linewidth = 0.5,
               arrow = arrow(length = unit(1.4, "mm"), type = "closed"), show.legend = FALSE) +
  geom_segment(data = r[r$cut_lo, ], aes(x = XL[1] + 0.6, xend = XL[1], yend = y), linewidth = 0.5,
               arrow = arrow(length = unit(1.4, "mm"), type = "closed"), show.legend = FALSE) +
  geom_point(data = r[abs(r$est) <= XL[2], ], aes(x = est), size = 1.8) +
  scale_colour_manual(values = sp_col, labels = sp_lab, breaks = species, name = NULL) +
  scale_x_continuous(limits = XL, breaks = seq(-6, 6, 3), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq_along(levels(r$cov)), labels = levels(r$cov),
                     limits = c(0.5, length(levels(r$cov)) + 0.5), expand = c(0, 0)) +
  facet_wrap(~ process, nrow = 1) +
  labs(x = "Standardised effect (logit scale)", y = NULL) +
  theme_fig + theme(legend.position = "bottom", legend.text = element_text(face = "italic"),
                    legend.key.width = unit(4, "mm"), panel.spacing.x = unit(6, "mm"),
                    panel.grid.major.y = element_line(colour = "grey96", linewidth = 0.3))
save_fig(p3, "Fig3_drivers", 140, 90)
print(r[r$cut_hi | r$cut_lo, c("species", "process", "covariate", "est", "lo", "hi")])

sessionInfo()
