###############################################################################
# 48_figures_dynamics.R
#
# Two main-text figures of the dynamic results:
#   Fig. 2  SPA effects on colonisation (recruitment) and extinction
#           (retention): odds ratios of SPA cover (cell fully inside vs fully
#           outside SPAs) under successive control sets (scripts 43-46)
#   Fig. 3  Covariates of colonisation and extinction in the primary models
#           (script 41): predicted annual probability against each selected
#           covariate (standardised) over its 5th-95th percentile range in the
#           cells where the species was recorded, other covariates at their mean, with
#           95 % intervals from 1,000 multivariate-normal coefficient draws
#
# Inputs:  results/core_spa/spa_in_dynamics.csv, spa_in_dynamics_checks.csv,
#          spa_neighbourhood.csv, spa_full_model.csv; results/{sp}_model_object.rds
# Outputs: results/manuscript_figures/submission/Fig2_spa_dynamics.{png,pdf,tif}
#          results/manuscript_figures/submission/Fig3_drivers.{png,pdf,tif}
###############################################################################

library(here)
library(unmarked)
library(MASS)
library(ggplot2)

set.seed(42)

species <- c("otitar", "tettet", "pteori", "ptealc")
sp_lab <- c(otitar = "O. tarda", tettet = "T. tetrax", pteori = "P. orientalis", ptealc = "P. alchata")
sp_col <- c(otitar = "#D55E00", tettet = "#CC79A7", pteori = "#009E73", ptealc = "#0072B2")
out_dir <- here("results", "manuscript_figures", "submission")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
theme_fig <- theme_minimal(base_size = 8) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        strip.text = element_text(face = "italic"), legend.position = "none",
        plot.title = element_text(face = "bold", size = 8))
save_fig <- function(p, name, w, h) {
  for (ext in c("png", "pdf")) ggsave(file.path(out_dir, paste0(name, ".", ext)), p, width = w, height = h, units = "mm", dpi = 600)
  ggsave(file.path(out_dir, paste0(name, ".tif")), p, width = w, height = h, units = "mm", dpi = 600, compression = "lzw")
}

# --- Fig. 2: SPA effects under successive control sets ----------------------
cd <- here("results", "core_spa")
a <- read.csv(file.path(cd, "spa_in_dynamics.csv")); a <- a[a$model == "spa", ]
a$set <- "Habitat and climate"
b <- read.csv(file.path(cd, "spa_in_dynamics_checks.csv"))
b <- b[b$check %in% c("effort", "det_spa", "coords"), ]
b$set <- c(effort = "+ survey effort", det_spa = "+ SPA in detection", coords = "+ coordinates")[b$check]
n <- read.csv(file.path(cd, "spa_neighbourhood.csv")); n <- n[n$neighbourhood == "both", ]
n$set <- "+ neighbourhood"
f <- read.csv(file.path(cd, "spa_full_model.csv")); f$set <- "All controls"
keep <- c("species", "term", "estimate", "lo", "hi", "set", "max_se")
d <- rbind(a[, keep], b[, keep], n[, keep], f[, keep])
# fits that are not identifiable (non-finite standard errors or above 5) are
# left out (P. alchata with survey effort)
d <- d[d$term %in% c("col(spa)", "ext(spa)") & is.finite(d$lo) & is.finite(d$max_se) & d$max_se <= 5, ]
d$process <- factor(ifelse(d$term == "col(spa)", "A  Colonisation (recruitment)", "B  Extinction (retention)"))
lev <- c("Habitat and climate", "+ survey effort", "+ SPA in detection", "+ coordinates", "+ neighbourhood", "All controls")
d$set <- factor(d$set, levels = rev(lev))
d$sp <- factor(sp_lab[d$species], levels = sp_lab)
LIM <- c(0.02, 200)
d$lo_c <- pmax(exp(d$lo), LIM[1]); d$hi_c <- pmin(exp(d$hi), LIM[2])
p2 <- ggplot(d, aes(y = set, colour = species)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  geom_errorbar(aes(xmin = lo_c, xmax = hi_c), width = 0.25, linewidth = 0.4, orientation = "y") +
  geom_point(aes(x = exp(estimate)), size = 1.6) +
  scale_x_log10(limits = LIM, breaks = c(0.1, 1, 10, 100), labels = c("0.1", "1", "10", "100")) +
  scale_colour_manual(values = sp_col) +
  facet_grid(process ~ sp) +
  labs(x = "Odds ratio, cell inside vs outside SPAs (log scale)", y = NULL) +
  theme_fig + theme(strip.text.y = element_text(face = "bold", angle = 0, hjust = 0))
save_fig(p2, "Fig2_spa_dynamics", 190, 110)

# --- Fig. 3: covariates of colonisation and extinction ----------------------
lab <- c(pseudo_steppe = "Pseudo-steppe", bio1 = "Mean temperature", bio2 = "Diurnal temperature range",
         tree_cover = "Tree cover", grass_cover = "Herbaceous cover", topo_aspect = "Aspect",
         Land_Cover_Type_1_Percent_Class_12 = "Cropland (MODIS)", Land_Cover_Type_1_Percent_Class_7 = "Open shrubland (MODIS)")
rows <- list()
for (sp in species) {
  m <- readRDS(here("results", paste0(sp, "_model_object.rds")))
  cf <- coef(m); dr <- mvrnorm(1000, cf, vcov(m))
  u <- m@data; y <- getY(u); rec <- rowSums(y == 1, na.rm = TRUE) > 0
  n_t <- u@numPrimary
  cov_of <- function(v) {                       # values in recorded cells
    sc <- siteCovs(u); ys <- yearlySiteCovs(u)
    if (v %in% names(sc)) sc[[v]][rec] else ys[[v]][rep(rec, each = n_t)]
  }
  for (sub in c("col", "ext")) {
    terms <- sub(paste0("^", sub, "\\((.*)\\)$"), "\\1", grep(paste0("^", sub, "\\("), names(cf), value = TRUE))
    terms <- setdiff(terms, "Int")
    for (v in terms) {
      xs <- seq(quantile(cov_of(v), 0.05, na.rm = TRUE), quantile(cov_of(v), 0.95, na.rm = TRUE), length.out = 60)
      eta <- outer(dr[, paste0(sub, "(Int)")], rep(1, length(xs))) + outer(dr[, paste0(sub, "(", v, ")")], xs)
      pr <- plogis(eta)
      rows[[paste(sp, sub, v)]] <- data.frame(species = sp, process = sub, covariate = v, x = xs,
                                              med = apply(pr, 2, median), lo = apply(pr, 2, quantile, 0.025),
                                              hi = apply(pr, 2, quantile, 0.975))
    }
  }
}
r <- do.call(rbind, rows)
r$process <- factor(ifelse(r$process == "col", "Colonisation", "Extinction"))
r$covariate <- factor(lab[r$covariate], levels = lab)
r$sp <- factor(sp_lab[r$species], levels = sp_lab)
r$panel <- factor(paste0(r$process, "\n", r$covariate))
r$panel <- factor(r$panel, levels = unique(r$panel[order(r$process, r$covariate)]))
p3 <- ggplot(r, aes(x, med, colour = species, fill = species)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.5) +
  scale_colour_manual(values = sp_col, labels = sp_lab, breaks = species, name = NULL) +
  scale_fill_manual(values = sp_col, labels = sp_lab, breaks = species, name = NULL) +
  facet_wrap(~ panel, scales = "free", ncol = 5) +
  labs(x = "Covariate (standard deviations from the mean; 5th to 95th percentile in recorded cells)",
       y = "Annual probability") +
  theme_fig + theme(legend.position = "bottom", legend.text = element_text(face = "italic"),
                    strip.text = element_text(face = "plain"))
save_fig(p3, "Fig3_drivers", 190, 120)
print(table(r$species, r$covariate, r$process)[, , 1])

sessionInfo()
