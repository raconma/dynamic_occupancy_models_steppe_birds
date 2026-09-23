###############################################################################
# 29_manuscript_figures_tables.R
#
# Assembles the main-text figures and the supplementary tables of the
# Biological Conservation draft from the outputs of scripts 21-28. No new
# analyses: every number comes from those outputs.
#
# Figures (results/manuscript_figures/):
#   fig1_coverage.png          modelled cells, cells with detections and SPAs
#   fig2_naive_vs_corrected.png copy of script 26 output (primary model)
#   fig3_atlas_status.png      copy of script 22 output (status vs II Atlas)
#   fig4a_habitat.png, fig4b_spa.png  core vs other habitat and SPA vs
#                              matched cells under the three model specifications
# Tables (results/manuscript_tables/), one CSV each:
#   tableS_detection.csv       detection coefficients (primary models)
#   detection_summary.csv      mean detection per checklist and per occupied
#                              surveyed cell-year, cells with detections
#   tableS_first_transition.csv anchored rates with and without 2017->18
#   tableS_unanchored.csv      events in cells with and without detections
#   tableS_habitat_spa_all.csv class rates, all thresholds and specifications
#   tableS_spa_balance.csv     covariate balance of the SPA matching
#   tableS_spa_sites.csv       SPAs holding the recent net change
#   tableS_atlas_effort.csv    II Atlas comparison by effort threshold
###############################################################################

library(here)
library(unmarked)
library(sf)
library(ggplot2)

source(here("R", "spatial_layers.R"))
set.seed(42)  # no stochastic steps

species <- c("otitar", "ptealc", "pteori", "tettet")
species_names <- c(otitar = "O. tarda", ptealc = "P. alchata",
                   pteori = "P. orientalis", tettet = "T. tetrax")
variants <- c(primary = "Primary model", eps_null = "Extinction constant",
              steprep_v3 = "stepRep detection")
fig_dir <- here("results", "manuscript_figures")
tab_dir <- here("results", "manuscript_tables")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)
res <- function(...) here("results", ...)
variant_dir <- function(base, v) if (v == "primary") res(base) else res(base, v)

###############################################################################
# FIGURE 1: coverage
###############################################################################

cells <- do.call(rbind, lapply(species, function(sp) {
  d <- read.csv(res("net_change_maps", paste0(sp, "_cells.csv")))
  data.frame(species = factor(species_names[sp], levels = species_names),
             x = d$x, y = d$y, detected = d$ever_detected)
}))
spa_plot <- st_simplify(load_spa_union(), dTolerance = 500)
spain <- if (requireNamespace("rnaturalearth", quietly = TRUE)) {
  st_transform(rnaturalearth::ne_countries(scale = 50, country = "Spain",
                                           returnclass = "sf"), UTM30)
} else NULL
p1 <- ggplot()
if (!is.null(spain)) p1 <- p1 + geom_sf(data = spain, fill = "grey97", colour = "grey60", linewidth = 0.2)
p1 <- p1 +
  geom_tile(data = cells[!cells$detected, ], aes(x, y), fill = "grey80",
            width = 4500, height = 4500) +
  geom_sf(data = spa_plot, fill = NA, colour = "#3D6B45", linewidth = 0.15) +
  geom_tile(data = cells[cells$detected, ], aes(x, y), fill = "#A63446",
            width = 9000, height = 9000) +
  coord_sf(crs = UTM30, xlim = c(-60000, 1130000), ylim = c(3980000, 4880000),
           expand = FALSE) +
  facet_wrap(~ species) +
  labs(x = NULL, y = NULL,
       caption = paste("Grey: surveyed cells (>= 3 eBird checklists, 2017-2023).",
                       "Red: cells with >= 1 detection (enlarged for visibility).",
                       "Green outlines: Special Protection Areas.")) +
  theme_bw(base_size = 9) +
  theme(strip.text = element_text(face = "italic"), axis.text = element_blank(),
        axis.ticks = element_blank(), panel.grid = element_blank())
ggsave(file.path(fig_dir, "fig1_coverage.png"), p1, width = 8, height = 7, dpi = 300)

file.copy(res("naive_vs_corrected_anchored", "fig_naive_vs_corrected.png"),
          file.path(fig_dir, "fig2_naive_vs_corrected.png"), overwrite = TRUE)
file.copy(res("atlas10", "fig_status_vs_atlas2.png"),
          file.path(fig_dir, "fig3_atlas_status.png"), overwrite = TRUE)

###############################################################################
# FIGURE 4: extinction contrasts across specifications
###############################################################################

rates <- do.call(rbind, lapply(names(variants), function(v) {
  r <- read.csv(file.path(variant_dir("core_spa", v), "class_rates.csv"))
  r$variant <- v
  r[r$period == "both", ]
}))
matched <- do.call(rbind, lapply(names(variants), function(v) {
  m <- read.csv(file.path(variant_dir("core_spa", v), "spa_matched.csv"))
  m$variant <- v
  m[m$period == "both" & m$spec == "base", ]
}))
keep_sp <- c("otitar", "pteori", "tettet")   # P. alchata: extinction not identifiable
hab <- rates[rates$species %in% keep_sp &
             rates$class %in% c("core_strict0.5", "rest_strict0.5",
                                "core_strict0.7", "rest_strict0.7"), ]
hab$threshold <- ifelse(grepl("0.5", hab$class), ">= 50 % pseudo-steppe", ">= 70 % pseudo-steppe")
hab$habitat <- ifelse(grepl("^core", hab$class), "Core", "Other")
hab$species <- factor(species_names[hab$species], levels = species_names[keep_sp])
hab$variant <- factor(variants[hab$variant], levels = variants)
p4a <- ggplot(hab, aes(x = ext_rate_med, y = variant, colour = habitat)) +
  geom_pointrange(aes(xmin = ext_rate_lo, xmax = ext_rate_hi),
                  position = position_dodge(width = 0.55), size = 0.2) +
  facet_grid(species ~ threshold) +
  scale_colour_manual(values = c(Core = "#3D6B45", Other = "#A63446"), name = NULL) +
  labs(x = "Extinction rate per occupied cell-year, 2018-2023", y = NULL,
       title = "(a) Core vs other habitat") +
  theme_bw(base_size = 9) +
  theme(strip.text.y = element_text(face = "italic"), legend.position = "bottom")

mt <- matched[matched$species %in% keep_sp, ]
mt$species <- factor(species_names[mt$species], levels = species_names[keep_sp])
mt$variant <- factor(variants[mt$variant], levels = variants)
p4b <- ggplot(mt, aes(x = diff_ext_rate_med, y = variant)) +
  geom_vline(xintercept = 0, colour = "grey50", linetype = "dashed") +
  geom_pointrange(aes(xmin = diff_lo, xmax = diff_hi,
                      shape = ifelse(point_only, "Cell resampling only", "Full uncertainty")),
                  size = 0.25) +
  scale_shape_manual(values = c("Full uncertainty" = 16, "Cell resampling only" = 1),
                     name = NULL) +
  facet_grid(species ~ .) +
  labs(x = "Extinction rate, SPA cells minus matched unprotected cells", y = NULL,
       title = "(b) Special Protection Areas") +
  theme_bw(base_size = 9) +
  theme(strip.text.y = element_text(face = "italic"), legend.position = "bottom")
ggsave(file.path(fig_dir, "fig4a_habitat.png"), p4a, width = 7, height = 5.5, dpi = 300)
ggsave(file.path(fig_dir, "fig4b_spa.png"), p4b, width = 5, height = 5.5, dpi = 300)

###############################################################################
# SUPPLEMENTARY TABLES
###############################################################################

det <- do.call(rbind, lapply(species, function(sp) {
  m <- readRDS(res(paste0(sp, "_model_object.rds")))
  cf <- coef(m); se <- sqrt(diag(vcov(m)))
  i <- grep("^p\\(", names(cf))
  data.frame(species = sp, term = sub("^p\\((.*)\\)$", "\\1", names(cf)[i]),
             estimate = round(cf[i], 3), se = round(se[i], 3), row.names = NULL)
}))
write.csv(det, file.path(tab_dir, "tableS_detection.csv"), row.names = FALSE)

# Detection summaries over realised visits: mean per-checklist detection and
# probability of at least one detection in an occupied surveyed cell-year
source(here("R", "colext_utils.R"))
det_sum <- do.call(rbind, lapply(species, function(sp) {
  m <- readRDS(res(paste0(sp, "_model_object.rds")))
  D <- colext_design(m)
  p <- colext_rates(coef(m), D)$p
  ok <- !is.na(D$y)
  by_year <- lapply(1:D$n_t, function(t) {
    cols <- ((t - 1) * D$J + 1):(t * D$J)
    list(visits = rowSums(ok[, cols, drop = FALSE]),
         p_star = 1 - apply(1 - ifelse(ok[, cols], p[, cols], 0), 1, prod))
  })
  visits <- sapply(by_year, `[[`, "visits"); p_star <- sapply(by_year, `[[`, "p_star")
  sv <- visits > 0
  data.frame(species = sp, mean_p_checklist = mean(p[ok]),
             median_visits_surveyed_cell_year = median(visits[sv]),
             mean_p_detect_occupied_cell_year = mean(p_star[sv]),
             n_cells = D$M, n_cells_detected = sum(rowSums(D$y, na.rm = TRUE) > 0))
}))
write.csv(det_sum, file.path(tab_dir, "detection_summary.csv"), row.names = FALSE)

nv <- read.csv(res("naive_vs_corrected_anchored", "summary.csv"))
nv <- nv[nv$metric %in% c("gamma", "eps", "ratio"),
         c("species", "period", "metric", "naive", "corrected_med", "corrected_lo",
           "corrected_hi", "n_col_opp", "n_ext_opp")]
write.csv(nv, file.path(tab_dir, "tableS_first_transition.csv"), row.names = FALSE)

cls <- read.csv(res("net_change_maps", "summary_by_class.csv"))
cls <- cls[cls$class %in% c("all", "det", "never_det"), ]
share <- read.csv(res("net_change_maps", "validation.csv"))[, c("species", "share_ext_never_detected")]
write.csv(merge(cls, share, by = "species"), file.path(tab_dir, "tableS_unanchored.csv"),
          row.names = FALSE)

write.csv(rates, file.path(tab_dir, "tableS_habitat_spa_all.csv"), row.names = FALSE)
bal <- do.call(rbind, lapply(names(variants), function(v) {
  b <- read.csv(file.path(variant_dir("core_spa", v), "spa_balance.csv")); b$variant <- v; b
}))
write.csv(bal[bal$variant == "primary", ], file.path(tab_dir, "tableS_spa_balance.csv"),
          row.names = FALSE)
file.copy(res("core_spa", "spa_losses_by_site.csv"),
          file.path(tab_dir, "tableS_spa_sites.csv"), overwrite = TRUE)
file.copy(res("atlas10", "status_vs_atlas2.csv"),
          file.path(tab_dir, "tableS_atlas_effort.csv"), overwrite = TRUE)

print(list.files(fig_dir)); print(list.files(tab_dir))
sessionInfo()
