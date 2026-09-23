###############################################################################
# 32_submission_figures.R
#
# Final figures for the Biological Conservation submission, drawn at print
# size (190 mm, double column) from the outputs of scripts 21, 23, 26 and 29.
# No new analyses. Script 29 keeps the working figures; this script only
# changes presentation.
#
#   Fig. 1  naive vs detection-corrected dynamics on anchored transitions
#           (extinction rate, colonisation rate, extinctions per colonisation,
#           expected net change)
#   Fig. 2  great bustard: (a) map of expected net change in recorded cells
#           with Special Protection Areas, (b) partition of the net change
#           by protection and habitat class
#   Fig. 3  (a) marginal/core extinction ratio, (b) SPA minus matched
#           extinction, three model specifications
#   Fig. S1 surveyed cells and cells with detections
#   Fig. S2 II Atlas squares classified in 2023
#   Fig. S3 expected net change maps, little bustard and black-bellied
#           sandgrouse
#   Graphical abstract
#
# P. alchata is omitted wherever the plotted quantity depends on its
# non-identifiable extinction rate.
#
# Outputs: results/manuscript_figures/submission/*.pdf, *.tif (600 dpi),
#          *.png (300 dpi, for embedding in the manuscript file)
###############################################################################

library(here)
library(sf)
library(ggplot2)

source(here("R", "spatial_layers.R"))
set.seed(42)  # no stochastic steps

out_dir <- here("results", "manuscript_figures", "submission")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
res <- function(...) here("results", ...)
variant_dir <- function(base, v) if (v == "primary") res(base) else res(base, v)

# Species order follows the text; labels are plotmath strings
sp_all <- c("otitar", "tettet", "pteori", "ptealc")
sp_est <- c("otitar", "tettet", "pteori")          # extinction identifiable
sp_common <- c(otitar = "Great bustard", tettet = "Little bustard",
               pteori = "Black-bellied sandgrouse", ptealc = "Pin-tailed sandgrouse")
sp_latin <- c(otitar = "Otis tarda", tettet = "Tetrax tetrax",
              pteori = "Pterocles orientalis", ptealc = "Pterocles alchata")
# Common names are quoted so that plotmath keeps hyphens ("Black-bellied")
sp_expr <- setNames(sprintf('"%s"~(italic("%s"))', sp_common, sp_latin), names(sp_common))
sp_expr2 <- setNames(sprintf('atop("%s", italic("%s"))', sp_common, sp_latin), names(sp_common))
variants <- c(primary = "Primary model", eps_null = "Constant extinction",
              steprep_v3 = "Steppe-rep. detection")

col_loss <- "#B2182B"; col_gain <- "#2166AC"; col_effort <- "#2C6E91"
theme_sub <- theme_bw(base_size = 8) +
  theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey95"),
        legend.key.size = unit(3, "mm"), plot.title = element_text(size = 8, face = "bold"))

save_fig <- function(p, name, w, h) {
  ggsave(file.path(out_dir, paste0(name, ".pdf")), p, width = w, height = h, units = "mm",
         device = cairo_pdf, bg = "white")
  ggsave(file.path(out_dir, paste0(name, ".tif")), p, width = w, height = h, units = "mm",
         dpi = 600, compression = "lzw", bg = "white")
  ggsave(file.path(out_dir, paste0(name, ".png")), p, width = w, height = h, units = "mm",
         dpi = 300, bg = "white")
}

# Mainland Spain outline (largest polygon) and terrestrial SPAs
spain <- st_transform(rnaturalearth::ne_countries(scale = 50, country = "Spain",
                                                  returnclass = "sf"), UTM30)
spain_parts <- st_cast(st_geometry(spain), "POLYGON")
spain_main <- st_sf(geometry = spain_parts[which.max(st_area(spain_parts))])
spa_sites <- load_spa_sites()
spa_land <- suppressWarnings(st_intersection(st_simplify(spa_sites, dTolerance = 300), spain_main))
spa_union <- suppressWarnings(st_intersection(st_simplify(load_spa_union(), dTolerance = 300), spain_main))
bb <- st_bbox(spain_main)
map_coord <- coord_sf(crs = UTM30, xlim = c(bb["xmin"], bb["xmax"]), ylim = c(bb["ymin"], bb["ymax"]),
                      expand = TRUE)
map_theme <- theme_void(base_size = 8) +
  theme(legend.position = "bottom", strip.text = element_text(size = 8, margin = margin(2, 0, 2, 0)),
        plot.background = element_rect(fill = "white", colour = NA))

hab <- read.csv(res("core_spa", "cell_habitat.csv"))
atlas2 <- load_atlas2()
atlas_col <- c(otitar = "OTITAR", ptealc = "PTEALC", pteori = "PTEORI", tettet = "TETTET")
atlas_range <- function(sp) atlas2[as.numeric(st_drop_geometry(atlas2)[[atlas_col[sp]]]) > 0, ]

###############################################################################
# FIG. 1: naive vs detection-corrected
###############################################################################

nv <- read.csv(res("naive_vs_corrected_anchored", "summary.csv"))
nv <- nv[nv$period == "recent" & nv$species %in% sp_est, ]
nv$sp <- factor(nv$species, levels = rev(sp_est))
cr <- read.csv(res("core_spa", "class_rates.csv"))
net <- cr[cr$variant == "primary" & cr$class == "det_all" & cr$period == "both" &
            cr$species %in% sp_est, ]
net$sp <- factor(net$species, levels = rev(sp_est))
ylab_sp <- scale_y_discrete(labels = function(x) parse(text = sp_expr2[x]))

rate_panel <- function(metric, title, log = FALSE, ref = NULL, show_y = FALSE) {
  d <- nv[nv$metric == metric, ]
  p <- ggplot(d, aes(y = sp)) +
    { if (!is.null(ref)) geom_vline(xintercept = ref, colour = "grey50", linetype = "dashed", linewidth = 0.3) } +
    geom_segment(aes(x = naive, xend = corrected_med, yend = sp), colour = "grey70", linewidth = 0.3) +
    geom_linerange(aes(xmin = corrected_lo, xmax = corrected_hi), linewidth = 0.5) +
    geom_point(aes(x = corrected_med, shape = "Detection-corrected"), size = 1.8) +
    geom_point(aes(x = naive, shape = "Naive"), size = 1.8) +
    scale_shape_manual(values = c("Detection-corrected" = 16, "Naive" = 4), name = NULL) +
    ylab_sp + labs(x = NULL, y = NULL, title = title) + theme_sub +
    theme(legend.position = "none")
  if (log) p <- p + scale_x_log10(breaks = c(0.5, 1, 2, 4))
  if (!show_y) p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  p
}
p1a <- rate_panel("eps", "a  Extinction rate\n", show_y = TRUE)
p1b <- rate_panel("gamma", "b  Colonisation rate\n")
p1c <- rate_panel("ratio", "c  Extinctions per\n    colonisation (log)", log = TRUE, ref = 1)
p1d <- ggplot(net, aes(y = sp)) +
  geom_vline(xintercept = 0, colour = "grey50", linetype = "dashed", linewidth = 0.3) +
  geom_linerange(aes(xmin = net_lo, xmax = net_hi), linewidth = 0.5) +
  geom_point(aes(x = net_med), size = 1.8, shape = 16) +
  ylab_sp + labs(x = NULL, y = NULL, title = "d  Expected net\n    change (cells)") + theme_sub +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
leg1 <- cowplot::get_legend(rate_panel("eps", "")+ theme(legend.position = "bottom"))
row1 <- cowplot::plot_grid(p1a, p1b, p1c, p1d, nrow = 1, align = "h", axis = "tb",
                           rel_widths = c(1.55, 1, 1, 1))
fig1 <- cowplot::plot_grid(row1, leg1, ncol = 1, rel_heights = c(1, 0.12))
save_fig(fig1, "Fig1_naive_vs_corrected", 190, 62)

###############################################################################
# FIG. 2: great bustard map and partition of the net change
###############################################################################

ot <- read.csv(res("net_change_maps", "otitar_cells.csv"))
ot <- ot[ot$ever_detected, ]
ot$spa_frac <- hab$spa_frac[match(ot$cells, hab$cells)]
ot$abs_net <- abs(ot$net_both)
ot$dir <- ifelse(ot$net_both < 0, "Loss", "Gain")
small <- ot[ot$abs_net < 0.05, ]; big <- ot[ot$abs_net >= 0.05, ]
sites <- read.csv(res("core_spa", "spa_losses_by_site.csv"))
top2 <- head(sites[sites$species == "otitar", ][order(sites$net_both[sites$species == "otitar"]), ], 2)
lab_pts <- st_point_on_surface(spa_sites[match(top2$spa_code, spa_sites$localId), ])
lab_xy <- data.frame(st_coordinates(lab_pts), lab = c("1", "2"))
p2a <- ggplot() +
  geom_sf(data = spain_main, fill = "white", colour = "grey55", linewidth = 0.2) +
  geom_sf(data = atlas_range("otitar"), fill = "grey88", colour = NA) +
  geom_sf(data = spa_union, fill = "#CFE3D1", colour = "#5E8C66", linewidth = 0.1, alpha = 0.8) +
  geom_point(data = small, aes(x, y), shape = 21, size = 0.7, fill = "white", colour = "grey40",
             stroke = 0.25) +
  geom_point(data = big[order(-big$abs_net), ], aes(x, y, size = abs_net, fill = dir), shape = 21,
             colour = "grey20", stroke = 0.25) +
  geom_text(data = lab_xy, aes(X, Y, label = lab), size = 2.6, fontface = "bold",
            nudge_x = 22000, nudge_y = 18000) +
  scale_size_area(max_size = 4.5, breaks = c(0.25, 0.5, 1), name = "Expected net change\n(cells, absolute)") +
  scale_fill_manual(values = c(Loss = col_loss, Gain = col_gain), name = NULL) +
  map_coord + map_theme +
  guides(fill = guide_legend(override.aes = list(size = 2.5), order = 1),
         size = guide_legend(order = 2)) +
  labs(title = "a") + theme(plot.title = element_text(face = "bold", size = 9))
cls <- c(det_all = "All recorded cells", spa_in = "SPA cells", spa_out = "Unprotected cells",
         core_strict0.5 = "Core (≥ 50 %)", rest_strict0.5 = "Marginal (< 50 %)",
         core_strict0.7 = "Core (≥ 70 %)", rest_strict0.7 = "Marginal (< 70 %)")
grp <- c(det_all = "All", spa_in = "Protection", spa_out = "Protection",
         core_strict0.5 = "Habitat, 50 %", rest_strict0.5 = "Habitat, 50 %",
         core_strict0.7 = "Habitat, 70 %", rest_strict0.7 = "Habitat, 70 %")
part <- cr[cr$species == "otitar" & cr$variant == "primary" & cr$period == "both" &
             cr$class %in% names(cls), ]
part$lab <- factor(sprintf("%s (%d)", cls[part$class], part$n_cells),
                   levels = rev(sprintf("%s (%d)", cls[names(cls)],
                                        part$n_cells[match(names(cls), part$class)])))
part$grp <- factor(grp[part$class], levels = unique(grp))
p2b <- ggplot(part, aes(y = lab)) +
  geom_vline(xintercept = 0, colour = "grey50", linetype = "dashed", linewidth = 0.3) +
  geom_linerange(aes(xmin = net_lo, xmax = net_hi), linewidth = 0.5) +
  geom_point(aes(x = net_med), size = 1.8, colour = col_loss) +
  facet_grid(grp ~ ., scales = "free_y", space = "free_y") +
  labs(x = "Expected net change (cells)", y = NULL, title = "b") + theme_sub +
  theme(strip.text.y = element_text(angle = 0, hjust = 0), plot.title = element_text(face = "bold", size = 9))
fig2 <- cowplot::plot_grid(p2a, p2b, nrow = 1, rel_widths = c(1.75, 1))
save_fig(fig2, "Fig2_great_bustard", 190, 105)

###############################################################################
# FIG. 3: habitat and SPA contrasts
###############################################################################

vlev <- rev(names(variants))
ctr <- do.call(rbind, lapply(names(variants), function(v) {
  x <- read.csv(file.path(variant_dir("core_spa", v), "core_contrast.csv")); x$variant <- v
  x[x$period == "both" & x$species %in% sp_est & x$definition %in% c("strict0.5", "strict0.7"), ]
}))
mt <- do.call(rbind, lapply(names(variants), function(v) {
  x <- read.csv(file.path(variant_dir("core_spa", v), "spa_matched.csv")); x$variant <- v
  x[x$period == "both" & x$species %in% sp_est & x$spec %in% c("base", "effort"), ]
}))
prep <- function(d) {
  d$sp <- factor(d$species, levels = sp_est, labels = sp_expr[sp_est])
  d$var <- factor(d$variant, levels = vlev, labels = variants[vlev])   # primary drawn on top
  d$primary <- d$variant == "primary"
  d$fix <- ifelse(d$point_only, "Coefficients fixed at point estimates", "Full uncertainty")
  d
}
ctr <- prep(ctr); mt <- prep(mt)
RL <- c(0.03, 200)
ctr$thr <- factor(ifelse(ctr$definition == "strict0.5", "≥ 50 % pseudo-steppe", "≥ 70 % pseudo-steppe"))
ctr$lo_p <- pmax(ctr$ratio_lo, RL[1]); ctr$hi_p <- pmin(ctr$ratio_hi, RL[2])
ctr$lo_cut <- ctr$ratio_lo < RL[1]; ctr$hi_cut <- ctr$ratio_hi > RL[2]
pd <- position_dodge(width = 0.6)
band <- data.frame(var = factor(variants["primary"], levels = variants[vlev]))
p3a <- ggplot(ctr, aes(y = var, colour = thr, group = thr)) +
  geom_tile(data = band, aes(y = var, x = 1), width = Inf, height = 1, fill = "grey93",
            colour = NA, inherit.aes = FALSE) +
  geom_vline(xintercept = 1, colour = "grey50", linetype = "dashed", linewidth = 0.3) +
  geom_linerange(aes(xmin = lo_p, xmax = hi_p), position = pd, linewidth = 0.5) +
  geom_point(aes(x = ratio_med, shape = fix), position = pd, size = 1.7, fill = "white") +
  geom_point(data = ctr[ctr$lo_cut, ], aes(x = RL[1]), shape = "<", size = 2.2, position = pd,
             show.legend = FALSE) +
  geom_point(data = ctr[ctr$hi_cut, ], aes(x = RL[2]), shape = ">", size = 2.2, position = pd,
             show.legend = FALSE) +
  scale_x_log10(limits = RL, breaks = c(0.1, 1, 10, 100), labels = c("0.1", "1", "10", "100")) +
  scale_colour_manual(values = c("grey55", "black"), name = NULL) +
  scale_shape_manual(values = c("Full uncertainty" = 16, "Coefficients fixed at point estimates" = 21),
                     name = NULL) +
  scale_y_discrete(limits = unname(variants[c("steprep_v3", "eps_null", "primary")])) +
  facet_grid(. ~ sp, labeller = label_parsed) +
  labs(x = "Extinction ratio, marginal / core habitat (log scale)\n← lower in marginal      higher in marginal →",
       y = NULL, title = "a") + theme_sub +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold", size = 9)) +
  guides(shape = "none")
mt$match <- factor(ifelse(mt$spec == "base", "Matched on habitat", "Matched on habitat and survey effort"))
p3b <- ggplot(mt, aes(y = var, colour = match, group = match)) +
  geom_tile(data = band, aes(y = var, x = 0), width = Inf, height = 1, fill = "grey93",
            colour = NA, inherit.aes = FALSE) +
  geom_vline(xintercept = 0, colour = "grey50", linetype = "dashed", linewidth = 0.3) +
  geom_linerange(aes(xmin = diff_lo, xmax = diff_hi), position = pd, linewidth = 0.5) +
  geom_point(aes(x = diff_ext_rate_med, shape = fix), position = pd, size = 1.7, fill = "white") +
  scale_colour_manual(values = c("black", col_effort), name = NULL) +
  scale_shape_manual(values = c("Full uncertainty" = 16, "Coefficients fixed at point estimates" = 21),
                     name = NULL) +
  scale_y_discrete(limits = unname(variants[c("steprep_v3", "eps_null", "primary")])) +
  facet_grid(. ~ sp, labeller = label_parsed) +
  labs(x = "Extinction rate, SPA cells minus matched unprotected cells\n← lower in SPAs      higher in SPAs →",
       y = NULL, title = "b") + theme_sub +
  theme(legend.position = "bottom", legend.box = "vertical", legend.spacing.x = unit(2, "mm"),
        legend.text = element_text(margin = margin(r = 6)), plot.title = element_text(face = "bold", size = 9))
fig3 <- cowplot::plot_grid(p3a, p3b, ncol = 1, align = "v", axis = "lr", rel_heights = c(1, 1.12))
save_fig(fig3, "Fig3_contrasts", 190, 130)

###############################################################################
# FIG. S1: coverage
###############################################################################

cov <- do.call(rbind, lapply(sp_all, function(sp) {
  d <- read.csv(res("net_change_maps", paste0(sp, "_cells.csv")))
  data.frame(sp = sp, x = d$x, y = d$y,
             status = ifelse(d$ever_detected, "Cell with ≥ 1 detection", "Surveyed cell"))
}))
cov$sp <- factor(cov$sp, levels = sp_all, labels = sp_expr[sp_all])
pS1 <- ggplot() +
  geom_sf(data = spain_main, fill = "white", colour = "grey55", linewidth = 0.2) +
  geom_tile(data = cov[cov$status == "Surveyed cell", ], aes(x, y, fill = status), width = 4500, height = 4500) +
  geom_sf(data = spa_union, fill = NA, colour = "#5E8C66", linewidth = 0.12) +
  geom_tile(data = cov[cov$status != "Surveyed cell", ], aes(x, y, fill = status), width = 9000, height = 9000) +
  scale_fill_manual(values = c("Surveyed cell" = "grey80", "Cell with ≥ 1 detection" = col_loss), name = NULL) +
  facet_wrap(~ sp, labeller = label_parsed) + map_coord + map_theme
save_fig(pS1, "FigS1_coverage", 190, 150)

###############################################################################
# FIG. S2: II Atlas squares classified in 2023 (P. alchata not classified)
###############################################################################

atlas2$sq <- seq_len(nrow(atlas2))
st10 <- do.call(rbind, lapply(sp_est, function(sp) {
  q <- read.csv(res("atlas10", paste0(sp, "_squares.csv")))
  q <- q[q$atlas2 & q$status_eff10 != "insufficient effort", c("sq_id", "status_eff10")]
  q$sp <- sp; q
}))
st_map <- merge(atlas2[, "sq"], st10, by.x = "sq", by.y = "sq_id")
st_map$sp <- factor(st_map$sp, levels = sp_est, labels = sp_expr[sp_est])
st_map$status <- factor(st_map$status_eff10, levels = c("retained", "uncertain", "lost"),
                        labels = c("Retained", "Uncertain", "Lost"))
rng <- do.call(rbind, lapply(sp_est, function(sp) {
  r <- atlas_range(sp)[, "geometry"]; r$sp <- factor(sp, levels = sp_est, labels = sp_expr[sp_est]); r
}))
pS2 <- ggplot() +
  geom_sf(data = spain_main, fill = "white", colour = "grey55", linewidth = 0.2) +
  geom_sf(data = rng, fill = "grey88", colour = NA) +
  geom_sf(data = st_centroid(st_map), aes(fill = status), shape = 22, size = 1.6, colour = "grey20",
          stroke = 0.2) +
  scale_fill_manual(values = c(Retained = col_gain, Uncertain = "#F2C14E", Lost = col_loss),
                    name = "Status in 2023 (II Atlas squares, ≥ 10 surveyed cell-years)") +
  facet_wrap(~ sp, nrow = 1, labeller = label_parsed) + map_coord + map_theme
save_fig(pS2, "FigS2_atlas_status", 190, 80)

###############################################################################
# FIG. S3: net change maps, little bustard and black-bellied sandgrouse
###############################################################################

oth <- do.call(rbind, lapply(c("tettet", "pteori"), function(sp) {
  d <- read.csv(res("net_change_maps", paste0(sp, "_cells.csv")))
  d <- d[d$ever_detected, ]
  data.frame(sp = sp, x = d$x, y = d$y, net = d$net_both)
}))
oth$sp <- factor(oth$sp, levels = c("tettet", "pteori"), labels = sp_expr[c("tettet", "pteori")])
oth$abs_net <- abs(oth$net); oth$dir <- ifelse(oth$net < 0, "Loss", "Gain")
rng2 <- do.call(rbind, lapply(c("tettet", "pteori"), function(sp) {
  r <- atlas_range(sp)[, "geometry"]
  r$sp <- factor(sp, levels = c("tettet", "pteori"), labels = sp_expr[c("tettet", "pteori")]); r
}))
pS3 <- ggplot() +
  geom_sf(data = spain_main, fill = "white", colour = "grey55", linewidth = 0.2) +
  geom_sf(data = rng2, fill = "grey88", colour = NA) +
  geom_sf(data = spa_union, fill = "#CFE3D1", colour = "#5E8C66", linewidth = 0.1, alpha = 0.8) +
  geom_point(data = oth[oth$abs_net < 0.05, ], aes(x, y), shape = 21, size = 0.7, fill = "white",
             colour = "grey40", stroke = 0.25) +
  geom_point(data = oth[oth$abs_net >= 0.05, ], aes(x, y, size = abs_net, fill = dir), shape = 21,
             colour = "grey20", stroke = 0.25) +
  scale_size_area(max_size = 4.5, breaks = c(0.25, 0.5, 1), name = "Expected net change\n(cells, absolute)") +
  scale_fill_manual(values = c(Loss = col_loss, Gain = col_gain), name = NULL) +
  facet_wrap(~ sp, nrow = 1, labeller = label_parsed) + map_coord + map_theme
save_fig(pS3, "FigS3_net_change_other_species", 190, 90)

###############################################################################
# GRAPHICAL ABSTRACT (numbers as in the abstract)
###############################################################################

ga1 <- ggplot(data.frame(k = factor(c("Raw records", "Detection-corrected"),
                                    levels = c("Raw records", "Detection-corrected")),
                         m = c(nv$naive[nv$species == "otitar" & nv$metric == "ratio"],
                               nv$corrected_med[nv$species == "otitar" & nv$metric == "ratio"]),
                         lo = c(NA, nv$corrected_lo[nv$species == "otitar" & nv$metric == "ratio"]),
                         hi = c(NA, nv$corrected_hi[nv$species == "otitar" & nv$metric == "ratio"])),
              aes(x = k, y = m)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_linerange(aes(ymin = lo, ymax = hi), na.rm = TRUE, linewidth = 0.6) +
  geom_point(size = 2.5, colour = col_loss) +
  labs(x = NULL, y = "Extinctions per colonisation", title = "Great bustard: about three losses\nper colonisation") +
  theme_sub
eff <- mt[mt$species == "otitar" & mt$variant == "primary" & mt$spec == "effort", ]
ga3 <- ggplot(eff, aes(y = 1, x = ratio_med)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_linerange(aes(xmin = ratio_lo, xmax = ratio_hi), linewidth = 0.6) +
  geom_point(size = 2.5, colour = col_loss) +
  scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4)) +
  labs(x = "Extinction ratio, SPA / matched cells", y = NULL,
       title = "No detectable reduction\ninside SPAs") +
  theme_sub + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
ga2 <- p2a + labs(title = "About half of the losses\ninside SPAs") + guides(size = "none", fill = "none") +
  theme(plot.title = element_text(size = 8, face = "bold"))
ga <- cowplot::plot_grid(ga1, ga2, ga3, nrow = 1, rel_widths = c(1, 1.3, 1))
save_fig(ga, "Graphical_abstract", 200, 80)

print(list.files(out_dir))
sessionInfo()
