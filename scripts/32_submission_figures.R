###############################################################################
# 32_submission_figures.R
#
# Final figures for the Biological Conservation submission, drawn at print
# size (190 mm, double column) from the outputs of scripts 21, 23, 26 and 29,
# in the style of the earlier GCB figures (one Okabe-Ito colour per species,
# open axes, bold panel letters). No new analyses. Script 29 keeps the
# working figures; this script only changes presentation.
#
#   Fig. 1  naive vs detection-corrected dynamics on anchored transitions,
#           four species (P. alchata extinction marked as not identifiable)
#           (extinction rate, colonisation rate, extinctions per colonisation,
#           expected net change)
#   Fig. 2  maps of expected net change in recorded cells for the three
#           species with identifiable extinction, with Special Protection
#           Areas, and partition of the net change by protection and habitat
#   Fig. 3  (a) marginal/core extinction ratio, (b) SPA minus matched
#           extinction, three model specifications
#   Fig. S1 surveyed cells and cells with detections
#   Fig. S2 II Atlas squares classified in 2023
#   Graphical abstract
#
# P. alchata is omitted from Figs. 2, 3 and S2, whose quantities depend on
# its non-identifiable extinction rate; Fig. 1 shows it with that flagged.
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

# Species order follows the text. One colour per species throughout, as in
# the GCB figures (Okabe-Ito)
sp_all <- c("otitar", "tettet", "pteori", "ptealc")
sp_est <- c("otitar", "tettet", "pteori")          # extinction identifiable
sp_short <- c(otitar = "O. tarda", tettet = "T. tetrax", pteori = "P. orientalis",
              ptealc = "P. alchata")
sp_col <- c(otitar = "#D55E00", tettet = "#CC79A7", pteori = "#009E73", ptealc = "#0072B2")
sp_col_lab <- setNames(sp_col, sp_short[names(sp_col)])
variants <- c(primary = "Primary model", eps_null = "Constant extinction",
              steprep_v3 = "Steppe-rep. detection")
vlev <- c("steprep_v3", "eps_null", "primary")     # bottom to top

theme_gcb <- theme_minimal(base_size = 8) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
        axis.line = element_line(colour = "grey35", linewidth = 0.3),
        axis.ticks = element_line(colour = "grey35", linewidth = 0.3),
        axis.text = element_text(colour = "grey20"),
        plot.title = element_text(size = 8.5, margin = margin(b = 3)),
        legend.position = "bottom", legend.key.size = unit(3, "mm"),
        legend.text = element_text(margin = margin(r = 8)),
        strip.text = element_text(face = "italic", size = 8),
        panel.spacing.x = unit(4, "mm"),
        plot.background = element_rect(fill = "white", colour = NA))
italic_y <- theme(axis.text.y = element_text(face = "italic", colour = "grey20"))
ptitle <- function(letter, text) bquote(bold(.(letter)) ~~ .(text))
errbar <- function(...) geom_errorbar(..., orientation = "y", width = 0.25, linewidth = 0.45)

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
spa_union <- suppressWarnings(st_intersection(st_simplify(load_spa_union(), dTolerance = 300), spain_main))
bb <- st_bbox(spain_main)
map_coord <- coord_sf(crs = UTM30, xlim = c(bb["xmin"], bb["xmax"]), ylim = c(bb["ymin"], bb["ymax"]),
                      expand = TRUE)
map_theme <- theme_void(base_size = 8) +
  theme(legend.position = "bottom", plot.title = element_text(size = 8.5, margin = margin(b = 2)),
        strip.text = element_text(face = "italic", size = 8.5, margin = margin(2, 0, 2, 0)),
        plot.background = element_rect(fill = "white", colour = NA))
col_spa_fill <- "#DDEBDD"; col_spa_line <- "#7FA67F"

atlas2 <- load_atlas2()
atlas_col <- c(otitar = "OTITAR", ptealc = "PTEALC", pteori = "PTEORI", tettet = "TETTET")
atlas_range <- function(sp) atlas2[as.numeric(st_drop_geometry(atlas2)[[atlas_col[sp]]]) > 0, ]
cr <- read.csv(res("core_spa", "class_rates.csv"))
cr <- cr[cr$variant == "primary" & !cr$point_only & cr$period == "both" & cr$species %in% sp_est, ]
cr$sp <- factor(sp_short[cr$species], levels = sp_short[sp_est])

###############################################################################
# FIG. 1: naive vs detection-corrected
###############################################################################

# All four species. P. alchata extinction is not identifiable (posterior
# spanning several orders of magnitude, median at the boundary), so its
# corrected extinction and ratio are drawn as a grey interval without a
# median, and its net change is marked "n.i."
nv <- read.csv(res("naive_vs_corrected_anchored", "summary.csv"))
nv <- nv[nv$period == "recent" & nv$species %in% sp_all, ]
nv$sp <- factor(sp_short[nv$species], levels = rev(sp_short[sp_all]))
nv <- transform(nv, med = corrected_med, lo = corrected_lo, hi = corrected_hi)
nv$ni <- nv$species == "ptealc" & nv$metric %in% c("eps", "ratio")
net <- cr[cr$class == "det_all", ]
net$sp <- factor(as.character(net$sp), levels = rev(sp_short[sp_all]))
net <- transform(net, med = net_med, lo = net_lo, hi = net_hi, ni = FALSE)
RL1 <- c(0.3, 6)                                   # x limits of the ratio panel

rate_panel <- function(d, letter, text, xlab, log = FALSE, ref = NULL, naive = TRUE, show_y = FALSE,
                       ni_label = NULL, sp_levels = sp_all) {
  if (log) d <- transform(d, lo_cut = lo < RL1[1], lo = pmax(lo, RL1[1]), hi = pmin(hi, RL1[2]))
  est <- d[!d$ni, ]; nid <- d[d$ni, ]
  p <- ggplot(est, aes(y = sp, colour = sp))
  if (!is.null(ref)) p <- p + geom_vline(xintercept = ref, colour = "grey50", linetype = "dashed",
                                         linewidth = 0.3)
  if (naive) p <- p + geom_segment(aes(x = naive, xend = med, yend = sp), colour = "grey65",
                                   linetype = "dotted", linewidth = 0.35)
  p <- p + errbar(aes(xmin = lo, xmax = hi)) +
    geom_point(aes(x = med, shape = "Detection-corrected"), size = 2.4)
  if (nrow(nid) > 0) {
    p <- p + errbar(data = nid, aes(xmin = lo, xmax = hi), colour = "grey75") +
      geom_label(data = nid, aes(x = sqrt_or_mid(lo, hi, log), label = "not identifiable"),
                 colour = "grey40", fill = "white", size = 2.2, vjust = -0.5, linewidth = 0,
                 label.padding = unit(0.3, "mm"))
    if (log) p <- p + geom_point(data = nid[nid$lo_cut, ], aes(x = RL1[1]), shape = "<", size = 2.2,
                                 colour = "grey60")
  }
  if (naive) p <- p + geom_point(data = d, aes(x = naive, shape = "Naive"), size = 2.4, stroke = 0.9)
  if (!is.null(ni_label)) p <- p + geom_label(data = ni_label, aes(x = x, y = sp, label = "n.i."),
                                              colour = "grey40", fill = "white", size = 2.4,
                                              linewidth = 0, inherit.aes = FALSE)
  p <- p + scale_colour_manual(values = sp_col_lab, guide = "none") +
    scale_shape_manual(values = c("Naive" = 4, "Detection-corrected" = 16), name = NULL) +
    scale_y_discrete(limits = rev(unname(sp_short[sp_levels]))) +
    labs(x = xlab, y = NULL, title = ptitle(letter, text)) + theme_gcb + italic_y +
    theme(legend.position = "none")
  if (log) p <- p + scale_x_log10(limits = RL1, breaks = c(0.5, 1, 2, 4))
  if (!show_y) p <- p + theme(axis.text.y = element_blank())
  p
}
# label position: middle of the interval (geometric middle on a log axis)
sqrt_or_mid <- function(lo, hi, log) if (log) sqrt(lo * hi) else (lo + hi) / 2
p1a <- rate_panel(nv[nv$metric == "eps", ], "A", "Extinction", "Per occupied cell-year", show_y = TRUE)
p1b <- rate_panel(nv[nv$metric == "gamma", ], "B", "Colonisation", "Per unoccupied cell-year")
p1c <- rate_panel(nv[nv$metric == "ratio", ], "C", "Extinctions per colonisation", "Ratio (log scale)",
                  log = TRUE, ref = 1)
p1d <- rate_panel(net, "D", "Net change", "Expected change (cells)", ref = 0, naive = FALSE,
                  ni_label = data.frame(x = 0, sp = sp_short[["ptealc"]]))
leg1 <- cowplot::get_plot_component(
  rate_panel(nv[nv$metric == "gamma", ], "A", "", "") +
    guides(shape = guide_legend(override.aes = list(colour = "grey30"))) +
    theme(legend.position = "bottom"), "guide-box-bottom")
row1 <- cowplot::plot_grid(p1a, p1b, p1c, p1d, nrow = 1, align = "h", axis = "tb",
                           rel_widths = c(1.3, 1, 1, 1))
fig1 <- cowplot::plot_grid(row1, leg1, ncol = 1, rel_heights = c(1, 0.1))
save_fig(fig1, "Fig1_naive_vs_corrected", 190, 72)

###############################################################################
# FIG. 2: where occupied cells are lost, three species, and partition of the
# net change by protection and habitat
###############################################################################

# The two SPAs with the largest great bustard losses, labelled on panel A
sites <- read.csv(res("core_spa", "spa_losses_by_site.csv"))
ot_sites <- sites[sites$species == "otitar", ]
top2 <- head(ot_sites[order(ot_sites$net_both), ], 2)
lab_pts <- st_point_on_surface(spa_sites[match(top2$spa_code, spa_sites$localId), ])
lab_xy <- data.frame(st_coordinates(lab_pts), lab = c("1", "2"))

map_sp <- function(sp, letter) {
  d <- read.csv(res("net_change_maps", paste0(sp, "_cells.csv")))
  d <- d[d$ever_detected, ]
  d$abs_net <- abs(d$net_both)
  d$dir <- factor(ifelse(d$net_both < 0, "Loss", "Gain"), levels = c("Loss", "Gain"))
  small <- d[d$abs_net < 0.05, ]
  big <- d[d$abs_net >= 0.05, ]
  big <- big[order(-big$abs_net), ]
  ggplot() +
    geom_sf(data = spain_main, fill = "white", colour = "grey45", linewidth = 0.2) +
    geom_sf(data = atlas_range(sp), fill = "grey90", colour = NA) +
    geom_sf(data = spa_union, fill = col_spa_fill, colour = col_spa_line, linewidth = 0.08) +
    geom_point(data = small, aes(x, y), shape = 21, size = 0.6, fill = "white", colour = "grey45",
               stroke = 0.2) +
    geom_point(data = big, aes(x, y, size = abs_net, shape = dir), colour = sp_col[sp],
               stroke = 0.7) +
    scale_shape_manual(values = c(Loss = 16, Gain = 1), drop = FALSE, name = NULL) +
    scale_size_area(max_size = 4, limits = c(0, 1.1), breaks = c(0.25, 0.5, 1),
                    name = "Absolute net change (cells)") +
    labs(title = bquote(bold(.(letter)) ~~ italic(.(unname(sp_short[sp]))))) +
    map_coord + map_theme
}
m2 <- lapply(seq_along(sp_est), function(i) map_sp(sp_est[i], LETTERS[i]))
ga_map <- m2[[1]] + theme(legend.position = "none")   # unlabelled copy for the graphical abstract
m2[[1]] <- m2[[1]] + geom_text(data = lab_xy, aes(X, Y, label = lab), size = 2.6, fontface = "bold",
                               nudge_x = 24000, nudge_y = 20000)
leg_map <- cowplot::get_plot_component(
  m2[[1]] + guides(shape = guide_legend(override.aes = list(size = 2.5, colour = "grey30"), order = 1),
                   size = guide_legend(override.aes = list(colour = "grey30"), order = 2)),
  "guide-box-bottom")
m2 <- lapply(m2, function(p) p + theme(legend.position = "none"))

cls <- c(det_all = "All recorded cells", spa_in = "SPA cells", spa_out = "Unprotected cells",
         core_strict0.5 = "Core (≥ 50 %)", rest_strict0.5 = "Marginal (< 50 %)",
         core_strict0.7 = "Core (≥ 70 %)", rest_strict0.7 = "Marginal (< 70 %)")
grp <- c(det_all = "All", spa_in = "Protection", spa_out = "Protection",
         core_strict0.5 = "Habitat", rest_strict0.5 = "Habitat",
         core_strict0.7 = "Habitat", rest_strict0.7 = "Habitat")
part <- cr[cr$class %in% names(cls), ]
part$lab <- factor(cls[part$class], levels = rev(cls))
part$grp <- factor(grp[part$class], levels = unique(grp))
part$sp <- factor(as.character(part$sp), levels = rev(sp_short[sp_est]))   # O. tarda on top
pd2 <- position_dodge(width = 0.7)
p2d <- ggplot(part, aes(y = lab, colour = sp, group = sp)) +
  geom_vline(xintercept = 0, colour = "grey50", linetype = "dashed", linewidth = 0.3) +
  geom_errorbar(aes(xmin = net_lo, xmax = net_hi), orientation = "y", width = 0.35, linewidth = 0.4,
                position = pd2) +
  geom_point(aes(x = net_med), size = 1.8, position = pd2) +
  scale_colour_manual(values = sp_col_lab, breaks = unname(sp_short[sp_est]), name = NULL) +
  facet_grid(grp ~ ., scales = "free_y", space = "free_y") +
  labs(x = "Expected net change (cells)", y = NULL, title = ptitle("D", "Net change by class of cell")) +
  theme_gcb +
  theme(strip.text.y = element_blank(), legend.text = element_text(face = "italic", margin = margin(r = 8)),
        panel.spacing.y = unit(2, "mm"), legend.margin = margin(0, 0, 0, 0))
fig2 <- cowplot::plot_grid(
  cowplot::plot_grid(m2[[1]], m2[[2]], m2[[3]], p2d, nrow = 2),
  leg_map, ncol = 1, rel_heights = c(1, 0.06))
save_fig(fig2, "Fig2_net_change_maps", 190, 180)

###############################################################################
# FIG. 3: habitat and SPA contrasts
###############################################################################

ctr <- do.call(rbind, lapply(names(variants), function(v) {
  x <- read.csv(file.path(variant_dir("core_spa", v), "core_contrast.csv")); x$variant <- v
  x[x$period == "both" & x$species %in% sp_est & x$definition %in% c("strict0.5", "strict0.7"), ]
}))
mt <- do.call(rbind, lapply(names(variants), function(v) {
  x <- read.csv(file.path(variant_dir("core_spa", v), "spa_matched.csv")); x$variant <- v
  x[x$period == "both" & x$species %in% sp_est & x$spec %in% c("base", "effort"), ]
}))
prep <- function(d) {
  d$sp <- factor(sp_short[d$species], levels = sp_short[sp_est])
  d$var <- factor(variants[d$variant], levels = unname(variants[vlev]))
  # open symbols: coefficients fixed at their point estimates
  d$pt_fill <- ifelse(d$point_only, "white", sp_col[d$species])
  d
}
ctr <- prep(ctr); mt <- prep(mt)
RL <- c(0.03, 200)
ctr$grp <- factor(ifelse(ctr$definition == "strict0.5", "Core ≥ 50 % pseudo-steppe",
                         "Core ≥ 70 % pseudo-steppe"))
ctr$lo_p <- pmax(ctr$ratio_lo, RL[1]); ctr$hi_p <- pmin(ctr$ratio_hi, RL[2])
mt$grp <- factor(ifelse(mt$spec == "base", "Matched on habitat and initial occupancy",
                        "Also matched on survey effort"),
                 levels = c("Matched on habitat and initial occupancy", "Also matched on survey effort"))
band <- data.frame(var = factor(variants["primary"], levels = unname(variants[vlev])))
pd <- position_dodge(width = 0.55)
contrast_panel <- function(d, x, lo, hi, ref, letter, text, xlab, shapes) {
  ggplot(d, aes(y = var, colour = sp, group = grp)) +
    geom_tile(data = band, aes(y = var, x = ref), width = Inf, height = 1, fill = "grey95",
              colour = NA, inherit.aes = FALSE) +
    geom_vline(xintercept = ref, colour = "grey50", linetype = "dashed", linewidth = 0.3) +
    geom_errorbar(aes(xmin = .data[[lo]], xmax = .data[[hi]]), orientation = "y", width = 0.3,
                  linewidth = 0.45, position = pd) +
    geom_point(aes(x = .data[[x]], shape = grp, fill = pt_fill), size = 2, stroke = 0.6, position = pd) +
    scale_colour_manual(values = sp_col_lab, guide = "none") +
    scale_fill_identity() +
    scale_shape_manual(values = shapes, name = NULL) +
    guides(shape = guide_legend(override.aes = list(colour = "grey30", fill = "grey30"))) +
    scale_y_discrete(limits = unname(variants[vlev])) +
    facet_grid(. ~ sp) +
    labs(x = xlab, y = NULL, title = ptitle(letter, text)) + theme_gcb
}
p3a <- contrast_panel(ctr, "ratio_med", "lo_p", "hi_p", 1, "A", "Marginal vs core habitat",
                      "Extinction ratio, marginal / core (log scale)", c(21, 24)) +
  geom_point(data = ctr[ctr$ratio_lo < RL[1], ], aes(x = RL[1]), shape = "<", size = 2.2,
             position = pd, show.legend = FALSE) +
  geom_point(data = ctr[ctr$ratio_hi > RL[2], ], aes(x = RL[2]), shape = ">", size = 2.2,
             position = pd, show.legend = FALSE) +
  scale_x_log10(limits = RL, breaks = c(0.1, 1, 10, 100), labels = c("0.1", "1", "10", "100"))
p3b <- contrast_panel(mt, "diff_ext_rate_med", "diff_lo", "diff_hi", 0, "B",
                      "Estimated SPA\u2013control difference in local extinction",
                      "Extinction rate, SPA minus matched unprotected cells", c(21, 22))
fig3 <- cowplot::plot_grid(p3a, p3b, ncol = 1, align = "v", axis = "lr")
save_fig(fig3, "Fig3_contrasts", 190, 130)

###############################################################################
# FIG. S1: coverage
###############################################################################

cov <- do.call(rbind, lapply(sp_all, function(sp) {
  d <- read.csv(res("net_change_maps", paste0(sp, "_cells.csv")))
  data.frame(species = sp, x = d$x, y = d$y, detected = d$ever_detected)
}))
cov$sp <- factor(sp_short[cov$species], levels = sp_short[sp_all])
pS1 <- ggplot() +
  geom_sf(data = spain_main, fill = "white", colour = "grey45", linewidth = 0.2) +
  geom_tile(data = cov[!cov$detected, ], aes(x, y), fill = "grey82", width = 4500, height = 4500) +
  geom_sf(data = spa_union, fill = NA, colour = col_spa_line, linewidth = 0.12) +
  geom_tile(data = cov[cov$detected, ], aes(x, y, fill = sp), width = 9000, height = 9000) +
  scale_fill_manual(values = sp_col_lab, guide = "none") +
  facet_wrap(~ sp) + map_coord + map_theme
save_fig(pS1, "FigS1_coverage", 190, 150)

###############################################################################
# FIG. S2: II Atlas squares classified in 2023 (P. alchata not classified)
###############################################################################

atlas2$sq <- seq_len(nrow(atlas2))
st10 <- do.call(rbind, lapply(sp_est, function(sp) {
  q <- read.csv(res("atlas10", paste0(sp, "_squares.csv")))
  q <- q[q$atlas2 & q$status_eff10 != "insufficient effort", c("sq_id", "status_eff10")]
  q$species <- sp; q
}))
st_map <- merge(atlas2[, "sq"], st10, by.x = "sq", by.y = "sq_id")
st_map$sp <- factor(sp_short[st_map$species], levels = sp_short[sp_est])
st_map$status <- factor(st_map$status_eff10, levels = c("retained", "uncertain", "lost"),
                        labels = c("Retained", "Uncertain", "Lost"))
rng <- do.call(rbind, lapply(sp_est, function(sp) {
  r <- atlas_range(sp)[, "geometry"]; r$sp <- factor(sp_short[sp], levels = sp_short[sp_est]); r
}))
pS2 <- ggplot() +
  geom_sf(data = spain_main, fill = "white", colour = "grey45", linewidth = 0.2) +
  geom_sf(data = rng, fill = "grey90", colour = NA) +
  geom_sf(data = st_centroid(st_map), aes(fill = status), shape = 22, size = 1.6, colour = "grey20",
          stroke = 0.2) +
  scale_fill_manual(values = c(Retained = "#56B4E9", Uncertain = "#F0E442", Lost = "#D55E00"),
                    name = "Status in 2023 (II Atlas squares, ≥ 10 surveyed cell-years)") +
  facet_wrap(~ sp, nrow = 1) + map_coord + map_theme
save_fig(pS2, "FigS2_atlas_status", 190, 80)

###############################################################################
# GRAPHICAL ABSTRACT (numbers as in the abstract). The bird silhouettes of the
# GCB conceptual figure were added by hand and are not drawn here.
###############################################################################

ga1 <- rate_panel(nv[nv$metric == "ratio" & nv$species %in% sp_est, ], "", "",
                  "Extinctions per colonisation (log)", log = TRUE, ref = 1, show_y = TRUE,
                  sp_levels = sp_est) +
  labs(title = "Correcting for detection\nchanges the trajectories")
ga2 <- ga_map + labs(title = "Great bustard losses: about\nhalf of them inside SPAs")
eff <- mt[mt$variant == "primary" & mt$spec == "effort", ]
eff$sp <- factor(as.character(eff$sp), levels = rev(sp_short[sp_est]))
ga3 <- ggplot(eff, aes(y = sp, colour = sp)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  errbar(aes(xmin = diff_lo, xmax = diff_hi)) +
  geom_point(aes(x = diff_ext_rate_med), size = 2.4) +
  scale_colour_manual(values = sp_col_lab, guide = "none") +
  labs(x = "Extinction, SPA minus matched cells", y = NULL,
       title = "SPAs: no lower extinction for bustards;\nsandgrouse reduction not robust") +
  theme_gcb + italic_y
ga <- cowplot::plot_grid(ga1, ga2, ga3, nrow = 1, rel_widths = c(1, 1.05, 1))
save_fig(ga, "Graphical_abstract", 200, 80)

print(list.files(out_dir))
sessionInfo()
