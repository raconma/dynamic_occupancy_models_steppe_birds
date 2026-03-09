###############################################################################
# fig_forest_heatmap.R
#
# 4-panel figure: A (psi forest), B (gamma forest), C (epsilon forest),
#                 D (heatmap summary of gamma/epsilon effects)
###############################################################################

library(unmarked)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)

# ── Species metadata ─────────────────────────────────────────────────────────
species <- c("otitar", "ptealc", "pteori", "tettet")
sp_labels <- c(otitar = "O. tarda", ptealc = "P. alchata",
               pteori = "P. orientalis", tettet = "T. tetrax")
sp_pch     <- c(otitar = 21, ptealc = 22, pteori = 23, tettet = 24)
sp_colours <- c(otitar = "#D55E00", ptealc = "#0072B2",
                pteori = "#009E73", tettet = "#CC79A7")

# ── Covariate short labels ──────────────────────────────────────────────────
cov_labels <- c(
  "Land_Cover_Type_1_Percent_Class_6"  = "LC6",
  "Land_Cover_Type_1_Percent_Class_7"  = "LC7",
  "Land_Cover_Type_1_Percent_Class_12" = "LC12",
  "Land_Cover_Type_1_Percent_Class_13" = "LC13",

"NDVI" = "NDVI", "pr" = "Precipitation", "tmmn" = "Min. temp.",
  "tmmx" = "Max. temp.", "bio1" = "Mean annual temp.",
  "bio2" = "Diurnal range", "tree_cover" = "Tree cover",
  "grass_cover" = "Grass cover", "topo_elev" = "Elevation",
  "topo_aspect" = "Aspect"
)

# ── Driver category colours ─────────────────────────────────────────────────
driver_cat <- c(
  bio1 = "static", bio2 = "static", tree_cover = "static",
  grass_cover = "static", topo_elev = "static", topo_aspect = "static",
  pr = "climate", tmmn = "climate", tmmx = "climate",
  Land_Cover_Type_1_Percent_Class_6  = "landcover",
  Land_Cover_Type_1_Percent_Class_7  = "landcover",
  Land_Cover_Type_1_Percent_Class_12 = "landcover",
  Land_Cover_Type_1_Percent_Class_13 = "landcover",
  NDVI = "ndvi"
)

cat_colours <- c(static = "grey60", climate = "#2166AC",
                 landcover = "#E07B39", ndvi = "#35978F")

# ── Model formulas (which covariates per species x process) ──────────────────
model_covs <- list(
  otitar  = list(psi = c("bio1","bio2","tree_cover","grass_cover","topo_elev"),
                 col = c("NDVI","pr","tmmn","tmmx"),
                 ext = c("Land_Cover_Type_1_Percent_Class_6",
                          "Land_Cover_Type_1_Percent_Class_13","tmmx")),
  ptealc  = list(psi = c("bio1","bio2","tree_cover","grass_cover","topo_aspect"),
                 col = c("pr"),
                 ext = c("pr","tmmx")),
  pteori  = list(psi = c("bio2","tree_cover","grass_cover"),
                 col = c("Land_Cover_Type_1_Percent_Class_7","NDVI","tmmn","tmmx"),
                 ext = c("Land_Cover_Type_1_Percent_Class_12","NDVI","pr")),
  tettet  = list(psi = c("bio2","tree_cover","grass_cover","topo_elev"),
                 col = c("Land_Cover_Type_1_Percent_Class_12"),
                 ext = c("Land_Cover_Type_1_Percent_Class_12"))
)

# ── Extract coefficients ─────────────────────────────────────────────────────
all_rows <- list()

for (sp in species) {
  mod <- readRDS(file.path("results", paste0(sp, "_model_object.rds")))
  b   <- coef(mod)
  se  <- sqrt(diag(vcov(mod)))

  for (proc in c("psi", "col", "ext")) {
    covs <- model_covs[[sp]][[proc]]
    for (cv in covs) {
      nm <- paste0(proc, "(", cv, ")")
      est <- b[nm]; s <- se[nm]
      z   <- abs(est / s)
      # flag huge SE in psi as non-significant
      sig <- ifelse(proc == "psi" & s > 5, FALSE, z >= 1.96)
      all_rows[[length(all_rows) + 1]] <- data.frame(
        species  = sp,
        process  = proc,
        covar    = cv,
        label    = unname(cov_labels[cv]),
        estimate = est,
        se       = s,
        sig      = sig,
        category = unname(driver_cat[cv]),
        stringsAsFactors = FALSE
      )
    }
  }
}

df <- bind_rows(all_rows)
rownames(df) <- NULL

# ── Sort covariates by mean |effect| (largest at top) per panel ──────────────
sort_covs <- function(d) {
  ord <- d %>%
    group_by(label) %>%
    summarise(mean_abs = mean(abs(estimate)), .groups = "drop") %>%
    arrange(mean_abs)
  d$label <- factor(d$label, levels = ord$label)
  d
}

df_psi <- sort_covs(filter(df, process == "psi"))
df_col <- sort_covs(filter(df, process == "col"))
df_ext <- sort_covs(filter(df, process == "ext"))

# ── Helper: forest plot ──────────────────────────────────────────────────────
make_forest <- function(d, title, annotate_ptealc_gamma = FALSE) {

  dodge <- position_dodge(width = 0.6)

  p <- ggplot(d, aes(x = estimate, y = label, group = species)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    # non-sig error bars in grey
    geom_errorbarh(data = filter(d, !sig),
                   aes(xmin = estimate - 1.96 * se, xmax = estimate + 1.96 * se),
                   height = 0, colour = "grey70", linewidth = 0.4,
                   position = dodge) +
    # sig error bars in dark
    geom_errorbarh(data = filter(d, sig),
                   aes(xmin = estimate - 1.96 * se, xmax = estimate + 1.96 * se),
                   height = 0, colour = "grey20", linewidth = 0.4,
                   position = dodge) +
    geom_point(aes(shape = species, fill = species,
                   alpha = sig), colour = "grey20",
               size = 3, stroke = 0.4,
               position = dodge) +
    scale_shape_manual(
      values = sp_pch, labels = sp_labels, name = "Species"
    ) +
    scale_fill_manual(
      values = sp_colours, labels = sp_labels, name = "Species"
    ) +
    scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.3), guide = "none") +
    labs(x = "Standardised coefficient", y = NULL, title = title) +
    theme_bw(base_size = 9) +
    theme(plot.title = element_text(face = "bold", size = 10),
          panel.grid.minor = element_blank(),
          legend.position = "bottom",
          legend.box = "horizontal",
          legend.key.size = unit(3, "mm"))

  if (annotate_ptealc_gamma) {
    p <- p + annotate("text", x = 0, y = 0.6, hjust = 0.5,
                      label = expression(italic("P. alchata")~gamma~"excluded (n = 16 events)"),
                      colour = "grey50", size = 2.5)
  }

  p
}

# ── Build panels A, B, C ────────────────────────────────────────────────────
pA <- make_forest(df_psi, expression(bold(A)~~psi[1]~"(initial occupancy)"))
pB <- make_forest(df_col, expression(bold(B)~~gamma~"(colonisation)"),
                  annotate_ptealc_gamma = TRUE)
pC <- make_forest(df_ext, expression(bold(C)~~epsilon~"(extinction)"))

# ── Panel D: heatmap ────────────────────────────────────────────────────────
# Rows = covariates from gamma/epsilon; cols = species x process
df_ge <- filter(df, process %in% c("col", "ext"))

# build full grid
all_covs_ge <- unique(df_ge$label)
proc_labels <- c(col = "gam", ext = "eps")

grid <- expand.grid(
  label   = all_covs_ge,
  species = species,
  process = c("col", "ext"),
  stringsAsFactors = FALSE
)

hm <- left_join(grid, df_ge, by = c("label", "species", "process"))

hm$col_label <- paste0(sp_labels[hm$species], " ", proc_labels[hm$process])
# order columns
col_order <- as.vector(outer(species, c("col", "ext"),
                             function(s, p) paste0(sp_labels[s], " ", proc_labels[p])))
hm$col_label <- factor(hm$col_label, levels = col_order)

# cell colour categories
hm <- hm %>%
  mutate(
    cell_cat = case_when(
      is.na(estimate)    ~ "not_in_model",
      sig & estimate > 0 ~ "sig_pos",
      sig & estimate < 0 ~ "sig_neg",
      !sig & estimate > 0 ~ "ns_pos",
      !sig & estimate <= 0 ~ "ns_neg"
    ),
    beta_text = ifelse(is.na(estimate), "", sprintf("%.2f", estimate))
  )

cell_fills <- c(sig_pos = "#B2182B", sig_neg = "#2166AC",
                ns_pos = "#FDDBC7", ns_neg = "#D1E5F0",
                not_in_model = "white")

# sort rows same as the col/ext forest panels: by mean |effect|
row_ord <- df_ge %>%
  group_by(label) %>%
  summarise(mean_abs = mean(abs(estimate)), .groups = "drop") %>%
  arrange(mean_abs)
hm$label <- factor(hm$label, levels = row_ord$label)

pD <- ggplot(hm, aes(x = col_label, y = label)) +
  geom_tile(aes(fill = cell_cat), colour = "grey80", linewidth = 0.3) +
  geom_text(aes(label = beta_text), size = 2.5) +
  scale_fill_manual(values = cell_fills,
                    labels = c(sig_pos = "Sig. +", sig_neg = "Sig. -",
                               ns_pos = "NS +", ns_neg = "NS -",
                               not_in_model = "Not in model"),
                    name = "Effect") +
  labs(x = NULL, y = NULL,
       title = expression(bold(D)~~"Coefficient summary")) +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        plot.title = element_text(face = "bold", size = 10),
        panel.grid = element_blank(),
        legend.position = "bottom",
        legend.key.size = unit(3, "mm"))

# ── Assemble ─────────────────────────────────────────────────────────────────
fig_combined <- (pA | pB | pC | pD) +
  plot_layout(widths = c(1, 1, 1, 1.2), guides = "collect") &
  theme(legend.position = "bottom")

dir.create("figs", showWarnings = FALSE)

# 1. Combined
ggsave("figs/pub_fig4_forest_heatmap_combined.png", fig_combined,
       width = 240, height = 180, units = "mm", dpi = 300)
ggsave("figs/pub_fig4_forest_heatmap_combined.pdf", fig_combined,
       width = 240, height = 180, units = "mm")

# 2. Forest only (A+B+C)
fig_forest <- (pA | pB | pC) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave("figs/pub_fig4a_forest_only.png", fig_forest,
       width = 180, height = 180, units = "mm", dpi = 300)
ggsave("figs/pub_fig4a_forest_only.pdf", fig_forest,
       width = 180, height = 180, units = "mm")

# 3. Alias for forest-only (wider)
ggsave("figs/pub_fig2_forest_3submodels.png", fig_forest,
       width = 220, height = 180, units = "mm", dpi = 300)
ggsave("figs/pub_fig2_forest_3submodels.pdf", fig_forest,
       width = 220, height = 180, units = "mm")

# 4. Heatmap only
ggsave("figs/pub_fig4b_heatmap_only.png", pD,
       width = 180, height = 180, units = "mm", dpi = 300)
ggsave("figs/pub_fig4b_heatmap_only.pdf", pD,
       width = 180, height = 180, units = "mm")

cat("Done – all figures saved to figs/\n")
