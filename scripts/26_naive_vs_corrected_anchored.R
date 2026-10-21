###############################################################################
# 26_naive_vs_corrected_anchored.R
#
# Naive vs detection-corrected colonisation and extinction on the same
# data-anchored transitions: cells with >= 1 detection, transitions with
# surveys in both years (script 21). Replaces the intercept-based comparison
# of script 13 (Fig. 2 of the GCB submission).
#
#   naive      from detection histories: gamma = colonisations / unoccupied
#              opportunities, epsilon = extinctions / occupied opportunities
#   corrected  posterior (forward-backward) expected events on the same
#              transitions: gamma = sum(col) / sum(P(unoccupied at t)),
#              epsilon = sum(ext) / sum(P(occupied at t))
#   event ratio  extinctions / colonisations (> 1 = net loss)
#   net        colonisations - extinctions (expected cells)
# Periods: all transitions (2017-2023), recent (2018-2023) and, as a
# sensitivity to thinly surveyed cell-years, recent transitions whose two
# years both have >= 2 checklists (recent_min2): with p ~ 0.2-0.35 per
# checklist, a single checklist misses most presences; and recent transitions
# excluding those that involve 2020 (recent_no2020), when COVID-19
# restrictions shifted checklists in recorded cells later in the breeding
# season (script 37), which could confound late-season movements with local
# extinction because the detection model has no seasonal term.
# With env STABLE=1 the script instead compares recent transitions with the
# subset whose two years differ less than twofold in checklists
# (recent_stable; robustness plan E) and writes summary_stable_effort.csv,
# leaving summary.csv untouched.
# Uncertainty for the corrected values: 1,000 multivariate-normal draws
# (point_only = TRUE when the fit's variances are not finite; see
# draw_coefs() in R/colext_utils.R).
#
# Inputs:  results/{sp}_model_object.rds
# Outputs: results/naive_vs_corrected_anchored/summary.csv
#          results/naive_vs_corrected_anchored/fig_naive_vs_corrected.png
###############################################################################

library(here)
library(unmarked)
library(MASS)
library(ggplot2)

source(here("R", "colext_utils.R"))
source(here("R", "spa_matching.R"))    # block_sample()
set.seed(42)

# Model variant: "primary" (default) or a robustness refit (scripts 28 and 41),
# e.g. Rscript scripts/23_core_spa_analysis.R trend
variant <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(variant)) variant <- "primary"

species <- c("otitar", "ptealc", "pteori", "tettet")
if (nzchar(Sys.getenv("SPECIES"))) species <- strsplit(Sys.getenv("SPECIES"), ",")[[1]]   # restricted run
species_names <- c(otitar = "O. tarda", ptealc = "P. alchata",
                   pteori = "P. orientalis", tettet = "T. tetrax")
N_DRAWS <- 1000
STABLE <- Sys.getenv("STABLE") == "1"
hab <- read.csv(here("results", "core_spa", "cell_habitat.csv"))
out_dir <- here("results", "naive_vs_corrected_anchored")
if (variant != "primary") out_dir <- file.path(out_dir, variant)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
q <- function(x, p) as.numeric(quantile(x, p, na.rm = TRUE))

rows <- list()
for (sp in species) {
  message("=== ", sp, " ===")
  mod <- readRDS(colext_model_path(sp, variant))
  D <- colext_design(mod)
  n_t <- D$n_t
  naive <- sapply(1:n_t, function(t) {
    yt <- D$y[, ((t - 1) * D$J + 1):(t * D$J), drop = FALSE]
    ifelse(rowSums(!is.na(yt)) > 0, as.integer(rowSums(yt, na.rm = TRUE) > 0),
           NA_integer_)
  })
  det <- rowSums(naive, na.rm = TRUE) > 0
  z0 <- naive[, 1:(n_t - 1)]; z1 <- naive[, 2:n_t]
  both <- !is.na(z0) & !is.na(z1)
  anchored <- both & det            # recycled row-wise over transitions
  periods <- list(all = 1:(n_t - 1), recent = 2:(n_t - 1), recent_min2 = 2:(n_t - 1),
                  recent_no2020 = setdiff(2:(n_t - 1), which(2017:2022 %in% 2019:2020)))
  stopifnot(n_t == 7)   # transitions indexed 2017->18 ... 2022->23
  # Checklists per cell-year and transitions with >= 2 in both years
  visits <- sapply(1:n_t, function(t) rowSums(!is.na(D$y[, ((t - 1) * D$J + 1):(t * D$J), drop = FALSE])))
  min2 <- visits[, 1:(n_t - 1)] >= 2 & visits[, 2:n_t] >= 2
  v0 <- visits[, 1:(n_t - 1)]; v1 <- visits[, 2:n_t]
  stable <- pmax(v0, v1) < 2 * pmin(v0, v1)      # less than twofold change in checklists
  if (STABLE) periods <- list(recent = 2:(n_t - 1), recent_stable = 2:(n_t - 1))
  mask_of <- function(p) {
    w <- anchored; w[, -periods[[p]]] <- FALSE
    if (p == "recent_min2") w <- w & min2
    if (p == "recent_stable") w <- w & stable
    w
  }

  # 50 x 50 km blocks of recorded cells, resampled with each draw so that
  # intervals also carry sampling variability among cells
  occ_tab <- read.csv(colext_input_path(sp, variant))[D$keep, ]
  h <- hab[match(occ_tab$cells, hab$cells), ]
  stopifnot(!anyNA(h$cells))
  blk <- paste(floor(h$x / 5e4), floor(h$y / 5e4))
  rec <- which(det & rowSums(anchored) > 0)   # cells with anchored transitions
  post_blk <- array(NA_real_, c(N_DRAWS, length(periods), 2),
                    dimnames = list(NULL, names(periods), c("ratio", "net")))

  cf <- coef(mod)
  draws <- draw_coefs(mod, N_DRAWS)
  point_only <- attr(draws, "point_only")
  post <- array(NA_real_, c(N_DRAWS, length(periods), 4),
                dimnames = list(NULL, names(periods), c("gamma", "eps", "ratio", "net")))

  for (k in seq_len(N_DRAWS)) {
    fb <- forward_backward(colext_rates(draws[k, ], D), D)
    occ_t <- fb$z1[, 1:(n_t - 1)]
    for (p in names(periods)) {
      w <- mask_of(p)
      E <- sum(fb$ext[w]); C <- sum(fb$col[w])
      post[k, p, ] <- c(C / sum(1 - occ_t[w]), E / sum(occ_t[w]), E / C, C - E)
      b <- rec[block_sample(blk[rec])]
      Eb <- sum((fb$ext * w)[b, ]); Cb <- sum((fb$col * w)[b, ])
      post_blk[k, p, ] <- c(Eb / Cb, Cb - Eb)
    }
  }

  for (p in names(periods)) {
    w <- mask_of(p)
    n_col <- sum(w & z0 == 0 & z1 == 1); n_col_opp <- sum(w & z0 == 0)
    n_ext <- sum(w & z0 == 1 & z1 == 0); n_ext_opp <- sum(w & z0 == 1)
    naive_v <- c(gamma = n_col / n_col_opp, eps = n_ext / n_ext_opp,
                 ratio = n_ext / n_col, net = n_col - n_ext)
    for (m in c("gamma", "eps", "ratio", "net")) {
      d <- post[, p, m]
      rows[[paste(sp, p, m)]] <- data.frame(
        species = sp, variant = variant, period = p, metric = m,
        naive = naive_v[m], corrected_med = median(d),
        corrected_lo = q(d, .025), corrected_hi = q(d, .975),
        p_corrected_gt_naive = mean(d > naive_v[m]),
        blk_lo = if (m %in% c("ratio", "net")) q(post_blk[, p, m], .025) else NA,
        blk_hi = if (m %in% c("ratio", "net")) q(post_blk[, p, m], .975) else NA,
        n_col_events = n_col, n_col_opp = n_col_opp,
        n_ext_events = n_ext, n_ext_opp = n_ext_opp, n_transitions = sum(w),
        point_only = point_only, row.names = NULL)
    }
  }
}

res <- do.call(rbind, rows); rownames(res) <- NULL
write.csv(res, file.path(out_dir, if (STABLE) "summary_stable_effort.csv" else "summary.csv"), row.names = FALSE)
if (STABLE) { print(res[res$metric %in% c("ratio", "net"), ], digits = 3); quit(save = "no") }
options(width = 200)
print(res, digits = 3)

fd <- res[res$period == "recent" & res$metric %in% c("gamma", "eps"), ]
fd$species <- factor(species_names[fd$species], levels = rev(species_names))
fd$metric <- factor(fd$metric, levels = c("gamma", "eps"),
                    labels = c("Colonisation (gamma)", "Extinction (epsilon)"))
p <- ggplot(fd, aes(y = species)) +
  geom_pointrange(aes(x = corrected_med, xmin = corrected_lo, xmax = corrected_hi,
                      shape = "Detection-corrected"), size = 0.3) +
  geom_point(aes(x = naive, shape = "Naive"), size = 2.4) +
  scale_shape_manual(values = c("Detection-corrected" = 16, "Naive" = 4), name = NULL) +
  scale_x_log10() +
  facet_wrap(~ metric, scales = "free_x") +
  labs(x = "Annual rate, 2018-2023 (log scale)", y = NULL,
       caption = paste("Same transitions for both: cells with >= 1 detection,",
                       "surveyed in both years. Bars: 95% intervals, 1,000 draws.")) +
  theme_bw(base_size = 9) +
  theme(axis.text.y = element_text(face = "italic"), legend.position = "bottom")
ggsave(file.path(out_dir, "fig_naive_vs_corrected.png"), p, width = 7, height = 3.6, dpi = 300)

sessionInfo()
