###############################################################################
# 33_spa_controls_natura2000.R
#
# Are the "unprotected" control cells of the SPA contrast (script 23) really
# outside the protected-area network? Controls are cells with < 10 % SPA
# cover, but they can lie in Natura 2000 sites designated under the Habitats
# Directive (SCIs/SACs). This script
#   (1) computes, for every model cell, the exact fraction inside any Natura
#       2000 site and inside Habitats Directive sites (SCI + SAC);
#   (2) reproduces the script 23 matching (primary model; base and effort
#       specifications; cells with anchored recent transitions) and reports
#       how many controls lie in Natura 2000;
#   (3) repeats the matched contrast with the control pool restricted to
#       cells with < 10 % cover of any Natura 2000 site, using the same
#       coefficient draws and block resamples as the original pool, so the
#       two contrasts differ only in the controls.
#
# Inputs:  results/{sp}_model_object.rds, colext inputs (colext_input_path()),
#          results/core_spa/cell_habitat.csv (script 23), Natura 2000
#          (R/spatial_layers.R)
# Outputs: results/core_spa/cell_natura2000.csv
#          results/core_spa/spa_controls_natura2000.csv
#          results/core_spa/spa_matched_outside_n2k.csv
###############################################################################

library(here)
library(unmarked)
library(MASS)
library(sf)

source(here("R", "colext_utils.R"))
source(here("R", "spatial_layers.R"))
source(here("R", "spa_matching.R"))
set.seed(42)

species <- c("otitar", "tettet", "pteori", "ptealc")   # extinction identifiable in all four
N_DRAWS <- 1000
N2K_MAX <- 0.1                               # control cells: < 10 % Natura 2000
hab_dir <- here("results", "core_spa")

###############################################################################
# 1. NATURA 2000 COVER PER CELL (cached; exact polygon intersection)
###############################################################################

hab <- read.csv(file.path(hab_dir, "cell_habitat.csv"))
n2k_path <- file.path(hab_dir, "cell_natura2000.csv")
if (!file.exists(n2k_path)) {
  polys <- wc_cell_polygons(hab$cells)
  n2k <- load_n2k_sites()
  write.csv(data.frame(cells = hab$cells, n2k_frac = cover_by_polygons(polys, n2k),
                       hd_frac = cover_by_polygons(polys, n2k[!grepl("specialProtectionArea",
                                                                     n2k$desig), ])),
            n2k_path, row.names = FALSE)
}
n2k_cov <- read.csv(n2k_path)
hab$n2k_frac <- n2k_cov$n2k_frac[match(hab$cells, n2k_cov$cells)]
hab$hd_frac <- n2k_cov$hd_frac[match(hab$cells, n2k_cov$cells)]
stopifnot(!anyNA(hab$n2k_frac))

###############################################################################
# 2-3. CONTROL COVER AND CONTRAST WITH CONTROLS OUTSIDE NATURA 2000
###############################################################################

desc_rows <- list(); contrast_rows <- list()
for (sp in species) {
  message("=== ", sp)
  s <- spa_setup(sp, "primary", hab)
  pools <- list(all = s$ctl, outside_n2k = s$ctl[s$h$n2k_frac[s$ctl] < N2K_MAX])
  matches <- list()
  for (pl in names(pools)) {
    mm <- spa_matches(s, specs = c("base", "effort"), ctl = pools[[pl]])
    for (nm in names(mm)) {
      m <- mm[[nm]]; ctl <- pools[[pl]]
      matches[[paste(pl, nm)]] <- m
      desc_rows[[paste(sp, pl, nm)]] <- data.frame(
        species = sp, pool = pl, spec = nm, n_treated = length(s$tr),
        n_pool = length(ctl), n_pool_in_n2k = sum(s$h$n2k_frac[ctl] >= N2K_MAX),
        n_pool_in_hd = sum(s$h$hd_frac[ctl] >= N2K_MAX),
        n_unique_controls = length(unique(m)),
        n_unique_controls_in_n2k = sum(s$h$n2k_frac[unique(m)] >= N2K_MAX),
        share_pairs_control_in_n2k = mean(s$h$n2k_frac[m] >= N2K_MAX),
        smd_after_steppe = smd(s$h$steppe_strict[s$tr], s$h$steppe_strict[m]),
        smd_after_lpsi = smd(s$lpsi[s$tr], s$lpsi[m]),
        smd_after_visits = smd(s$visits[s$tr], s$visits[m]))
    }
  }
  draws <- draw_coefs(s$mod, N_DRAWS)
  point_only <- attr(draws, "point_only")
  diff_d <- ratio_d <- matrix(NA_real_, N_DRAWS, length(matches),
                              dimnames = list(NULL, names(matches)))
  for (k in seq_len(N_DRAWS)) {
    fb <- forward_backward(colext_rates(draws[k, ], s$D), s$D)
    E <- rowSums((fb$ext * s$both)[, s$recent])
    O <- rowSums((fb$z1[, 1:(s$n_t - 1)] * s$both)[, s$recent])
    b <- block_sample(s$blk[s$tr])
    for (mm in names(matches)) {
      t_i <- s$tr[b]; c_i <- matches[[mm]][b]
      r_t <- sum(E[t_i]) / sum(O[t_i]); r_c <- sum(E[c_i]) / sum(O[c_i])
      diff_d[k, mm] <- r_t - r_c; ratio_d[k, mm] <- r_t / r_c
    }
  }
  for (mm in names(matches)) {
    sd <- summ_diff(diff_d[, mm]); sr <- summ_ratio(ratio_d[, mm])
    contrast_rows[[paste(sp, mm)]] <- data.frame(
      species = sp, pool = sub(" .*", "", mm), spec = sub(".* ", "", mm),
      point_only = point_only, n_undefined = sd[["n_undefined"]],
      diff_med = sd[["med"]], diff_lo = sd[["lo"]], diff_hi = sd[["hi"]],
      p_spa_higher = sd[["p_pos"]],
      ratio_med = sr[["med"]], ratio_lo = sr[["lo"]], ratio_hi = sr[["hi"]])
  }
}

desc_df <- do.call(rbind, desc_rows); rownames(desc_df) <- NULL
contrast_df <- do.call(rbind, contrast_rows); rownames(contrast_df) <- NULL
write.csv(desc_df, file.path(hab_dir, "spa_controls_natura2000.csv"), row.names = FALSE)
write.csv(contrast_df, file.path(hab_dir, "spa_matched_outside_n2k.csv"), row.names = FALSE)
print(desc_df, digits = 3)
print(contrast_df, digits = 3)

sessionInfo()
