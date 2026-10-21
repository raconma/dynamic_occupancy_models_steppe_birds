###############################################################################
# 42_spa_target_species.R
#
# SPA contrast restricted to SPAs designated for each species. The treated
# cells of script 23 are cells with >= 50 % cover of any SPA, which includes
# SPAs designated for other birds (e.g. wetlands or forest raptors). Here a
# cell is treated for a species when >= 50 % of it lies in SPAs whose Natura
# 2000 Standard Data Form lists that species with a significant population
# (site assessment A, B or C; D = non-significant). Controls are the same as
# in script 23 (cells with < 10 % cover of any SPA). The matched contrast is
# repeated with the same matching, draws and block resampling as script 23,
# and the SPA designation dates of the treated cells are summarised.
#
# Inputs:  results/{sp}_model_object.rds, colext inputs, cell_habitat.csv
#          (script 23), SPA polygons (R/spatial_layers.R), Natura 2000
#          tabular data, end 2024 rev. 1 (EEA; data-raw/data/
#          natura2000_tabular/csv/: SPECIES and NATURA2000SITES tables)
# Outputs: data/derived/natura2000_target_species_es.csv
#          results/core_spa/cell_spa_target.csv
#          results/core_spa/spa_matched_target_species.csv
###############################################################################

library(here)
library(unmarked)
library(MASS)
library(sf)

source(here("R", "colext_utils.R"))
source(here("R", "spatial_layers.R"))
source(here("R", "spa_matching.R"))
set.seed(42)

species <- c("otitar", "tettet", "pteori", "ptealc")
sci <- c(otitar = "Otis tarda", tettet = "Tetrax tetrax",
         pteori = "Pterocles orientalis", ptealc = "Pterocles alchata")
N_DRAWS <- 1000
hab_dir <- here("results", "core_spa")
tab_dir <- here("data-raw", "data", "natura2000_tabular", "csv")

# --- 1. SPAs listing each species with a significant population -------------
spp <- read.csv(file.path(tab_dir, "Natura2000_end2024_rev1_SPECIES.csv"), fileEncoding = "UTF-8-BOM")
sites <- read.csv(file.path(tab_dir, "Natura2000_end2024_rev1_NATURA2000SITES.csv"), fileEncoding = "UTF-8-BOM")
spp <- spp[spp$COUNTRY_CODE == "ES" & spp$SPECIESNAME %in% sci, ]
target <- merge(spp[, c("SITECODE", "SPECIESNAME", "POPULATION_TYPE", "POPULATION")],
                sites[sites$COUNTRY_CODE == "ES", c("SITECODE", "SITETYPE", "DATE_SPA")],
                by = "SITECODE")
dir.create(here("data", "derived"), showWarnings = FALSE)
write.csv(target, here("data", "derived", "natura2000_target_species_es.csv"), row.names = FALSE)

spa <- load_spa_sites()
hab <- read.csv(file.path(hab_dir, "cell_habitat.csv"))
tgt_path <- file.path(hab_dir, "cell_spa_target.csv")
if (!file.exists(tgt_path)) {
  polys <- wc_cell_polygons(hab$cells)
  cov <- sapply(species, function(sp) {
    codes <- unique(target$SITECODE[target$SPECIESNAME == sci[[sp]] &
                                      target$POPULATION %in% c("A", "B", "C")])
    # MITECO codes carry a "_SPA" suffix (ES0000204_SPA); EEA codes do not
    cover_by_polygons(polys, spa[sub("_SPA$", "", spa$localId) %in% codes, ])
  })
  write.csv(data.frame(cells = hab$cells, cov), tgt_path, row.names = FALSE)
}
tgt <- read.csv(tgt_path)

# --- 2. Matched contrast with treated cells in SPAs designated for the species
rows <- list()
for (sp in species) {
  message("=== ", sp)
  s <- spa_setup(sp, "primary", hab)
  frac <- tgt[[sp]][match(s$h$cells, tgt$cells)]
  stopifnot(!anyNA(frac))
  tr_all <- s$tr
  s$tr <- which(s$anch & frac >= 0.5)
  mm <- spa_matches(s, specs = c("base", "effort"))
  # designation dates of the SPAs covering the treated cells (largest share)
  codes <- unique(target$SITECODE[target$SPECIESNAME == sci[[sp]] & target$POPULATION %in% c("A", "B", "C")])
  dates <- target$DATE_SPA[match(codes, target$SITECODE)]
  draws <- draw_coefs(s$mod, N_DRAWS)
  diff_d <- ratio_d <- matrix(NA_real_, N_DRAWS, length(mm), dimnames = list(NULL, names(mm)))
  for (k in seq_len(N_DRAWS)) {
    fb <- forward_backward(colext_rates(draws[k, ], s$D), s$D)
    E <- rowSums((fb$ext * s$both)[, s$recent])
    O <- rowSums((fb$z1[, 1:(s$n_t - 1)] * s$both)[, s$recent])
    b <- block_sample(s$blk[s$tr])
    for (nm in names(mm)) {
      t_i <- s$tr[b]; c_i <- mm[[nm]][b]
      r_t <- sum(E[t_i]) / sum(O[t_i]); r_c <- sum(E[c_i]) / sum(O[c_i])
      diff_d[k, nm] <- r_t - r_c; ratio_d[k, nm] <- r_t / r_c
    }
  }
  for (nm in names(mm)) {
    sd <- summ_diff(diff_d[, nm]); sr <- summ_ratio(ratio_d[, nm])
    rows[[paste(sp, nm)]] <- data.frame(
      species = sp, spec = nm, n_treated_any_spa = length(tr_all), n_treated = length(s$tr),
      n_unique_controls = length(unique(mm[[nm]])), n_blocks_treated = length(unique(s$blk[s$tr])),
      n_target_spas = length(codes), median_date_spa = as.character(median(as.Date(dates), na.rm = TRUE)),
      smd_after_steppe = smd(s$h$steppe_strict[s$tr], s$h$steppe_strict[mm[[nm]]]),
      smd_after_lpsi = smd(s$lpsi[s$tr], s$lpsi[mm[[nm]]]),
      diff_med = sd[["med"]], diff_lo = sd[["lo"]], diff_hi = sd[["hi"]], p_spa_higher = sd[["p_pos"]],
      ratio_med = sr[["med"]], ratio_lo = sr[["lo"]], ratio_hi = sr[["hi"]])
  }
}
out <- do.call(rbind, rows); rownames(out) <- NULL
write.csv(out, file.path(hab_dir, "spa_matched_target_species.csv"), row.names = FALSE)
print(out, digits = 3)

sessionInfo()
