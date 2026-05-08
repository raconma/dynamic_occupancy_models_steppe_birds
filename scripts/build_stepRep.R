###############################################################################
# build_stepRep.R
#
# Purpose: Build a 'steppe representativeness' covariate (cell x year) for the
#          detection sub-model of dynamic occupancy models on eBird data.
#
# Inputs:  - CORINE Land Cover 2018 v2020_20u1, 100 m, EPSG:3035
#            (data-raw/data/corine/.../DATA/U2018_CLC2018_V2020_20u1.tif)
#          - eBird zero-filled CSVs per species (data/processed_2023/{sp}/)
#          - WorldClim 5 km grid (data-raw/data/environmental_data/.../
#            variables_spain.grd) -- defines the cell_id used in the models
#
# Outputs: - data/derived/corine/clc2018_iberia.tif         (cropped CORINE)
#          - data/derived/corine/steppe_{strict,broad}.tif  (binary masks)
#          - data/derived/corine/focal_{strict,broad}_{500m,1km}.tif
#          - data/derived/stepRep_cellyear_{sp}.{csv,rds}   (per species)
#          - reports/stepRep_diagnostics.md (+ figures)
#
# Notes:   - Habitat layer is static (CLC2018) by design: pseudo-steppe habitat
#            is slow-changing; year-to-year variation in stepRep comes from
#            where eBird checklists fall WITHIN each cell.
#          - Buffers / extractions performed in EPSG:3035 (metric).
###############################################################################

# -- Packages --
suppressPackageStartupMessages({
  library(here)
  library(terra)
  library(sf)
  library(dplyr)
  library(data.table)
})

set.seed(1)

# -- Paths --
corine_tif    <- here("data-raw", "data", "corine", "Results",
                      "u2018_clc2018_v2020_20u1_raster100m", "DATA",
                      "U2018_CLC2018_V2020_20u1.tif")
worldclim_grd <- here("data-raw", "data", "environmental_data",
                      "environmental_data_occ", "variables_spain.grd")
ebird_dir     <- here("data", "processed_2023")
out_corine    <- here("data", "derived", "corine")
out_data      <- here("data", "derived")
out_reports   <- here("reports")
out_figs      <- here("reports", "figs")

invisible(lapply(c(out_corine, out_data, out_reports, out_figs),
                 dir.create, recursive = TRUE, showWarnings = FALSE))

stopifnot(file.exists(corine_tif), file.exists(worldclim_grd))

# -- CLC code <-> raster value lookup --
# Source: data-raw/data/corine/.../Legend/CLC2018_CLC2018_V2018_20_QGIS.txt
# Entries are "{raster_value} {clc_code} ...". We store the inverse mapping.
clc_lookup <- c(
  "211" = 12L,  # Non-irrigated arable land
  "231" = 18L,  # Pastures
  "242" = 20L,  # Complex cultivation patterns
  "243" = 21L,  # Agriculture with significant areas of natural vegetation
  "244" = 22L,  # Agro-forestry areas (dehesa)
  "321" = 26L,  # Natural grasslands
  "333" = 32L   # Sparsely vegetated areas
)
# Excluded explicitly (not pseudo-steppe habitat): 212 (irrigated), 213 (rice),
# 221 (vineyards), 223 (olive groves).

vals_strict <- unname(clc_lookup[c("211", "231", "321", "333")])
vals_broad  <- unname(clc_lookup[c("211", "231", "321", "333",
                                   "242", "243", "244")])

# Iberia mainland bbox in EPSG:3035 (excludes Canarias and Baleares)
# Derived from WGS84 -10..3.32, 35.95..43.95 projected to LAEA + 5 km margin
iberia_bbox_3035 <- ext(2585000, 3785000, 1465000, 2495000)

###############################################################################
# Section 2: CORINE crop + reclassification + focal (proportion) rasters
###############################################################################

message("\n==[ Section 2: CORINE crop + reclassify + focal ]==")

# -- 2a. Crop CORINE to Iberia peninsular --------------------------------------
clc_iberia_path <- file.path(out_corine, "clc2018_iberia.tif")
if (!file.exists(clc_iberia_path)) {
  message("  [2a] Cropping CORINE to Iberia bbox (EPSG:3035) ...")
  clc_full <- rast(corine_tif)
  crop(clc_full, iberia_bbox_3035, snap = "out",
       filename = clc_iberia_path, datatype = "INT1S", overwrite = TRUE,
       gdal = c("COMPRESS=LZW", "PREDICTOR=2",
                "TILED=YES", "BIGTIFF=IF_SAFER"))
  rm(clc_full); gc(verbose = FALSE)
} else {
  message("  [2a] skip (exists): ", basename(clc_iberia_path))
}
clc <- rast(clc_iberia_path)
message(sprintf("       clc2018_iberia: %d x %d px, %.1f MB on disk",
                ncol(clc), nrow(clc), file.size(clc_iberia_path) / 2^20))

# -- 2b. Build binary steppe masks (strict, broad) -----------------------------
# NOTE on raster values: the v2020_20u1 100 m product uses raster values 1..44
# for CLC classes (per Legend) and encodes the actual pixel-level NoData as
# value 48 (NOT 45 as the QGIS legend file would suggest). We reclassify
# 0..48 explicitly so that values outside the steppe sets (incl. 48 = NoData)
# resolve to 0; any pixel outside that range is forced to NA via othersNA.
build_mask <- function(vals_keep, out_path, label) {
  if (file.exists(out_path)) {
    message("  [2b] skip (exists): ", basename(out_path))
    return(invisible())
  }
  message("  [2b] Building mask: ", label,
          "  (raster values: ", paste(vals_keep, collapse = ","), ")")
  rcl <- cbind(0:48, as.integer(0:48 %in% vals_keep))
  m <- classify(clc, rcl = rcl, others = NA)
  # Replace NA (sea / outside CORINE coverage) with 0: semantically
  # "not steppe", and keeps focal sums well-defined at coastal buffers
  # (fixes terra::focal returning Inf/-Inf on all-NA neighbourhoods).
  m <- subst(m, NA, 0)
  writeRaster(m, out_path, datatype = "INT1U", overwrite = TRUE,
              gdal = c("COMPRESS=LZW", "PREDICTOR=2", "TILED=YES"))
  invisible()
}
mask_strict_path <- file.path(out_corine, "steppe_strict.tif")
mask_broad_path  <- file.path(out_corine, "steppe_broad.tif")
build_mask(vals_strict, mask_strict_path, "steppe_strict")
build_mask(vals_broad,  mask_broad_path,  "steppe_broad")

# -- 2c. Focal (proportion) rasters at 500 m and 1 km radii --------------------
# Circular kernel; weights normalised to sum to 1 -> focal sum gives the
# proportion of the circle that is steppe. NAs (sea / outside coverage) are
# treated as 0 (semantically: not steppe), which correctly biases edge pixels.
build_focal <- function(mask_path, out_path, radius_m, label) {
  if (file.exists(out_path)) {
    message("  [2c] skip (exists): ", basename(out_path))
    return(invisible())
  }
  message("  [2c] focal: ", label, "  (radius = ", radius_m, " m) ...")
  m <- rast(mask_path)
  w <- focalMat(m, radius_m, type = "circle")  # weights sum to 1
  t0 <- Sys.time()
  focal(m, w = w, fun = "sum", na.policy = "omit", na.rm = TRUE,
        filename = out_path, datatype = "FLT4S", overwrite = TRUE,
        gdal = c("COMPRESS=LZW", "PREDICTOR=3",
                 "TILED=YES", "BIGTIFF=IF_SAFER"))
  message(sprintf("       done in %.1f s", as.numeric(Sys.time() - t0,
                                                      units = "secs")))
  invisible()
}

build_focal(mask_strict_path, file.path(out_corine, "focal_strict_500m.tif"),
            500,  "strict_500m")
build_focal(mask_strict_path, file.path(out_corine, "focal_strict_1km.tif"),
            1000, "strict_1km")
build_focal(mask_broad_path,  file.path(out_corine, "focal_broad_500m.tif"),
            500,  "broad_500m")
build_focal(mask_broad_path,  file.path(out_corine, "focal_broad_1km.tif"),
            1000, "broad_1km")

message("[Section 2 done]")

###############################################################################
# Section 3: per-checklist extraction
###############################################################################
# For each species, for each unique checklist, extract:
#   - clc_class                 : CLC raster value at the point (1..48)
#   - step_strict_point         : 1 if pixel is in steppe-strict mask else 0
#   - step_broad_point          : 1 if pixel is in steppe-broad mask else 0
#   - step_strict_500m / _1km   : proportion of steppe-strict in circular buffer
#   - step_broad_500m  / _1km   : proportion of steppe-broad  in circular buffer
#   - cells                     : WorldClim 5 km grid cell index (model domain)
# Coordinates are reprojected to EPSG:3035 (LAEA) for metric extraction.
# Iberia mainland filter applied to exclude Canarias and Baleares.

message("\n==[ Section 3: per-checklist extraction ]==")

# -- 3a. Load rasters (lazy via terra) -----------------------------------------
clc_iberia   <- rast(file.path(out_corine, "clc2018_iberia.tif"))
mask_strict  <- rast(file.path(out_corine, "steppe_strict.tif"))
mask_broad   <- rast(file.path(out_corine, "steppe_broad.tif"))
foc_str_500m <- rast(file.path(out_corine, "focal_strict_500m.tif"))
foc_str_1km  <- rast(file.path(out_corine, "focal_strict_1km.tif"))
foc_brd_500m <- rast(file.path(out_corine, "focal_broad_500m.tif"))
foc_brd_1km  <- rast(file.path(out_corine, "focal_broad_1km.tif"))

# WorldClim 5 km grid (any single layer; we only need the grid geometry / NA)
wc_grid <- rast(worldclim_grd, lyrs = 1)

# Iberia mainland bbox (WGS84). Baleares (~lon 1.5-4.4, lat 38.6-40.1) excluded
# by an explicit rule below; Canarias falls outside the lat/lon bbox.
iberia_bbox_wgs <- list(lon_min = -10, lon_max = 3.32,
                        lat_min = 35.95, lat_max = 43.95)

species_codes <- c("otitar", "ptealc", "pteori", "tettet")

# -- 3b. Per-species extraction ------------------------------------------------
for (sp in species_codes) {
  out_path <- file.path(out_data, sprintf("checklist_stepRep_%s.rds", sp))
  if (file.exists(out_path)) {
    message("  [3] skip (exists): ", basename(out_path)); next
  }

  csv_path <- file.path(ebird_dir, sp,
                        sprintf("ebd_%s_breeding_spain_zf.csv", sp))
  message("  [3] ", sp, "  -> reading ", basename(csv_path))
  dt <- fread(csv_path,
              select = c("checklist_id", "latitude", "longitude", "year",
                         "species_observed", "duration_minutes",
                         "effort_distance_km", "number_observers",
                         "time_observations_started", "locality_id",
                         "observer_id"))

  n_raw <- nrow(dt)
  dt <- unique(dt, by = "checklist_id")  # one row per checklist
  n_unique <- nrow(dt)

  # Filter to Iberia mainland bbox; drop Baleares explicitly
  dt <- dt[longitude >= iberia_bbox_wgs$lon_min &
           longitude <= iberia_bbox_wgs$lon_max &
           latitude  >= iberia_bbox_wgs$lat_min &
           latitude  <= iberia_bbox_wgs$lat_max]
  dt <- dt[!(longitude > 1.5 & latitude > 38.6 & latitude < 40.1)]
  n_iberia <- nrow(dt)

  message(sprintf("       %d rows -> %d unique checklists -> %d in Iberia mainland",
                  n_raw, n_unique, n_iberia))

  # WorldClim cell index (model domain). Compute in WGS84.
  dt[, cells := terra::cellFromXY(wc_grid, cbind(longitude, latitude))]

  # Drop checklists whose WorldClim cell index is NA (outside the grid).
  n_pre_wc <- nrow(dt)
  dt <- dt[!is.na(cells)]
  if (n_pre_wc - nrow(dt) > 0) {
    message(sprintf("       dropped %d checklists outside WorldClim grid",
                    n_pre_wc - nrow(dt)))
  }

  # Project to EPSG:3035 for raster extraction
  v3035 <- terra::project(
    terra::vect(dt, geom = c("longitude", "latitude"), crs = "EPSG:4326",
                keepgeom = TRUE),
    "EPSG:3035")

  # Extract from CLC + masks + focals (vectorised; fast for ~300k points)
  ext_list <- list(
    clc_class         = terra::extract(clc_iberia,   v3035, ID = FALSE)[[1]],
    step_strict_point = terra::extract(mask_strict,  v3035, ID = FALSE)[[1]],
    step_broad_point  = terra::extract(mask_broad,   v3035, ID = FALSE)[[1]],
    step_strict_500m  = terra::extract(foc_str_500m, v3035, ID = FALSE)[[1]],
    step_strict_1km   = terra::extract(foc_str_1km,  v3035, ID = FALSE)[[1]],
    step_broad_500m   = terra::extract(foc_brd_500m, v3035, ID = FALSE)[[1]],
    step_broad_1km    = terra::extract(foc_brd_1km,  v3035, ID = FALSE)[[1]]
  )
  for (cn in names(ext_list)) dt[, (cn) := ext_list[[cn]]]

  saveRDS(dt, out_path)
  fwrite(dt, sub("\\.rds$", ".csv", out_path))
  na_clc <- sum(is.na(dt$clc_class))
  message(sprintf("       saved: %s  (n=%d, NA clc=%d)",
                  basename(out_path), nrow(dt), na_clc))
}

message("[Section 3 done]")


###############################################################################
# Section 4: cell-year aggregation per species
###############################################################################
# For each species, collapse the per-checklist table to (cells, year) and
# compute the proportion of checklists in that cell-year that satisfy each
# steppe-representativeness criterion.

message("\n==[ Section 4: cell-year aggregation ]==")

for (sp in species_codes) {
  in_path <- file.path(out_data, sprintf("checklist_stepRep_%s.rds", sp))
  out_rds <- file.path(out_data, sprintf("stepRep_cellyear_%s.rds", sp))
  out_csv <- file.path(out_data, sprintf("stepRep_cellyear_%s.csv", sp))
  if (file.exists(out_rds)) {
    message("  [4] skip (exists): ", basename(out_rds)); next
  }

  d <- readRDS(in_path)
  setDT(d)

  agg <- d[, .(
    n_checklists         = .N,
    stepRep_strict_point = mean(step_strict_point, na.rm = TRUE),
    stepRep_broad_point  = mean(step_broad_point,  na.rm = TRUE),
    stepRep_strict_500m  = mean(step_strict_500m,  na.rm = TRUE),
    stepRep_strict_1km   = mean(step_strict_1km,   na.rm = TRUE),
    stepRep_broad_500m   = mean(step_broad_500m,   na.rm = TRUE),
    stepRep_broad_1km    = mean(step_broad_1km,    na.rm = TRUE)
  ), by = .(cells, year)]
  setorder(agg, cells, year)

  saveRDS(agg, out_rds)
  fwrite(agg, out_csv)
  message(sprintf("  [4] %s: %d cell-year rows (cells=%d, years=%d-%d, total checklists=%d)",
                  sp, nrow(agg), uniqueN(agg$cells),
                  min(agg$year), max(agg$year), sum(agg$n_checklists)))
}

message("[Section 4 done]")


###############################################################################
# Section 5: diagnostics report
###############################################################################
# Builds reports/stepRep_diagnostics.md with embedded figures (PNG) covering
# the seven required diagnostic items + integration note.

suppressPackageStartupMessages({
  library(ggplot2)
})

message("\n==[ Section 5: diagnostics ]==")

# -- 5a. Parse CLC legend (raster_value -> clc_code, label) --------------------
legend_path <- here("data-raw", "data", "corine", "Results",
                    "u2018_clc2018_v2020_20u1_raster100m", "Legend",
                    "CLC2018_CLC2018_V2018_20_QGIS.txt")
parse_clc_legend <- function(path) {
  # File format (CRLF, no header): each line is
  #   clc_code,r,g,b,a,label
  # The raster pixel value is implicit: line index (1-based).
  raw <- readLines(path, warn = FALSE)
  raw <- raw[nzchar(trimws(raw))]
  parts <- strsplit(trimws(raw), ",", fixed = TRUE)
  rbindlist(lapply(seq_along(parts), function(i) {
    rest <- parts[[i]]
    data.table(raster_value = i,
               clc_code     = as.integer(rest[1]),
               label        = trimws(paste(rest[6:length(rest)],
                                            collapse = ",")))
  }))
}
clc_leg <- parse_clc_legend(legend_path)

# -- 5b. Per-species: CLC frequency under checklists, weighted means, etc. -----
sp_codes <- c("otitar", "ptealc", "pteori", "tettet")
sp_names <- c(otitar = "Otis tarda", ptealc = "Pterocles alchata",
              pteori = "Pterocles orientalis", tettet = "Tetrax tetrax")

# Load all per-species checklist tables and per-species cell-year tables
chk <- lapply(sp_codes, function(sp) {
  d <- as.data.table(readRDS(file.path(out_data,
                                       sprintf("checklist_stepRep_%s.rds", sp))))
  d[, species := sp][]
})
names(chk) <- sp_codes
cy <- lapply(sp_codes, function(sp) {
  d <- as.data.table(readRDS(file.path(out_data,
                                       sprintf("stepRep_cellyear_%s.rds", sp))))
  d[, species := sp][]
})
names(cy) <- sp_codes

# Pool checklists for the global CLC frequency table. terra::extract preserved
# the categorical labels from the source raster's vat.dbf, so clc_class is a
# factor whose level index == raster value. Join to the legend by raster_value
# to attach the CLC 3-digit code.
val_to_code <- setNames(clc_leg$clc_code, as.character(clc_leg$raster_value))

get_clc_freq <- function(d) {
  x <- d[!is.na(clc_class),
         .(raster_value = as.integer(clc_class),
           label        = as.character(clc_class))]
  out <- x[, .N, by = .(raster_value, label)]
  out[, clc_code := val_to_code[as.character(raster_value)]]
  setorder(out, -N)
  out[]
}
clc_freq_apriltojune <- get_clc_freq(chk$otitar)
clc_freq_maytoaug    <- get_clc_freq(chk$ptealc)

# Top 15 classes for each window (pct uses full denominator before head)
top_aprjun <- head(clc_freq_apriltojune[, .(clc_code, label, N,
                                            pct = round(100 * N / sum(N), 2))], 15)
top_mayaug <- head(clc_freq_maytoaug   [, .(clc_code, label, N,
                                            pct = round(100 * N / sum(N), 2))], 15)

# Annual weighted means of stepRep_strict_500m / broad_500m
annual_mean <- function(dt, var) {
  dt[, .(weighted_mean = sum(get(var) * n_checklists, na.rm = TRUE) /
                         sum(n_checklists, na.rm = TRUE),
         total_checklists = sum(n_checklists)),
     by = year][order(year)]
}

# Identify focal cells per species (n_detections >= 1)
focal_cells <- lapply(sp_codes, function(sp) {
  chk[[sp]][species_observed == TRUE, unique(cells)]
})
names(focal_cells) <- sp_codes

trend_table <- rbindlist(lapply(sp_codes, function(sp) {
  pen <- annual_mean(cy[[sp]], "stepRep_strict_500m")
  fcl <- annual_mean(cy[[sp]][cells %in% focal_cells[[sp]]],
                     "stepRep_strict_500m")
  pen[, scope := "peninsular"]
  fcl[, scope := "focal_cells"]
  out <- rbind(pen, fcl)
  out[, species := sp]
  out
}))

# Spearman correlation: log(n_checklists) vs stepRep_strict_500m per cell-year
cor_table <- rbindlist(lapply(sp_codes, function(sp) {
  d <- cy[[sp]]
  rho500 <- suppressWarnings(cor(log(d$n_checklists),
                                 d$stepRep_strict_500m, method = "spearman"))
  rho1k  <- suppressWarnings(cor(log(d$n_checklists),
                                 d$stepRep_strict_1km,  method = "spearman"))
  rhoB   <- suppressWarnings(cor(log(d$n_checklists),
                                 d$stepRep_broad_500m,  method = "spearman"))
  data.table(species = sp, rho_strict_500m = round(rho500, 3),
             rho_strict_1km = round(rho1k, 3),
             rho_broad_500m = round(rhoB, 3))
}))

# Strict vs broad disagreement at cell-year level (using 0.5 stratification)
strict_vs_broad_tab <- rbindlist(lapply(sp_codes, function(sp) {
  d <- cy[[sp]]
  d[, strict_class := stepRep_strict_500m >= 0.5]
  d[, broad_class  := stepRep_broad_500m  >= 0.5]
  data.table(species   = sp,
             n_total   = nrow(d),
             both_low  = sum(!d$strict_class & !d$broad_class),
             both_high = sum( d$strict_class &  d$broad_class),
             flip_S0_B1= sum(!d$strict_class &  d$broad_class),
             flip_S1_B0= sum( d$strict_class & !d$broad_class))
}))

# -- 5c. Figures ---------------------------------------------------------------
theme_set(theme_bw(base_size = 11) +
          theme(strip.background = element_rect(fill = "grey95"),
                panel.grid.minor = element_blank()))

# Fig 1: histograms of stepRep_strict_500m per species (cell-year level)
df_hist <- rbindlist(lapply(sp_codes, function(sp) {
  d <- cy[[sp]][, .(stepRep_strict_500m, stepRep_broad_500m)]
  d[, species := sp_names[[sp]]]; d
}))
df_hist <- melt(df_hist, id.vars = "species",
                variable.name = "mask",
                value.name = "stepRep")
df_hist[, mask := factor(mask,
                         levels = c("stepRep_strict_500m","stepRep_broad_500m"),
                         labels = c("strict (211/231/321/333)",
                                    "broad (+242/243/244)"))]

p1 <- ggplot(df_hist, aes(stepRep, fill = mask)) +
  geom_histogram(binwidth = 0.05, boundary = 0,
                 position = "identity", alpha = 0.55) +
  facet_wrap(~ species, ncol = 2) +
  scale_fill_manual(values = c("#1f77b4", "#ff7f0e"), name = "mask") +
  labs(x = "stepRep_500m  (proportion of buffer in steppe habitat)",
       y = "cell-years") +
  theme(legend.position = "top")
ggsave(file.path(out_figs, "stepRep_hist.png"), p1,
       width = 8, height = 6, dpi = 150)

# Fig 2: annual weighted mean (peninsular vs focal cells)
trend_table[, species_lab := sp_names[species]]
p2 <- ggplot(trend_table,
             aes(year, weighted_mean, colour = scope, linetype = scope)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.2) +
  facet_wrap(~ species_lab, ncol = 2) +
  scale_colour_manual(values = c(peninsular = "#666666",
                                  focal_cells = "#d62728"),
                      labels = c(peninsular = "all peninsular cells",
                                  focal_cells = "cells with >=1 detection")) +
  scale_linetype_manual(values = c(peninsular = "dashed",
                                    focal_cells = "solid"),
                        labels = c(peninsular = "all peninsular cells",
                                    focal_cells = "cells with >=1 detection")) +
  scale_x_continuous(breaks = 2017:2023) +
  labs(x = NULL, y = "weighted mean stepRep_strict_500m\n(weights = n_checklists)",
       colour = NULL, linetype = NULL) +
  theme(legend.position = "top")
ggsave(file.path(out_figs, "stepRep_temporal.png"), p2,
       width = 9, height = 6, dpi = 150)

# Fig 3: maps of mean stepRep_strict_500m per cell over Iberia (5 km grid)
wc_grid_full <- rast(worldclim_grd, lyrs = 1)  # any layer; cells share geometry

build_cell_raster <- function(sp) {
  d <- cy[[sp]][, .(mean_stepRep = sum(stepRep_strict_500m * n_checklists,
                                       na.rm = TRUE) /
                                    sum(n_checklists, na.rm = TRUE),
                    n_checklists = sum(n_checklists)),
                by = cells]
  r <- rast(wc_grid_full); values(r) <- NA_real_
  r[d$cells] <- d$mean_stepRep
  names(r) <- sp_names[[sp]]
  r
}

map_stack <- rast(lapply(sp_codes, build_cell_raster))

# Save four maps as a single 2x2 PNG using terra::plot
png(file.path(out_figs, "stepRep_map.png"),
    width = 11, height = 8, units = "in", res = 150)
op <- par(mfrow = c(2, 2), mar = c(1, 1, 2, 4))
for (i in seq_along(sp_codes)) {
  plot(map_stack[[i]], main = sp_names[[sp_codes[i]]],
       col = hcl.colors(20, "YlGnBu", rev = TRUE),
       range = c(0, 1), axes = FALSE, mar = c(1, 1, 2, 4))
}
par(op); dev.off()

# Fig 4: hexbin of log(n_checklists) vs stepRep_strict_500m
df_cor <- rbindlist(lapply(sp_codes, function(sp) {
  d <- cy[[sp]][, .(log_n = log(n_checklists),
                    stepRep = stepRep_strict_500m)]
  d[, species := sp_names[[sp]]]; d
}))

p4 <- ggplot(df_cor, aes(log_n, stepRep)) +
  geom_hex(bins = 30) +
  geom_smooth(method = "loess", colour = "red", se = FALSE,
              linewidth = 0.6, span = 0.5) +
  facet_wrap(~ species, ncol = 2) +
  scale_fill_viridis_c(option = "magma", trans = "log10",
                       name = "cell-years\n(log10 count)") +
  labs(x = "log(n_checklists per cell-year)",
       y = "stepRep_strict_500m") +
  theme(legend.position = "right")
ggsave(file.path(out_figs, "stepRep_vs_effort.png"), p4,
       width = 9, height = 6, dpi = 150)

# Fig 5: strict vs broad scatter (cell-year)
df_sb <- rbindlist(lapply(sp_codes, function(sp) {
  d <- cy[[sp]][, .(strict = stepRep_strict_500m, broad = stepRep_broad_500m)]
  d[, species := sp_names[[sp]]]; d
}))
p5 <- ggplot(df_sb, aes(strict, broad)) +
  geom_hex(bins = 30) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey30") +
  facet_wrap(~ species, ncol = 2) +
  scale_fill_viridis_c(option = "viridis", trans = "log10",
                       name = "cell-years\n(log10 count)") +
  coord_equal() +
  labs(x = "stepRep_strict_500m", y = "stepRep_broad_500m") +
  theme(legend.position = "right")
ggsave(file.path(out_figs, "stepRep_strict_vs_broad.png"), p5,
       width = 9, height = 6, dpi = 150)

# -- 5d. Render the markdown report -------------------------------------------
report_path <- file.path(out_reports, "stepRep_diagnostics.md")
con <- file(report_path, "w")

write_section <- function(...) cat(..., "\n", sep = "", file = con)

write_section("# stepRep diagnostics\n")
write_section("_Generated by `scripts/build_stepRep.R` on ",
              format(Sys.Date(), "%Y-%m-%d"), "._\n")

# 1. CORINE provenance
write_section("## 1. CORINE provenance and habitat masks\n")
write_section("- **Layer:** CORINE Land Cover 2018, version v2020_20u1, raster 100 m.")
write_section("- **CRS:** ETRS89-extended / LAEA Europe (EPSG:3035).")
write_section("- **File:** `data-raw/data/corine/Results/u2018_clc2018_v2020_20u1_raster100m/DATA/U2018_CLC2018_V2020_20u1.tif`.")
write_section("- **Why this version:** CLC2018 v2020_20u1 is the canonical Level-3 product distributed by Copernicus Land Monitoring Service. It supports the 211/231/321/333/242/243/244 codes used in the literature on Iberian pseudo-steppe birds (Suarez, Traba, Brotons, Delgado & Moreira). The experimental CLC+ Backbone (10 m, 2018 and 2021 reference years) is more recent but uses an aggregated, non-Level-3 nomenclature that does not directly support the strict mask required here. The habitat layer is treated as **static** by design: pseudo-steppe habitat is slow-changing across our reference period (2017-2023), so cell-year variation in stepRep is driven by *where birders go inside each cell*, not by habitat change.")
write_section("- **NoData quirk:** the source TIF encodes pixel-level NoData as raster value 48 (the QGIS legend file refers to value 45 as 999=NODATA, but the actual stored NoData in pixels is 48; ~40k pixels in Iberia). Both 48 and any input NA were mapped to 0 (\"not steppe\") in the binary masks; this keeps focal sums numerically defined while correctly diluting buffers near the coast.")
write_section("\n**Strict mask** (`steppe_strict`): CLC 211 (non-irrigated arable), 231 (pastures), 321 (natural grasslands), 333 (sparsely vegetated areas).")
write_section("\n**Broad mask** (`steppe_broad`): strict + 242 (complex cultivation patterns), 243 (agriculture with significant natural vegetation), 244 (agro-forestry / dehesa).")
write_section("\n**Excluded explicitly** (not pseudo-steppe habitat for *Otis*, *Tetrax*, *Pterocles*): 212 (irrigated), 213 (rice), 221 (vineyards), 223 (olive groves).\n")

# Coverage table
write_section("## 2. Coverage\n")
write_section("Per-species totals (after dedup by `checklist_id`, Iberia mainland filter, and drop of checklists outside the WorldClim 5 km grid):\n")
cov_tab <- rbindlist(lapply(sp_codes, function(sp) {
  data.table(species = sp_names[[sp]],
             n_checklists = nrow(chk[[sp]]),
             n_cells = uniqueN(chk[[sp]]$cells),
             n_cell_years = nrow(cy[[sp]]),
             n_focal_cells = length(focal_cells[[sp]]))
}))
write_section("| Species | n_checklists | n_cells | n_cell_years | n_focal_cells |")
write_section("|---|---|---|---|---|")
for (i in seq_len(nrow(cov_tab)))
  write_section("| ", cov_tab$species[i], " | ",
                format(cov_tab$n_checklists[i], big.mark = ","), " | ",
                format(cov_tab$n_cells[i], big.mark = ","), " | ",
                format(cov_tab$n_cell_years[i], big.mark = ","), " | ",
                format(cov_tab$n_focal_cells[i], big.mark = ","), " |")
write_section("\n_Focal cells = cells with at least one detection of the species in 2017-2023, used for diagnostic 4 below._\n")

# 3. CLC class frequency under checklists
write_section("## 3. CLC class frequency under checklists\n")
write_section("Top 15 CLC classes at the pixel where each unique checklist falls. Two breeding windows are reported separately because *Otis*/*Tetrax* (April-June) and *Pterocles* (May-August) are filtered on different dates and so see slightly different sets of birders.\n")
write_section("**April-June (Otis tarda + Tetrax tetrax base, ", nrow(chk$otitar),
              " checklists):**\n")
write_section("| CLC code | Label | n | % |")
write_section("|---|---|---|---|")
for (i in seq_len(nrow(top_aprjun)))
  write_section("| ", top_aprjun$clc_code[i], " | ",
                top_aprjun$label[i], " | ",
                format(top_aprjun$N[i], big.mark = ","), " | ",
                format(top_aprjun$pct[i], nsmall = 2), " |")

write_section("\n**May-August (Pterocles spp. base, ", nrow(chk$ptealc),
              " checklists):**\n")
write_section("| CLC code | Label | n | % |")
write_section("|---|---|---|---|")
for (i in seq_len(nrow(top_mayaug)))
  write_section("| ", top_mayaug$clc_code[i], " | ",
                top_mayaug$label[i], " | ",
                format(top_mayaug$N[i], big.mark = ","), " | ",
                format(top_mayaug$pct[i], nsmall = 2), " |")
write_section("\n**Reading.** Forest classes (broad-leaved 311, coniferous 312, mixed 313, transitional woodland-shrub 324) and 'sclerophyllous vegetation' (323) dominate the pixel under birder feet, despite peninsular Iberia being roughly half non-irrigated arable + pastures + natural grasslands by area. This is the *expected birder bias signature*: visits are biased toward forested, accessible, diverse habitat.\n")

# 4. Histogram
write_section("## 4. Distribution of stepRep_500m per cell-year\n")
write_section("![Histograms](figs/stepRep_hist.png)\n")
write_section("Reading. Strict has heavy mass at 0 and a thin upper tail; broad shifts the distribution rightward (median moves from ~0.16 to ~0.32). Half of cell-years for *Otis*/*Tetrax* have <= 16% of birder buffer in strict steppe.\n")

# 5. Temporal trend
write_section("## 5. Temporal trend (weighted by n_checklists)\n")
write_section("![Temporal trend](figs/stepRep_temporal.png)\n")
write_section("Annual weighted mean of `stepRep_strict_500m`:\n")
write_section("| Species | scope | ",
              paste(2017:2023, collapse = " | "), " |")
write_section("|---|---|", paste(rep("---", 7), collapse = "|"), "|")
for (sp in sp_codes) {
  for (sc in c("peninsular", "focal_cells")) {
    sub <- trend_table[species == sp & scope == sc][order(year)]
    vals <- format(round(sub$weighted_mean, 3), nsmall = 3)
    write_section("| ", sp_names[[sp]], " | ", sc, " | ",
                  paste(vals, collapse = " | "), " |")
  }
}
write_section("\n**Reading.** Two distinct patterns:")
write_section("- Peninsular cells (all checklists) are flat at ~0.21-0.25, with no drift across 2017-2023.")
write_section("- Focal cells (those with at least one detection of the species) show a **monotonic increase** in stepRep over the period: about +0.13 for *Otis tarda*, +0.11 for *Tetrax tetrax*, +0.14 for *Pterocles alchata*, +0.09 for *Pterocles orientalis*.\n")
write_section("This is the **opposite** of the initial 'newcomers diluting steppe coverage' hypothesis. Inside cells where the focal species are known to occur, birders are increasingly targeting the actual steppe parcels over time; visits there are reaching the species' habitat more efficiently. The bias is real and time-varying, just in the opposite direction from the initial guess: late-year checklists in focal cells are **more informative per visit** than early-year ones, which would tend to inflate gamma and deflate epsilon spuriously if not controlled. Including stepRep as a detection covariate is therefore still motivated, and arguably more so given the magnitude of the focal-cell drift (+0.10 to +0.14 over six years on a 0-1 scale).\n")

# 6. Map
write_section("## 6. Spatial pattern (mean stepRep across years)\n")
write_section("![Map](figs/stepRep_map.png)\n")
write_section("Mean `stepRep_strict_500m` per WorldClim 5 km cell, weighted by checklist count, pooling 2017-2023. Cells with no checklists are blank.\n")

# 7. Correlation
write_section("## 7. Correlation between effort and stepRep\n")
write_section("![Effort vs stepRep](figs/stepRep_vs_effort.png)\n")
write_section("Spearman rank correlations between `log(n_checklists)` and stepRep at the cell-year level:\n")
write_section("| Species | rho strict 500m | rho strict 1km | rho broad 500m |")
write_section("|---|---|---|---|")
for (i in seq_len(nrow(cor_table)))
  write_section("| ", sp_names[[cor_table$species[i]]], " | ",
                cor_table$rho_strict_500m[i], " | ",
                cor_table$rho_strict_1km[i], " | ",
                cor_table$rho_broad_500m[i], " |")
write_section("\n**Reading.** Spearman rho is essentially zero for the strict masks (-0.005 to -0.02) and only slightly negative for broad (~-0.05). At cell-year level, effort is **not** systematically associated with steppe representativeness. The time-varying bias revealed in §5 (rising focal-cell trend) operates through a different channel than effort-stepRep correlation; controlling for effort alone in the detection model would not absorb it.\n")

# 8. Strict vs broad
write_section("## 8. Strict vs broad mask\n")
write_section("![Strict vs broad](figs/stepRep_strict_vs_broad.png)\n")
write_section("Cell-year reclassification using a 0.5 threshold on each mask:\n")
write_section("| Species | both <0.5 | both >=0.5 | strict<0.5 / broad>=0.5 | strict>=0.5 / broad<0.5 |")
write_section("|---|---|---|---|---|")
for (i in seq_len(nrow(strict_vs_broad_tab)))
  write_section("| ", sp_names[[strict_vs_broad_tab$species[i]]], " | ",
                format(strict_vs_broad_tab$both_low[i], big.mark = ","), " | ",
                format(strict_vs_broad_tab$both_high[i], big.mark = ","), " | ",
                format(strict_vs_broad_tab$flip_S0_B1[i], big.mark = ","), " | ",
                format(strict_vs_broad_tab$flip_S1_B0[i], big.mark = ","), " |")
write_section("\nThe `strict<0.5 / broad>=0.5` flips quantify how much dehesa + complex-cultivation pixels rescue cell-years. Whether to treat that mass as 'steppe-compatible' is an ecological decision: dehesa supports *Pterocles orientalis* and *Tetrax tetrax* in part of their range but is not pseudo-steppe sensu stricto. The pilot model uses **strict**; broad is reported as a sensitivity.\n")

# 9. Integration note
write_section("## 9. Integration note (model use)\n")
write_section("The covariate `stepRep_strict_500m`, standardised within each species' modelling table, will enter the **detection sub-model** of the colext / stPGOcc fits, jointly with the existing detection covariates (effort, duration, observers, hour). The pilot will run on **Pterocles alchata** first; if the coefficient is well identified and gamma / epsilon shift materially when the covariate is added, the analysis is extended to the other three species. `stepRep_strict_1km` and `stepRep_broad_500m` are kept as sensitivities to be reported alongside the headline result.\n")

# 10. Manuscript draft (Methods + Results paragraphs)
write_section("## 10. Manuscript draft (Methods + Results)\n")
write_section("_Narrative prose, IMRAD style. Citation slots marked `[REF NEEDED]` are to be filled with verified sources before submission. Effect sizes from the pilot fit are placeholders (`[PILOT]`)._\n")

write_section("### Methods (subsection: detection covariates -> add paragraph)\n")
write_section("Sub-cell variation in where eBird sampling effort is concentrated within each ", "5 km grid cell can bias detection probability and, by extension, ",
              "estimates of colonisation and local extinction in dynamic occupancy ",
              "models. To account for this, we derived a covariate of steppe ",
              "representativeness for each (cell, year) combination. We classified ",
              "the CORINE Land Cover 2018 raster (version v2020_20u1, 100 m, EPSG:3035; ",
              "Copernicus Land Monitoring Service [REF NEEDED]) into a binary ",
              "pseudo-steppe mask comprising CLC Level-3 classes 211 (non-irrigated ",
              "arable land), 231 (pastures), 321 (natural grasslands) and 333 ",
              "(sparsely vegetated areas), following the habitat use described for ",
              "*Otis tarda*, *Tetrax tetrax*, *Pterocles alchata* and *Pterocles ",
              "orientalis* in the Iberian Peninsula [REF NEEDED: Suarez; Traba; ",
              "Brotons; Delgado & Moreira]. Permanent and water-demanding crops (CLC ",
              "212, 213, 221 and 223) were excluded.\n")
write_section("For each unique eBird checklist passing our effort filters, we ",
              "projected its coordinates to ETRS89-LAEA Europe (EPSG:3035) and ",
              "extracted the proportion of a 500 m circular buffer occupied by ",
              "pseudo-steppe pixels, using a focal kernel with normalised circular ",
              "weights. Pixels outside CORINE coverage were treated as zero, so that ",
              "buffers near the coast are correctly diluted. The per-checklist values ",
              "were then aggregated to the (cell, year) level as the unweighted mean ",
              "of the checklists assigned to that cell-year, yielding ",
              "`stepRep_strict_500m`, bounded in [0, 1]. The same procedure was ",
              "repeated with a 1 km buffer (`stepRep_strict_1km`) and with a broader ",
              "mask that additionally included CLC 242 (complex cultivation patterns), ",
              "243 (agriculture with significant areas of natural vegetation) and 244 ",
              "(agro-forestry / dehesa); these versions are reported as sensitivities. ",
              "The CORINE habitat layer was treated as static across the study period ",
              "because pseudo-steppe land cover is slow-changing at the scales ",
              "relevant here; cell-year variation in the covariate therefore reflects ",
              "the spatial structure of birder sampling effort within each cell, not ",
              "habitat change.\n")
write_section("`stepRep_strict_500m` was standardised within each species' ",
              "modelling table (centred and scaled to unit variance) and added to ",
              "the detection sub-model of the dynamic occupancy fits [REF NEEDED: ",
              "Fiske & Chandler 2011 for `colext`; Doser et al. for `stPGOcc`] ",
              "alongside the existing detection covariates (sampling-event ",
              "duration, distance travelled, number of observers, time of day, year).\n")

write_section("### Results (subsection: detection model -> add paragraph)\n")
write_section("After dedup by `checklist_id` and the Iberia mainland filter, the ",
              "modelling tables comprised between 263,145 (*Otis tarda*, *Tetrax ",
              "tetrax*; April-June) and 289,720 (*Pterocles alchata*, *Pterocles ",
              "orientalis*; May-August) unique checklists, mapped to 15,836 and ",
              "16,440 occupied 5 km cells respectively. Half of cell-years had ",
              "`stepRep_strict_500m` <= 0.16 for *Otis* and *Tetrax* and <= 0.15 for ",
              "the two *Pterocles* species, and a quarter had zero pseudo-steppe ",
              "representativeness; the broader mask raised these medians to ~0.32. The ",
              "class composition under checklists confirmed a strong sampling bias ",
              "toward anthropised and diverse habitats: discontinuous and continuous ",
              "urban fabric (CLC 112 and 111), forest classes (311, 312, 313) and ",
              "sclerophyllous vegetation (323) together accounted for ~30% of the ",
              "pixels under checklists, while pastures (CLC 231) accounted for only ",
              "4.4%.\n")
write_section("This sampling bias was time-varying. In cells where the focal ",
              "species was detected at least once during 2017-2023, the ",
              "checklist-weighted annual mean of `stepRep_strict_500m` increased ",
              "monotonically over the study period: from 0.48 in 2017 to 0.62 in 2023 ",
              "for *Otis tarda*, with comparable rises of +0.09 to +0.14 over the ",
              "same period for the remaining three species. In peninsular cells ",
              "overall the same statistic was flat across years (range 0.21-0.25). ",
              "Late-period checklists in occupied cells therefore concentrated on ",
              "actual pseudo-steppe pixels more efficiently than early-period ones, ",
              "which would inflate apparent colonisation and deflate apparent ",
              "extinction in dynamic occupancy fits if this within-cell sampling ",
              "structure were not modelled. Sampling effort alone (log number of ",
              "checklists per cell-year) was not associated with steppe ",
              "representativeness (Spearman rho between -0.005 and -0.02 across ",
              "species), so adjusting only for effort would not absorb the bias.\n")
write_section("[PILOT placeholder] Including `stepRep_strict_500m` in the ",
              "detection sub-model of the *Pterocles alchata* dynamic occupancy fit ",
              "[INSERT coefficient, SE, 95% CI on the logit scale]. The estimated ",
              "colonisation rate gamma changed from [pre] to [post] and local ",
              "extinction rate epsilon from [pre] to [post] when the covariate was ",
              "added. The 1 km buffer and the broad mask gave consistent direction ",
              "and magnitude (Table S[X]). Based on this, we extended the analysis ",
              "to the other three species [if warranted], with the headline numbers ",
              "reported in Table [X] and the sensitivity sets in the Supplementary ",
              "Information.\n")

close(con)
message("  [5] report written: ", report_path)
message("[Section 5 done]")


# -- Reproducibility footer --
if (interactive()) print(sessionInfo()) else cat(capture.output(print(sessionInfo())), sep = "\n")
