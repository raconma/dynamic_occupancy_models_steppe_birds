###############################################################################
# spatial_layers.R
#
# Shared spatial layers for scripts 21-27 (all returned in ETRS89 / UTM 30N):
#   load_spa_sites()   Natura 2000 Special Protection Areas (MITECO 2023), one
#                      feature per SPA with code and name
#   load_spa_union()   the same SPAs dissolved into one geometry
#   load_n2k_sites()   all Natura 2000 sites (SPAs, SCIs, SACs), optionally
#                      filtered by designation
#   load_atlas2()      II Breeding Bird Atlas 10x10 km grid (aves_spain.shp;
#                      matches MITECO INB 2007 in GBIF, fieldwork 1998-2002)
#   wc_cell_polygons() polygons of the WorldClim 2.5 arc-min cells that define
#                      the model sites (the "cells" id in occ_wide_dynamic.csv)
#   cover_by_polygons(), cover_by_raster()  exact cover fractions per cell
###############################################################################

UTM30 <- 25830

load_spa_union <- function() {
  sf::st_union(load_spa_sites())
}

# Natura 2000 sites (MITECO 2023), one feature per site with code, name and
# designation, geometry repaired. `desig` is a pattern on the INSPIRE
# designation (specialProtectionArea, siteOfCommunityImportance,
# specialAreaOfConservation); NULL keeps all sites.
load_n2k_sites <- function(desig = NULL) {
  n2k <- sf::st_read(here::here("data-raw", "data", "threats", "rn2000",
                                "PS_NATURA2000_2023_PB.shp"), quiet = TRUE)
  if (!is.null(desig)) n2k <- n2k[grepl(desig, n2k$desig), ]
  # One SPA ("Estepas cerealistas de los rios Jarama y Henares", a key great
  # bustard site) holds a one-vertex hole that makes st_make_valid() fail.
  # Drop rings with < 4 vertices (degenerate by definition) before repairing.
  drop_degenerate_rings <- function(g) {
    sf::st_multipolygon(lapply(g, function(poly)
      Filter(function(r) nrow(r) >= 4, poly)))
  }
  sf::st_geometry(n2k) <- sf::st_sfc(lapply(sf::st_geometry(n2k),
                                            drop_degenerate_rings),
                                     crs = sf::st_crs(n2k))
  n2k <- sf::st_make_valid(sf::st_transform(n2k, UTM30))
  n2k[, c("localId", "SOName", "desig")]
}

# Individual SPAs (code and name)
load_spa_sites <- function() {
  load_n2k_sites("specialProtectionArea")[, c("localId", "SOName")]
}

load_atlas2 <- function() {
  sf::st_transform(sf::st_read(here::here(
    "data-raw", "data", "validation", "atlas_biodiversidad",
    "atlas_biodiversidad", "aves_spain.shp"), quiet = TRUE), UTM30)
}

wc_grid <- function() {
  terra::rast(here::here("data-raw", "data", "environmental_data",
                         "environmental_data_occ", "variables_spain.grd"),
              lyrs = 1)
}

# Square lon/lat polygons of the given WorldClim cell ids, in UTM 30N
wc_cell_polygons <- function(cell_ids, grid = wc_grid()) {
  xy <- terra::xyFromCell(grid, cell_ids)
  h <- terra::res(grid) / 2
  polys <- lapply(seq_len(nrow(xy)), function(i) {
    x <- xy[i, 1]; y <- xy[i, 2]
    sf::st_polygon(list(rbind(c(x - h[1], y - h[2]), c(x + h[1], y - h[2]),
                              c(x + h[1], y + h[2]), c(x - h[1], y + h[2]),
                              c(x - h[1], y - h[2]))))
  })
  sf::st_transform(sf::st_sf(cells = cell_ids,
                             geometry = sf::st_sfc(polys, crs = 4326)), UTM30)
}

# Exact cover fractions of model cells. terra::project(method = "average")
# onto the WorldClim grid shifts the averaging window and misstates cover by
# up to ~0.1, so cover is computed per cell polygon instead.
#   cover_by_polygons(): fraction of each cell inside a set of polygons
#                        (exact area intersection in UTM 30N)
#   cover_by_raster():   fraction of each cell whose raster values fall in
#                        each class set: every raster pixel is assigned to the
#                        cell containing its centre (terra::rasterize) and
#                        class membership is averaged per cell (terra::zonal),
#                        one pass over the cropped raster; `na_values` are
#                        treated as missing, not as absence
cover_by_polygons <- function(polys, sites) {
  u <- sf::st_union(sf::st_transform(sites, UTM30))
  ov <- suppressWarnings(sf::st_intersection(polys, sf::st_sf(geometry = u)))
  a <- tapply(as.numeric(sf::st_area(ov)), ov$cells, sum)
  out <- setNames(rep(0, nrow(polys)), polys$cells)
  out[names(a)] <- a
  as.numeric(out) / as.numeric(sf::st_area(polys))
}

cover_by_raster <- function(polys, r, class_sets, na_values = NULL) {
  v <- terra::vect(sf::st_transform(polys, terra::crs(r)))
  rc <- terra::crop(r, terra::ext(v) + 1000)
  id <- terra::rasterize(v, rc, field = seq_len(nrow(polys)))
  out <- as.data.frame(lapply(class_sets, function(cls) {
    m <- rc %in% cls
    if (!is.null(na_values)) m <- terra::ifel(rc %in% na_values, NA, m)
    z <- terra::zonal(m, id, fun = "mean", na.rm = TRUE)
    z[[2]][match(seq_len(nrow(polys)), z[[1]])]
  }))
  out
}
