###############################################################################
# spatial_layers.R
#
# Shared spatial layers for scripts 21-23 (all returned in ETRS89 / UTM 30N):
#   load_spa_union()   Natura 2000 Special Protection Areas (MITECO 2023)
#   load_atlas2()      II Breeding Bird Atlas 10x10 km grid (aves_spain.shp;
#                      matches MITECO INB 2007 in GBIF, fieldwork 1998-2002)
#   wc_cell_polygons() polygons of the WorldClim 2.5 arc-min cells that define
#                      the model sites (the "cells" id in occ_wide_dynamic.csv)
###############################################################################

UTM30 <- 25830

load_spa_union <- function() {
  spa <- sf::st_read(here::here("data-raw", "data", "threats", "rn2000",
                                "PS_NATURA2000_2023_PB.shp"), quiet = TRUE)
  spa <- spa[grepl("specialProtectionArea", spa$desig), ]
  # One SPA ("Estepas cerealistas de los rios Jarama y Henares", a key great
  # bustard site) holds a one-vertex hole that makes st_make_valid() fail.
  # Drop rings with < 4 vertices (degenerate by definition) before repairing.
  drop_degenerate_rings <- function(g) {
    sf::st_multipolygon(lapply(g, function(poly)
      Filter(function(r) nrow(r) >= 4, poly)))
  }
  sf::st_geometry(spa) <- sf::st_sfc(lapply(sf::st_geometry(spa),
                                            drop_degenerate_rings),
                                     crs = sf::st_crs(spa))
  sf::st_union(sf::st_make_valid(sf::st_transform(spa, UTM30)))
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
