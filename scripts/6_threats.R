knitr::opts_chunk$set(echo = TRUE, message=FALSE, comment = '', fig.width = 6, fig.height = 6)
#clear environment
rm(list=ls())


library(terra)
library(rnaturalearth)
library(sf)

###################################################################################################
#####################                     Richness                     ############################
###################################################################################################


richness <- rast("../data/richness.tif")

spain_admin <- ne_states(country = "Spain", returnclass = "sf")
spain_peninsula <- spain_admin[!spain_admin$name %in% c("Canary Islands", "Ceuta", "Melilla"), ]
spain_peninsula <- st_transform(spain_peninsula, "EPSG:25830")

###################################################################################################
#####################                       Solar                      ############################
###################################################################################################

# https://signa.ign.es/signa/
photovoltaic_plants <- read.csv("../data/threats/solar_plants.csv")
photovoltaic_sf <- st_as_sf(photovoltaic_plants, coords = c("x", "y"), crs = 4326)
photovoltaic_sf <- st_transform(photovoltaic_sf, "EPSG:25830")
plot(st_geometry(photovoltaic_sf), col = "black", pch = 19, cex = 0.5)
plot(richness)
plot(st_geometry(spain_peninsula), add = TRUE, border = "black", lwd = 1)
coords <- st_coordinates(photovoltaic_sf)
points(coords[, 1], coords[, 2], col = "black", pch = 19, cex = 0.5)

# How many solar plants are built in areas with at least one steppe bird
photovoltaic_values <- extract(richness, st_coordinates(photovoltaic_sf))

total_rows <- nrow(photovoltaic_values)

for(i in 1:6){
  rows_not_zero <- sum(photovoltaic_values[["sum"]] >= i, na.rm = TRUE)
  percentage_not_zero <- (rows_not_zero / total_rows) * 100
  print(i)
  print(percentage_not_zero)
}
