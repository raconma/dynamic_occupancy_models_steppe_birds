knitr::opts_chunk$set(echo = TRUE, message=FALSE, comment = '', fig.width = 6, fig.height = 6)
#clear environment
rm(list=ls())

library(terra)
library(rnaturalearth)
library(sf)
library(raster)

sdm_threshold <- function(sdm, occs, quantile = .25, binary = FALSE){
  occs_col <- occs[, c("longitude", "latitude")] 
  occPredVals <- terra::extract(sdm, occs_col)
  thresh <- quantile(na.omit(occPredVals[,2]), c(quantile))
  sdm_thresh <- sdm
  sdm_thresh[sdm_thresh < thresh] <- 0
  sdm_thresh[sdm_thresh >= thresh] <- 1
  return(sdm_thresh)
}

species <- c("cirpyg", "falnau", "otitar", "ptealc", "pteori", "tettet")

raster_list <- list()

predicted_occupancy_ref <- rast(paste0("../data/", species[1], "/", species[1], "_OccuMap.tif"))

for (i in species) {
  predicted_occupancy <- rast(paste0("../data/", i, "/", i, "_OccuMap.tif"))
  initial_occurrence <- read.csv(paste0("../data/", i, "/", i, "_occ_species_observed.csv"))
  
  predicted_occupancy <- resample(predicted_occupancy, predicted_occupancy_ref)
  
  plot(predicted_occupancy)
  points(initial_occurrence$longitude, initial_occurrence$latitude, col = "black", pch = 19, cex = 0.5)
  
  threshold_25 <- sdm_threshold(predicted_occupancy, initial_occurrence, quantile = .25)
  plot(threshold_25)
  
  raster_list[[i]] <- threshold_25
}

raster_stack <- rast(raster_list)

sum_raster <- sum(raster_stack)

plot(sum_raster)

writeRaster(sum_raster, "../data/richness.tif", overwrite = TRUE)
