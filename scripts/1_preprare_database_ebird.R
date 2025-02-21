knitr::opts_chunk$set(echo = TRUE, message=FALSE, comment = '', fig.width = 6, fig.height = 6)
#clear environment
rm(list=ls())

# load packages
library(raster)
library(ggplot2)
library(ggmap)
#library(ggsn) # no longer available
library(viridisLite)
library(auk)
library(lubridate)
library(sf)
library(gridExtra)
library(tidyverse)
library(grid)
#library(velox) # no longer available
#library(MODIS) # no longer available
library(dggridR)
library(unmarked)
library(ebirdst)
library(MuMIn)
library(AICcmodavg)
library(fields)
library(purrr)
library(dplyr)
library(knitr)
library(pander)

# https://cornelllabofornithology.github.io/ebird-best-practices/intro.html

# Resolve namespace conflicts
select <- dplyr::select
filter <- dplyr::filter

# Create a new folder: "data" within your working directory where we'll store data and save files
data_dir <- "../../data/"
if (!dir.exists(data_dir)) {
  dir.create(data_dir)
}

#auk_set_awk_path()
ebd <- auk_ebd("../raw_data/ebird_raw_mar2024/ebd_ES_smp_relMar-2024.txt")
# Read the raw data downloaded from eBird
ebd <- auk_ebd("../raw_data/ebird_raw_mar2024/ebd_ES_smp_relMar-2024.txt", 
               file_sampling = "../raw_data/ebird_raw_mar2024/ebd_sampling_relMar-2024/ebd_sampling_relMar-2024.txt")

# Declare the studied species
species <- c(c("Otis tarda", "otitar"),
             c("Pterocles alchata", "ptealc"),
             c("Pterocles orientalis", "pteori"),
             c("Tetrax tetrax", "tettet")
             )

# Loop through every species
for (i in seq(1, length(species), by = 2)) {
  # Create a new folder: "data/especies_name_abbreviated" for each species
  data_dir <- paste0("../data/",species[i + 1])
  if (!dir.exists(data_dir)) {
    dir.create(data_dir)
  }
  
  # Breeding season for each species
  if (species[i] == "Otis tarda" | species[i] == "Tetrax tetrax") {
    date_range <- c("*-04-01", "*-06-30")
  } else if (species[i] == "Pterocles alchata" | species[i] == "Pterocles orientalis") {
    date_range <- c("*-05-01", "*-08-31")
  }
  print((species[i]))
  print(date_range)
  
  # Filter the data by species, country code, observation dates and observation protocols
  ebd_filters <- ebd %>% 
    auk_species(species[i]) %>% 
    auk_country("ES") %>% 
    auk_date(date = date_range) %>% 
    auk_protocol(protocol = c("Stationary", "Traveling")) %>% 
    auk_complete()
  
  # Save the data as a text file
  ####################### DELETE FILES FIRST #######################
  f_ebd <- file.path(data_dir, paste0("ebd_",species[i + 1],"_breeding.txt"))
  print(f_ebd)
  f_sampling <- file.path(data_dir, paste0("ebd_checklists_breeding_",species[i + 1],".txt"))
  if (!file.exists(f_ebd)) {
    auk_filter(ebd_filters, file = f_ebd, file_sampling = f_sampling)
  }
  
  # Add the sampling event data (SED) so we also have the non-detection data (sampled places with no detection of our species)
  ebd_zf <- auk_zerofill(f_ebd, f_sampling, collapse = TRUE)
  
  time_to_decimal <- function(x) {
    x <- hms(x, quiet = TRUE)
    hour(x) + minute(x) / 60 + second(x) / 3600
  }
  # Rename
  names(ebd_zf)[names(ebd_zf) == "observation_count.x"] <- "observation_count"
  names(ebd_zf)[names(ebd_zf) == "scientific_name.x"] <- "scientific_name"
  
  ebd_zf$scientific_name.y <- NULL
  ebd_zf$observation_count.y <- NULL
  
  # Format
  ebd_zf <- ebd_zf %>% 
    mutate(
      observation_count = if_else(observation_count == "X", NA_character_, observation_count),
      observation_count = as.integer(observation_count),
      effort_distance_km = if_else(protocol_type != "Traveling", 0, effort_distance_km),
      time_observations_started = time_to_decimal(time_observations_started),
      year = year(observation_date),
      day_of_year = yday(observation_date)
    )
  
  # eBird recommends 6 hours and 10 km as the maximum duration and distance for a checklist
  # but for now we'll use 5 hours and 5 km
  ebd_zf_filtered <- ebd_zf %>% 
    filter(
      duration_minutes <= 5 * 60,
      effort_distance_km <= 5,
      year >= 2017,
      year <= 2022, # landcover variable only goes up to 2022-01-01T00:00:00Z
      number_observers <= 10
    )
  
  # Select the columns we want to keep in our final csv file
  ebird <- ebd_zf_filtered %>% 
    select(
      checklist_id, observer_id, sampling_event_identifier,
      scientific_name,
      observation_count, species_observed, 
      state_code, locality_id, latitude, longitude,
      protocol_type, all_species_reported,
      observation_date, year, day_of_year,
      time_observations_started, 
      duration_minutes, effort_distance_km,
      number_observers
    )
  
  # Store the data as .csv file
  write.csv(ebird, paste0("../data/", species[i + 1], "/ebd_", species[i + 1], "_breeding_spain_zf.csv"), na = "", row.names = FALSE)
}
