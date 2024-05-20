#load packages
library(raster)
library(ggplot2)
library(ggmap)
#library(ggsn)
library(viridisLite)
library(auk)
library(lubridate)
library(sf)
library(gridExtra)
library(tidyverse)
library(grid)
#library(velox) #mal
#library(MODIS) #mal
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

# resolve namespace conflicts
select <- dplyr::select
filter <- dplyr::filter

# create a new folder: "Data" within your working directory where we'll store data and save files
data_dir <- "./data/steppe_species"
if (!dir.exists(data_dir)) {
  dir.create(data_dir)
}

auk_set_awk_path()
ebd <- auk_ebd("./data/ebird_raw_mar2024/ebd_ES_smp_relMar-2024.txt")
#read in eBird data downloaded from eBird as a text (.txt) file
ebd <- auk_ebd("data/ebird_raw_mar2024/ebd_ES_smp_relMar-2024.txt", 
               file_sampling = "./data/ebird_raw_mar2024/ebd_sampling_relMar-2024/ebd_sampling_relMar-2024.txt")


especies <- list(
  c("Circus pygargus", "cirpyg"),
  c("Chersophilus duponti", "chedup"),
  c("Otis tarda", "otttar"),
  c("Falco naumanni", "falnau"),
  c("Pterocles alchata", "ptealc"),
  c("Pterocles orientalis", "pteori"),
  c("Tetrax tetrax", "tettet")
)

for (i in seq_along(especies)) {
  ebd_filters <- ebd %>% 
    auk_species(especies[[i]][1]) %>% 
    auk_country("ES") %>% 
    auk_date(date = c("*-04-01", "*-07-30")) %>% 
    auk_protocol(protocol = c("Stationary", "Traveling")) %>% 
    auk_complete()
  
  f_ebd <- file.path(data_dir, paste0("ebd_", especies[[i]][2], "_breeding.txt"))
  f_sampling <- file.path(data_dir, paste0("ebd_checklists_breeding_", especies[[i]][2], ".txt"))
  
  if (!file.exists(f_ebd)) {
    auk_filter(ebd_filters, file = f_ebd, file_sampling = f_sampling)
  }
  
  ebd_zf <- auk_zerofill(f_ebd, f_sampling, collapse = TRUE)
  
  time_to_decimal <- function(x) {
    x <- hms(x, quiet = TRUE)
    hour(x) + minute(x) / 60 + second(x) / 3600
  }
  
  names(ebd_zf)[names(ebd_zf) == "observation_count.x"] <- "observation_count"
  names(ebd_zf)[names(ebd_zf) == "scientific_name.x"] <- "scientific_name"
  
  ebd_zf$scientific_name.y <- NULL
  ebd_zf$observation_count.y <- NULL
  
  ebd_zf <- ebd_zf %>% 
    mutate(
      observation_count = if_else(observation_count == "X", NA_character_, observation_count),
      observation_count = as.integer(observation_count),
      effort_distance_km = if_else(protocol_type != "Traveling", 0, effort_distance_km),
      time_observations_started = time_to_decimal(time_observations_started),
      year = year(observation_date),
      day_of_year = yday(observation_date)
    )
  
  ebd_zf_filtered <- ebd_zf %>% 
    filter(
      duration_minutes <= 5 * 60,
      effort_distance_km <= 5,
      year >= 2010,
      number_observers <= 10
    )
  
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
  
  #csv2 separa por ; y csv por ,
  write_csv2(ebird, paste0("./data/ebd_", especies[[i]][2], "_breeding_spain_zf.csv"), na = "")
}
