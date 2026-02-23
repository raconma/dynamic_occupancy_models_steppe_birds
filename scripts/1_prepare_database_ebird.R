###############################################################################
# 1_prepare_database_ebird.R
#
# Purpose: Filter raw eBird data for 4 steppe bird species in Spain.
#          Produces zero-filled detection histories per species.
#
# Inputs:  data/raw/ebird_raw_mar2024/ebd_ES_smp_relMar-2024.txt
#          data/raw/ebird_raw_mar2024/ebd_sampling_relMar-2024.txt
#
# Outputs: data/processed/{species}/ebd_{species}_breeding_spain_zf.csv
#
# Note:    eBird data must be downloaded manually from https://ebird.org/data/download
#          See data-raw/get_data.R for download instructions.
###############################################################################

# -- Load packages --
library(here)
library(auk)
library(dplyr)
library(lubridate)

# Resolve namespace conflicts
select <- dplyr::select
filter <- dplyr::filter

# -- Read the raw eBird data --
ebd <- auk_ebd(
  here("data", "raw", "ebird_raw_mar2024", "ebd_ES_smp_relMar-2024.txt"),
  file_sampling = here("data", "raw", "ebird_raw_mar2024", "ebd_sampling_relMar-2024.txt")
)

# -- Declare the studied species --
# Named list: scientific name -> abbreviated code
species_list <- list(
  "Otis tarda"           = "otitar",
  "Pterocles alchata"    = "ptealc",
  "Pterocles orientalis" = "pteori",
  "Tetrax tetrax"        = "tettet"
)

# -- Helper: convert time string to decimal hours --
time_to_decimal <- function(x) {
  x <- hms(x, quiet = TRUE)
  hour(x) + minute(x) / 60 + second(x) / 3600
}

# -- Loop through every species --
for (sp_name in names(species_list)) {
  sp_code <- species_list[[sp_name]]
  message("Processing: ", sp_name, " (", sp_code, ")")

  # Create output directory
  sp_dir <- here("data", "processed", sp_code)
  if (!dir.exists(sp_dir)) dir.create(sp_dir, recursive = TRUE)

  # Define breeding season dates per species
  if (sp_name %in% c("Otis tarda", "Tetrax tetrax")) {
    date_range <- c("*-04-01", "*-06-30")
  } else {
    date_range <- c("*-05-01", "*-08-31")
  }
  message("  Breeding season: ", date_range[1], " to ", date_range[2])

  # Filter the raw eBird data
  ebd_filters <- ebd %>%
    auk_species(sp_name) %>%
    auk_country("ES") %>%
    auk_date(date = date_range) %>%
    auk_protocol(protocol = c("Stationary", "Traveling")) %>%
    auk_complete()

  # Write filtered data to intermediate text files
  f_ebd      <- file.path(sp_dir, paste0("ebd_", sp_code, "_breeding.txt"))
  f_sampling <- file.path(sp_dir, paste0("ebd_checklists_breeding_", sp_code, ".txt"))

  if (!file.exists(f_ebd)) {
    auk_filter(ebd_filters, file = f_ebd, file_sampling = f_sampling)
  }

  # Zero-fill: create detection/non-detection dataset
  ebd_zf <- auk_zerofill(f_ebd, f_sampling, collapse = TRUE)

  # Fix column names from zero-fill merge
  names(ebd_zf)[names(ebd_zf) == "observation_count.x"] <- "observation_count"
  names(ebd_zf)[names(ebd_zf) == "scientific_name.x"]   <- "scientific_name"
  ebd_zf$scientific_name.y   <- NULL
  ebd_zf$observation_count.y <- NULL

  # Format variables
  ebd_zf <- ebd_zf %>%
    mutate(
      observation_count        = if_else(observation_count == "X", NA_character_, observation_count),
      observation_count        = as.integer(observation_count),
      effort_distance_km       = if_else(protocol_type != "Traveling", 0, effort_distance_km),
      time_observations_started = time_to_decimal(time_observations_started),
      year        = year(observation_date),
      day_of_year = yday(observation_date)
    )

  # Apply effort filters (eBird best practices: max 5h, 5km, 10 observers)
  ebd_zf_filtered <- ebd_zf %>%
    filter(
      duration_minutes   <= 5 * 60,
      effort_distance_km <= 5,
      year >= 2017,
      year <= 2022,          # landcover data available up to 2022
      number_observers   <= 10
    )

  # Select final columns
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

  # Save output
  out_path <- here("data", "processed", sp_code,
                   paste0("ebd_", sp_code, "_breeding_spain_zf.csv"))
  write.csv(ebird, out_path, na = "", row.names = FALSE)
  message("  Saved: ", out_path, " (", nrow(ebird), " rows)")
}

message("Step 1 complete: eBird databases prepared for all species.")
