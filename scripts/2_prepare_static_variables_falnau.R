knitr::opts_chunk$set(echo = TRUE, message=FALSE, comment = '', fig.width = 6, fig.height = 6)
#clear environment
rm(list=ls())

library(unmarked)
library(here)
library(auk)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(ggplot2)
library(readr)
library(dplyr)
library(purrr)
library(knitr)
library(raster)
library(tidyr)
library(psych)
library(HH)
library(pander)
library(MuMIn)
library(stringr)
library(AICcmodavg)
library(ncdf4)
library(terra)
# resolve namespace conflicts
select <- dplyr::select
filter <- dplyr::filter

#set global options

theme_set(theme_bw())

set.seed(1)

occ_raw <- read.csv("../data/falnau/ebd_falnau_breeding_spain_zf.csv")

occ_raw <- occ_raw %>% 
  filter(year >= 2017,
         year <= 2022)

# Step 1. Check the data frame structure
str(occ_raw)

#format data
occ <- filter_repeat_visits(occ_raw, 
                            min_obs = 3, max_obs = 10,
                            annual_closure = TRUE,
                            date_var = "observation_date",
                            site_vars = c("locality_id", "observer_id"))

# Step 2. Check how many distinct sites are in the dataset.

n_distinct(occ$site)

# Step 3. Get a summary.  
# Checking that, e.g., missing values where we're not expecting them 
# and that the number of records for each level of factor variables
# and summary statistics for continuous variables makes senses.
summary(occ) 

# Step 4. Tabulations of relevant variables 
with(occ, table(species_observed , observation_count)) 
with(occ, table(species_observed , year)) 
# no indication of extreme counts, etc.

occ_species_observed <- occ %>% filter(species_observed == TRUE)

write.csv(occ_species_observed, "../data/falnau/falnau_occ_species_observed.csv", row.names = FALSE)


################################################################################
##################         LOAD BIOCLIMATIC VARIABLES         ##################
################################################################################

# Load the environmental variables
variables <- stack("../data/environmental_data/environmental_data_occ/variables_spain.grd")
names(variables)

# Extract the environmental variables at the site locations
occ_var <- occ %>% 
  cbind(as.data.frame(terra::extract(variables,
                                     occ[, c("longitude", "latitude")],
                                     cellnumbers=T )))

summary(occ_var) #check for NA in the environmental variables

# Number of sites  
n_distinct(occ_var$site)

occ_var <- occ_var %>% # Drop NA for the environmental variables
  drop_na(bio12,tree_cover) 

n_distinct(occ_var$site)

names(occ_var)

# Check for correlation between variables
var_cor <- cor(occ_var[, c(24:44)])
var_dist <- abs(as.dist(var_cor))
var_clust <- hclust(1-var_dist) #Cluster variables 
# Plot correlation cluster 
plot(var_clust)  

names(occ_var)

# Choose the variables to include in the model taking the correlation in account 
occ_var_std <- occ_var %>% mutate_at(c("bio1", "tree_cover", "bio2", "grass_cover", "bio12"), ~(scale(.) %>% as.vector))

# Check the mean and standard deviation of the variables
mean(occ_var_std$tree_cover)
sd(occ_var_std$tree_cover)


################################################################################
##################      CONVERT FROM LONG TO WIDE FORMAT      ##################
################################################################################

# Convert the data from long to wide format year by year
# It is useful to do it right now so we can append the dynamic variables more easily

occ_2017 <- filter(occ_var_std, year==2017)
occ_wide.2017 <- format_unmarked_occu(occ_2017, 
                                      site_id = "site", 
                                      response = "species_observed",
                                      site_covs = c("locality_id", "n_observations", "cells",
                                                    "latitude", "longitude", "bio1",
                                                    "tree_cover", "bio2",
                                                    "grass_cover", "bio12"),
                                      obs_covs = c("time_observations_started", 
                                                   "duration_minutes", 
                                                   "effort_distance_km", 
                                                   "number_observers", 
                                                   "protocol_type")) %>%
  rename_with(~ paste0(., ".2017"), -locality_id)

occ_2018 <- filter(occ_var_std, year==2018)
occ_wide.2018 <- format_unmarked_occu(occ_2018, 
                                      site_id = "site", 
                                      response = "species_observed",
                                      site_covs = c("locality_id", "n_observations", "cells",
                                                    "latitude", "longitude", "bio1",
                                                    "tree_cover", "bio2",
                                                    "grass_cover", "bio12"),
                                      obs_covs = c("time_observations_started", 
                                                   "duration_minutes", 
                                                   "effort_distance_km", 
                                                   "number_observers", 
                                                   "protocol_type")) %>%
  rename_with(~ paste0(., ".2018"), -locality_id)

occ_2019 <- filter(occ_var_std, year==2019)
occ_wide.2019 <- format_unmarked_occu(occ_2019, 
                                      site_id = "site", 
                                      response = "species_observed",
                                      site_covs = c("locality_id", "n_observations", "cells",
                                                    "latitude", "longitude", "bio1",
                                                    "tree_cover", "bio2",
                                                    "grass_cover", "bio12"),
                                      obs_covs = c("time_observations_started", 
                                                   "duration_minutes", 
                                                   "effort_distance_km", 
                                                   "number_observers", 
                                                   "protocol_type")) %>%
  rename_with(~ paste0(., ".2019"), -locality_id)

occ_2020 <- filter(occ_var_std, year==2020)
occ_wide.2020 <- format_unmarked_occu(occ_2020, 
                                      site_id = "site", 
                                      response = "species_observed",
                                      site_covs = c("locality_id", "n_observations", "cells",
                                                    "latitude", "longitude", "bio1",
                                                    "tree_cover", "bio2",
                                                    "grass_cover", "bio12"),
                                      obs_covs = c("time_observations_started", 
                                                   "duration_minutes", 
                                                   "effort_distance_km", 
                                                   "number_observers", 
                                                   "protocol_type")) %>%
  rename_with(~ paste0(., ".2020"), -locality_id)


occ_2021 <- filter(occ_var_std, year==2021)
occ_wide.2021 <- format_unmarked_occu(occ_2021, 
                                      site_id = "site", 
                                      response = "species_observed",
                                      site_covs = c("locality_id", "n_observations", "cells",
                                                    "latitude", "longitude", "bio1",
                                                    "tree_cover", "bio2",
                                                    "grass_cover", "bio12"),
                                      obs_covs = c("time_observations_started", 
                                                   "duration_minutes", 
                                                   "effort_distance_km", 
                                                   "number_observers", 
                                                   "protocol_type")) %>%
  rename_with(~ paste0(., ".2021"), -locality_id)

occ_2022 <- filter(occ_var_std, year==2022)
occ_wide.2022 <- format_unmarked_occu(occ_2022, 
                                      site_id = "site", 
                                      response = "species_observed",
                                      site_covs = c("locality_id", "n_observations", "cells",
                                                    "latitude", "longitude", "bio1",
                                                    "tree_cover", "bio2",
                                                    "grass_cover", "bio12"),
                                      obs_covs = c("time_observations_started", 
                                                   "duration_minutes", 
                                                   "effort_distance_km", 
                                                   "number_observers", 
                                                   "protocol_type")) %>%
  rename_with(~ paste0(., ".2022"), -locality_id)

# Append the data year by year
occ.total.1 <- occ_wide.2017 %>% 
  full_join(occ_wide.2018, by = c("locality_id"="locality_id"), multiple="all")
occ.total.2 <- occ.total.1 %>% 
  full_join(occ_wide.2019, by = c("locality_id"="locality_id"), multiple="all") 
occ.total.3 <- occ.total.2 %>% 
  full_join(occ_wide.2020, by = c("locality_id"="locality_id"), multiple="all")
occ.total.4 <- occ.total.3 %>% 
  full_join(occ_wide.2021, by = c("locality_id"="locality_id"), multiple="all")
occ_wide <- occ.total.4 %>% 
  full_join(occ_wide.2022, by = c("locality_id"="locality_id"), multiple="all") 

dim(occ_wide)
# Remove duplicated columns
occ_wide <- occ_wide %>% distinct(locality_id, .keep_all = TRUE)
dim(occ_wide)

#head(occ_wide)
# Convert the detection histories in 1 = presence / 0 = absence instead of TRUE/FALSE
cols <- sapply(occ_wide, is.logical) #Select columns that are logical (TRUE/FALSE)
occ_wide[,cols] <- lapply(occ_wide[,cols], as.numeric) # Transform to numeric; 1 = TRUE, 0 = FALSE.

occ_wide <- occ_wide %>% 
  rename( bio1 = bio1.2017, tree_cover = tree_cover.2017, bio2 = bio2.2017,
          grass_cover = grass_cover.2017, bio12 = bio12.2017, cells = cells.2017,
          latitude = latitude.2017, longitude = longitude.2017)

occ_wide <- occ_wide %>% 
  mutate(bio1 = coalesce(bio1.2018, bio1.2019, bio1.2020, bio1.2021, bio1.2022),
         tree_cover = coalesce(tree_cover.2018, tree_cover.2019, tree_cover.2020, tree_cover.2021, tree_cover.2022),
         bio2 = coalesce(bio2.2018, bio2.2019, bio2.2020, bio2.2021, bio2.2022),
         grass_cover = coalesce(grass_cover.2018, grass_cover.2019, grass_cover.2020, grass_cover.2021, grass_cover.2022),
         bio12 = coalesce(bio12.2018, bio12.2019, bio12.2020, bio12.2021, bio12.2022),
         cells = coalesce(cells.2018, cells.2019, cells.2020, cells.2021, cells.2022),
         latitude = coalesce(latitude.2018, latitude.2019, latitude.2020, latitude.2021, latitude.2022),
         longitude = coalesce(longitude.2018, longitude.2019, longitude.2020, longitude.2021, longitude.2022),
  )


# Check for duplicates
duplicated(occ_wide$cells)
# Only retain non-duplicated cells:
occ_wide_clean <- occ_wide[!duplicated(occ_wide$cells),]

occ_wide_clean <- occ_wide_clean %>% # Drop NA
  drop_na(cells) 

dim(occ_wide_clean)

# calculate the percent decrease in the number of sites
1 - nrow(occ_wide_clean) / nrow(occ_wide)

# Store the sites in wide format for downloading the dynamic variables with Google Earth Engine
occ_wide_latlong <- occ_wide_clean %>% 
  select(cells, latitude, longitude)
write.csv(occ_wide_latlong, "../data/falnau/falnau_occ_wide_latlong.csv", row.names = FALSE)

write.csv(occ_wide_clean, "../data/falnau/falnau_occ_wide_static.csv", row.names = FALSE)