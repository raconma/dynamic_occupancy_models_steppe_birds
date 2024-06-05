rm(list=ls())

## Packages ---------------------------------------------------------------
install.load.package <- function(x) {
  if (!require(x, character.only = TRUE))
    install.packages(x, repos='http://cran.us.r-project.org')
  require(x, character.only = TRUE)
}
package_vec <- c(
  "readr","dplyr","tidyr", "raster", "rnaturalearth", "rnaturalearthdata", "readxl", "sf", "viridisLite",
  "HH", "terra", "dggridR", "AICcmodavg", "unmarked", "MuMIn", "auk" # names of the packages required placed here as character objects
)

sapply(package_vec, install.load.package)

## Functionality ----------------------------------------------------------
`%nin%` <- Negate(`%in%`) # a function for negation of %in% function

# resolve namespace conflicts
select <- dplyr::select
filter <- dplyr::filter

# load dataset
occ_var <- read.csv("../../data/tettet/ebd_tettet_breeding_spain_variables.csv")

occ_var <- filter(occ_var) %>%
  drop_na(EVI,
          Land_Cover_Type_1_Percent_Class_0,
          Land_Cover_Type_1_Percent_Class_13,
          Majority_Land_Cover_Type_1,
          NDVI,
          pr,
          tmmn,
          tmmx
  ) # Eliminate NAs from variables

# Create base table for wide format, creating "site" variable (each of our instances)
occ_filtered <- filter_repeat_visits(occ_var,
                            min_obs = 3, max_obs = 6,
                            annual_closure = TRUE,
                            date_var = "observation_date",
                            site_vars = c("locality_id", "observer_id"))

# Create wide format table for each year

occ_2010 <- filter(occ_filtered, year==2010)
occ_wide.2010 <- format_unmarked_occu(occ_2010, 
                                      site_id = "site", 
                                      response = "species_observed",
                                      site_covs = c("EVI",
                                                    "Land_Cover_Type_1_Percent_Class_0",
                                                    "Land_Cover_Type_1_Percent_Class_13",
                                                    "Majority_Land_Cover_Type_1",
                                                    "NDVI",
                                                    "latitude",
                                                    "locality_id",
                                                    "longitude",
                                                    "pr",
                                                    "state_code",
                                                    "tmmn",
                                                    "tmmx",
                                                    "year"
                                      ),
                                      obs_covs = c("duration_minutes", 
                                                   "effort_distance_km", 
                                                   "number_observers",
                                                   "observation_date",
                                                   "time_observations_started")) %>%
  rename_with(~ paste0(., ".2010"), -locality_id)

occ_2011 <- filter(occ_filtered, year==2011)
occ_wide.2011 <- format_unmarked_occu(occ_2011, 
                                      site_id = "site", 
                                      response = "species_observed",
                                      site_covs = c("EVI",
                                                    "Land_Cover_Type_1_Percent_Class_0",
                                                    "Land_Cover_Type_1_Percent_Class_13",
                                                    "Majority_Land_Cover_Type_1",
                                                    "NDVI",
                                                    "latitude",
                                                    "locality_id",
                                                    "longitude",
                                                    "pr",
                                                    "state_code",
                                                    "tmmn",
                                                    "tmmx",
                                                    "year"
                                      ),
                                      obs_covs = c("duration_minutes", 
                                                   "effort_distance_km", 
                                                   "number_observers",
                                                   "observation_date",
                                                   "time_observations_started")) %>%
  rename_with(~ paste0(., ".2011"), -locality_id)

occ_2012 <- filter(occ_filtered, year==2012)
occ_wide.2012 <- format_unmarked_occu(occ_2012, 
                                      site_id = "site", 
                                      response = "species_observed",
                                      site_covs = c("EVI",
                                                    "Land_Cover_Type_1_Percent_Class_0",
                                                    "Land_Cover_Type_1_Percent_Class_13",
                                                    "Majority_Land_Cover_Type_1",
                                                    "NDVI",
                                                    "latitude",
                                                    "locality_id",
                                                    "longitude",
                                                    "pr",
                                                    "state_code",
                                                    "tmmn",
                                                    "tmmx",
                                                    "year"
                                      ),
                                      obs_covs = c("duration_minutes", 
                                                   "effort_distance_km", 
                                                   "number_observers",
                                                   "observation_date",
                                                   "time_observations_started")) %>%
  rename_with(~ paste0(., ".2012"), -locality_id)

occ_2013 <- filter(occ_filtered, year==2013)
occ_wide.2013 <- format_unmarked_occu(occ_2013, 
                                      site_id = "site", 
                                      response = "species_observed",
                                      site_covs = c("EVI",
                                                    "Land_Cover_Type_1_Percent_Class_0",
                                                    "Land_Cover_Type_1_Percent_Class_13",
                                                    "Majority_Land_Cover_Type_1",
                                                    "NDVI",
                                                    "latitude",
                                                    "locality_id",
                                                    "longitude",
                                                    "pr",
                                                    "state_code",
                                                    "tmmn",
                                                    "tmmx",
                                                    "year"
                                      ),
                                      obs_covs = c("duration_minutes", 
                                                   "effort_distance_km", 
                                                   "number_observers",
                                                   "observation_date",
                                                   "time_observations_started")) %>%
  rename_with(~ paste0(., ".2013"), -locality_id)

occ_2014 <- filter(occ_filtered, year==2014)
occ_wide.2014 <- format_unmarked_occu(occ_2014, 
                                      site_id = "site", 
                                      response = "species_observed",
                                      site_covs = c("EVI",
                                                    "Land_Cover_Type_1_Percent_Class_0",
                                                    "Land_Cover_Type_1_Percent_Class_13",
                                                    "Majority_Land_Cover_Type_1",
                                                    "NDVI",
                                                    "latitude",
                                                    "locality_id",
                                                    "longitude",
                                                    "pr",
                                                    "state_code",
                                                    "tmmn",
                                                    "tmmx",
                                                    "year"
                                      ),
                                      obs_covs = c("duration_minutes", 
                                                   "effort_distance_km", 
                                                   "number_observers",
                                                   "observation_date",
                                                   "time_observations_started")) %>%
  rename_with(~ paste0(., ".2014"), -locality_id)

occ_2015 <- filter(occ_filtered, year==2015)
occ_wide.2015 <- format_unmarked_occu(occ_2015, 
                                      site_id = "site", 
                                      response = "species_observed",
                                      site_covs = c("EVI",
                                                    "Land_Cover_Type_1_Percent_Class_0",
                                                    "Land_Cover_Type_1_Percent_Class_13",
                                                    "Majority_Land_Cover_Type_1",
                                                    "NDVI",
                                                    "latitude",
                                                    "locality_id",
                                                    "longitude",
                                                    "pr",
                                                    "state_code",
                                                    "tmmn",
                                                    "tmmx",
                                                    "year"
                                      ),
                                      obs_covs = c("duration_minutes", 
                                                   "effort_distance_km", 
                                                   "number_observers",
                                                   "observation_date",
                                                   "time_observations_started")) %>%
  rename_with(~ paste0(., ".2015"), -locality_id)

occ_2016 <- filter(occ_filtered, year==2016)
occ_wide.2016 <- format_unmarked_occu(occ_2016, 
                                      site_id = "site", 
                                      response = "species_observed",
                                      site_covs = c("EVI",
                                                    "Land_Cover_Type_1_Percent_Class_0",
                                                    "Land_Cover_Type_1_Percent_Class_13",
                                                    "Majority_Land_Cover_Type_1",
                                                    "NDVI",
                                                    "latitude",
                                                    "locality_id",
                                                    "longitude",
                                                    "pr",
                                                    "state_code",
                                                    "tmmn",
                                                    "tmmx",
                                                    "year"
                                      ),
                                      obs_covs = c("duration_minutes", 
                                                   "effort_distance_km", 
                                                   "number_observers",
                                                   "observation_date",
                                                   "time_observations_started")) %>%
  rename_with(~ paste0(., ".2016"), -locality_id)

occ_2017 <- filter(occ_filtered, year==2017)
occ_wide.2017 <- format_unmarked_occu(occ_2017, 
                                      site_id = "site", 
                                      response = "species_observed",
                                      site_covs = c("EVI",
                                                    "Land_Cover_Type_1_Percent_Class_0",
                                                    "Land_Cover_Type_1_Percent_Class_13",
                                                    "Majority_Land_Cover_Type_1",
                                                    "NDVI",
                                                    "latitude",
                                                    "locality_id",
                                                    "longitude",
                                                    "pr",
                                                    "state_code",
                                                    "tmmn",
                                                    "tmmx",
                                                    "year"
                                      ),
                                      obs_covs = c("duration_minutes", 
                                                   "effort_distance_km", 
                                                   "number_observers",
                                                   "observation_date",
                                                   "time_observations_started")) %>%
  rename_with(~ paste0(., ".2017"), -locality_id)

occ_2018 <- filter(occ_filtered, year==2018)
occ_wide.2018 <- format_unmarked_occu(occ_2018, 
                                      site_id = "site", 
                                      response = "species_observed",
                                      site_covs = c("EVI",
                                                    "Land_Cover_Type_1_Percent_Class_0",
                                                    "Land_Cover_Type_1_Percent_Class_13",
                                                    "Majority_Land_Cover_Type_1",
                                                    "NDVI",
                                                    "latitude",
                                                    "locality_id",
                                                    "longitude",
                                                    "pr",
                                                    "state_code",
                                                    "tmmn",
                                                    "tmmx",
                                                    "year"
                                      ),
                                      obs_covs = c("duration_minutes", 
                                                   "effort_distance_km", 
                                                   "number_observers",
                                                   "observation_date",
                                                   "time_observations_started")) %>%
  rename_with(~ paste0(., ".2018"), -locality_id)

occ_2019 <- filter(occ_filtered, year==2019)
occ_wide.2019 <- format_unmarked_occu(occ_2019, 
                                      site_id = "site", 
                                      response = "species_observed",
                                      site_covs = c("EVI",
                                                    "Land_Cover_Type_1_Percent_Class_0",
                                                    "Land_Cover_Type_1_Percent_Class_13",
                                                    "Majority_Land_Cover_Type_1",
                                                    "NDVI",
                                                    "latitude",
                                                    "locality_id",
                                                    "longitude",
                                                    "pr",
                                                    "state_code",
                                                    "tmmn",
                                                    "tmmx",
                                                    "year"
                                      ),
                                      obs_covs = c("duration_minutes", 
                                                   "effort_distance_km", 
                                                   "number_observers",
                                                   "observation_date",
                                                   "time_observations_started")) %>%
  rename_with(~ paste0(., ".2019"), -locality_id)

occ_2020 <- filter(occ_filtered, year==2020)
occ_wide.2020 <- format_unmarked_occu(occ_2020, 
                                      site_id = "site", 
                                      response = "species_observed",
                                      site_covs = c("EVI",
                                                    "Land_Cover_Type_1_Percent_Class_0",
                                                    "Land_Cover_Type_1_Percent_Class_13",
                                                    "Majority_Land_Cover_Type_1",
                                                    "NDVI",
                                                    "latitude",
                                                    "locality_id",
                                                    "longitude",
                                                    "pr",
                                                    "state_code",
                                                    "tmmn",
                                                    "tmmx",
                                                    "year"
                                      ),
                                      obs_covs = c("duration_minutes", 
                                                   "effort_distance_km", 
                                                   "number_observers",
                                                   "observation_date",
                                                   "time_observations_started")) %>%
  rename_with(~ paste0(., ".2020"), -locality_id)

occ_2021 <- filter(occ_filtered, year==2021)
occ_wide.2021 <- format_unmarked_occu(occ_2021, 
                                      site_id = "site", 
                                      response = "species_observed",
                                      site_covs = c("EVI",
                                                    "Land_Cover_Type_1_Percent_Class_0",
                                                    "Land_Cover_Type_1_Percent_Class_13",
                                                    "Majority_Land_Cover_Type_1",
                                                    "NDVI",
                                                    "latitude",
                                                    "locality_id",
                                                    "longitude",
                                                    "pr",
                                                    "state_code",
                                                    "tmmn",
                                                    "tmmx",
                                                    "year"
                                      ),
                                      obs_covs = c("duration_minutes", 
                                                   "effort_distance_km", 
                                                   "number_observers",
                                                   "observation_date",
                                                   "time_observations_started")) %>%
  rename_with(~ paste0(., ".2021"), -locality_id)

occ.total <- occ_wide.2010 %>% 
  full_join(occ_wide.2011, by = c("locality_id"="locality_id"), multiple="all")
occ.total.1 <- occ.total %>% 
  full_join(occ_wide.2012, by = c("locality_id"="locality_id"), multiple="all")
occ.total.2 <- occ.total.1 %>% 
  full_join(occ_wide.2013, by = c("locality_id"="locality_id"), multiple="all") 
occ.total.3 <- occ.total.2 %>% 
  full_join(occ_wide.2014, by = c("locality_id"="locality_id"), multiple="all") 
occ.total.4 <- occ.total.3 %>% 
  full_join(occ_wide.2015, by = c("locality_id"="locality_id"), multiple="all") 
occ.total.5 <- occ.total.4 %>% 
  full_join(occ_wide.2016, by = c("locality_id"="locality_id"), multiple="all") 
occ.total.6 <- occ.total.5 %>% 
  full_join(occ_wide.2017, by = c("locality_id"="locality_id"), multiple="all") 
occ.total.7 <- occ.total.6 %>% 
  full_join(occ_wide.2018, by = c("locality_id"="locality_id"), multiple="all") 
occ.total.8 <- occ.total.7 %>% 
  full_join(occ_wide.2019, by = c("locality_id"="locality_id"), multiple="all")
occ.total.9 <- occ.total.8 %>% 
  full_join(occ_wide.2020, by = c("locality_id"="locality_id"), multiple="all") 
occ_wide_total <- occ.total.9 %>% 
  full_join(occ_wide.2021, by = c("locality_id"="locality_id"), multiple="all") 

dim(occ_wide_total)

# Remove duplicated columns
occ_wide_total <- occ_wide_total %>% distinct(locality_id, .keep_all = TRUE)
dim(occ_wide_total)

####### ¿Es necesario hacerlo?
occ_wide_total <- occ_wide_total %>% 
  mutate(EVI = coalesce(EVI.2010, EVI.2011, EVI.2012, EVI.2013, EVI.2014, EVI.2015, EVI.2016, EVI.2017, EVI.2018, EVI.2019, EVI.2020, EVI.2021),
         Land_Cover_Type_1_Percent_Class_0 = coalesce(Land_Cover_Type_1_Percent_Class_0.2010, Land_Cover_Type_1_Percent_Class_0.2011, Land_Cover_Type_1_Percent_Class_0.2012, Land_Cover_Type_1_Percent_Class_0.2013, Land_Cover_Type_1_Percent_Class_0.2014, Land_Cover_Type_1_Percent_Class_0.2015, Land_Cover_Type_1_Percent_Class_0.2016, Land_Cover_Type_1_Percent_Class_0.2017, Land_Cover_Type_1_Percent_Class_0.2018, Land_Cover_Type_1_Percent_Class_0.2019, Land_Cover_Type_1_Percent_Class_0.2020, Land_Cover_Type_1_Percent_Class_0.2021),
         Land_Cover_Type_1_Percent_Class_13 = coalesce(Land_Cover_Type_1_Percent_Class_13.2010, Land_Cover_Type_1_Percent_Class_13.2011, Land_Cover_Type_1_Percent_Class_13.2012, Land_Cover_Type_1_Percent_Class_13.2013, Land_Cover_Type_1_Percent_Class_13.2014, Land_Cover_Type_1_Percent_Class_13.2015, Land_Cover_Type_1_Percent_Class_13.2016, Land_Cover_Type_1_Percent_Class_13.2017, Land_Cover_Type_1_Percent_Class_13.2018, Land_Cover_Type_1_Percent_Class_13.2019, Land_Cover_Type_1_Percent_Class_13.2020, Land_Cover_Type_1_Percent_Class_13.2021),
         Majority_Land_Cover_Type_1 = coalesce(Majority_Land_Cover_Type_1.2010, Majority_Land_Cover_Type_1.2011, Majority_Land_Cover_Type_1.2012, Majority_Land_Cover_Type_1.2013, Majority_Land_Cover_Type_1.2014, Majority_Land_Cover_Type_1.2015, Majority_Land_Cover_Type_1.2016, Majority_Land_Cover_Type_1.2017, Majority_Land_Cover_Type_1.2018, Majority_Land_Cover_Type_1.2019, Majority_Land_Cover_Type_1.2020, Majority_Land_Cover_Type_1.2021),
         NDVI = coalesce(NDVI.2010, NDVI.2011, NDVI.2012, NDVI.2013, NDVI.2014, NDVI.2015, NDVI.2016, NDVI.2017, NDVI.2018, NDVI.2019, NDVI.2020, NDVI.2021),
         pr = coalesce(pr.2010, pr.2011, pr.2012, pr.2013, pr.2014, pr.2015, pr.2016, pr.2017, pr.2018, pr.2019, pr.2020, pr.2021),
         tmmn = coalesce(tmmn.2010, tmmn.2011, tmmn.2012, tmmn.2013, tmmn.2014, tmmn.2015, tmmn.2016, tmmn.2017, tmmn.2018, tmmn.2019, tmmn.2020, tmmn.2021),
         tmmx = coalesce(tmmx.2010, tmmx.2011, tmmx.2012, tmmx.2013, tmmx.2014, tmmx.2015, tmmx.2016, tmmx.2017, tmmx.2018, tmmx.2019, tmmx.2020, tmmx.2021)
         )

dim(occ_wide_total)

# Convert the detection histories in 1= presence/ 0= absence instead of TRUE/FALSE
cols <- sapply(occ_wide_total, is.logical) #Select columns that are logical (TRUe/FALSE)
occ_wide_total[,cols] <- lapply(occ_wide_total[,cols], as.numeric) # Transform to numeric; 1= TRUE, 0= FALSE.

# Export the table as csv
write.table(occ_wide_total, "../../data/tettet/occ_tettet_dynamic_wide_table.csv", sep = ",", row.names = FALSE)
