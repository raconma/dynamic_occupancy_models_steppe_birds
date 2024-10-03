knitr::opts_chunk$set(echo = TRUE, message=FALSE, comment = '', fig.width = 6, fig.height = 6)
#clear environment
rm(list=ls())

library(readr)
library(dplyr)
library(tidyr)

# resolve namespace conflicts
select <- dplyr::select
filter <- dplyr::filter


################################################################################
##################            LOAD DYNAMIC VARIABLES          ##################
################################################################################

# Load the data in wide format with the static variables
occ_wide_clean <- read.csv("../data/tettet/tettet_occ_wide_static.csv")

# Load the dynamic variables
dyn_var <- read.csv("../data/tettet/tettet_dynamic_variables.csv")

dyn_var <- dyn_var[, !names(dyn_var) %in% c("system.index", "cells", "latitude", "longitude", ".geo")]

# Merge the dynamic variables with the data
occ_wide_clean <- cbind(occ_wide_clean, dyn_var)

occ_wide_clean <- occ_wide_clean %>%
  drop_na(NDVI_2017, NDVI_2018, NDVI_2019, NDVI_2020, NDVI_2021, NDVI_2022)

write.csv(occ_wide_clean, "../data/tettet/tettet_occ_wide_dynamic.csv", row.names = FALSE)
