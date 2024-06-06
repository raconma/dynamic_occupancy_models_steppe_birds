rm(list=ls())

## Packages ---------------------------------------------------------------
install.load.package <- function(x) {
  if (!require(x, character.only = TRUE))
    install.packages(x, repos='http://cran.us.r-project.org')
  require(x, character.only = TRUE)
}
package_vec <- c(
  "readr","dplyr","tidyr", "raster", "rnaturalearth", "rnaturalearthdata", "readxl", "sf", "viridisLite",
  "HH", "terra", "dggridR", "AICcmodavg", "unmarked", "MuMIn", "auk", "psych", "regclass" # names of the packages required placed here as character objects
)

sapply(package_vec, install.load.package)

## Functionality ----------------------------------------------------------
`%nin%` <- Negate(`%in%`) # a function for negation of %in% function

# resolve namespace conflicts
select <- dplyr::select
filter <- dplyr::filter

# load dataset
occ <- read.csv("../../data/falnau/occ_falnau_dynamic_wide_table.csv")

names(occ)

selected_columns <- occ[, c("EVI",
                            "Land_Cover_Type_1_Percent_Class_0", 
                            "Land_Cover_Type_1_Percent_Class_13", 
                            "Majority_Land_Cover_Type_1", 
                            "NDVI", 
                            "pr", 
                            "tmmn", 
                            "tmmx")]

colnames(selected_columns) <- c("EVI", 
                                "LC_0", 
                                "LC_13", 
                                "Maj_LC", 
                                "NDVI", 
                                "pr", 
                                "tmmn", 
                                "tmmx")

var_cor <- cor(selected_columns)
var_dist <- abs(as.dist(var_cor))
var_clust <- hclust(1 - var_dist) # Cluster variables
#Plot correlation cluster
plot(var_clust, hang = -1)


pairs.panels(select(occ, c("Majority_Land_Cover_Type_1",
                               "Land_Cover_Type_1_Percent_Class_13",
                               "tmmn",
                               "NDVI",
                               "pr"
                               )), main = "Variables") # drop irrelevant columns

# VIF of predictor variables
testVIF <- select(occ, c("Majority_Land_Cover_Type_1",
                         "Land_Cover_Type_1_Percent_Class_13",
                         "tmmn",
                         "NDVI",
                         "pr"
))
HH::vif(testVIF)

# Standardize (rescaling) variables so they match magnitude order
occ_std <- occ %>% mutate_at(c("Majority_Land_Cover_Type_1",
                                       "Land_Cover_Type_1_Percent_Class_13",
                                       "tmmn",
                                       "NDVI",
                                       "pr"), ~(scale(.) %>% as.vector))

# Check if each variable follows a normal distribution (mean = 0, sd = 1)
mean(occ_std$Majority_Land_Cover_Type_1)
sd(occ_std$Majority_Land_Cover_Type_1)

mean(occ_std$Land_Cover_Type_1_Percent_Class_13)
sd(occ_std$Land_Cover_Type_1_Percent_Class_13)

mean(occ_std$tmmn)
sd(occ_std$tmmn)

mean(occ_std$NDVI)
sd(occ_std$NDVI)

mean(occ_std$pr)
sd(occ_std$pr)

# Check for duplicates, spatial subsampling
dup <- duplicated(occ_std$cells)
if(length(dup) > 0){
  # Only retain non-duplicated cells:
  occ_wide_clean <- occ_std[!duplicated(occ_std$cells),]
  
  # calculate the percent decrease in the number of sites
  1 - nrow(occ_wide_clean) / nrow(occ_std)
} else {
  occ_wide_clean <- occ_std
}


# Save results
write.csv(occ_wide_clean, "../../data/falnau/occ_falnau_dynamic_wide_table_clean.csv", row.names = FALSE)
