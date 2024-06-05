rm(list=ls())

## Packages ---------------------------------------------------------------
install.load.package <- function(x) {
  if (!require(x, character.only = TRUE))
    install.packages(x, repos='http://cran.us.r-project.org')
  require(x, character.only = TRUE)
}
package_vec <- c(
  "readr","dplyr","tidyr", "raster", "rnaturalearth", "rnaturalearthdata", "readxl", "sf", "viridisLite",
  "HH", "terra", "dggridR", "AICcmodavg", "unmarked", "MuMIn", "auk", "psych", "regclass", "lubridate" # names of the packages required placed here as character objects
)

sapply(package_vec, install.load.package)

## Functionality ----------------------------------------------------------
`%nin%` <- Negate(`%in%`) # a function for negation of %in% function

# resolve namespace conflicts
select <- dplyr::select
filter <- dplyr::filter

# load dataset
occ_wide_clean <- read.csv("../../data/tettet/occ_tettet_dynamic_wide_table_clean.csv")

names(occ_wide_clean)

# Survey covariates
duration <- as.matrix(occ_wide_clean[, c(21:26, 70:75, 119:124, 168:173, 217:222, 266:271, 315:320, 364:369, 413:418, 462:467, 511:516, 560:565)])
head(duration) 
effort <- as.matrix(occ_wide_clean[, c(27:32, 76:81, 125:130, 174:179, 223:228, 272:277, 321:326, 370:375, 419:424, 468:473, 517:522, 566:571)])
head(effort)
observers <- as.matrix(occ_wide_clean[, c(33:38, 82:87, 131:136, 180:185, 229:234, 278:283, 327:332, 376:381, 425:430, 474:479, 523:528, 572:577)])
head(observers)
time <- occ_wide_clean[, c(45:50, 94:99, 143:148, 192:197, 241:246, 290:295, 339:344, 388:393, 437:442, 486:491, 535:540, 584:589)]
head(time)
date <- as.matrix(occ_wide_clean[, c(39:44, 88:93, 137:142, 186:191, 235:240, 284:289, 333:338, 382:387, 431:436, 480:485, 529:534, 578:583)])
head(date)
days_of_year <- as.matrix(apply(date, 2, function(col) yday(ymd(col)))) # Transform date (yyyy-mm-dd) into numbre (day of the year)
head(days_of_year)

# Observations / Detection histories
detections <- occ_wide_clean[, c(2:7, 52:57, 101:106, 150:155, 199:204, 248:253, 297:302, 346:351, 395:400, 444:449, 493:498, 542:547)]
names(detections)

# convert all_star column to 1s and 0s
detections <- detections %>% mutate_if(is.logical, as.integer)
y.cross <- as.matrix(detections)  
head(y.cross) 

# Set all observations for which no survey covariate is available to NA. 
is.na(y.cross)
y.cross[is.na(time) != is.na(y.cross)] <- NA



# Yearly site covariates
Majority_Land_Cover_Type_1 <- as.matrix(occ_wide_clean[, c(11, 61, 110, 159, 208, 257, 306, 355, 404, 453, 502, 551)])
head(Majority_Land_Cover_Type_1)

Land_Cover_Type_1_Percent_Class_13 <- as.matrix(occ_wide_clean[, c(10, 60, 109, 158, 207, 256, 305, 354, 403, 452, 501, 550)])
head(Land_Cover_Type_1_Percent_Class_13)

tmmn <- as.matrix(occ_wide_clean[, c(18, 68, 117, 166, 215, 264, 313, 362, 411, 460, 509, 558)])
head(tmmn)

NDVI <- as.matrix(occ_wide_clean[, c(12, 62, 111, 160, 209, 258, 307, 356, 405, 454, 503, 552)])
head(NDVI)

pr <- as.matrix(occ_wide_clean[, c(16, 66, 115, 164, 213, 262, 311, 360, 409, 458, 507, 556)])
head(pr)

# Yearly site covariates
n <- 6720   # number of sites
T <- 12    # number of primary periods
J <- 6    # number of secondary periods

years <- data.frame(matrix(rep(2010:2021, each=n), n, T))
years <- data.frame(lapply(years, as.factor))

# Standardise survey covariates
time
time <- scale(time)
duration <- scale(duration)
effort <- scale(effort)
observers <- scale(observers)
days_of_year <- scale(days_of_year)

# Standardise site covariates
Majority_Land_Cover_Type_1 <- scale(Majority_Land_Cover_Type_1)
Land_Cover_Type_1_Percent_Class_13 <- scale(Land_Cover_Type_1_Percent_Class_13)
tmmn <- scale(tmmn)
NDVI <- scale(NDVI)
pr <- scale(pr)

# Verify matrix dimensions
dim(time)
dim(duration)
dim(effort)
dim(observers)
dim(days_of_year)
dim(y.cross)
dim(years)

dim(Majority_Land_Cover_Type_1)
dim(Land_Cover_Type_1_Percent_Class_13)
dim(tmmn)
dim(NDVI)
dim(pr)

print(unique(Majority_Land_Cover_Type_1))
print(unique(Land_Cover_Type_1_Percent_Class_13))

# Make unmarkedMultFrame for dynamic occupancy model
occ_umf <- unmarkedMultFrame(y = y.cross, # detection histories
                             yearlySiteCovs = list(Majority_Land_Cover_Type_1 = Majority_Land_Cover_Type_1, 
                                                   Land_Cover_Type_1_Percent_Class_13 = Land_Cover_Type_1_Percent_Class_13,
                                                   tmmn = tmmn, NDVI = NDVI, pr = pr), # list of site covariates
                             obsCovs = list(duration = duration, effort = effort,
                                            observers = observers, time = time, days_of_year = days_of_year), # list of survey covariates
                             numPrimary = 12) # number of primary time periods (here, number of years)