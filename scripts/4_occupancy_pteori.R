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

################################################################################
##################       PREPARE THE DATA FOR THE MODELS      ##################
################################################################################

occ_wide_clean <- read.csv("../data/pteori/pteori_occ_wide_dynamic.csv")

names(occ_wide_clean)

# Survey covariates
duration <- as.matrix(occ_wide_clean[, c(32:41, 102:111, 172:181, 242:251, 312:321, 382:391)])
head(duration) 
effort <- as.matrix(occ_wide_clean[, c(42:51, 112:121, 182:191, 252:261, 322:331, 392:401)])
head(effort)
observers <- as.matrix(occ_wide_clean[, c(52:61, 122:131, 192:201, 262:271, 332:341, 402:411)])
head(observers)
time <- occ_wide_clean[, c(22:31, 92:101, 162:171, 232:241, 302:311, 372:381)]
head(time)

# Observations / Detection histories
detections <- occ_wide_clean[, c(2:11, 73:82, 143:152, 213:222, 283:292, 353:362)]
names(detections)

# convert all_star column to 1s and 0s
detections <- detections %>% mutate_if(is.logical, as.integer)
y.cross <- as.matrix(detections)  
head(y.cross) 

# Set all observations for which no survey covariate is available to NA. 
is.na(y.cross)
y.cross[is.na(time) != is.na(y.cross)] <- NA

# Site covariates
siteCovs <- occ_wide_clean[, c(17:21)]
head(siteCovs)

# Yearly site covariates
EVI <- occ_wide_clean[, c(422:427)]
head(EVI)
Land_Cover_Type_1_Percent_Class_0 <- occ_wide_clean[, c(428:433)]
head(Land_Cover_Type_1_Percent_Class_0)
Land_Cover_Type_1_Percent_Class_10 <- occ_wide_clean[, c(434:439)]
head(Land_Cover_Type_1_Percent_Class_10)
Land_Cover_Type_1_Percent_Class_12 <- occ_wide_clean[, c(440:445)]
head(Land_Cover_Type_1_Percent_Class_12)
Land_Cover_Type_1_Percent_Class_13 <- occ_wide_clean[, c(446:451)]
head(Land_Cover_Type_1_Percent_Class_13)
Land_Cover_Type_1_Percent_Class_14 <- occ_wide_clean[, c(452:457)]
head(Land_Cover_Type_1_Percent_Class_14)
Land_Cover_Type_1_Percent_Class_6 <- occ_wide_clean[, c(458:463)]
head(Land_Cover_Type_1_Percent_Class_6)
Land_Cover_Type_1_Percent_Class_7 <- occ_wide_clean[, c(464:469)]
head(Land_Cover_Type_1_Percent_Class_7)
NDVI <- occ_wide_clean[, c(476:481)]
head(NDVI)
pr <- occ_wide_clean[, c(482:487)]
head(pr)
tmmn <- occ_wide_clean[, c(488:493)]
head(tmmn)
tmmx <- occ_wide_clean[, c(494:499)]
head(tmmx)

# Standardise survey covariates
time <- scale(time)
duration <- scale(duration)
effort <- scale(effort)
observers <- scale(observers)

# Standardise site covariates
siteCovs <- scale(siteCovs)

# Standardise yearly site covariates
EVI <- scale(EVI)
Land_Cover_Type_1_Percent_Class_0 <- scale(Land_Cover_Type_1_Percent_Class_0)
Land_Cover_Type_1_Percent_Class_13 <- scale(Land_Cover_Type_1_Percent_Class_13)
NDVI <- scale(NDVI)
pr <- scale(pr)
tmmn <- scale(tmmn)
tmmx <- scale(tmmx)

# Verify matrix dimensions
dim(time)
dim(duration)
dim(effort)
dim(observers)
dim(y.cross)

dim(siteCovs)
dim(NDVI)

# Yearly site covariates
n <- nrow(occ_wide_clean)   # number of sites
T <- 6     # number of primary periods
J <- 10    # number of secondary periods

years <- data.frame(matrix(rep(2017:2022, each=n), n, T))
years <- data.frame(lapply(years, as.factor))

dim(years)

# Make unmarkedMultFrame for dynamic occupancy model
occ_umf <- unmarkedMultFrame(y = y.cross, # detection histories
                             siteCovs = data.frame(siteCovs),  # (static) site covariates
                             yearlySiteCovs = list(years = years,
                                                   EVI = EVI,
                                                   Land_Cover_Type_1_Percent_Class_0 = Land_Cover_Type_1_Percent_Class_0,
                                                   Land_Cover_Type_1_Percent_Class_6 = Land_Cover_Type_1_Percent_Class_6,
                                                   Land_Cover_Type_1_Percent_Class_7 = Land_Cover_Type_1_Percent_Class_7,
                                                   Land_Cover_Type_1_Percent_Class_10 = Land_Cover_Type_1_Percent_Class_10,
                                                   Land_Cover_Type_1_Percent_Class_12 = Land_Cover_Type_1_Percent_Class_12,
                                                   Land_Cover_Type_1_Percent_Class_13 = Land_Cover_Type_1_Percent_Class_13,
                                                   Land_Cover_Type_1_Percent_Class_14 = Land_Cover_Type_1_Percent_Class_14,
                                                   NDVI = NDVI,
                                                   pr = pr,
                                                   tmmn = tmmn,
                                                   tmmx = tmmx), # list of yearly (dynamic) site covariates
                             obsCovs = list(duration = duration, effort = effort,
                                            observers = observers, time = time), # list of survey covariates
                             numPrimary = 6) # number of primary time periods (here, number of years)

################################################################################
##################                    MODELS                  ##################
################################################################################

#psi occup
#gamma colon
#epsilon extin
#p detect

Mod.final <- colext(psiformula = ~ bio2 + tree_cover + grass_cover, 
                    gammaformula = ~ Land_Cover_Type_1_Percent_Class_0 + Land_Cover_Type_1_Percent_Class_13 + NDVI + pr + tmmx, 
                    epsilonformula = ~ Land_Cover_Type_1_Percent_Class_6 + Land_Cover_Type_1_Percent_Class_7 + Land_Cover_Type_1_Percent_Class_12 + Land_Cover_Type_1_Percent_Class_13 + tmmn, 
                    pformula = ~ time + duration + effort + observers, 
                    data = occ_umf)

# Inspect the fitted model:
summary(Mod.final)

model_occ <- Mod.final

# Check the model fit
GOF <- parboot(Mod.final, # global (most parameterised) model
               nsim = 10) # number of simulations / better with 1000

GOF # view output

occ_gof <- mb.gof.test(Mod.final, nsim = 3, plot.hist = FALSE)


################################################################################
##################              RESPONSE VARIABLES            ##################
################################################################################

# Occupancy
occformulaList<-c(
  ~bio2,
  ~tree_cover,
  ~grass_cover)
#get obsCovs as data.frame
siteCovs.df<-as.data.frame(occ_umf@siteCovs)

#now fit models

#create empty list to save predicted detection Probs
occCovariatePredOcc<-list()

#for loop to iterate over detection covariates and produce predicted estimates for detection prob
for(i in 1:length(occformulaList)){
  print(occformulaList[[i]])
  new.occ_model <- colext(psiformula = occformulaList[[i]],
                          gammaformula = ~ 1,
                          epsilonformula = ~ 1, 
                          pformula = ~ 1, 
                          data = occ_umf)
  
  #get det covariate name
  varName<-gsub("~","",as.character(occformulaList[[i]]))[2]
  
  #get range of values from data to simulate new data
  varData<-siteCovs.df[,names(siteCovs.df)==varName]
  
  #simulate new data
  sim.data<-data.frame(seq(min(varData,na.rm=TRUE), max(varData,na.rm=TRUE),length.out=100))
  colnames(sim.data)<-c(varName)
  #get predicted values using sim.data
  pred.df<-predict(new.occ_model, type="psi", newdata=sim.data, appendData=TRUE)
  colnames(pred.df)[5]<-"occ_cov"
  
  #add column with covariate name
  pred.df$CovariateName<-varName
  
  occCovariatePredOcc<-rbind(occCovariatePredOcc, pred.df)
  
}

#make CovariateName a factor
occCovariatePredOcc$CovariateName<-as.factor(occCovariatePredOcc$CovariateName)

# Plot it
occPlot<-ggplot(data=occCovariatePredOcc)+
  geom_ribbon(aes(x=occ_cov, ymin=lower, ymax=upper), fill="gray80")+
  geom_line(aes(x=occ_cov, y=Predicted),color="seagreen4")+
  labs(x="Covariate Values", y="Probability of Site Occupancy")+
  #ylim(0,1)+
  theme(panel.border=element_rect(color="black",fill="transparent"), panel.background = element_rect(fill="white"))

occPlotFacet<-occPlot+facet_wrap(.~CovariateName, scales="free",ncol=3)
occPlotFacet

################################################################################

# Colonization

colformulaList<-c(
  ~Land_Cover_Type_1_Percent_Class_0,
  ~Land_Cover_Type_1_Percent_Class_13,
  ~NDVI,
  ~pr,
  ~tmmx
)

#get obsCovs as data.frame
siteCovs.df<-as.data.frame(occ_umf@yearlySiteCovs)

#now fit models

#create empty list to save predicted detection Probs
colCovariatePred<-list()

#for loop to iterate over detection covariates and produce predicted estimates for detection prob
for(i in 1:length(colformulaList)){
  new.occ_model <- colext(psiformula = ~ 1,
                          gammaformula =  colformulaList[[i]],
                          epsilonformula = ~ 1, 
                          pformula = ~ 1,
                          data = occ_umf)
  
  #get col covariate name
  varName<-gsub("~","",as.character(colformulaList[[i]]))[2]
  
  #get range of values from data to simulate new data
  varData<-siteCovs.df[,names(siteCovs.df)==varName]
  
  #simulate new data
  sim.data<-data.frame(seq(min(varData,na.rm=TRUE), max(varData,na.rm=TRUE),length.out=100))
  colnames(sim.data)<-c(varName)
  #get predicted values using sim.data
  pred.df<-predict(new.occ_model, type="col", newdata=sim.data, appendData=TRUE)
  colnames(pred.df)[5]<-"col_cov"
  
  #add column with covariate name
  pred.df$CovariateName<-varName
  
  colCovariatePred<-rbind(colCovariatePred, pred.df)
  
}

#make CovariateName a factor
colCovariatePred$CovariateName<-as.factor(colCovariatePred$CovariateName)

# Plot it
colPlot<-ggplot(data=colCovariatePred)+
  geom_ribbon(aes(x=col_cov, ymin=lower, ymax=upper), fill="gray80")+
  geom_line(aes(x=col_cov, y=Predicted),color="royalblue3")+
  labs(x="Covariate Values", y="Colonization Probability")+
  #ylim(0,1)+
  theme(panel.border=element_rect(color="black",fill="transparent"), panel.background = element_rect(fill="white"))

colPlotFacet<-colPlot+facet_wrap(.~CovariateName, scales="free",ncol=3)
colPlotFacet

################################################################################

# Extinction

extformulaList<-c(
  ~Land_Cover_Type_1_Percent_Class_6,
  ~Land_Cover_Type_1_Percent_Class_7,
  ~Land_Cover_Type_1_Percent_Class_12,
  ~Land_Cover_Type_1_Percent_Class_13,
  ~tmmn
)

#get obsCovs as data.frame
siteCovs.df<-as.data.frame(occ_umf@yearlySiteCovs)

#now fit models

#create empty list to save predicted detection Probs
extCovariatePred<-list()

#for loop to iterate over detection covariates and produce predicted estimates for detection prob
for(i in 1:length(extformulaList)){
  new.occ_model <- colext(psiformula = ~ 1,
                          gammaformula = ~ 1,
                          epsilonformula = extformulaList[[i]], 
                          pformula = ~ 1,
                          data = occ_umf)
  
  #get col covariate name
  varName<-gsub("~","",as.character(extformulaList[[i]]))[2]
  
  #get range of values from data to simulate new data
  varData<-siteCovs.df[,names(siteCovs.df)==varName]
  
  #simulate new data
  sim.data<-data.frame(seq(min(varData,na.rm=TRUE), max(varData,na.rm=TRUE),length.out=100))
  colnames(sim.data)<-c(varName)
  #get predicted values using sim.data
  pred.df<-predict(new.occ_model, type="ext", newdata=sim.data, appendData=TRUE)
  colnames(pred.df)[5]<-"ext_cov"
  
  #add column with covariate name
  pred.df$CovariateName<-varName
  
  extCovariatePred<-rbind(extCovariatePred, pred.df)
  
}

#make CovariateName a factor
extCovariatePred$CovariateName<-as.factor(extCovariatePred$CovariateName)

# Plot it
extPlot<-ggplot(data=extCovariatePred)+
  geom_ribbon(aes(x=ext_cov, ymin=lower, ymax=upper), fill="gray80")+
  geom_line(aes(x=ext_cov, y=Predicted),color="royalblue3")+
  labs(x="Covariate Values", y="Extinction Probability")+
  #ylim(0,1)+
  theme(panel.border=element_rect(color="black",fill="transparent"), panel.background = element_rect(fill="white"))

extPlotFacet<-extPlot+facet_wrap(.~CovariateName, scales="free",ncol=3)
extPlotFacet

################################################################################

# Detection 
detformulaList<-c(
  ~time,
  ~duration,
  ~effort,
  ~observers
)


#get obsCovs as data.frame
obsCovs.df<-as.data.frame(occ_umf@obsCovs)
#now fit models

#create empty list to save predicted detection Probs
DetCovariatePred<-list()

#for loop to iterate over detection covariates and produce predicted estimates for detection prob
for(i in 1:length(detformulaList)){
  print(detformulaList[[i]])
  new.occ_model <- colext(psiformula = ~ 1,
                          gammaformula = ~ 1,
                          epsilonformula = ~ 1, 
                          pformula = detformulaList[[i]],
                          data = occ_umf)
  
  #get det covariate name
  varName<-gsub("~","",as.character(detformulaList[[i]]))[2]
  
  #get range of values from data to simulate new data
  varData<-obsCovs.df[,names(obsCovs.df)==varName]
  
  #simulate new data
  sim.data<-data.frame(seq(min(varData,na.rm=TRUE), max(varData,na.rm=TRUE),length.out=100))
  colnames(sim.data)<-c(varName)
  #get predicted values using sim.data
  pred.df<-predict(new.occ_model, type="det", newdata=sim.data, appendData=TRUE)
  colnames(pred.df)[5]<-"det_cov"
  
  #add column with covariate name
  pred.df$CovariateName<-varName
  
  DetCovariatePred<-rbind(DetCovariatePred, pred.df)
  
}

#make CovariateName a factor
DetCovariatePred$CovariateName<-as.factor(DetCovariatePred$CovariateName)

# Plot it
detPlot<-ggplot(data=DetCovariatePred)+
  geom_ribbon(aes(x=det_cov, ymin=lower, ymax=upper), fill="gray80")+
  geom_line(aes(x=det_cov, y=Predicted),color="royalblue3")+
  labs(x="Covariate Values", y="Detection Probability")+
  #ylim(0,1)+
  theme(panel.border=element_rect(color="black",fill="transparent"), panel.background = element_rect(fill="white"))

detPlotFacet<-detPlot+facet_wrap(.~CovariateName, scales="free",ncol=3)
detPlotFacet

################################################################################
##################               Occupation map               ##################
################################################################################

# Load variables
variables <- brick(stack(here("../data/environmental_data/environmental_data_occ/variables_spain.grd")))
# Select the variables that we used to calibrate our occupancy models
variables_selection <- c("bio2", "tree_cover", "grass_cover")
variables.sel <- variables[[variables_selection]]


# Transform in a dataset that we can use to predict our model
p_variables <- data.frame(rasterToPoints(variables.sel) )
# prediction surface
pred_surface <- p_variables
# Note that should have the same names of the variables
names(pred_surface)[names(pred_surface) == 'x'] <- 'longitude'
names(pred_surface)[names(pred_surface) == 'y'] <- 'latitude'
# Drop NA from the environmental variables. 
pred_surface <- pred_surface %>% 
  drop_na(tree_cover, bio2)


# We need to apply the same re-scaling of our covariates as we did before, so that the new data are at the same scale as our coefficient estimates.

pred_surface_std <- pred_surface %>%
  mutate_at(c("bio2", "tree_cover", "grass_cover"), ~(scale(.) %>% as.vector))

occ_pred <- predict(model_occ, 
                    newdata = as.data.frame(pred_surface_std[,c("bio2", "tree_cover", "grass_cover")]), type = "psi")


pred_occ <- bind_cols(pred_surface_std, occ_prob = occ_pred$Predicted,occ_se = occ_pred$SE) %>% 
  dplyr::select(longitude, latitude, occ_prob, occ_se)

#save result
write.csv(pred_occ,here("../data/pteori/occ_pteori_prediction.csv"), row.names = FALSE)

pred_occ <- read_csv(here("../data/pteori/occ_pteori_prediction.csv"))
#We'll want to set the map projection so every time we plot data (vector or raster), everything will line up correctly
map_proj <- "+proj=longlat +datum=WGS84 +no_defs"

r_pred <- pred_occ %>% 
  # convert to spatial features
  st_as_sf(coords = c("longitude", "latitude"), crs = map_proj) %>% 
  st_transform(crs = raster::projection(variables.sel[[1]])) %>% 
  # rasterize
  rasterize(variables.sel[[1]])

# Select the prediction with the standard error
r_pred <- r_pred[[c("occ_prob", "occ_se")]]


spain <- ne_countries(country = 'spain', scale = "medium", returnclass = "sf")
spain_crop <- rmapshaper::ms_filter_islands(spain,
                                            min_area = 100000000000,
                                            drop_null_geometries=TRUE)


# project predictions
r_pred_proj <- projectRaster(r_pred, crs = map_proj, method = "ngb")
r_pred_proj_crop<- crop(r_pred_proj, spain_crop)
r_pred_proj_DE <-mask(r_pred_proj_crop,spain_crop)

#First, with occupancy
r_pred_occu_DE <- subset(r_pred_proj_DE, 1, drop=TRUE)
# Convert the landscape data RasterLayer objects to data frames for ggplot
r_pred_occu_DE.df <- as.data.frame(r_pred_occu_DE, xy = TRUE, na.rm = TRUE)

# Now, we can add the Standard error in the predictions
r_pred_occuSE_DE<-subset(r_pred_proj_DE, 2, drop=TRUE)
# Convert the landscape data RasterLayer objects to data frames for ggplot
r_pred_occuSE_DE.df <- as.data.frame(r_pred_occuSE_DE, xy = TRUE, na.rm = TRUE)

#add Statistic Colum and standardize column names
r_pred_occu_DE.df$Statistic<-"Mean p(Occupancy)"
colnames(r_pred_occu_DE.df)[3]<-"Probability"
r_pred_occuSE_DE.df$Statistic<-"SE"
colnames(r_pred_occuSE_DE.df)[3]<-"Probability"

#combine data.frames
r_pred_comb<-rbind(r_pred_occu_DE.df, r_pred_occuSE_DE.df)
#make Statistic a factor
r_pred_comb$Statistic<-as.factor(r_pred_comb$Statistic)

# Plot the maps
PrOccuMap<-ggplot()+
  geom_tile(data=r_pred_comb, aes(x=x, y=y, fill=Probability, color=Probability))+
  #geom_map(data=nh.poly, map=nh.poly, aes(y=lat, x=long, map_id=id),color="black",fill="transparent",alpha=0.8, inherit.aes = FALSE)+
  scale_fill_viridis_c(name="p(occupancy)")+
  scale_color_viridis_c(name="p(occupancy)")+
  theme(panel.border=element_rect(color="black",fill="transparent"))+
  labs(x="Longitude",y="Latitude")+
  guides(alpha="none", color="none")
PrOccuMap.facet <- PrOccuMap+facet_grid(~Statistic)
PrOccuMap.facet

################################################################################

raster_initial <- rasterFromXYZ(r_pred_comb[, c("x", "y", "Probability")])
crs(raster_initial) <- CRS("+init=EPSG:25830")
writeRaster(raster_initial, "../data/pteori/pteori_OccuMap.tif", format = "GTiff", overwrite = TRUE)

################################################################################

# Simulating occupancy # We use calibrated data for the simulations instead of the rasters, because we do not have maps for the dynamic year covariates
# Predict initial occupancy probability per site
summary(occ_umf)
model_prediction <- Mod.final
summary(model_prediction)
#Run simulations
#S<-nrow(pred_surface_std) # Sites

predict_data <- occ_wide_clean %>% 
  drop_na("NDVI_2017", "NDVI_2018"   ,"NDVI_2019" , "NDVI_2020" , "NDVI_2021" , "NDVI_2022")


S<-nrow(predict_data) # Sites
T<-occ_umf@numPrimary # times (i.e. primary sample periods)


psi1<-predict(model_prediction, type="psi",newdata=as.data.frame(predict_data[,c("bio2", 
                                                                                 "tree_cover",
                                                                                 "grass_cover")]))[1] # initial occ

#Get colonisation & extinction predictions using yearly + site covariates
col<-data.frame(matrix(ncol=T,nrow=S))
ext<-data.frame(matrix(ncol=T,nrow=S))

# for col & ext need to send the each year separately
for(i in 1:T){
  pp <- c(2017:2022)
  #spXdets_sdf<-detGGA_sdf[detGGA_sdf$CommonName ==speciesName,]
  Land_Cover_Type_1_Percent_Class_0 <- as.data.frame(occ_wide_clean[,c(paste0('Land_Cover_Type_1_Percent_Class_0_',(pp[i])))])
  Land_Cover_Type_1_Percent_Class_13 <- as.data.frame(occ_wide_clean[,c(paste0('Land_Cover_Type_1_Percent_Class_13_',(pp[i])))])
  NDVI <- as.data.frame(occ_wide_clean[,c(paste0('NDVI_',(pp[i])))])
  pr <- as.data.frame(occ_wide_clean[,c(paste0('pr_',(pp[i])))])
  tmmx <- as.data.frame(occ_wide_clean[,c(paste0('tmmx_',(pp[i])))])
  new.dat<-cbind(Land_Cover_Type_1_Percent_Class_0,
                 Land_Cover_Type_1_Percent_Class_13,NDVI,pr,tmmx) %>% 
    drop_na() %>% 
    scale(.) %>% 
    as.data.frame(.)
  
  names(new.dat) <- c("Land_Cover_Type_1_Percent_Class_0",
                      "Land_Cover_Type_1_Percent_Class_13",
                      "NDVI",
                      "pr",
                      "tmmx"
  )
  if(nrow(new.dat) > nrow(predict_data)) {
    new.dat <- new.dat[c(1:nrow(predict_data)),]
  } else {
    new.dat <- new.dat
  }
  col[,i] <- predict(model_prediction, type = "col", newdata=new.dat)[1]
}

# for col & ext need to send the each year separately
for(i in 1:T){
  pp <- c(2017:2022)
  #spXdets_sdf<-detGGA_sdf[detGGA_sdf$CommonName ==speciesName,]
  Land_Cover_Type_1_Percent_Class_6 <- as.data.frame(occ_wide_clean[,c(paste0('Land_Cover_Type_1_Percent_Class_6_',(pp[i])))])
  Land_Cover_Type_1_Percent_Class_7 <- as.data.frame(occ_wide_clean[,c(paste0('Land_Cover_Type_1_Percent_Class_7_',(pp[i])))])  
  Land_Cover_Type_1_Percent_Class_12 <- as.data.frame(occ_wide_clean[,c(paste0('Land_Cover_Type_1_Percent_Class_12_',(pp[i])))])  
  Land_Cover_Type_1_Percent_Class_13 <- as.data.frame(occ_wide_clean[,c(paste0('Land_Cover_Type_1_Percent_Class_13_',(pp[i])))])  
  tmmn <- as.data.frame(occ_wide_clean[,c(paste0('tmmn_',(pp[i])))])  
  new.dat<-cbind(Land_Cover_Type_1_Percent_Class_6,Land_Cover_Type_1_Percent_Class_7,
                 Land_Cover_Type_1_Percent_Class_12,Land_Cover_Type_1_Percent_Class_13,tmmn) %>% 
    drop_na() %>% 
    scale(.) %>% 
    as.data.frame(.)
  
  names(new.dat) <- c("Land_Cover_Type_1_Percent_Class_6",
                      "Land_Cover_Type_1_Percent_Class_7",
                      "Land_Cover_Type_1_Percent_Class_12", 
                      "Land_Cover_Type_1_Percent_Class_13",
                      "tmmn")
  if(nrow(new.dat) > nrow(predict_data)) {
    new.dat <- new.dat[c(1:nrow(predict_data)),]
  } else {
    new.dat <- new.dat
  }
  ext[,i]<- predict(model_prediction, type = "ext", newdata=new.dat)[1]
}

## run stochastic sims - entire time period. 
nsim<-5000
Zs<-list()  
prev<-rep(NA,nsim)
prevT<-matrix(NA,nsim,T)
for (nn in 1:nsim){
  Z<-matrix(NA,S,T)
  Z[,1]<-(runif(S)<psi1)*1  #simulate starting state
  for (ii in 2:T){
    tmp <- Z[,ii-1]*(1-ext[,ii])+(1-Z[,ii-1])*col[,ii]
    Z[,ii] <- (runif(S)<tmp)*1
  }
  Zs[[nn]]<-Z
  prev[nn]<-mean(Z[,T])  #calculate prevalence at T, as a way to summarize results
  prevT[nn,]<-colMeans(Z)
}

hist(prev,xlim=c(0,0.2),col=rgb(1,0,0,0.3))  #plot hist of prevalences across sims
abline(v=mean(prev),col="red",lw=3)  #and their mean
#abline(v=psis[T],col="blue",lty=2,lw=2)  #and the probability obtained directly

##Also plot prev over time (+ var)
#pdf(file=paste0(dir.out,speciesName,"/",speciesName,"_Prevalence.pdf"))
#par(mfrow=c(3,1))
plot(colMeans(prevT), ylab="Prevalence", xlab="Year", ylim=c(min(prevT),max(prevT)))
plot(colMeans(col),ylab="Colonisation", xlab="Year",ylim=c(min(col),max(col)))
plot(colMeans(ext),ylab="Extinction", xlab="Year",ylim=c(min(ext),max(ext)))
dev.off()

# Simulation results

# Summarise the mean occupancy (prevalence) over time:
mean_prev <- data.frame(year=as.character(2017:2022), mean_prev=colMeans(prevT), sd_prev=apply(prevT,2,sd))
col_prev <- data.frame(year=as.character(2017:2022), mean_prev=colMeans(col), sd_prev=apply(col,2,sd))
ext_prev <- data.frame(year=as.character(2017:2022), mean_prev=colMeans(ext), sd_prev=apply(ext,2,sd))
# Plot mean occupancy (prevalence) and sd:
library(ggplot2)
ggplot(mean_prev, aes(x=year, y=mean_prev, group=1)) + 
  geom_errorbar(aes(ymin=mean_prev-sd_prev, ymax=mean_prev+sd_prev), width=.1) +
  geom_line() +
  geom_point() +
  xlab("Year") + ylab("Predicted mean occupancy")

ggplot(ext_prev, aes(x=year, y=mean_prev, group=1)) + 
  geom_errorbar(aes(ymin=mean_prev-sd_prev, ymax=mean_prev+sd_prev), width=.1) +
  geom_line() +
  geom_point() +
  xlab("Year") + ylab("Predicted mean extinction")

ggplot(col_prev, aes(x=year, y=mean_prev, group=1)) + 
  geom_errorbar(aes(ymin=mean_prev-sd_prev, ymax=mean_prev+sd_prev), width=.1) +
  geom_line() +
  geom_point() +
  xlab("Year") + ylab("Predicted mean colonization")
