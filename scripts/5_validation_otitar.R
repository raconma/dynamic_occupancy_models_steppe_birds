library(sf)
library(terra)
library(dplyr)
library(ggplot2)
library(caret)
library(pROC)
library(Metrics)
library(spdep)
library(blockCV)

# Spanish biodiversity atlas data for validation
atlas <- st_read("../data/validation/atlas_biodiversidad/atlas_biodiversidad/aves_spain.shp")

# Our model predictions
pred_occ <- read.csv("../data/otitar/occ_otitar_prediction.csv")
names(pred_occ) <- tolower(names(pred_occ))
stopifnot(all(c("longitude","latitude","occ_prob") %in% names(pred_occ)))

################################################################################
##################            PREDICTIONS TO RASTER           ##################
################################################################################

# Create spatial raster from predictions
r_pred <- terra::rast(
  data.frame(
    x = pred_occ$longitude,
    y = pred_occ$latitude,
    occ_prob = pred_occ$occ_prob
  ),
  type = "xyz"
)

# Assign CRS manually
crs(r_pred) <- "+proj=longlat"

# Reproject continuous probability raster
r_pred_proj <- project(r_pred, crs(atlas))

# Extract median probability per atlas polygon
ext_prob <- terra::extract(r_pred_proj, atlas, fun = median, na.rm = TRUE)
atlas$pred_prob_mean <- ext_prob$occ_prob

# Data for validation
atlas_valid <- atlas %>% filter(!is.na(pred_prob_mean))

################################################################################
##################      SPATIAL CROSS-VALIDATION (BLOCKCV)    ##################
################################################################################

set.seed(123)

sf_cv <- atlas_valid

sb <- spatialBlock(
  speciesData = sf_cv,
  species = "OTITAR",
  theRange = 50000,
  k = 5,
  selection = "random",
  iteration = 100,
  showBlocks = FALSE
)

sf_cv$fold <- sb$foldID

results_cv <- list()

for (i in 1:max(sf_cv$fold)) {
  
  train <- sf_cv[sf_cv$fold != i, ]
  test  <- sf_cv[sf_cv$fold == i, ]
  
  # Threshold from training ROC
  roc_train <- roc(train$OTITAR, train$pred_prob_mean, quiet = TRUE)
  th <- coords(roc_train, "best", ret = "threshold") %>% as.numeric()
  
  # Binary map for test fold
  test$bin <- ifelse(test$pred_prob_mean >= th, 1, 0)
  
  # Metrics
  auc_i <- pROC::auc(pROC::roc(test$OTITAR, test$pred_prob_mean, quiet = TRUE))
  rmse_i <- rmse(test$OTITAR, test$pred_prob_mean)
  
  cm <- table(test$OTITAR, test$bin)
  TP <- cm["1","1"]; TN <- cm["0","0"]
  FP <- cm["0","1"]; FN <- cm["1","0"]
  
  sens_i <- TP/(TP+FN)
  spec_i <- TN/(TN+FP)
  TSS_i  <- sens_i + spec_i - 1
  
  results_cv[[i]] <- data.frame(
    fold=i, threshold=th, AUC=auc_i, RMSE=rmse_i,
    Sens=sens_i, Spec=spec_i, TSS=TSS_i
  )
}

df_cv <- do.call(rbind, results_cv)
print(df_cv)
summary(df_cv)

################################################################################
##################        OPTIMAL THRESHOLD                   ##################
################################################################################

opt_th <- mean(df_cv$threshold)
cat("Spatially robust threshold =", opt_th, "\n")

################################################################################
##################        APPLY THRESHOLD TO FINAL MAP        ##################
################################################################################

# Binary raster
r_pred_bin <- r_pred >= opt_th
names(r_pred_bin) <- "occ_bin"

# Reproject
r_pred_bin_proj <- project(r_pred_bin, crs(atlas), method = "bilinear")

# Extract binary presence/absence
ext_bin <- terra::extract(r_pred_bin_proj, atlas, fun = max, na.rm = TRUE)

atlas$pred_bin_final <- ifelse(
  is.na(atlas$pred_prob_mean), NA,
  ifelse(ext_bin$occ_bin > 0, 1, 0)
)

################################################################################
##################            CONTINUOUS VALIDATION            #################
################################################################################

atlas_deciles <- atlas_valid %>%
  mutate(prob_bin = ntile(pred_prob_mean, 10)) %>%
  group_by(prob_bin) %>%
  summarise(mean_prob = mean(pred_prob_mean),
            obs_prev  = mean(as.numeric(OTITAR)),
            n = n())

spearman_corr <- cor.test(atlas_deciles$mean_prob,
                          atlas_deciles$obs_prev,
                          method="spearman")

r2 <- summary(lm(obs_prev ~ mean_prob, data=atlas_deciles))$r.squared

cat("Spearman rho:", round(spearman_corr$estimate,3),
    "| p =", signif(spearman_corr$p.value,3), "\n")
cat("R^2:", round(r2,3), "\n")

ggplot(atlas_deciles, aes(x=mean_prob, y=obs_prev)) +
  geom_point(aes(size=n)) +
  geom_abline(slope=1, intercept=0, linetype="dashed") +
  geom_smooth(method="lm", se=TRUE, color="blue") +
  theme_minimal() +
  labs(x="Mean predicted probability",
       y="Observed prevalence",
       title="Continuous validation (calibration plot)")

################################################################################
##################           SPATIAL VALIDATION (RESIDUALS)    #################
################################################################################

atlas_resid <- atlas_valid %>%
  mutate(resid = pred_prob_mean - as.numeric(OTITAR))

ggplot(atlas_resid) +
  geom_sf(aes(fill = resid), color=NA) +
  scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0) +
  ggtitle("Residuals (Predicted - Observed)") +
  theme_minimal()

# Moran's I
sf_valid <- atlas_resid %>% filter(!is.na(resid))
sp_valid <- sf::as_Spatial(sf_valid)

nb <- poly2nb(sp_valid, queen=TRUE)
lw <- nb2listw(nb, style="W", zero.policy=TRUE)

moran_res <- moran.test(sf_valid$resid, lw, zero.policy=TRUE)
print(moran_res)

################################################################################
##################                VISUALIZATION               ##################
################################################################################

map_pred <- ggplot(atlas) +
  geom_sf(aes(fill = factor(pred_bin_final))) +
  scale_fill_manual(values=c("white","darkgreen"),
                    name="Prediction",
                    labels=c("Absence (0)", "Presence (1)")) +
  ggtitle("Model prediction (Optimal threshold)") +
  theme_bw()

map_obs <- ggplot(atlas) +
  geom_sf(aes(fill = factor(OTITAR))) +
  scale_fill_manual(values=c("grey80","red"),
                    name="Atlas OTITAR",
                    labels=c("Absence (0)", "Presence (1)")) +
  ggtitle("Biodiversity atlas") +
  theme_bw()

print(map_pred)
print(map_obs)
