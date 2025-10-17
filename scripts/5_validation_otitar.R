library(sf)
library(terra)
library(dplyr)
library(ggplot2)
library(caret)
library(pROC)
library(Metrics)
library(spdep)

# Spanish biodiversity atlas data for validation
atlas <- st_read("../data/validation/atlas_biodiversidad/atlas_biodiversidad/aves_spain.shp")
names(atlas)
head(atlas)
table(atlas$OTITAR, useNA="ifany")

# Our model predictions
pred_occ <- read.csv("../data/otitar/occ_otitar_prediction.csv")

################################################################################
##################            PREDICTIONS TO RASTER           ##################
################################################################################

# Check prediction data
names(pred_occ) <- tolower(names(pred_occ))
stopifnot(all(c("longitude", "latitude", "occ_prob") %in% names(pred_occ)))

# Create spatial raster from predictions
r_pred <- terra::rast(
  data.frame(x = pred_occ$longitude,
             y = pred_occ$latitude,
             occ_prob = pred_occ$occ_prob),
  type = "xyz"
)

# Asign CRS manually to avoid problems
crs(r_pred) <- rast(crs = "+proj=longlat")

plot(r_pred, main = "Occupancy probability map (OTITAR)")

# Reproject continuous probability raster
r_pred_proj <- project(r_pred, crs(atlas))

# Extract median probability per atlas polygon
ext_prob <- terra::extract(r_pred_proj, atlas, fun = median, na.rm = TRUE)
atlas$pred_prob_mean <- ext_prob$occ_prob

################################################################################
##################       OPTIMAL THRESHOLD -> BINARY MAP      ##################
################################################################################

# Create a common subset for both maps
atlas_valid <- atlas %>% filter(!is.na(pred_prob_mean))
cat("Common subset cells:", nrow(atlas_valid), "\n")

# Function to generate threshold based on sensitivity and specificity (70/30 and 80/20)
run_split_validation <- function(data, split = 0.7, seed = 123) {
  set.seed(seed)
  n <- nrow(data)
  idx_train <- sample(seq_len(n), size = round(split * n))
  train <- data[idx_train, ]
  test  <- data[-idx_train, ]
  
  roc_train <- pROC::roc(train$OTITAR, train$pred_prob_mean, quiet = TRUE)
  coords_best <- pROC::coords(roc_train, "best", best.method = "youden",
                              ret = c("threshold","sensitivity","specificity"),
                              transpose = FALSE)
  thresh <- as.numeric(coords_best["threshold"])
  
  test$pred_bin_opt <- ifelse(test$pred_prob_mean >= thresh, 1, 0)
  cm_test <- caret::confusionMatrix(as.factor(test$pred_bin_opt),
                                    as.factor(test$OTITAR),
                                    positive="1")
  auc_test <- pROC::auc(pROC::roc(test$OTITAR, test$pred_prob_mean))
  rmse_test <- Metrics::rmse(test$OTITAR, test$pred_prob_mean)
  
  cm <- table(test$OTITAR, test$pred_bin_opt)
  TP <- cm["1","1"]; TN <- cm["0","0"]; FP <- cm["0","1"]; FN <- cm["1","0"]
  sens <- TP/(TP+FN); spec <- TN/(TN+FP)
  TSS <- sens + spec - 1
  
  list(split = split, threshold = thresh, auc = auc_test, rmse = rmse_test,
       sens = sens, spec = spec, TSS = TSS, confusion = cm_test,
       train = train, test = test)  # >>> CAMBIO: añadimos train y test
}

# Compare 70/30 vs 80/20 splits
res70 <- run_split_validation(atlas_valid, 0.7)
res80 <- run_split_validation(atlas_valid, 0.8)

print(res70$confusion)
cat(paste("Threshold:", round(res70$threshold,3),
          "| AUC:", round(res70$auc,3),
          "| RMSE:", round(res70$rmse,3),
          "| TSS:", round(res70$TSS,3), "\n"))

print(res80$confusion)
cat(paste("Threshold:", round(res80$threshold,3),
          "| AUC:", round(res80$auc,3),
          "| RMSE:", round(res80$rmse,3),
          "| TSS:", round(res80$TSS,3), "\n"))

################################################################################
##################        APPLY OPTIMAL THRESHOLD TO MAP      ##################
################################################################################

opt_th <- res80$threshold
cat("Optimal threshold =", round(opt_th, 4), "\n")

# Create binary raster
r_pred_bin_opt <- r_pred >= opt_th
names(r_pred_bin_opt) <- "occ_bin"

# Plot binario óptimo
plot(r_pred_bin_opt, main = paste("Binary map (optimal threshold =", round(opt_th, 3), ")"))

# Reproject predictions to match atlas CRS
r_pred_bin_proj_opt <- project(r_pred_bin_opt, crs(atlas), method = "bilinear")
ext_opt <- terra::extract(r_pred_bin_proj_opt, atlas, fun = max, na.rm = TRUE)

atlas$pred_bin_opt <- ifelse(is.na(atlas$pred_prob_mean), NA,
                             ifelse(ext_opt$occ_bin > 0, 1, 0))

################################################################################
##################          VALIDATION WITH TEST              ##################
################################################################################

# Just test set from 80/20 split
atlas_valid_opt <- res80$test
cat("Test cells (20%):", nrow(atlas_valid_opt), "\n")

# Confusion matrix
cm_opt <- table(Observed = atlas_valid_opt$OTITAR,
                Predicted = atlas_valid_opt$pred_bin_opt)
print("Confusion Matrix:")
print(cm_opt)

cm_opt_caret <- caret::confusionMatrix(as.factor(atlas_valid_opt$pred_bin_opt),
                                       as.factor(atlas_valid_opt$OTITAR),
                                       positive = "1")
print(cm_opt_caret)

# AUC check
roc_obj_opt <- pROC::roc(atlas_valid_opt$OTITAR, atlas_valid_opt$pred_prob_mean)
auc_value_opt <- pROC::auc(roc_obj_opt)
cat(paste("AUC =", round(auc_value_opt, 3), "\n"))

# RMSE check
rmse_value_opt <- Metrics::rmse(atlas_valid_opt$OTITAR, atlas_valid_opt$pred_prob_mean)
cat(paste("RMSE =", round(rmse_value_opt, 3), "\n"))

# TSS check
cm_opt <- table(atlas_valid_opt$OTITAR, atlas_valid_opt$pred_bin_opt)
TP <- cm_opt[2,2]; TN <- cm_opt[1,1]; FP <- cm_opt[1,2]; FN <- cm_opt[2,1]  # >>> CORREGIDO
sensitivity_opt <- TP / (TP + FN)
specificity_opt <- TN / (TN + FP)
TSS_opt <- sensitivity_opt + specificity_opt - 1
cat(paste("TSS =", round(TSS_opt, 3), "\n"))

# ROC curve check
plot(roc_obj_opt, col="darkred", lwd=2,
     main=paste("ROC curve (AUC =", round(auc_value_opt,3), ")"))
abline(a=0, b=1, lty=2)

cat(paste("Cells:", nrow(atlas_valid_opt),
          "\nAUC:", round(auc_value_opt,3),
          "\nRMSE:", round(rmse_value_opt,3),
          "\nTSS:", round(TSS_opt,3),
          "\nBalanced Accuracy:", round(cm_opt_caret$byClass["Balanced Accuracy"],3),
          "\nKappa:", round(cm_opt_caret$overall["Kappa"],3), "\n"))

################################################################################
##################            CONTINUOUS VALIDATION            #################
################################################################################

atlas_deciles <- atlas_valid %>%
  mutate(prob_bin = ntile(pred_prob_mean, 10)) %>%
  group_by(prob_bin) %>%
  summarise(mean_prob = mean(pred_prob_mean),
            obs_prev = mean(as.numeric(OTITAR)),
            n = n())

spearman_corr <- cor.test(atlas_deciles$mean_prob, atlas_deciles$obs_prev, method="spearman")
r2 <- summary(lm(obs_prev ~ mean_prob, data=atlas_deciles))$r.squared

cat("Spearman rho:", round(spearman_corr$estimate,3), "(p =", signif(spearman_corr$p.value,3), ")\n")
cat("R²:", round(r2,3), "\n")

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

atlas_resid <- atlas_valid %>% mutate(resid = pred_prob_mean - as.numeric(OTITAR))

ggplot(atlas_resid) +
  geom_sf(aes(fill = resid), color=NA) +
  scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0) +
  ggtitle("Residuals (Predicted - Observed)") +
  theme_minimal()

sf_valid <- atlas_resid %>% filter(!is.na(resid))
sp_valid <- sf::as_Spatial(sf_valid)
nb <- spdep::poly2nb(sp_valid, queen=TRUE)
lw <- spdep::nb2listw(nb, style="W", zero.policy=TRUE)
moran_res <- spdep::moran.test(sf_valid$resid, lw, zero.policy=TRUE)
print(moran_res)

################################################################################
##################                VISUALIZATION               ##################
################################################################################

map_pred <- ggplot(atlas) +
  geom_sf(aes(fill = factor(pred_bin_opt))) +
  scale_fill_manual(values=c("white", "darkgreen"),
                    name="Prediction",
                    labels=c("Absence (0)", "Presence (1)")) +
  ggtitle("Model prediction (Optimal threshold)") +
  theme_bw()

map_obs <- ggplot(atlas) +
  geom_sf(aes(fill = factor(OTITAR))) +
  scale_fill_manual(values=c("grey80", "red"),
                    name="Atlas OTITAR",
                    labels=c("Absence (0)", "Presence (1)")) +
  ggtitle("Biodiversity atlas") +
  theme_bw()

print(map_pred)
print(map_obs)
