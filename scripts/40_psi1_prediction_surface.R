###############################################################################
# 40_psi1_prediction_surface.R
#
# Initial-occupancy (psi1) prediction over the WorldClim 2.5' grid of mainland
# Spain from the current primary models, for the atlas validation of
# scripts/5_validation.R (Table S10). Replaces the prediction step of
# scripts/4_occupancy_models.R, whose output predates the current fits.
#
# Covariate transformation: the fitted models use scale() of the input-table
# columns, which script 2 had already standardised with scaling_params, so a
# raw raster value r enters the model as (r - m) / s with
#   m = center + scale * mean(table column),  s = scale * sd(table column).
# The surface is checked against psi1 of the fitted model at its own sites.
#
# Inputs:  results/{sp}_model_object.rds, data/processed_2023/{sp}/
#          ({sp}_occ_wide_dynamic.csv, {sp}_scaling_params.rds), environmental
#          and topography rasters (as in script 4)
# Outputs: data/processed_2023/{sp}/occ_{sp}_prediction.csv (previous file
#          kept as occ_{sp}_prediction_previous.csv)
###############################################################################

library(here)
library(unmarked)
library(terra)

source(here("R", "colext_utils.R"))
set.seed(42)  # no stochastic steps

species <- c("otitar", "tettet", "pteori", "ptealc")
vars <- rast(here("data-raw", "data", "environmental_data", "environmental_data_occ",
                  "variables_spain.grd"))
topo <- c(resample(rast(here("data-raw", "data", "topology_data", "topo_aspect_masked.tif")),
                   vars, method = "bilinear"),
          resample(rast(here("data-raw", "data", "topology_data", "topo_elev_masked.tif")),
                   vars, method = "bilinear"))
names(topo) <- c("topo_aspect", "topo_elev")
vars <- c(vars, topo)

chk <- list()
for (sp in species) {
  mod <- readRDS(colext_model_path(sp))
  psi_vars <- all.vars(mod@formlist$psi)
  tab <- read.csv(here("data", "processed_2023", sp, paste0(sp, "_occ_wide_dynamic.csv")))
  sp_par <- readRDS(here("data", "processed_2023", sp, paste0(sp, "_scaling_params.rds")))
  surf <- as.data.frame(vars[[psi_vars]], xy = TRUE, na.rm = TRUE)
  nd <- surf[, psi_vars, drop = FALSE]
  for (v in psi_vars) {
    m <- sp_par[[v]]$center + sp_par[[v]]$scale * mean(tab[[v]], na.rm = TRUE)
    s <- sp_par[[v]]$scale * sd(tab[[v]], na.rm = TRUE)
    nd[[v]] <- (nd[[v]] - m) / s
  }
  pr <- predict(mod, type = "psi", newdata = nd)
  out <- data.frame(longitude = surf$x, latitude = surf$y, occ_prob = pr$Predicted, occ_se = pr$SE)
  # Check: surface value at the model's cells against the fitted psi1
  D <- colext_design(mod)
  cells <- tab$cells[D$keep]
  sv <- out$occ_prob[match(cells, cellFromXY(vars, cbind(out$longitude, out$latitude)))]
  chk[[sp]] <- data.frame(species = sp, n = sum(!is.na(sv)),
                          cor_logit = cor(qlogis(sv), qlogis(colext_rates(coef(mod), D)$psi1),
                                          use = "complete.obs"))
  f <- here("data", "processed_2023", sp, paste0("occ_", sp, "_prediction.csv"))
  bk <- sub("\\.csv$", "_previous.csv", f)
  if (file.exists(f) && !file.exists(bk)) file.copy(f, bk)
  write.csv(out, f, row.names = FALSE)
}
print(do.call(rbind, chk))

sessionInfo()
