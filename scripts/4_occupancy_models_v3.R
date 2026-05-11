###############################################################################
# 4_occupancy_models_v3.R
#
# Same data prep and detection sub-model as v2 (+ stepRep_obs) but with
# gamma and epsilon stripped down to intercept-only. Tests whether the
# 2-4 covariates currently in gamma / epsilon are worth their parameter
# cost given the data's poor identifiability of those sub-models.
###############################################################################

suppressPackageStartupMessages({
  library(here); library(unmarked); library(dplyr); library(tidyr); library(data.table)
})

source(here("R", "model_configs.R"))

YEARS   <- 2017:2023
T_years <- length(YEARS)
J_reps  <- 10
OUT_SUFFIX <- "_v3"
STEPREP_VARIANT <- "stepRep_strict_500m"
species_codes <- c("otitar","ptealc","pteori","tettet")

get_obs_cols <- function(df, prefix, years, reps = 1:10) {
  cols <- c()
  for (yr in years) for (k in reps) {
    cn <- paste0(prefix, ".", k, ".", yr)
    if (cn %in% names(df)) cols <- c(cols, cn)
  }
  cols
}
expand_matrix <- function(mat, J) {
  e <- kronecker(mat, matrix(1, 1, J))
  colnames(e) <- as.vector(sapply(colnames(mat),
                                   function(n) paste0(n, ".", 1:J)))
  e
}

results <- list()

for (sp in species_codes) {
  tryCatch({
    message("\n== ", sp, " ==")
    cfg <- get_model_config(sp)
    # Detection: v2 formula (with stepRep_obs)
    p_form <- update(cfg$p_formula, ~ . + stepRep_obs)
    # gamma and epsilon: INTERCEPT-ONLY
    g_form <- ~ 1
    e_form <- ~ 1

    d <- read.csv(here("data","processed_2023", sp,
                        paste0(sp, "_occ_wide_dynamic.csv")))

    duration  <- as.matrix(d[, get_obs_cols(d, "duration_minutes", YEARS)])
    effort    <- as.matrix(d[, get_obs_cols(d, "effort_distance_km", YEARS)])
    observers <- as.matrix(d[, get_obs_cols(d, "number_observers", YEARS)])
    timev     <- d[, get_obs_cols(d, "time_observations_started", YEARS)]
    detn      <- d[, get_obs_cols(d, "y", YEARS)]
    detn      <- detn %>% mutate(across(where(is.logical), as.integer))
    y.cross   <- as.matrix(detn)
    y.cross[is.na(timev) != is.na(y.cross)] <- NA

    siteCovs <- d[, c("bio1","bio2","tree_cover","grass_cover",
                       "topo_aspect","topo_elev")]

    EVI  <- d[, paste0("EVI_",  YEARS)]
    NDVI <- d[, paste0("NDVI_", YEARS)]
    pr   <- d[, paste0("pr_",   YEARS)]
    tmmn <- d[, paste0("tmmn_", YEARS)]
    tmmx <- d[, paste0("tmmx_", YEARS)]
    step <- d[, paste0(STEPREP_VARIANT, "_", YEARS)]

    siteCovs <- scale(siteCovs)
    duration <- scale(duration); effort <- scale(effort)
    observers<- scale(observers); timev <- scale(timev)
    EVI<-scale(EVI); NDVI<-scale(NDVI); pr<-scale(pr); tmmn<-scale(tmmn); tmmx<-scale(tmmx)
    step <- scale(step)

    NDVI_obs <- expand_matrix(NDVI, J_reps)
    pr_obs   <- expand_matrix(pr,   J_reps)
    topo_aspect_obs <- expand_matrix(siteCovs[,"topo_aspect",drop=FALSE], J_reps*T_years)
    topo_elev_obs   <- expand_matrix(siteCovs[,"topo_elev",  drop=FALSE], J_reps*T_years)
    step_obs <- expand_matrix(step, J_reps)

    n <- nrow(d)
    years_df <- data.frame(matrix(rep(YEARS, each=n), n, T_years))
    years_df <- data.frame(lapply(years_df, as.factor))

    occ_umf <- unmarkedMultFrame(
      y = y.cross,
      siteCovs = data.frame(siteCovs),
      yearlySiteCovs = list(years = years_df, EVI = EVI, NDVI = NDVI,
                             pr = pr, tmmn = tmmn, tmmx = tmmx, stepRep = step),
      obsCovs = list(duration = duration, effort = effort, observers = observers,
                      time = timev, NDVI_obs = NDVI_obs, pr_obs = pr_obs,
                      topo_aspect_obs = topo_aspect_obs,
                      topo_elev_obs = topo_elev_obs,
                      stepRep_obs = step_obs),
      numPrimary = T_years
    )

    cat("  fitting v3: psi", deparse(cfg$psi_formula),
        " | gamma ~ 1 | epsilon ~ 1 | p", deparse(p_form), "\n")
    fit <- colext(psiformula = cfg$psi_formula,
                  gammaformula = g_form,
                  epsilonformula = e_form,
                  pformula = p_form,
                  data = occ_umf)
    cat("  AIC =", round(fit@AIC, 2), "\n")
    cat("  gamma intercept   :", coef(fit, type="col"),
        " SE:", sqrt(diag(vcov(fit, type="col"))), "\n")
    cat("  epsilon intercept :", coef(fit, type="ext"),
        " SE:", sqrt(diag(vcov(fit, type="ext"))), "\n")
    cat("  stepRep_obs coef  :", coef(fit, type="det")["p(stepRep_obs)"], "\n")

    out_path <- here("results", paste0(sp, OUT_SUFFIX, "_model_object.rds"))
    saveRDS(fit, out_path)
    results[[sp]] <- fit
  }, error = function(e) {
    message("  [ERROR] ", sp, ": ", conditionMessage(e))
  })
}

# Summary table v1 / v2 / v3
v1_aics <- sapply(species_codes, function(sp) {
  m <- readRDS(here("/Users/gfandos/Documents/GitHub/dynamic_occupancy_models_steppe_birds/results",
                     paste0(sp, "_model_object.rds")))
  m@AIC
})
v2_aics <- sapply(species_codes, function(sp) {
  m <- readRDS(here("/Users/gfandos/Documents/GitHub/dynamic_occupancy_models_steppe_birds/results/stepRep_v2_run",
                     paste0(sp, "_v2_model_object.rds")))
  m@AIC
})
v3_aics <- sapply(species_codes, function(sp) {
  if (!is.null(results[[sp]])) results[[sp]]@AIC else NA
})

tab <- data.table(species = species_codes,
                  AIC_v1 = round(v1_aics, 2),
                  AIC_v2 = round(v2_aics, 2),
                  AIC_v3 = round(v3_aics, 2))
tab[, d_v2_v1 := AIC_v2 - AIC_v1]
tab[, d_v3_v2 := AIC_v3 - AIC_v2]
tab[, d_v3_v1 := AIC_v3 - AIC_v1]
print(tab)
fwrite(tab, here("results", "v1_v2_v3_aic.csv"))
cat("\nv3 fits saved with _v3 suffix.\n")
