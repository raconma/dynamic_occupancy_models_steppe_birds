###############################################################################
# 4_occupancy_models_v4.R
#
# Detection sub-model as v2 (+ stepRep_obs), gamma reduced to intercept-only
# (colonisation at identifiability floor for these species), and epsilon with
# ONE covariate per species, chosen from the strongest |z| of v1, avoiding
# NDVI for pteori (per docs/decisions_gcb_v4.md to keep climate-vs-land-use
# attribution clean):
#
#   otitar : eps ~ Land_Cover_Type_1_Percent_Class_13   (v1 z = 3.36)
#   ptealc : eps ~ 1                                    (v1: neither cov sig)
#   pteori : eps ~ Land_Cover_Type_1_Percent_Class_12   (v1 z = -3.16; LC12 over NDVI)
#   tettet : eps ~ Land_Cover_Type_1_Percent_Class_12   (v1 only cov)
###############################################################################

suppressPackageStartupMessages({
  library(here); library(unmarked); library(dplyr); library(tidyr); library(data.table)
})
source(here("R", "model_configs.R"))

YEARS   <- 2017:2023
T_years <- length(YEARS)
J_reps  <- 10
OUT_SUFFIX <- "_v4"
STEPREP_VARIANT <- "stepRep_strict_500m"
species_codes <- c("otitar","ptealc","pteori","tettet")

# v4 epsilon formula per species
eps_formula_v4 <- list(
  otitar = ~ Land_Cover_Type_1_Percent_Class_13,
  ptealc = ~ 1,
  pteori = ~ Land_Cover_Type_1_Percent_Class_12,
  tettet = ~ Land_Cover_Type_1_Percent_Class_12
)

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
    p_form <- update(cfg$p_formula, ~ . + stepRep_obs)
    g_form <- ~ 1
    e_form <- eps_formula_v4[[sp]]

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
    LC0  <- d[, paste0("Land_Cover_Type_1_Percent_Class_0_",  YEARS)]
    LC6  <- d[, paste0("Land_Cover_Type_1_Percent_Class_6_",  YEARS)]
    LC7  <- d[, paste0("Land_Cover_Type_1_Percent_Class_7_",  YEARS)]
    LC10 <- d[, paste0("Land_Cover_Type_1_Percent_Class_10_", YEARS)]
    LC12 <- d[, paste0("Land_Cover_Type_1_Percent_Class_12_", YEARS)]
    LC13 <- d[, paste0("Land_Cover_Type_1_Percent_Class_13_", YEARS)]
    LC14 <- d[, paste0("Land_Cover_Type_1_Percent_Class_14_", YEARS)]
    step <- d[, paste0(STEPREP_VARIANT, "_", YEARS)]

    siteCovs <- scale(siteCovs)
    duration<-scale(duration); effort<-scale(effort)
    observers<-scale(observers); timev<-scale(timev)
    EVI<-scale(EVI); NDVI<-scale(NDVI); pr<-scale(pr); tmmn<-scale(tmmn); tmmx<-scale(tmmx)
    LC0<-scale(LC0); LC6<-scale(LC6); LC7<-scale(LC7); LC10<-scale(LC10)
    LC12<-scale(LC12); LC13<-scale(LC13); LC14<-scale(LC14)
    step<-scale(step)

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
                             pr = pr, tmmn = tmmn, tmmx = tmmx,
                             Land_Cover_Type_1_Percent_Class_0  = LC0,
                             Land_Cover_Type_1_Percent_Class_6  = LC6,
                             Land_Cover_Type_1_Percent_Class_7  = LC7,
                             Land_Cover_Type_1_Percent_Class_10 = LC10,
                             Land_Cover_Type_1_Percent_Class_12 = LC12,
                             Land_Cover_Type_1_Percent_Class_13 = LC13,
                             Land_Cover_Type_1_Percent_Class_14 = LC14,
                             stepRep = step),
      obsCovs = list(duration = duration, effort = effort, observers = observers,
                      time = timev, NDVI_obs = NDVI_obs, pr_obs = pr_obs,
                      topo_aspect_obs = topo_aspect_obs,
                      topo_elev_obs = topo_elev_obs,
                      stepRep_obs = step_obs),
      numPrimary = T_years
    )

    cat("  psi    :", deparse(cfg$psi_formula), "\n")
    cat("  gamma  : ~ 1\n")
    cat("  eps    :", deparse(e_form), "\n")
    cat("  det    :", deparse(p_form), "\n")
    fit <- colext(psiformula = cfg$psi_formula,
                  gammaformula = g_form,
                  epsilonformula = e_form,
                  pformula = p_form,
                  data = occ_umf)
    cat("  AIC =", round(fit@AIC, 2), "\n")
    cat("  gamma intercept   :", round(coef(fit, type="col"),3),
        "  SE:", round(sqrt(diag(vcov(fit, type="col"))),3),
        "  gamma_prob =", signif(plogis(coef(fit, type="col"))[1], 3), "\n")
    eps_co <- coef(fit, type="ext")
    eps_se <- sqrt(diag(vcov(fit, type="ext")))
    for (i in seq_along(eps_co))
      cat("  eps", names(eps_co)[i], ":", round(eps_co[i],3),
          "  SE:", round(eps_se[i],3), "\n")
    cat("  eps baseline (intercept-only on logit) -> eps_prob =",
        signif(plogis(eps_co[1]),3), "\n")
    cat("  stepRep_obs coef  :",
        round(coef(fit, type="det")["p(stepRep_obs)"],3), "\n")

    saveRDS(fit, here("results", paste0(sp, OUT_SUFFIX, "_model_object.rds")))
    sink(here("results", paste0(sp, OUT_SUFFIX, "_model_summary.txt")))
    cat("Species:", sp, "  v4: gamma~1, eps with 1 cov, p+stepRep_obs\n")
    cat("Date:", as.character(Sys.time()), "\n\n")
    print(summary(fit))
    sink()
    results[[sp]] <- fit
  }, error = function(e) {
    message("  [ERROR] ", sp, ": ", conditionMessage(e))
  })
}

# Build comparison table v1/v2/v3/v4
REPO <- "/Users/gfandos/Documents/GitHub/dynamic_occupancy_models_steppe_birds"
v1_aics <- sapply(species_codes, function(sp)
  readRDS(file.path(REPO,"results", paste0(sp, "_model_object.rds")))@AIC)
v2_aics <- sapply(species_codes, function(sp)
  readRDS(file.path(REPO,"results","stepRep_v2_run",
                     paste0(sp, "_v2_model_object.rds")))@AIC)
v3_aics <- sapply(species_codes, function(sp) {
  fp <- here("results", paste0(sp, "_v3_model_object.rds"))
  if (file.exists(fp)) readRDS(fp)@AIC else NA
})
v4_aics <- sapply(species_codes, function(sp) {
  if (!is.null(results[[sp]])) results[[sp]]@AIC else NA
})

tab <- data.table(species = species_codes,
                  AIC_v1 = round(v1_aics, 2),
                  AIC_v2 = round(v2_aics, 2),
                  AIC_v3 = round(v3_aics, 2),
                  AIC_v4 = round(v4_aics, 2))
tab[, d_v4_v2 := AIC_v4 - AIC_v2]
tab[, d_v4_v1 := AIC_v4 - AIC_v1]
print(tab)
fwrite(tab, here("results", "v1_v2_v3_v4_aic.csv"))

# epsilon baseline summary (intercept-only on logit -> probability)
eps_summary <- rbindlist(lapply(species_codes, function(sp) {
  if (is.null(results[[sp]])) return(NULL)
  fit <- results[[sp]]
  co <- coef(fit, type="ext")
  se <- sqrt(diag(vcov(fit, type="ext")))
  data.table(species = sp,
             intercept_logit = round(co[[1]], 3),
             intercept_SE    = round(se[[1]], 3),
             baseline_eps    = round(plogis(co[[1]]), 4),
             baseline_eps_lo = round(plogis(co[[1]] - 1.96*se[[1]]), 4),
             baseline_eps_hi = round(plogis(co[[1]] + 1.96*se[[1]]), 4),
             cov_in_eps      = paste(names(co)[-1], collapse=", "))
}))
print(eps_summary)
fwrite(eps_summary, here("results", "v4_eps_baseline.csv"))
cat("\nv4 fits saved.\n")
