###############################################################################
# 50_weather_anomalies_dynamics.R
#
# Do year-to-year weather and vegetation conditions change colonisation and
# extinction? The covariates of the primary models (script 41) describe where
# rates are higher (static or nearly static across 2017-2023); this script
# adds the within-cell deviation of the breeding-season covariates that vary
# between years. Covariates are those of the dynamic submodels of the
# earlier MSc analysis: TerraClimate minimum and maximum temperature and
# precipitation, and MODIS NDVI. MODIS land cover is not used because its
# within-cell variation between years is negligible (5-8 % of its total SD).
#
# Pre-specified (fixed before fitting, no selection):
#   anomaly     value in year t minus the cell mean over 2017-2023, divided by
#               the SD of the anomalies over all cell-years. Minimum and
#               maximum temperature anomalies are correlated (r = 0.87) and
#               are averaged into one temperature anomaly. Cells with a
#               missing value get an anomaly of 0 (the cell mean).
#   timing      conditions of the breeding season of year t enter the
#               transition from t to t + 1 (unmarked convention, as in the
#               MSc analysis); the alternative timing (year t + 1) is not fitted
#   model A     primary model + temperature, precipitation and NDVI anomalies
#               in colonisation and extinction
#   model B     as A, each anomaly split into its regional-year component
#               (mean anomaly of all cells in year t; identified by the six
#               transitions only and confounded with any other year-level
#               change such as survey effort) and its local component (the
#               deviation of the cell from that year mean)
# Reports coefficients with 95 % Wald intervals, the AIC difference from the
# primary model, convergence and the largest SE of the fit.
#
# Inputs:  results/{sp}_model_object.rds (script 41),
#          data/processed_2023/{sp}/{sp}_occ_wide_dynamic.csv
# Outputs: results/robustness/weather_anomalies_dynamics.csv
###############################################################################

library(here)
library(unmarked)

set.seed(42)

species <- c("otitar", "tettet", "pteori", "ptealc")
YEARS <- 2017:2023; T_years <- length(YEARS)
out_dir <- here("results", "robustness"); dir.create(out_dir, showWarnings = FALSE)

rhs <- function(f) attr(terms(f), "term.labels")
form <- function(x) reformulate(if (length(x)) x else "1")
yearly <- function(d, var) as.matrix(d[, paste0(var, "_", YEARS)])            # M x T
starts_from <- function(cf, sub, labels) {
  unname(sapply(c("(Intercept)", labels), function(l) {
    nm <- paste0(sub, "(", if (l == "(Intercept)") "Int" else l, ")")
    if (nm %in% names(cf)) cf[[nm]] else 0
  }))
}
# within-cell anomaly (M x T), NA set to the cell mean
anomaly <- function(x) { a <- x - rowMeans(x, na.rm = TRUE); a[is.na(a)] <- 0; a }
vec <- function(m) as.vector(t(m))                                            # site-major, as yearlySiteCovs

rows <- list()
for (sp in species) {
  message("=== ", sp)
  mod <- readRDS(here("results", paste0(sp, "_model_object.rds")))
  fl <- mod@formlist; u <- mod@data
  d <- read.csv(here("data", "processed_2023", sp, paste0(sp, "_occ_wide_dynamic.csv")), check.names = FALSE)
  stopifnot(nrow(d) == numSites(u))

  raw <- list(temp = (anomaly(yearly(d, "tmmn")) + anomaly(yearly(d, "tmmx"))) / 2,
              prec = anomaly(yearly(d, "pr")),
              ndvi = anomaly(yearly(d, "NDVI")))
  raw_sd <- sapply(raw, sd)                                                   # SD in original units
  ysc <- yearlySiteCovs(u)
  for (v in names(raw)) {
    a <- raw[[v]] / raw_sd[[v]]
    reg <- matrix(colMeans(a), nrow(a), T_years, byrow = TRUE)
    ysc[[paste0(v, "_anom")]] <- vec(a)
    ysc[[paste0(v, "_reg")]] <- vec(reg)
    ysc[[paste0(v, "_loc")]] <- vec(a - reg)
  }
  umf <- unmarkedMultFrame(y = getY(u), siteCovs = siteCovs(u), yearlySiteCovs = ysc,
                           obsCovs = obsCovs(u), numPrimary = T_years)

  add <- list(A = paste0(names(raw), "_anom"),
              B = c(paste0(names(raw), "_reg"), paste0(names(raw), "_loc")))
  cf0 <- coef(mod)
  for (mdl in names(add)) {
    col <- c(rhs(fl$col), add[[mdl]]); ext <- c(rhs(fl$ext), add[[mdl]])
    st <- c(starts_from(cf0, "psi", rhs(fl$psi)), starts_from(cf0, "col", col),
            starts_from(cf0, "ext", ext), starts_from(cf0, "p", rhs(fl$det)))
    f1 <- tryCatch(
      colext(psiformula = fl$psi, gammaformula = form(col), epsilonformula = form(ext),
             pformula = fl$det, data = umf, starts = st, control = list(maxit = 3000)),
      error = function(e) { message("  ", mdl, " failed: ", conditionMessage(e)); NULL })
    if (is.null(f1)) next
    cf <- coef(f1)
    se <- tryCatch(suppressWarnings(sqrt(diag(vcov(f1)))),
                   error = function(e) setNames(rep(NA_real_, length(cf)), names(cf)))
    for (proc in c("col", "ext")) for (term in add[[mdl]]) {
      nm <- paste0(proc, "(", term, ")")
      rows[[paste(sp, mdl, nm)]] <- data.frame(
        species = sp, model = mdl, process = proc, term = term,
        estimate = cf[[nm]], se = se[[nm]], lo = cf[[nm]] - 1.96 * se[[nm]], hi = cf[[nm]] + 1.96 * se[[nm]],
        delta_aic = f1@AIC - mod@AIC, converged = f1@opt$convergence == 0,
        max_se = max(se, na.rm = TRUE),
        anomaly_sd_raw = raw_sd[[sub("_(anom|reg|loc)$", "", term)]])
    }
    print(do.call(rbind, rows[grepl(paste0("^", sp, " ", mdl, " "), names(rows))]), digits = 3)
  }
}
out <- do.call(rbind, rows); rownames(out) <- NULL
write.csv(out, file.path(out_dir, "weather_anomalies_dynamics.csv"), row.names = FALSE)
print(out, digits = 3)

sessionInfo()
