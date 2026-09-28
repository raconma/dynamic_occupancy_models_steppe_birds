###############################################################################
# 38_pooled_histories_refit.R
#
# Sensitivity of the primary models to how detection histories were built.
# In the data preparation (scripts/2_prepare_static_variables.R) repeat-visit
# sites were defined as locality x observer x year with 3-10 checklists
# (auk::filter_repeat_visits, Johnston et al. 2021), and then a single site
# was kept per cell (distinct(locality_id) and !duplicated(cells)), so the
# visits of other qualifying sites in the same cell and year were discarded.
#
# Here the histories are rebuilt pooling the visits of ALL qualifying sites
# of each cell and year (random subsample of 10 when there are more, fixed
# seed), for the same cells, with the same site and yearly covariates, and
# the primary model of each species is refitted with the same formulas.
# Observation covariates are built and scaled as in the original fits
# (column-wise scale(); environmental observation covariates unchanged).
#
# Outputs:
#   results/robustness/{sp}_pooled_model_object.rds   (variant "pooled" for
#       colext_model_path(); run scripts 23 and 26 with argument "pooled")
#   results/robustness/pooled_histories_summary.csv   (information retained
#       by the original single-site histories vs the pooled histories)
###############################################################################

library(here)
library(unmarked)
library(terra)
library(auk)

source(here("R", "colext_utils.R"))
set.seed(42)

species <- c("otitar", "tettet", "pteori", "ptealc")
YEARS <- 2017:2023; J <- 10; T_years <- length(YEARS)
grid <- rast(here("data-raw", "data", "environmental_data", "environmental_data_occ",
                  "variables_spain.grd"), lyrs = 1)
dir.create(here("results", "robustness"), showWarnings = FALSE)

summ_rows <- list()
for (sp in species) {
  message("=== ", sp)
  mod <- readRDS(colext_model_path(sp))
  umf <- mod@data
  occ <- read.csv(here("data", "processed_2023", sp, paste0(sp, "_occ_wide_dynamic.csv")))
  stopifnot(nrow(occ) == numSites(umf))
  # CSV cell ids must be the WorldClim cells of the rows' coordinates
  stopifnot(mean(cellFromXY(grid, cbind(occ$longitude, occ$latitude)) == occ$cells) > 0.99)

  # --- Qualifying repeat-visit sites, as in script 2 -----------------------
  zf <- read.csv(here("data", "processed_2023", sp, paste0("ebd_", sp, "_breeding_spain_zf.csv")))
  zf <- zf[zf$year %in% YEARS, ]
  rv <- filter_repeat_visits(zf, min_obs = 3, max_obs = 10, annual_closure = TRUE,
                             date_var = "observation_date",
                             site_vars = c("locality_id", "observer_id"))
  rv$cells <- cellFromXY(grid, cbind(rv$longitude, rv$latitude))
  rv <- rv[rv$cells %in% occ$cells, ]

  # --- Pooled histories: all qualifying visits of a cell and year ----------
  M <- nrow(occ)
  y <- dur <- eff <- obs <- tim <- matrix(NA_real_, M, J * T_years)
  key <- split(seq_len(nrow(rv)), paste(rv$cells, rv$year))
  for (i in seq_len(M)) for (t in seq_len(T_years)) {
    v <- key[[paste(occ$cells[i], YEARS[t])]]
    if (is.null(v)) next
    if (length(v) > J) v <- sort(sample(v, J))
    v <- v[order(rv$observation_date[v])]
    cols <- (t - 1) * J + seq_along(v)
    y[i, cols] <- as.integer(rv$species_observed[v])
    dur[i, cols] <- rv$duration_minutes[v]
    eff[i, cols] <- ifelse(is.na(rv$effort_distance_km[v]), 0, rv$effort_distance_km[v])
    obs[i, cols] <- rv$number_observers[v]
    tim[i, cols] <- rv$time_observations_started[v]
  }
  # As in the original fits: no detection where covariates are missing
  y[is.na(tim) | is.na(dur) | is.na(obs)] <- NA
  sc <- function(m) { s <- scale(m); s[is.na(m)] <- NA; s }

  oc0 <- obsCovs(umf)
  asmat <- function(v) matrix(v, M, J * T_years, byrow = TRUE)
  new_umf <- unmarkedMultFrame(
    y = y, siteCovs = siteCovs(umf),
    yearlySiteCovs = yearlySiteCovs(umf),   # site-major data frame, M x T rows
    obsCovs = list(duration = sc(dur), effort = sc(eff), observers = sc(obs), time = sc(tim),
                   NDVI_obs = asmat(oc0$NDVI_obs), pr_obs = asmat(oc0$pr_obs),
                   topo_aspect_obs = asmat(oc0$topo_aspect_obs),
                   topo_elev_obs = asmat(oc0$topo_elev_obs)),
    numPrimary = T_years)

  fl <- mod@formlist
  fit <- colext(psiformula = fl$psi, gammaformula = fl$col, epsilonformula = fl$ext,
                pformula = fl$det, data = new_umf, starts = coef(mod),
                control = list(maxit = 2000))
  saveRDS(fit, here("results", "robustness", paste0(sp, "_pooled_model_object.rds")))

  # --- Information retained: original vs pooled -----------------------------
  yo <- getY(umf)
  cy <- function(Y) sapply(seq_len(T_years), function(t) {
    b <- Y[, (t - 1) * J + 1:J, drop = FALSE]
    ifelse(rowSums(!is.na(b)) > 0, as.integer(rowSums(b, na.rm = TRUE) > 0), NA)
  })
  no <- cy(yo); np <- cy(y)
  anch <- function(n) { d <- rowSums(n, na.rm = TRUE) > 0
    s <- !is.na(n); sum((s[, 2:6] & s[, 3:7])[d, ]) }   # recent transitions, both years surveyed
  vcov_ok <- all(is.finite(suppressWarnings(diag(vcov(fit))))) &&
             all(suppressWarnings(diag(vcov(fit))) > 0)
  summ_rows[[sp]] <- data.frame(
    species = sp,
    cell_years_surveyed_original = sum(!is.na(no)), cell_years_surveyed_pooled = sum(!is.na(np)),
    detection_cell_years_original = sum(no == 1, na.rm = TRUE),
    detection_cell_years_pooled = sum(np == 1, na.rm = TRUE),
    pooled_detections_missing_in_original = sum(np == 1 & (is.na(no) | no == 0), na.rm = TRUE),
    cells_recorded_original = sum(rowSums(no, na.rm = TRUE) > 0),
    cells_recorded_pooled = sum(rowSums(np, na.rm = TRUE) > 0),
    anchored_transitions_original = anch(no), anchored_transitions_pooled = anch(np),
    aic_pooled = fit@AIC, converged = fit@opt$convergence == 0, finite_se = vcov_ok)
  print(summ_rows[[sp]])
}
out <- do.call(rbind, summ_rows)
write.csv(out, here("results", "robustness", "pooled_histories_summary.csv"), row.names = FALSE)
print(out)

sessionInfo()
