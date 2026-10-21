###############################################################################
# 39_pooled_histories_primary.R
#
# Makes the pooled detection histories (script 38) the basis of the analysis.
#
# In scripts/2_prepare_static_variables.R repeat-visit sites were defined as
# locality x observer x year with 3-10 checklists (auk::filter_repeat_visits,
# Johnston et al. 2021), and a single site was then kept per cell
# (distinct(locality_id), !duplicated(cells)), discarding the visits of the
# other qualifying sites in the same cell and year. Here, for every cell and
# year, the visits of ALL qualifying sites are pooled, up to 30 per cell-year.
# A random subsample of 10 visits (script 38) changed the anchored dynamics
# between subsamples: a quarter of the cell-years of recorded cells have more
# than 10 visits, with p ~ 0.2 per checklist ten visits still miss ~11 % of
# presences, and anchored events are few. Cell-years with more than 30 visits
# (heavily visited sites, up to 430) are thinned deterministically to 30
# visits evenly spaced through the season (no randomness); with 30 visits the
# probability of missing a present species is ~0.1 %, so the cap has almost
# no bearing on the inferred state. Then
#   1. the model input tables are rewritten with the pooled visits (detection
#      and checklist-effort columns, J per year; all other columns unchanged):
#        data/processed_2023/{sp}/{sp}_occ_wide_dynamic.csv
#        reports/script_4_v2_inputs/{sp}_occ_wide_dynamic.csv (stepRep inputs)
#      The single-site originals are kept as *_singlesite.csv.
#      The tables also carry the day of year of each visit (day_of_year) and
#      the proportion of pseudo-steppe within 500 m of it (steprep_visit).
#   2. the primary model and the stepRep v2 model are refitted on the pooled
#      histories with the original formulas, site and yearly covariates, and
#      environmental observation covariates (constant within a year, expanded
#      to J visits); checklist covariates are scaled over all visits. These
#      refits are kept in results/_pooled_original_formulas/ for comparison;
#      the models used downstream come from the covariate selection on the
#      pooled histories (script 41). Single-site originals are kept in
#      results/_single_site_histories/.
#
# Outputs: the files above; results/pooled_refit_summary.csv
###############################################################################

library(here)
library(unmarked)
library(terra)
library(auk)

set.seed(42)

species <- c("otitar", "tettet", "pteori", "ptealc")
YEARS <- 2017:2023; T_years <- length(YEARS)
grid <- rast(here("data-raw", "data", "environmental_data", "environmental_data_occ",
                  "variables_spain.grd"), lyrs = 1)
foc500 <- rast(here("data", "derived", "corine", "focal_strict_500m.tif"))
bk_dir <- here("results", "_single_site_histories")
dir.create(bk_dir, showWarnings = FALSE)
# refits with the original formulas (superseded by the covariate selection
# of script 41, which writes the primary and stepRep models)
orig_dir <- here("results", "_pooled_original_formulas")
dir.create(orig_dir, showWarnings = FALSE)

# Keep the single-site original once; later runs read from the backup
backup_once <- function(path, bk) {
  if (!file.exists(bk)) stopifnot(file.copy(path, bk))
  bk
}
ck <- c(y = "y", duration = "duration_minutes", effort = "effort_distance_km",
        observers = "number_observers", time = "time_observations_started",
        doy = "day_of_year", steprep = "steprep_visit")
col_of <- function(var, j, yr) paste0(var, ".", j, ".", yr)

rows <- list()
for (sp in species) {
  message("=== ", sp)
  # --- Qualifying repeat-visit sites and one pooled selection per cell-year --
  zf <- read.csv(here("data", "processed_2023", sp, paste0("ebd_", sp, "_breeding_spain_zf.csv")))
  zf <- zf[zf$year %in% YEARS, ]
  rv <- filter_repeat_visits(zf, min_obs = 3, max_obs = 10, annual_closure = TRUE,
                             date_var = "observation_date",
                             site_vars = c("locality_id", "observer_id"))
  rv$cells <- cellFromXY(grid, cbind(rv$longitude, rv$latitude))
  # within-cell steppe representativeness of each checklist: proportion of
  # pseudo-steppe (CORINE 2018 codes 211, 231, 321, 333) within 500 m of its
  # location (raster from build_stepRep.R; outside CORINE coverage = 0)
  pts <- project(vect(cbind(rv$longitude, rv$latitude), crs = "EPSG:4326"), crs(foc500))
  rv$steprep_visit <- extract(foc500, pts, ID = FALSE)[[1]]
  rv$steprep_visit[is.na(rv$steprep_visit)] <- 0
  rv <- rv[order(rv$cells, rv$year, rv$observation_date, rv$checklist_id), ]
  groups <- split(seq_len(nrow(rv)), paste(rv$cells, rv$year))
  groups <- groups[order(names(groups))]
  J_CAP <- 30
  J <- min(max(lengths(groups)), J_CAP)
  message("  visits per cell-year: max ", max(lengths(groups)), "; cell-years thinned to ",
          J, ": ", sum(lengths(groups) > J))
  vis <- function(cell, yr) {
    v <- groups[[paste(cell, yr)]]
    if (is.null(v)) return(NULL)
    v <- v[order(rv$observation_date[v], rv$time_observations_started[v])]
    if (length(v) > J) v <- v[unique(round(seq(1, length(v), length.out = J)))]
    rv[v, ]
  }

  # --- 1. Rewrite the model input tables ------------------------------------
  csvs <- c(primary = here("data", "processed_2023", sp, paste0(sp, "_occ_wide_dynamic.csv")),
            steprep = here("reports", "script_4_v2_inputs", paste0(sp, "_occ_wide_dynamic.csv")))
  pooled_tabs <- list()
  for (k in names(csvs)) {
    src <- backup_once(csvs[k], sub("\\.csv$", "_singlesite.csv", csvs[k]))
    d <- read.csv(src, check.names = FALSE)
    # replace the visit columns by J per year, ordered year-major as in the
    # original tables (y.1.2017 ... y.J.2017, y.1.2018, ...)
    old <- grepl(paste0("^(", paste(ck, collapse = "|"), ")\\.[0-9]+\\.[0-9]{4}$"), names(d))
    d <- d[, !old, drop = FALSE]
    mats <- lapply(ck, function(v) matrix(NA, nrow(d), J * T_years))
    for (i in seq_len(nrow(d))) for (t in seq_len(T_years)) {
      x <- vis(d$cells[i], YEARS[t])
      if (is.null(x)) next
      cols <- (t - 1) * J + seq_len(nrow(x))
      mats$y[i, cols] <- x$species_observed
      mats$duration[i, cols] <- x$duration_minutes
      mats$effort[i, cols] <- ifelse(is.na(x$effort_distance_km), 0, x$effort_distance_km)
      mats$observers[i, cols] <- x$number_observers
      mats$time[i, cols] <- x$time_observations_started
      mats$doy[i, cols] <- x$day_of_year
      mats$steprep[i, cols] <- x$steprep_visit
    }
    # as in the original tables: no detection record without effort data
    mats$y[is.na(mats$time) | is.na(mats$duration) | is.na(mats$observers)] <- NA
    for (v in names(ck)) {
      cn <- unlist(lapply(YEARS, function(yr) col_of(ck[[v]], seq_len(J), yr)))
      d[cn] <- as.data.frame(mats[[v]])
    }
    pooled_tabs[[k]] <- d   # written after the refit, with its final y
  }

  # --- 2. Refit the primary and stepRep v2 models on the pooled histories ---
  models <- c(primary = here("results", paste0(sp, "_model_object.rds")),
              steprep = here("results", "stepRep_v2_run", paste0(sp, "_v2_model_object.rds")))
  for (k in names(models)) {
    src_mod <- readRDS(backup_once(models[k], file.path(bk_dir, basename(models[k]))))
    umf <- src_mod@data
    d <- pooled_tabs[[k]]
    M <- numSites(umf); stopifnot(nrow(d) == M)
    wide <- function(var) as.matrix(d[, unlist(lapply(YEARS, function(yr) col_of(var, seq_len(J), yr)))])
    y <- wide("y"); storage.mode(y) <- "numeric"
    dur <- wide("duration_minutes"); eff <- wide("effort_distance_km")
    obs <- wide("number_observers"); tim <- wide("time_observations_started")
    y[is.na(tim) | is.na(dur) | is.na(obs)] <- NA   # as in the original fits
    # one mean and SD per covariate over all visits (column-wise scale()
    # of the original fits gives NaN for sparse late-visit columns with
    # constant values once J is large)
    sc <- function(m) (m - mean(m, na.rm = TRUE)) / sd(m, na.rm = TRUE)
    oc0 <- obsCovs(umf)
    keep_oc <- setdiff(names(oc0), c("duration", "effort", "observers", "time"))
    # environmental observation covariates are constant within a year:
    # take the first visit of each year in the source frame and expand to J
    J0 <- ncol(getY(umf)) / T_years
    expand_J <- function(v) {
      m <- matrix(v, M, J0 * T_years, byrow = TRUE)[, (seq_len(T_years) - 1) * J0 + 1, drop = FALSE]
      m[, rep(seq_len(T_years), each = J), drop = FALSE]
    }
    oc <- c(list(duration = sc(dur), effort = sc(eff), observers = sc(obs), time = sc(tim)),
            lapply(oc0[keep_oc], expand_J))
    # visits with a missing detection covariate carry no detection record
    # (unmarked would drop them when fitting; the event code must agree)
    # The table keeps the effort-based mask only; script 41 masks the
    # visits with missing values of its own detection covariates
    write.csv(d, csvs[k], row.names = FALSE)
    dv <- intersect(all.vars(src_mod@formlist$det), names(oc))
    y[Reduce(`|`, lapply(oc[dv], is.na))] <- NA
    new_umf <- unmarkedMultFrame(y = y, siteCovs = siteCovs(umf),
                                 yearlySiteCovs = yearlySiteCovs(umf), obsCovs = oc,
                                 numPrimary = T_years)
    fl <- src_mod@formlist
    t0 <- Sys.time()
    fit <- colext(psiformula = fl$psi, gammaformula = fl$col, epsilonformula = fl$ext,
                  pformula = fl$det, data = new_umf, starts = coef(src_mod),
                  control = list(maxit = 2000))
    saveRDS(fit, file.path(orig_dir, basename(models[k])))
    se <- suppressWarnings(sqrt(diag(vcov(fit))))
    rows[[paste(sp, k)]] <- data.frame(
      species = sp, model = k, n_sites = M, visits_per_year_max = J,
      sites_removed = length(fit@sitesRemoved),
      detections_single_site = sum(getY(umf) == 1, na.rm = TRUE),
      detections_pooled = sum(y == 1, na.rm = TRUE),
      aic = fit@AIC, converged = fit@opt$convergence == 0, finite_se = all(is.finite(se)),
      fit_seconds = as.numeric(Sys.time() - t0, units = "secs"))
    print(rows[[paste(sp, k)]])
  }
}
out <- do.call(rbind, rows); rownames(out) <- NULL
write.csv(out, here("results", "pooled_refit_summary.csv"), row.names = FALSE)
print(out)

sessionInfo()
