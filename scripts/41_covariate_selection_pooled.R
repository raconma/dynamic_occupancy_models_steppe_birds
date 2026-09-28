###############################################################################
# 41_covariate_selection_pooled.R
#
# Covariate selection on the pooled detection histories (script 39) and fit
# of the primary models used downstream.
#
# The covariate structures of the earlier fits were chosen on single-site
# histories and the selection was not recorded in code. Here it is redone on
# the pooled histories with rules fixed before fitting:
#
#   detection    see below; colonisation and extinction at all
#       initial-occupancy candidates while detection is selected.
#   detection    forward selection by AIC among survey duration, distance
#       travelled, number of observers, start time (+ quadratic), day of year
#       (+ quadratic), NDVI and within-cell steppe representativeness of the
#       checklist (proportion of pseudo-steppe within 500 m; script 39); a
#       quadratic term can enter only after its linear term. Initial
#       occupancy, colonisation and extinction at all their candidates.
#   initial occupancy  forward selection by AIC among bio1, bio2, tree
#       cover, herbaceous cover and aspect (elevation is excluded: |r| > 0.8
#       with bio1), with the selected detection submodel.
#   colonisation, then extinction  forward selection by AIC, at most three
#       slopes each, among the initial-occupancy candidates and three land-use
#       covariates: pseudo-steppe fraction of the cell (CORINE 2018 codes 211,
#       231, 321, 333; script 23) and MODIS MCD12C1 cropland (class 12) and
#       open shrubland (class 7) percentages (yearly; broad pseudo-steppe and
#       MODIS grassland excluded, |r| > 0.7 with pseudo-steppe and cropland).
#       Colonisation is selected with extinction on the initial-occupancy
#       covariates, and extinction with the selected colonisation submodel.
#       (Constant colonisation and extinction were rejected by AIC by 18 to
#       95 units: with gamma estimated over the whole grid, ~0.002, the
#       constant prior explained most apparent colonisations in recorded
#       cells as missed detections.)
#   A term is added when it lowers AIC by at least 2. A candidate fit is
#   rejected when the optimiser does not converge or a standard error of the
#   submodel being selected is not finite or exceeds 3 (logit scale,
#   standardised covariates), which is how separation shows. (A first rule on
#   slope magnitude, |slope| > 5, was dropped before any selection ran: it
#   rejected the starting model of O. tarda for its tree-cover effect, which
#   is well estimated.) After selection, while any colonisation or
#   extinction parameter, intercept included, has a standard error above 3,
#   the slope with the largest absolute estimate in that submodel is removed
#   and the model refitted. (The first version of this rule looked at slopes
#   only; it was extended after the P. alchata fit showed a colonisation
#   intercept of -18 with SE 3.2 alongside slopes of 7.6 and -8.3.)
#
# All candidate fits use the same detection record: visits with a missing
# value of any candidate detection covariate carry no detection. Checklist,
# day-of-year and yearly covariates are standardised over all cell-years;
# site covariates over cells.
#
# Inputs:  data/processed_2023/{sp}/{sp}_occ_wide_dynamic.csv (script 39)
# Outputs: results/{sp}_model_object.rds                      (primary)
#          the input table above, with y set to NA where the model's
#          detection covariates are missing (rows as in the model)
#          results/model_selection/selection_log.csv, selected_models.csv
#
# Usage:   Rscript scripts/41_covariate_selection_pooled.R [species ...]
###############################################################################

library(here)
library(unmarked)

set.seed(42)

species <- c("otitar", "tettet", "pteori", "ptealc")
args <- commandArgs(trailingOnly = TRUE)
if (length(args)) species <- args

YEARS <- 2017:2023; T_years <- length(YEARS)
DELTA_AIC <- 2
MAX_SE <- 3
DET_CAND <- c("duration", "effort", "observers", "time", "I(time^2)",
              "doy", "I(doy^2)", "NDVI_obs", "steprep")
PSI_CAND <- c("bio1", "bio2", "tree_cover", "grass_cover", "topo_aspect")
DYN_CAND <- c(PSI_CAND, "pseudo_steppe", "Land_Cover_Type_1_Percent_Class_12",
              "Land_Cover_Type_1_Percent_Class_7")
MAX_DYN_SLOPES <- 3
NEEDS <- c("I(time^2)" = "time", "I(doy^2)" = "doy")
SITE_COVS <- c("bio1", "bio2", "tree_cover", "grass_cover", "topo_aspect", "topo_elev")
YEARLY_COVS <- c("NDVI", "EVI", "pr", "tmmn", "tmmx",
                 paste0("Land_Cover_Type_1_Percent_Class_", c(0, 6, 7, 10, 12, 13, 14)))
out_dir <- here("results", "model_selection")
dir.create(out_dir, showWarnings = FALSE)

sc <- function(m) (m - mean(m, na.rm = TRUE)) / sd(m, na.rm = TRUE)
visit_cols <- function(var, J) unlist(lapply(YEARS, function(yr) paste0(var, ".", seq_len(J), ".", yr)))
yearly <- function(d, var) as.matrix(d[, paste0(var, "_", YEARS)])            # M x T
expand_J <- function(m, J) m[, rep(seq_len(T_years), each = J), drop = FALSE]  # M x (J T)
form <- function(x) reformulate(if (length(x)) x else "1")

# --- Unmarked frame from the input table ---------------------------------------
build_frame <- function(csv) {
  d <- read.csv(csv, check.names = FALSE)
  J <- max(as.integer(sub("^y\\.([0-9]+)\\..*$", "\\1",
                          grep("^y\\.[0-9]+\\.[0-9]{4}$", names(d), value = TRUE))))
  y <- as.matrix(d[, visit_cols("y", J)]); storage.mode(y) <- "numeric"
  oc <- list(duration  = sc(as.matrix(d[, visit_cols("duration_minutes", J)])),
             effort    = sc(as.matrix(d[, visit_cols("effort_distance_km", J)])),
             observers = sc(as.matrix(d[, visit_cols("number_observers", J)])),
             time      = sc(as.matrix(d[, visit_cols("time_observations_started", J)])),
             doy       = sc(as.matrix(d[, visit_cols("day_of_year", J)])),
             steprep   = sc(as.matrix(d[, visit_cols("steprep_visit", J)])),
             NDVI_obs  = sc(expand_J(yearly(d, "NDVI"), J)))
  # one detection record for every candidate model
  n_before <- sum(!is.na(y))
  y[Reduce(`|`, lapply(oc, is.na))] <- NA
  message("  visits without a detection covariate set to NA: ", n_before - sum(!is.na(y)))
  ysc <- lapply(YEARLY_COVS, function(v) sc(yearly(d, v))); names(ysc) <- YEARLY_COVS
  hab <- read.csv(here("results", "core_spa", "cell_habitat.csv"))
  scv <- d[, SITE_COVS]
  scv$pseudo_steppe <- hab$steppe_strict[match(d$cells, hab$cells)]
  stopifnot(!anyNA(scv$pseudo_steppe))
  umf <- unmarkedMultFrame(y = y, siteCovs = as.data.frame(scale(scv)),
                           yearlySiteCovs = ysc,
                           obsCovs = lapply(oc, function(m) as.vector(t(m))) |> as.data.frame(),
                           numPrimary = T_years)
  list(d = d, J = J, y = y, umf = umf)
}

# --- Fit with the rejection rules -------------------------------------------------
# terms: list(psi, col, ext, det); `check` names the submodel whose standard
# errors are checked ("psi", "col", "ext", "p"), or "all"
fit_checked <- function(umf, terms, check = "all") {
  fit <- tryCatch(colext(psiformula = form(terms$psi), gammaformula = form(terms$col),
                         epsilonformula = form(terms$ext), pformula = form(terms$det),
                         data = umf, control = list(maxit = 2000)),
                  error = function(e) NULL)
  if (is.null(fit)) return(list(fit = NULL, status = "error", se = NULL))
  # store evaluated formulas in the call (parboot re-evaluates it)
  fit@call$psiformula <- form(terms$psi); fit@call$gammaformula <- form(terms$col)
  fit@call$epsilonformula <- form(terms$ext); fit@call$pformula <- form(terms$det)
  se <- suppressWarnings(sqrt(diag(vcov(fit))))
  sub <- if (identical(check, "all")) rep(TRUE, length(se)) else startsWith(names(se), paste0(check, "("))
  status <- if (fit@opt$convergence != 0) "not converged" else
    if (!all(is.finite(se[sub]))) "non-finite SE" else
    if (any(se[sub] > MAX_SE)) "SE > 3" else "ok"
  list(fit = fit, status = status, se = se)
}

# --- Forward selection for one submodel ---------------------------------------------
forward <- function(umf, cand, terms_for, check, submodel, sp, max_terms = Inf) {
  sel <- character(0)
  cur <- fit_checked(umf, terms_for(sel), check)
  if (cur$status != "ok") stop(sp, " ", submodel, ": starting model rejected (", cur$status, ")")
  log <- list(data.frame(species = sp, submodel = submodel, step = 0, term = "(start)",
                         aic = cur$fit@AIC, delta_vs_current = 0, status = "ok", added = TRUE))
  step <- 0
  repeat {
    avail <- setdiff(cand, sel)
    avail <- avail[is.na(NEEDS[avail]) | NEEDS[avail] %in% sel]
    if (!length(avail) || length(sel) >= max_terms) break
    step <- step + 1
    tries <- lapply(avail, function(v) fit_checked(umf, terms_for(c(sel, v)), check))
    aic <- sapply(tries, function(t) if (t$status == "ok") t$fit@AIC else NA)
    delta <- aic - cur$fit@AIC
    best <- if (all(is.na(aic))) NA else which.min(aic)
    add <- !is.na(best) && -delta[best] >= DELTA_AIC
    log[[length(log) + 1]] <- data.frame(
      species = sp, submodel = submodel, step = step, term = avail,
      aic = sapply(tries, function(t) if (is.null(t$fit)) NA else t$fit@AIC),
      delta_vs_current = delta, status = sapply(tries, `[[`, "status"),
      added = seq_along(avail) %in% if (add) best else integer(0))
    message("  ", submodel, " step ", step, ": ",
            if (add) paste0("+ ", avail[best], " (dAIC ", round(delta[best], 1), ")") else "stop")
    if (!add) break
    sel <- c(sel, avail[best]); cur <- tries[[best]]
  }
  list(sel = sel, fit = cur$fit, log = do.call(rbind, log))
}

# Detection record written back to the input table (NA where the model's
# detection covariates are missing), so that event code and tables agree
write_y <- function(fr, fit, csv) {
  y <- getY(fit@data)
  d <- fr$d
  d[visit_cols("y", fr$J)] <- as.data.frame(matrix(ifelse(is.na(y), NA, y == 1), nrow(y)))
  write.csv(d, csv, row.names = FALSE)
}

logs <- list(); chosen <- list()
for (sp in species) {
  message("=== ", sp)
  csv_p <- here("data", "processed_2023", sp, paste0(sp, "_occ_wide_dynamic.csv"))
  fr <- build_frame(csv_p)
  t0 <- Sys.time()

  # 1. detection, other submodels at all initial-occupancy candidates
  det <- forward(fr$umf, DET_CAND,
                 function(x) list(psi = PSI_CAND, col = PSI_CAND, ext = PSI_CAND, det = x),
                 "p", "detection", sp)
  # 2. initial occupancy, with the selected detection submodel
  psi <- forward(fr$umf, PSI_CAND,
                 function(x) list(psi = x, col = PSI_CAND, ext = PSI_CAND, det = det$sel),
                 "psi", "initial occupancy", sp)
  # 3. colonisation, then extinction, by forward selection (<= 3 slopes each)
  col <- forward(fr$umf, DYN_CAND,
                 function(x) list(psi = psi$sel, col = x, ext = psi$sel, det = det$sel),
                 "col", "colonisation", sp, MAX_DYN_SLOPES)
  ext <- forward(fr$umf, DYN_CAND,
                 function(x) list(psi = psi$sel, col = col$sel, ext = x, det = det$sel),
                 "ext", "extinction", sp, MAX_DYN_SLOPES)
  #    while a colonisation or extinction SE exceeds 3, drop the largest
  #    absolute slope of that submodel
  terms <- list(psi = psi$sel, col = col$sel, ext = ext$sel, det = det$sel)
  dropped <- character(0)
  repeat {
    f <- fit_checked(fr$umf, terms, "all")
    if (is.null(f$fit)) stop(sp, ": final model failed")
    se_ce <- f$se[grepl("^(col|ext)\\(", names(f$se))]
    bad <- se_ce[!is.finite(se_ce) | se_ce > MAX_SE]
    if (!length(bad)) break
    sub <- sub("\\(.*$", "", names(bad)[which.max(replace(bad, !is.finite(bad), Inf))])
    cf <- coef(f$fit)
    slopes <- cf[startsWith(names(cf), paste0(sub, "(")) & !grepl("Int\\)$", names(cf))]
    if (!length(slopes)) break
    worst <- names(slopes)[which.max(abs(slopes))]
    term <- sub("^[a-z]+\\((.*)\\)$", "\\1", worst)
    terms[[sub]] <- setdiff(terms[[sub]], term)
    dropped <- c(dropped, worst)
    message("  dropped ", worst, " (estimate ", round(slopes[worst], 1), ")")
  }
  fit <- f$fit
  saveRDS(fit, here("results", paste0(sp, "_model_object.rds")))
  write_y(fr, fit, csv_p)

  logs[[sp]] <- rbind(det$log, psi$log, col$log, ext$log)
  se <- f$se
  chosen[[sp]] <- data.frame(
    species = sp,
    psi = paste(deparse(fit@formlist$psi), collapse = ""),
    col = paste(deparse(fit@formlist$col), collapse = ""),
    ext = paste(deparse(fit@formlist$ext), collapse = ""),
    det = paste(deparse(fit@formlist$det), collapse = ""),
    dropped = paste(dropped, collapse = "; "),
    aic = fit@AIC, converged = fit@opt$convergence == 0, finite_se = all(is.finite(se)),
    max_se = max(se), status = f$status,
    minutes = as.numeric(Sys.time() - t0, units = "mins"))
  write.csv(logs[[sp]], file.path(out_dir, paste0(sp, "_selection_log.csv")), row.names = FALSE)
  write.csv(chosen[[sp]], file.path(out_dir, paste0(sp, "_selected_model.csv")), row.names = FALSE)
  print(chosen[[sp]])
}

# combined tables (all species, whichever run produced them)
cmb <- function(pat) {
  f <- list.files(out_dir, pattern = pat, full.names = TRUE)
  do.call(rbind, lapply(f, read.csv))
}
write.csv(cmb("_selection_log\\.csv$"), file.path(out_dir, "selection_log.csv"), row.names = FALSE)
write.csv(cmb("_selected_model\\.csv$"), file.path(out_dir, "selected_models.csv"), row.names = FALSE)

sessionInfo()
