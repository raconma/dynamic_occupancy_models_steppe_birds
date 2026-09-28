###############################################################################
# 41_covariate_selection_pooled.R
#
# Covariate selection on the pooled detection histories (script 39) and fit
# of the models used downstream.
#
# The covariate structures of the earlier fits were chosen on single-site
# histories and the selection was not recorded in code. Here it is redone on
# the pooled histories with rules fixed before fitting:
#
#   colonisation and extinction  constant (gamma ~ 1, epsilon ~ 1), so that
#       spatial contrasts in the anchored dynamics (core/marginal, SPA) arise
#       from the detection histories and not from covariates
#   detection    forward selection by AIC among survey duration, distance
#       travelled, number of observers, start time (+ quadratic), day of year
#       (+ quadratic) and NDVI; a quadratic term can enter only after its
#       linear term. Initial occupancy is held at all its candidates.
#   initial occupancy  forward selection by AIC among bio1, bio2, tree
#       cover, herbaceous cover and aspect (elevation is excluded: |r| > 0.8
#       with bio1), with the selected detection submodel.
#   A term is added when it lowers AIC by at least 2. A candidate fit is
#   rejected when the optimiser does not converge or a standard error is not
#   finite or exceeds 3 (logit scale, standardised covariates), which is how
#   separation shows. (A first rule on slope magnitude, |slope| > 5, was
#   dropped before any selection ran: it rejected the starting model of
#   O. tarda for its tree-cover effect, which is well estimated: -7.6,
#   SE 1.1, in the final model.)
#
# All candidate fits use the same detection record: visits with a missing
# value of any candidate detection covariate carry no detection. Checklist,
# day-of-year and yearly covariates are standardised over all cell-years
# (the earlier fits standardised yearly covariates within each year, which
# removes their interannual variation); site covariates over cells.
#
# The stepRep sensitivity model (Appendix S2) adds within-cell steppe
# representativeness to the selected detection submodel.
#
# Inputs:  data/processed_2023/{sp}/{sp}_occ_wide_dynamic.csv and
#          reports/script_4_v2_inputs/{sp}_occ_wide_dynamic.csv (script 39)
# Outputs: results/{sp}_model_object.rds                      (primary)
#          results/robustness/{sp}_steprep_v3_model_object.rds (stepRep)
#          the input tables above, with y set to NA where the model's
#          detection covariates are missing (rows as in the models)
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
              "doy", "I(doy^2)", "NDVI_obs")
PSI_CAND <- c("bio1", "bio2", "tree_cover", "grass_cover", "topo_aspect")
NEEDS <- c("I(time^2)" = "time", "I(doy^2)" = "doy")
SITE_COVS <- c("bio1", "bio2", "tree_cover", "grass_cover", "topo_aspect", "topo_elev")
YEARLY_COVS <- c("NDVI", "EVI", "pr", "tmmn", "tmmx",
                 paste0("Land_Cover_Type_1_Percent_Class_", c(0, 6, 7, 10, 12, 13, 14)))
out_dir <- here("results", "model_selection")
dir.create(out_dir, showWarnings = FALSE)
dir.create(here("results", "robustness"), showWarnings = FALSE)

sc <- function(m) (m - mean(m, na.rm = TRUE)) / sd(m, na.rm = TRUE)
visit_cols <- function(var, J) unlist(lapply(YEARS, function(yr) paste0(var, ".", seq_len(J), ".", yr)))
yearly <- function(d, var) as.matrix(d[, paste0(var, "_", YEARS)])            # M x T
expand_J <- function(m, J) m[, rep(seq_len(T_years), each = J), drop = FALSE]  # M x (J T)

# --- Unmarked frame from an input table ----------------------------------------
build_frame <- function(csv, extra_det = NULL) {
  d <- read.csv(csv, check.names = FALSE)
  J <- max(as.integer(sub("^y\\.([0-9]+)\\..*$", "\\1",
                          grep("^y\\.[0-9]+\\.[0-9]{4}$", names(d), value = TRUE))))
  y <- as.matrix(d[, visit_cols("y", J)]); storage.mode(y) <- "numeric"
  oc <- list(duration  = sc(as.matrix(d[, visit_cols("duration_minutes", J)])),
             effort    = sc(as.matrix(d[, visit_cols("effort_distance_km", J)])),
             observers = sc(as.matrix(d[, visit_cols("number_observers", J)])),
             time      = sc(as.matrix(d[, visit_cols("time_observations_started", J)])),
             doy       = sc(as.matrix(d[, visit_cols("day_of_year", J)])),
             NDVI_obs  = sc(expand_J(yearly(d, "NDVI"), J)))
  for (v in names(extra_det)) oc[[v]] <- sc(expand_J(yearly(d, extra_det[[v]]), J))
  # one detection record for every candidate model
  n_before <- sum(!is.na(y))
  y[Reduce(`|`, lapply(oc, is.na))] <- NA
  message("  visits without a detection covariate set to NA: ", n_before - sum(!is.na(y)))
  ysc <- lapply(YEARLY_COVS, function(v) sc(yearly(d, v))); names(ysc) <- YEARLY_COVS
  umf <- unmarkedMultFrame(y = y, siteCovs = as.data.frame(scale(d[, SITE_COVS])),
                           yearlySiteCovs = ysc,
                           obsCovs = lapply(oc, function(m) as.vector(t(m))) |> as.data.frame(),
                           numPrimary = T_years)
  list(d = d, J = J, y = y, umf = umf)
}

# --- Fit with the rejection rules -------------------------------------------------
fit_checked <- function(umf, psi_terms, det_terms) {
  form <- function(x) reformulate(if (length(x)) x else "1")
  fit <- tryCatch(colext(psiformula = form(psi_terms), gammaformula = ~ 1,
                         epsilonformula = ~ 1, pformula = form(det_terms),
                         data = umf, control = list(maxit = 2000)),
                  error = function(e) NULL)
  if (is.null(fit)) return(list(fit = NULL, status = "error"))
  se <- suppressWarnings(sqrt(diag(vcov(fit))))
  status <- if (fit@opt$convergence != 0) "not converged" else
    if (!all(is.finite(se))) "non-finite SE" else
    if (any(se > MAX_SE)) "SE > 3" else "ok"
  list(fit = fit, status = status)
}

# --- Forward selection for one submodel ---------------------------------------------
forward <- function(umf, cand, terms_for, submodel, sp) {
  sel <- character(0)
  cur <- do.call(fit_checked, c(list(umf), terms_for(sel)))
  if (cur$status != "ok") stop(sp, " ", submodel, ": starting model rejected (", cur$status, ")")
  log <- list(data.frame(species = sp, submodel = submodel, step = 0, term = "(start)",
                         aic = cur$fit@AIC, delta_vs_current = 0, status = "ok", added = TRUE))
  step <- 0
  repeat {
    avail <- setdiff(cand, sel)
    avail <- avail[is.na(NEEDS[avail]) | NEEDS[avail] %in% sel]
    if (!length(avail)) break
    step <- step + 1
    tries <- lapply(avail, function(v) do.call(fit_checked, c(list(umf), terms_for(c(sel, v)))))
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
  csv_s <- here("reports", "script_4_v2_inputs", paste0(sp, "_occ_wide_dynamic.csv"))
  fr <- build_frame(csv_p)

  # 1. detection, initial occupancy at all candidates
  t0 <- Sys.time()
  det <- forward(fr$umf, DET_CAND, function(x) list(PSI_CAND, x), "detection", sp)
  # 2. initial occupancy, with the selected detection submodel
  psi <- forward(fr$umf, PSI_CAND, function(x) list(x, det$sel), "initial occupancy", sp)
  fit <- psi$fit
  saveRDS(fit, here("results", paste0(sp, "_model_object.rds")))
  write_y(fr, fit, csv_p)

  # 3. stepRep sensitivity: selected structure + steppe representativeness
  frs <- build_frame(csv_s, extra_det = list(stepRep_obs = "stepRep_strict_500m"))
  fs <- fit_checked(frs$umf, psi$sel, c(det$sel, "stepRep_obs"))
  message("  stepRep model: ", fs$status)
  if (!is.null(fs$fit)) {
    saveRDS(fs$fit, here("results", "robustness", paste0(sp, "_steprep_v3_model_object.rds")))
    write_y(frs, fs$fit, csv_s)
  }

  logs[[sp]] <- rbind(det$log, psi$log)
  se <- suppressWarnings(sqrt(diag(vcov(fit))))
  chosen[[sp]] <- data.frame(
    species = sp,
    psi = paste(deparse(fit@formlist$psi), collapse = ""),
    det = paste(deparse(fit@formlist$det), collapse = ""),
    aic = fit@AIC, converged = fit@opt$convergence == 0, finite_se = all(is.finite(se)),
    steprep_status = fs$status, steprep_aic = if (is.null(fs$fit)) NA else fs$fit@AIC,
    gamma = plogis(coef(fit)["col(Int)"]), epsilon = plogis(coef(fit)["ext(Int)"]),
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
