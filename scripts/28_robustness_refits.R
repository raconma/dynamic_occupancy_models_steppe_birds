###############################################################################
# 28_robustness_refits.R
#
# Refits the robustness variant of the colext models in which colonisation
# and extinction change linearly over the study period:
#
#   trend       primary model (results/{sp}_model_object.rds, script 41) with
#               gamma ~ t and epsilon ~ t, t = year - 2020 (no habitat
#               covariates, so spatial contrasts still come from the
#               detection histories)
#
# The stepRep variant (steprep_v3) is fitted by script 41. A variant with a
# separate gamma and epsilon for each year was not identifiable (years with
# almost no events sent their effects to the boundary; singular Hessian).
# Starting values: psi and p from the primary fit, the primary intercepts
# for gamma and epsilon and zero slopes. Fits whose Hessian gives non-finite standard errors are flagged:
# they cannot be propagated with multivariate-normal draws.
#
# Outputs: results/robustness/{sp}_trend_model_object.rds
#          results/robustness/refit_summary.csv
###############################################################################

library(here)
library(unmarked)

set.seed(42)

species <- c("otitar", "ptealc", "pteori", "tettet")
YEARS <- 2017:2023
out_dir <- here("results", "robustness")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

with_trend <- function(umf) {
  M <- numSites(umf); n_t <- umf@numPrimary
  ysc <- yearlySiteCovs(umf)
  ysc$t <- rep(YEARS, M) - 2020
  unmarkedMultFrame(y = getY(umf), siteCovs = siteCovs(umf), yearlySiteCovs = ysc,
                    obsCovs = obsCovs(umf), numPrimary = n_t)
}

rows <- list()
for (sp in species) {
  primary <- readRDS(here("results", paste0(sp, "_model_object.rds")))
  fl <- primary@formlist
  cf <- coef(primary)
  starts <- unname(c(cf[grep("^psi\\(", names(cf))], cf["col(Int)"], 0,
                     cf["ext(Int)"], 0, cf[grep("^p\\(", names(cf))]))
  message(sp, " / trend ...")
  t0 <- Sys.time()
  fit <- tryCatch(colext(psiformula = fl$psi, gammaformula = ~ t, epsilonformula = ~ t,
                         pformula = fl$det, data = with_trend(primary@data), starts = starts,
                         control = list(maxit = 2000)),
                  error = function(e) { message("  failed: ", conditionMessage(e)); NULL })
  if (is.null(fit)) {
    rows[[sp]] <- data.frame(species = sp, variant = "trend", converged = FALSE,
                             finite_se = NA, max_se = NA, aic = NA, aic_primary = primary@AIC,
                             fit_seconds = NA)
    next
  }
  saveRDS(fit, file.path(out_dir, paste0(sp, "_trend_model_object.rds")))
  se <- tryCatch(suppressWarnings(sqrt(diag(vcov(fit)))),
                 error = function(e) rep(NA_real_, length(coef(fit))))   # singular Hessian
  rows[[sp]] <- data.frame(
    species = sp, variant = "trend", converged = fit@opt$convergence == 0,
    finite_se = all(is.finite(se)), max_se = max(se, na.rm = TRUE), aic = fit@AIC, aic_primary = primary@AIC,
    fit_seconds = as.numeric(Sys.time() - t0, units = "secs"))
  print(rows[[sp]])
  print(round(cbind(est = coef(fit), se = se)[grep("^(col|ext)", names(coef(fit))), ], 2))
}

res <- do.call(rbind, rows); rownames(res) <- NULL
write.csv(res, file.path(out_dir, "refit_summary.csv"), row.names = FALSE)
print(res, digits = 3)

sessionInfo()
