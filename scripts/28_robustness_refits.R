###############################################################################
# 28_robustness_refits.R
#
# Refits three robustness variants of the primary colext models (script 41):
#
#   constant    colonisation and extinction constant (gamma ~ 1, epsilon ~ 1)
#   trend       colonisation and extinction on the primary habitat covariates
#               plus a linear trend, t = year - 2020
#   no_steprep  detection without within-cell steppe representativeness
#               (identical to the primary model where it was not selected)
#
# A variant with a separate gamma and epsilon for each year was not
# identifiable (years with almost no events sent their effects to the
# boundary; singular Hessian). Starting values come from the primary fit
# for the shared parameters and zero for the new ones. Fits whose Hessian
# gives non-finite standard errors are flagged: they cannot be propagated
# with multivariate-normal draws.
#
# Outputs: results/robustness/{sp}_{variant}_model_object.rds
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

rhs <- function(f) attr(terms(f), "term.labels")
form <- function(x) reformulate(if (length(x)) x else "1")
# starting values: primary estimates for parameters kept, 0 for new slopes
starts_from <- function(cf, sub, labels) {
  unname(sapply(c("(Intercept)", labels), function(l) {
    nm <- paste0(sub, "(", if (l == "(Intercept)") "Int" else l, ")")
    if (nm %in% names(cf)) cf[[nm]] else 0
  }))
}

rows <- list()
for (sp in species) {
  primary <- readRDS(here("results", paste0(sp, "_model_object.rds")))
  fl <- primary@formlist
  cf <- coef(primary)
  variants <- list(
    constant   = list(col = character(0), ext = character(0), det = rhs(fl$det)),
    trend      = list(col = c(rhs(fl$col), "t"), ext = c(rhs(fl$ext), "t"), det = rhs(fl$det)),
    no_steprep = list(col = rhs(fl$col), ext = rhs(fl$ext), det = setdiff(rhs(fl$det), "steprep")))
  for (v in names(variants)) {
    s <- variants[[v]]
    message(sp, " / ", v, " ...")
    starts <- c(starts_from(cf, "psi", rhs(fl$psi)), starts_from(cf, "col", s$col),
                starts_from(cf, "ext", s$ext), starts_from(cf, "p", s$det))
    t0 <- Sys.time()
    fit <- tryCatch(colext(psiformula = fl$psi, gammaformula = form(s$col), epsilonformula = form(s$ext),
                           pformula = form(s$det), data = with_trend(primary@data), starts = starts,
                           control = list(maxit = 2000)),
                    error = function(e) { message("  failed: ", conditionMessage(e)); NULL })
    if (is.null(fit)) {
      rows[[paste(sp, v)]] <- data.frame(species = sp, variant = v, converged = FALSE,
                                         finite_se = NA, max_se = NA, aic = NA,
                                         aic_primary = primary@AIC, fit_seconds = NA)
      next
    }
    saveRDS(fit, file.path(out_dir, paste0(sp, "_", v, "_model_object.rds")))
    se <- tryCatch(suppressWarnings(sqrt(diag(vcov(fit)))),
                   error = function(e) rep(NA_real_, length(coef(fit))))   # singular Hessian
    rows[[paste(sp, v)]] <- data.frame(
      species = sp, variant = v, converged = fit@opt$convergence == 0,
      finite_se = all(is.finite(se)), max_se = max(se, na.rm = TRUE), aic = fit@AIC,
      aic_primary = primary@AIC, fit_seconds = as.numeric(Sys.time() - t0, units = "secs"))
    print(rows[[paste(sp, v)]])
  }
}

res <- do.call(rbind, rows); rownames(res) <- NULL
write.csv(res, file.path(out_dir, "refit_summary.csv"), row.names = FALSE)
print(res, digits = 3)

sessionInfo()
