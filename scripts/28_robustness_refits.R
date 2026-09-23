###############################################################################
# 28_robustness_refits.R
#
# Refits two robustness variants of the colext models so that the spatial
# contrasts of scripts 23 and 26 can be repeated with extinction constant
# across cells (any core / SPA difference must then come from the data, not
# from the extinction covariates):
#
#   eps_null    primary model (results/{sp}_model_object.rds) with epsilon ~ 1
#   steprep_v3  stepRep detection model (results/stepRep_v2_run/
#               {sp}_v2_model_object.rds, detection + stepRep_obs) with
#               gamma ~ 1 and epsilon ~ 1 (the reduced specification of
#               Appendix S19; here epsilon ~ 1 also for P. orientalis)
#
# Starting values come from the source fit for the shared parameters.
# Fits whose Hessian gives non-finite standard errors are flagged: they
# cannot be propagated with multivariate-normal draws.
#
# Outputs: results/robustness/{sp}_{variant}_model_object.rds
#          results/robustness/refit_summary.csv
###############################################################################

library(here)
library(unmarked)

set.seed(42)

species <- c("otitar", "ptealc", "pteori", "tettet")
out_dir <- here("results", "robustness")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Starting values: source estimates for the submodels that are kept; for
# submodels reduced to an intercept, fixed plausible values (the source
# intercepts can be degenerate, e.g. ext(Int) = 52 in the otitar v2 fit)
START_COL <- qlogis(0.002)
START_EXT <- qlogis(0.08)
starts_for <- function(src, keep_col) {
  cf <- coef(src)
  pick <- function(pat) cf[grep(pat, names(cf))]
  col0 <- if (keep_col) pick("^col\\(") else START_COL
  unname(c(pick("^psi\\("), col0, START_EXT, pick("^p\\(")))
}

refit <- function(src, col_form, ext_form, starts) {
  fl <- src@formlist
  t0 <- Sys.time()
  fit <- colext(psiformula = fl$psi, gammaformula = col_form,
                epsilonformula = ext_form, pformula = fl$det,
                data = src@data, starts = starts)
  attr(fit, "fit_seconds") <- as.numeric(Sys.time() - t0, units = "secs")
  fit
}

rows <- list()
for (sp in species) {
  primary <- readRDS(here("results", paste0(sp, "_model_object.rds")))
  steprep <- readRDS(here("results", "stepRep_v2_run", paste0(sp, "_v2_model_object.rds")))

  variants <- list(
    eps_null = list(src = primary, col = primary@formlist$col, ext = ~ 1,
                    keep_col = TRUE),
    steprep_v3 = list(src = steprep, col = ~ 1, ext = ~ 1, keep_col = FALSE))

  for (v in names(variants)) {
    s <- variants[[v]]
    message(sp, " / ", v, " ...")
    fit <- tryCatch(refit(s$src, s$col, s$ext, starts_for(s$src, s$keep_col)),
                    error = function(e) { message("  failed: ", conditionMessage(e)); NULL })
    if (is.null(fit)) {
      rows[[paste(sp, v)]] <- data.frame(species = sp, variant = v, converged = FALSE,
                                         finite_se = NA, aic = NA, aic_source = s$src@AIC,
                                         gamma_int = NA, eps_int = NA, fit_seconds = NA)
      next
    }
    saveRDS(fit, file.path(out_dir, paste0(sp, "_", v, "_model_object.rds")))
    se <- suppressWarnings(sqrt(diag(vcov(fit))))
    cf <- coef(fit)
    rows[[paste(sp, v)]] <- data.frame(
      species = sp, variant = v, converged = fit@opt$convergence == 0,
      finite_se = all(is.finite(se)), aic = fit@AIC, aic_source = s$src@AIC,
      gamma_int = plogis(cf["col(Int)"]), eps_int = plogis(cf["ext(Int)"]),
      fit_seconds = attr(fit, "fit_seconds"))
  }
}

res <- do.call(rbind, rows); rownames(res) <- NULL
write.csv(res, file.path(out_dir, "refit_summary.csv"), row.names = FALSE)
print(res, digits = 3)

sessionInfo()
