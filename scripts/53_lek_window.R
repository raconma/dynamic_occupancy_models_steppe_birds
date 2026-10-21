###############################################################################
# 53_lek_window.R
#
# Robustness plan C (docs/robustness_plan_BiolCons.md): is the higher great
# bustard extinction in marginal habitat an artefact of birds leaving their
# leks within the April-June breeding window? Keeps only checklists made from
# 1 April to 15 May (lek period, before males move to post-breeding areas;
# Palacin et al. 2009), refits the primary formulas without reselection and
# saves the fit as the "lek" variant, so that scripts 26 and 23 can summarise
# it:
#   SPECIES=otitar Rscript scripts/26_naive_vs_corrected_anchored.R lek
#   SPECIES=otitar Rscript scripts/23_core_spa_analysis.R lek
#
# Inputs:  results/otitar_model_object.rds, data/processed_2023/otitar/
#          otitar_occ_wide_dynamic.csv (day_of_year of each visit)
# Outputs: results/robustness/otitar_lek_model_object.rds
###############################################################################

library(here)
library(unmarked)

set.seed(42)

SP <- "otitar"
YEARS <- 2017:2023
mod <- readRDS(here("results", paste0(SP, "_model_object.rds")))
fl <- mod@formlist; umf <- mod@data
d <- read.csv(here("data", "processed_2023", SP, paste0(SP, "_occ_wide_dynamic.csv")), check.names = FALSE)
y <- getY(umf); J <- ncol(y) / length(YEARS)
stopifnot(nrow(d) == nrow(y))

# last day of the lek window: 15 May (day 135, 136 in leap years)
doy <- as.matrix(d[, unlist(lapply(YEARS, function(yr) paste0("day_of_year.", seq_len(J), ".", yr)))])
last <- rep(ifelse(YEARS %% 4 == 0, 136, 135), each = J)
late <- !is.na(doy) & sweep(doy, 2, last, ">")
message("visits before: ", sum(!is.na(y)), "; dropped after 15 May: ", sum(late & !is.na(y)))
y[late] <- NA
umf@y <- y

fit <- colext(psiformula = fl$psi, gammaformula = fl$col, epsilonformula = fl$ext,
              pformula = fl$det, data = umf, starts = coef(mod), control = list(maxit = 3000))
message("converged: ", fit@opt$convergence == 0, "; largest SE: ",
        round(max(sqrt(diag(vcov(fit)))), 2))
dir.create(here("results", "robustness"), showWarnings = FALSE)
saveRDS(fit, here("results", "robustness", paste0(SP, "_lek_model_object.rds")))
print(cbind(primary = coef(mod), lek = coef(fit)), digits = 3)

sessionInfo()
