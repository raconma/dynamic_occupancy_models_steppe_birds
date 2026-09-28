###############################################################################
# 43_spa_in_dynamics.R
#
# Model-based view of SPA effectiveness for range dynamics. Instead of
# matching recorded cells after fitting (scripts 23, 36), the fraction of
# each cell inside SPAs enters the colonisation and extinction submodels of
# the primary model, so that the two functions of a protected area are
# estimated jointly with detection and all cells:
#   retention    effect of SPA cover on extinction (negative = SPAs retain)
#   recruitment  effect of SPA cover on colonisation (positive = SPAs gain)
# Both submodels are adjusted for pseudo-steppe fraction (added when not
# selected) besides the covariates selected in script 41.
#
# A second model splits SPA cover into SPAs whose Natura 2000 Standard Data
# Form declares a high-intensity negative agricultural pressure of
# intensification (threat codes A02 modification of cultivation practices,
# A07 biocides and fertilisers, A09 irrigation, A10 restructuring of
# agricultural land holdings) and the remaining SPAs. (Management-plan
# status could not be used: 139 of the 146 SPAs listing the four species
# declare an approved plan.)
#
# SPA cover terms are proportions of the cell (0 to 1, not standardised):
# coefficients are log-odds differences between a cell fully inside and a
# cell fully outside the SPAs. Observational: other differences between
# protected and unprotected cells can remain.
#
# Inputs:  results/{sp}_model_object.rds (script 41), cell_habitat.csv
#          (script 23), SPA polygons, Natura 2000 tabular IMPACT table
# Outputs: results/core_spa/cell_spa_pressure.csv
#          results/core_spa/spa_in_dynamics.csv
###############################################################################

library(here)
library(unmarked)
library(sf)

source(here("R", "spatial_layers.R"))
set.seed(42)

species <- c("otitar", "tettet", "pteori", "ptealc")
hab_dir <- here("results", "core_spa")
tab_dir <- here("data-raw", "data", "natura2000_tabular", "csv")
INTENS <- c("A02", "A07", "A09", "A10")

# --- 1. SPA cover split by declared agricultural intensification pressure ----
hab <- read.csv(file.path(hab_dir, "cell_habitat.csv"))
pr_path <- file.path(hab_dir, "cell_spa_pressure.csv")
if (!file.exists(pr_path)) {
  imp <- read.csv(file.path(tab_dir, "Natura2000_end2024_rev1_IMPACT.csv"), fileEncoding = "UTF-8-BOM")
  hi <- unique(imp$SITECODE[substr(imp$IMPACTCODE, 1, 3) %in% INTENS &
                              imp$INTENSITY == "HIGH" & imp$IMPACT_TYPE == "N"])
  spa <- load_spa_sites()
  code <- sub("_SPA$", "", spa$localId)
  polys <- wc_cell_polygons(hab$cells)
  write.csv(data.frame(cells = hab$cells,
                       spa_hi = cover_by_polygons(polys, spa[code %in% hi, ]),
                       spa_lo = cover_by_polygons(polys, spa[!code %in% hi, ])),
            pr_path, row.names = FALSE)
}
prs <- read.csv(pr_path)

rhs <- function(f) attr(terms(f), "term.labels")
form <- function(x) reformulate(if (length(x)) x else "1")
starts_from <- function(cf, sub, labels) {
  unname(sapply(c("(Intercept)", labels), function(l) {
    nm <- paste0(sub, "(", if (l == "(Intercept)") "Int" else l, ")")
    if (nm %in% names(cf)) cf[[nm]] else 0
  }))
}

rows <- list()
for (sp in species) {
  message("=== ", sp)
  mod <- readRDS(here("results", paste0(sp, "_model_object.rds")))
  fl <- mod@formlist
  tab <- read.csv(here("data", "processed_2023", sp, paste0(sp, "_occ_wide_dynamic.csv")))
  stopifnot(nrow(tab) == numSites(mod@data))
  sc <- siteCovs(mod@data)
  sc$spa <- hab$spa_frac[match(tab$cells, hab$cells)]
  sc$spa_hi <- prs$spa_hi[match(tab$cells, prs$cells)]
  sc$spa_lo <- prs$spa_lo[match(tab$cells, prs$cells)]
  stopifnot(!anyNA(sc$spa), !anyNA(sc$spa_hi))
  umf <- mod@data; siteCovs(umf) <- sc
  base_col <- union(rhs(fl$col), "pseudo_steppe"); base_ext <- union(rhs(fl$ext), "pseudo_steppe")
  specs <- list(adjusted = list(col = base_col, ext = base_ext),
                spa = list(col = c(base_col, "spa"), ext = c(base_ext, "spa")),
                pressure = list(col = c(base_col, "spa_hi", "spa_lo"), ext = c(base_ext, "spa_hi", "spa_lo")))
  cf0 <- coef(mod)
  for (nm in names(specs)) {
    s <- specs[[nm]]
    st <- c(starts_from(cf0, "psi", rhs(fl$psi)), starts_from(cf0, "col", s$col),
            starts_from(cf0, "ext", s$ext), starts_from(cf0, "p", rhs(fl$det)))
    fit <- tryCatch(colext(psiformula = fl$psi, gammaformula = form(s$col), epsilonformula = form(s$ext),
                           pformula = fl$det, data = umf, starts = st, control = list(maxit = 2000)),
                    error = function(e) { message("  ", nm, " failed: ", conditionMessage(e)); NULL })
    if (is.null(fit)) next
    cf <- coef(fit); se <- suppressWarnings(sqrt(diag(vcov(fit))))
    keep <- grepl("\\((spa|spa_hi|spa_lo)\\)$", names(cf))
    rows[[paste(sp, nm)]] <- data.frame(
      species = sp, model = nm, aic = fit@AIC, aic_primary = mod@AIC,
      converged = fit@opt$convergence == 0, max_se = max(se),
      term = if (any(keep)) names(cf)[keep] else NA,
      estimate = if (any(keep)) cf[keep] else NA, se = if (any(keep)) se[keep] else NA,
      lo = if (any(keep)) cf[keep] - 1.96 * se[keep] else NA,
      hi = if (any(keep)) cf[keep] + 1.96 * se[keep] else NA, row.names = NULL)
    print(rows[[paste(sp, nm)]], digits = 3)
  }
}
out <- do.call(rbind, rows); rownames(out) <- NULL
out$odds_ratio <- exp(out$estimate)
write.csv(out, file.path(hab_dir, "spa_in_dynamics.csv"), row.names = FALSE)
print(out, digits = 3)

sessionInfo()
