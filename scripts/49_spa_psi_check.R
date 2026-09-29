###############################################################################
# 49_spa_psi_check.R
#
# Is the higher colonisation of SPA cells (script 46) an artefact of how first
# detections are split between initial occupancy and colonisation? SPA and
# pseudo-steppe cover entered colonisation and extinction but not initial
# occupancy; if initial occupancy under-predicts 2017 occupancy in SPA cells,
# first detections there are forced into colonisation. Checks, each fitted
# with and without SPA cover in colonisation and extinction:
#   psi_spa       model of script 46 + SPA and pseudo-steppe cover in initial
#                 occupancy
#   psi_spa_t     as psi_spa, with survey effort of year t only (effort in
#                 t + 1 can respond to the outcome)
#   from2018      as psi_spa_t, fitted to 2018-2023 (2017 dropped: initial
#                 occupancy estimated for 2018)
# Also reports predicted annual colonisation probability for a typical cell
# (covariates at their mean) fully outside and fully inside SPAs.
#
# Inputs:  as scripts 45-46
# Outputs: results/core_spa/spa_psi_check.csv
###############################################################################

library(here)
library(unmarked)
library(terra)

set.seed(42)

species <- c("otitar", "tettet", "pteori", "ptealc")
YEARS <- 2017:2023; T_years <- length(YEARS)
hab <- read.csv(here("results", "core_spa", "cell_habitat.csv"))
grid <- rast(here("data-raw", "data", "environmental_data", "environmental_data_occ",
                  "variables_spain.grd"), lyrs = 1)
rhs <- function(f) attr(terms(f), "term.labels")
form <- function(x) reformulate(if (length(x)) x else "1")
sc01 <- function(v) (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE)

rows <- list()
for (sp in species) {
  message("=== ", sp)
  mod <- readRDS(here("results", paste0(sp, "_model_object.rds")))
  fl <- mod@formlist; u <- mod@data
  tab <- read.csv(here("data", "processed_2023", sp, paste0(sp, "_occ_wide_dynamic.csv")))
  h <- hab[match(tab$cells, hab$cells), ]
  y <- getY(u); M <- nrow(y); J <- ncol(y) / T_years

  # neighbourhood covariates as in script 45
  surf <- read.csv(here("data", "processed_2023", sp, paste0("occ_", sp, "_prediction.csv")))
  psi_cell <- rep(NA_real_, ncell(grid))
  psi_cell[cellFromXY(grid, cbind(surf$longitude, surf$latitude))] <- surf$occ_prob
  adj <- adjacent(grid, tab$cells, directions = "queen", pairs = TRUE)
  nb <- tapply(psi_cell[adj[, 2]], adj[, 1], mean, na.rm = TRUE)
  nb <- as.numeric(nb[as.character(tab$cells)]); nb[!is.finite(nb)] <- 0
  base_det <- rowSums(y[, seq_len(2 * J)] == 1, na.rm = TRUE) > 0
  xy <- cbind(h$x, h$y) / 1000; src <- which(base_det)
  dist_det <- sapply(seq_len(M), function(i) {
    d <- sqrt((xy[src, 1] - xy[i, 1])^2 + (xy[src, 2] - xy[i, 2])^2); min(d[src != i]) })
  nv <- sapply(seq_len(T_years), function(t) rowSums(!is.na(y[, (t - 1) * J + seq_len(J), drop = FALSE])))

  sc <- siteCovs(u)
  sc$spa <- h$spa_frac; sc$x_c <- sc01(h$x); sc$y_c <- sc01(h$y)
  sc$nb_psi1 <- sc01(qlogis(pmin(pmax(nb, 1e-4), 1 - 1e-4))); sc$dist_det <- sc01(log1p(dist_det))
  ysc <- yearlySiteCovs(u)
  ysc$effort_t <- sc01(log1p(as.vector(t(nv))))
  ysc$effort_next <- sc01(log1p(as.vector(t(cbind(nv[, -1], 0)))))
  oc <- obsCovs(u); oc$spa_det <- rep(h$spa_frac, each = ncol(y))
  umf7 <- unmarkedMultFrame(y = y, siteCovs = sc, yearlySiteCovs = ysc, obsCovs = oc, numPrimary = T_years)
  # 2018-2023: drop the first season (columns, yearly rows and visit rows)
  keep_col <- (J + 1):ncol(y)
  yr_keep <- rep(seq_len(T_years) > 1, M)
  oc6 <- oc[rep(seq_len(ncol(y)) > J, M), ]
  umf6 <- unmarkedMultFrame(y = y[, keep_col], siteCovs = sc, yearlySiteCovs = ysc[yr_keep, ],
                            obsCovs = oc6, numPrimary = T_years - 1)

  psi_terms <- union(rhs(fl$psi), c("pseudo_steppe", "spa"))
  det_terms <- c(rhs(fl$det), "spa_det")
  checks <- list(
    psi_spa   = list(umf = umf7, ctrl = c("effort_t", "effort_next", "x_c", "y_c", "nb_psi1", "dist_det")),
    psi_spa_t = list(umf = umf7, ctrl = c("effort_t", "x_c", "y_c", "nb_psi1", "dist_det")),
    from2018  = list(umf = umf6, ctrl = c("effort_t", "x_c", "y_c", "nb_psi1", "dist_det")))
  for (ck in names(checks)) {
    s <- checks[[ck]]
    base_col <- c(union(rhs(fl$col), "pseudo_steppe"), s$ctrl); base_ext <- union(rhs(fl$ext), "pseudo_steppe")
    fits <- list()
    for (with_spa in c(FALSE, TRUE)) {
      col <- c(base_col, if (with_spa) "spa"); ext <- c(base_ext, if (with_spa) "spa")
      fits[[as.character(with_spa)]] <- tryCatch(
        colext(psiformula = form(psi_terms), gammaformula = form(col), epsilonformula = form(ext),
               pformula = form(det_terms), data = s$umf, control = list(maxit = 3000)),
        error = function(e) { message("  ", ck, " failed: ", conditionMessage(e)); NULL })
    }
    f1 <- fits[["TRUE"]]; f0 <- fits[["FALSE"]]
    if (is.null(f1)) next
    cf <- coef(f1)
    se <- tryCatch(suppressWarnings(sqrt(diag(vcov(f1)))),
                   error = function(e) setNames(rep(NA_real_, length(cf)), names(cf)))
    g0 <- plogis(cf[["col(Int)"]]); g1 <- plogis(cf[["col(Int)"]] + cf[["col(spa)"]])
    for (term in c("psi(spa)", "col(spa)", "ext(spa)", "p(spa_det)")) {
      if (!term %in% names(cf)) next
      rows[[paste(sp, ck, term)]] <- data.frame(
        species = sp, check = ck, term = term, estimate = cf[[term]], se = se[[term]],
        lo = cf[[term]] - 1.96 * se[[term]], hi = cf[[term]] + 1.96 * se[[term]],
        delta_aic_spa = if (is.null(f0)) NA else f1@AIC - f0@AIC,
        converged = f1@opt$convergence == 0, max_se = max(se, na.rm = TRUE),
        gamma_outside = g0, gamma_inside = g1)
    }
    print(do.call(rbind, rows[grepl(paste0("^", sp, " ", ck, " "), names(rows))]), digits = 3)
  }
}
out <- do.call(rbind, rows); rownames(out) <- NULL
write.csv(out, here("results", "core_spa", "spa_psi_check.csv"), row.names = FALSE)
print(out, digits = 3)

sessionInfo()
