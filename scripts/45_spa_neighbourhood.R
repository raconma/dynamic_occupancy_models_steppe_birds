###############################################################################
# 45_spa_neighbourhood.R
#
# Is the higher colonisation of SPA cells (scripts 43-44) recruitment by the
# network or spatial contagion? Colonisation of these species is expected
# near existing populations, and SPAs were designated where the birds were.
# Two static neighbourhood covariates are added to colonisation:
#   nb_psi1   mean predicted initial occupancy of the 8 neighbouring grid
#             cells (psi1 surface of script 40, which covers unmodelled cells)
#   dist_det  distance to the nearest other cell with a detection in 2017 or
#             2018 (baseline years only; log(1 + km), standardised)
# Each is fitted with and without SPA cover in colonisation and extinction,
# on top of the SPA model of script 43 (primary covariates + pseudo-steppe),
# with SPA cover in detection in all fits (observers target SPAs; script 44).
# If the SPA effect on colonisation vanishes once neighbourhood is
# controlled, it reflects contagion rather than recruitment by the network.
#
# Inputs:  results/{sp}_model_object.rds, data/processed_2023/{sp}/
#          occ_{sp}_prediction.csv (script 40), cell_habitat.csv (script 23)
# Outputs: results/core_spa/spa_neighbourhood.csv
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
  fl <- mod@formlist; u <- mod@data
  tab <- read.csv(here("data", "processed_2023", sp, paste0(sp, "_occ_wide_dynamic.csv")))
  h <- hab[match(tab$cells, hab$cells), ]
  y <- getY(u); M <- nrow(y); J <- ncol(y) / T_years

  # neighbourhood 1: mean predicted psi1 of the 8 neighbours
  surf <- read.csv(here("data", "processed_2023", sp, paste0("occ_", sp, "_prediction.csv")))
  psi_cell <- rep(NA_real_, ncell(grid))
  psi_cell[cellFromXY(grid, cbind(surf$longitude, surf$latitude))] <- surf$occ_prob
  adj <- adjacent(grid, tab$cells, directions = "queen", pairs = TRUE)
  nb_psi1 <- tapply(psi_cell[adj[, 2]], adj[, 1], mean, na.rm = TRUE)
  nb_psi1 <- as.numeric(nb_psi1[as.character(tab$cells)])
  nb_psi1[!is.finite(nb_psi1)] <- 0

  # neighbourhood 2: distance to the nearest other cell detected in 2017-2018
  base_det <- rowSums(y[, seq_len(2 * J)] == 1, na.rm = TRUE) > 0
  xy <- cbind(h$x, h$y) / 1000
  src <- xy[base_det, , drop = FALSE]; src_id <- which(base_det)
  dist_det <- sapply(seq_len(M), function(i) {
    d <- sqrt((src[, 1] - xy[i, 1])^2 + (src[, 2] - xy[i, 2])^2)
    min(d[src_id != i])
  })

  sc <- siteCovs(u)
  sc$spa <- h$spa_frac
  sc$nb_psi1 <- sc01(qlogis(pmin(pmax(nb_psi1, 1e-4), 1 - 1e-4)))
  sc$dist_det <- sc01(log1p(dist_det))
  oc <- obsCovs(u); oc$spa_det <- rep(h$spa_frac, each = ncol(y))
  umf <- unmarkedMultFrame(y = y, siteCovs = sc, yearlySiteCovs = yearlySiteCovs(u),
                           obsCovs = oc, numPrimary = T_years)

  base_col <- union(rhs(fl$col), "pseudo_steppe"); base_ext <- union(rhs(fl$ext), "pseudo_steppe")
  det_terms <- c(rhs(fl$det), "spa_det")
  cf0 <- coef(mod)
  for (nb in c("none", "nb_psi1", "dist_det", "both")) {
    add <- switch(nb, none = character(0), both = c("nb_psi1", "dist_det"), nb)
    fits <- list()
    for (with_spa in c(FALSE, TRUE)) {
      col <- c(base_col, add, if (with_spa) "spa"); ext <- c(base_ext, if (with_spa) "spa")
      st <- c(starts_from(cf0, "psi", rhs(fl$psi)), starts_from(cf0, "col", col),
              starts_from(cf0, "ext", ext), starts_from(cf0, "p", det_terms))
      fits[[as.character(with_spa)]] <- tryCatch(
        colext(psiformula = fl$psi, gammaformula = form(col), epsilonformula = form(ext),
               pformula = form(det_terms), data = umf, starts = st, control = list(maxit = 2000)),
        error = function(e) { message("  ", nb, " failed: ", conditionMessage(e)); NULL })
    }
    f1 <- fits[["TRUE"]]; f0 <- fits[["FALSE"]]
    if (is.null(f1)) next
    cf <- coef(f1); se <- suppressWarnings(sqrt(diag(vcov(f1))))
    for (term in c("col(spa)", "ext(spa)", "col(nb_psi1)", "col(dist_det)")) {
      if (!term %in% names(cf)) next
      rows[[paste(sp, nb, term)]] <- data.frame(
        species = sp, neighbourhood = nb, term = term, estimate = cf[[term]], se = se[[term]],
        lo = cf[[term]] - 1.96 * se[[term]], hi = cf[[term]] + 1.96 * se[[term]],
        delta_aic_spa = if (is.null(f0)) NA else f1@AIC - f0@AIC, aic = f1@AIC,
        converged = f1@opt$convergence == 0, max_se = max(se))
    }
    print(do.call(rbind, rows[grepl(paste0("^", sp, " ", nb, " "), names(rows))]), digits = 3)
  }
}
out <- do.call(rbind, rows); rownames(out) <- NULL
write.csv(out, here("results", "core_spa", "spa_neighbourhood.csv"), row.names = FALSE)
print(out, digits = 3)

sessionInfo()
