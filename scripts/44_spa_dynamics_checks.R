###############################################################################
# 44_spa_dynamics_checks.R
#
# Robustness of the model-based SPA effects on colonisation (recruitment) and
# extinction (retention) of script 43. Each check refits the SPA model of
# script 43 (primary covariates + pseudo-steppe + SPA cover in colonisation
# and extinction) and the same model without SPA cover:
#   effort     checklists per cell-year in years t and t + 1 (log(1 + n),
#              standardised) added to colonisation and extinction: are
#              apparent colonisations inside SPAs a product of survey effort?
#   det_spa    SPA cover added to detection: do checklists inside SPAs
#              detect the species more often (observers targeting known
#              sites), which the model could otherwise read as colonisation?
#   coords     cell coordinates (standardised) added to colonisation and
#              extinction: regional differences in policy and in where SPAs
#              lie
#   recorded   fitted to cells with at least one detection only (the cells
#              of the anchored and matched analyses)
#   all        effort + det_spa + coords together
#
# Outputs: results/core_spa/spa_in_dynamics_checks.csv
###############################################################################

library(here)
library(unmarked)

set.seed(42)

species <- c("otitar", "tettet", "pteori", "ptealc")
hab <- read.csv(here("results", "core_spa", "cell_habitat.csv"))
YEARS <- 2017:2023; T_years <- length(YEARS)

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
  # checklists per cell-year and the same for the following year
  nv <- sapply(seq_len(T_years), function(t) rowSums(!is.na(y[, (t - 1) * J + seq_len(J), drop = FALSE])))
  nv_next <- cbind(nv[, -1], 0)
  sc <- siteCovs(u)
  sc$spa <- h$spa_frac; sc$x_c <- sc01(h$x); sc$y_c <- sc01(h$y)
  ysc <- yearlySiteCovs(u)
  ysc$effort_t <- sc01(log1p(as.vector(t(nv)))); ysc$effort_next <- sc01(log1p(as.vector(t(nv_next))))
  oc <- obsCovs(u); oc$spa_det <- rep(h$spa_frac, each = ncol(y))
  umf <- unmarkedMultFrame(y = y, siteCovs = sc, yearlySiteCovs = ysc, obsCovs = oc, numPrimary = T_years)
  det_cells <- which(rowSums(y == 1, na.rm = TRUE) > 0)

  base_col <- union(rhs(fl$col), "pseudo_steppe"); base_ext <- union(rhs(fl$ext), "pseudo_steppe")
  checks <- list(
    effort   = list(add = c("effort_t", "effort_next"), det = character(0), rows = NULL),
    det_spa  = list(add = character(0), det = "spa_det", rows = NULL),
    coords   = list(add = c("x_c", "y_c"), det = character(0), rows = NULL),
    recorded = list(add = character(0), det = character(0), rows = det_cells),
    all      = list(add = c("effort_t", "effort_next", "x_c", "y_c"), det = "spa_det", rows = NULL))
  cf0 <- coef(mod)
  for (ck in names(checks)) {
    s <- checks[[ck]]
    dat <- if (is.null(s$rows)) umf else umf[s$rows, ]
    det_terms <- c(rhs(fl$det), s$det)
    fits <- list()
    for (with_spa in c(FALSE, TRUE)) {
      col <- c(base_col, s$add, if (with_spa) "spa"); ext <- c(base_ext, s$add, if (with_spa) "spa")
      st <- c(starts_from(cf0, "psi", rhs(fl$psi)), starts_from(cf0, "col", col),
              starts_from(cf0, "ext", ext), starts_from(cf0, "p", det_terms))
      fits[[as.character(with_spa)]] <- tryCatch(
        colext(psiformula = fl$psi, gammaformula = form(col), epsilonformula = form(ext),
               pformula = form(det_terms), data = dat, starts = st, control = list(maxit = 2000)),
        error = function(e) { message("  ", ck, " failed: ", conditionMessage(e)); NULL })
    }
    f1 <- fits[["TRUE"]]; f0 <- fits[["FALSE"]]
    if (is.null(f1)) next
    cf <- coef(f1); se <- suppressWarnings(sqrt(diag(vcov(f1))))
    for (term in c("col(spa)", "ext(spa)", "p(spa_det)")) {
      if (!term %in% names(cf)) next
      rows[[paste(sp, ck, term)]] <- data.frame(
        species = sp, check = ck, term = term, estimate = cf[[term]], se = se[[term]],
        lo = cf[[term]] - 1.96 * se[[term]], hi = cf[[term]] + 1.96 * se[[term]],
        odds_ratio = exp(cf[[term]]), delta_aic_spa = if (is.null(f0)) NA else f1@AIC - f0@AIC,
        converged = f1@opt$convergence == 0, max_se = max(se), n_sites = numSites(dat))
    }
    print(do.call(rbind, rows[grepl(paste0("^", sp, " ", ck, " "), names(rows))]), digits = 3)
  }
}
out <- do.call(rbind, rows); rownames(out) <- NULL
write.csv(out, here("results", "core_spa", "spa_in_dynamics_checks.csv"), row.names = FALSE)
print(out, digits = 3)

sessionInfo()
