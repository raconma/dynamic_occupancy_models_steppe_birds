###############################################################################
# 54_rolling_validation.R
#
# Robustness plan A (docs/robustness_plan_BiolCons.md): rolling-origin
# validation of the dynamic model against static maps, including maps updated
# every year. Each test year (2021, 2022, 2023) is predicted from all earlier
# years; the primary formulas (script 41) are used throughout, without
# reselection within training windows.
#
# Models fitted to the training years:
#   dynamic   colext; occupancy in the test year = P(z = 1 | data) in the last
#             training year propagated with the fitted colonisation and
#             extinction (unconditional projection for cells without training
#             visits)
#   S1        fixed-state single-season model over all training years
#             (persistence map; the reference of script 51)
#   S2        single-season model fitted to the last training year only (a map
#             updated every year)
#   S3        fixed-state single-season model over the last three training
#             years (a recent atlas)
#   S4        habitat-only model, cell-years stacked as sites
# Static models give P(z = 1 | data) for cells with visits in their window and
# their habitat prediction otherwise. Each model predicts detection in the
# test year with its own detection submodel.
# Scores: held-out log-likelihood of the test-year detection histories per
# cell, with predicted occupancy bounded to [1e-6, 1 - 1e-6]; differences dynamic minus each static model, per test year and pooled
# over test years, with 95 % intervals from 1,000 resamples of 50 x 50 km
# blocks (cells keep all their test years). Cells: all surveyed in the test
# year, and those with a record in the training years.
#
# Inputs:  results/{sp}_model_object.rds, cell_habitat.csv
# Outputs: results/robustness/rolling_validation.csv
#          results/robustness/rolling_validation_scores.csv (AUC, Brier)
###############################################################################

library(here)
library(unmarked)

set.seed(42)

species <- c("otitar", "tettet", "pteori", "ptealc")
if (nzchar(Sys.getenv("SPECIES"))) species <- strsplit(Sys.getenv("SPECIES"), ",")[[1]]
YEARS <- 2017:2023; T_years <- length(YEARS)
TEST <- 5:7                                         # 2021, 2022, 2023
N_BOOT <- 1000; BLOCK_M <- 50000
hab <- read.csv(here("results", "core_spa", "cell_habitat.csv"))
out_dir <- here("results", "robustness"); dir.create(out_dir, showWarnings = FALSE)
STATIC <- c(S1 = "Persistence map (all training years)", S2 = "Map updated every year (last year)",
            S3 = "Recent atlas (last three years)", S4 = "Habitat-only model")

rhs <- function(f) attr(terms(f), "term.labels")
auc <- function(p, o) { r <- rank(p); n1 <- sum(o); n0 <- sum(!o); (sum(r[o]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
starts_from <- function(cf, sub, labels) {
  unname(sapply(c("(Intercept)", labels), function(l) {
    nm <- paste0(sub, "(", if (l == "(Intercept)") "Int" else l, ")")
    if (nm %in% names(cf)) cf[[nm]] else 0
  }))
}
# posterior where the window has data, model prediction elsewhere
post_or_prior <- function(post, prior, has) {
  stopifnot(length(post) %in% c(length(prior), sum(has)))
  out <- prior; out[has] <- if (length(post) == length(prior)) post[has] else post; out
}

cell_rows <- list(); score_rows <- list()
for (sp in species) {
  message("=== ", sp)
  mod <- readRDS(here("results", paste0(sp, "_model_object.rds")))
  fl <- mod@formlist; u <- mod@data; cf0 <- coef(mod)
  tab <- read.csv(here("data", "processed_2023", sp, paste0(sp, "_occ_wide_dynamic.csv")))
  h <- hab[match(tab$cells, hab$cells), ]
  y <- getY(u); M <- nrow(y); J <- ncol(y) / T_years
  sc <- siteCovs(u); ysc <- yearlySiteCovs(u); oc <- obsCovs(u)
  yr_row <- rep(seq_len(T_years), M); vis_yr <- rep(rep(seq_len(T_years), each = J), M)
  cols_of <- function(t) unlist(lapply(t, function(k) (k - 1) * J + seq_len(J)))
  occu_form <- as.formula(paste(deparse1(fl$det), "~", deparse1(fl$psi[[2]])))
  fixed_state <- function(win) {                   # one state per cell over the window
    yw <- y[, cols_of(win), drop = FALSE]
    f <- occu(occu_form, data = unmarkedFrameOccu(y = yw, siteCovs = sc, obsCovs = oc[vis_yr %in% win, , drop = FALSE]),
              control = list(maxit = 3000))
    has <- rowSums(!is.na(yw)) > 0
    list(fit = f, psi = post_or_prior(ranef(f)@post[, 2, 1], predict(f, type = "state", newdata = sc)$Predicted, has))
  }
  blk <- paste(floor(h$x / BLOCK_M), floor(h$y / BLOCK_M))

  for (T_te in TEST) {
    T_tr <- T_te - 1; tr <- seq_len(T_tr)
    message("  test ", YEARS[T_te])
    y_tr <- y[, cols_of(tr)]; y_te <- y[, cols_of(T_te)]
    oc_te <- oc[vis_yr == T_te, , drop = FALSE]
    surveyed <- rowSums(!is.na(y_te)) > 0
    det_tr <- rowSums(y_tr == 1, na.rm = TRUE) > 0
    has_tr <- rowSums(!is.na(y_tr)) > 0

    # dynamic
    umf_dyn <- unmarkedMultFrame(y = y_tr, siteCovs = sc, yearlySiteCovs = ysc[yr_row <= T_tr, , drop = FALSE],
                                 obsCovs = oc[vis_yr <= T_tr, , drop = FALSE], numPrimary = T_tr)
    st <- c(starts_from(cf0, "psi", rhs(fl$psi)), starts_from(cf0, "col", rhs(fl$col)),
            starts_from(cf0, "ext", rhs(fl$ext)), starts_from(cf0, "p", rhs(fl$det)))
    f_dyn <- colext(psiformula = fl$psi, gammaformula = fl$col, epsilonformula = fl$ext,
                    pformula = fl$det, data = umf_dyn, starts = st, control = list(maxit = 3000))
    psi_u <- predict(f_dyn, type = "psi", newdata = sc)$Predicted
    for (t in seq_len(T_tr - 1)) {
      nd <- cbind(sc, ysc[yr_row == t, , drop = FALSE])
      psi_u <- psi_u * (1 - predict(f_dyn, type = "ext", newdata = nd)$Predicted) +
        (1 - psi_u) * predict(f_dyn, type = "col", newdata = nd)$Predicted
    }
    z_last <- post_or_prior(ranef(f_dyn)@post[, 2, T_tr], psi_u, has_tr)
    nd_last <- cbind(sc, ysc[yr_row == T_tr, , drop = FALSE])
    psi_dyn <- z_last * (1 - predict(f_dyn, type = "ext", newdata = nd_last)$Predicted) +
      (1 - z_last) * predict(f_dyn, type = "col", newdata = nd_last)$Predicted

    # static references
    s1 <- fixed_state(tr); s2 <- fixed_state(T_tr); s3 <- fixed_state(max(1, T_tr - 2):T_tr)
    y_st <- do.call(rbind, lapply(tr, function(t) y[, cols_of(t), drop = FALSE]))
    oc_st <- do.call(rbind, lapply(tr, function(t) oc[vis_yr == t, , drop = FALSE]))
    f4 <- occu(occu_form, data = unmarkedFrameOccu(y = y_st, siteCovs = sc[rep(seq_len(M), T_tr), , drop = FALSE], obsCovs = oc_st),
               control = list(maxit = 3000))
    psi <- list(dynamic = psi_dyn, S1 = s1$psi, S2 = s2$psi, S3 = s3$psi,
                S4 = predict(f4, type = "state", newdata = sc)$Predicted)
    # predicted occupancy bounded away from 0 and 1 (1e-6), so that a model
    # that rules out an occupancy state the data then show gets a large but
    # finite penalty instead of minus infinity
    psi_raw <- psi
    psi <- lapply(psi, function(v) pmin(pmax(v, 1e-6), 1 - 1e-6))
    fits <- list(dynamic = f_dyn, S1 = s1$fit, S2 = s2$fit, S3 = s3$fit, S4 = f4)

    obs_det <- rowSums(y_te == 1, na.rm = TRUE) > 0
    at_bound <- sapply(psi_raw, function(v) sum(surveyed & ((obs_det & v < 1e-6) | (!obs_det & v > 1 - 1e-6))))
    if (any(at_bound > 0)) message("  cells where the bound was used: ", paste(names(at_bound), at_bound, collapse = ", "))
    ll <- sapply(names(psi), function(m) {
      p <- matrix(predict(fits[[m]], type = "det", newdata = oc_te)$Predicted, M, J, byrow = TRUE)
      p[is.na(y_te)] <- NA
      lik_occ <- apply(p^y_te * (1 - p)^(1 - y_te), 1, prod, na.rm = TRUE)
      pdet <- psi[[m]] * (1 - apply(1 - p, 1, prod, na.rm = TRUE))
      for (grp in c("all", "recorded")) {
        keep <- surveyed & (grp == "all" | det_tr)
        score_rows[[paste(sp, T_te, m, grp)]] <<- data.frame(
          species = sp, test_year = YEARS[T_te], model = m, cells = grp, n = sum(keep),
          auc = auc(pdet[keep], obs_det[keep]), brier = mean((pdet[keep] - obs_det[keep])^2),
          converged = fits[[m]]@opt$convergence == 0)
      }
      log(psi[[m]] * lik_occ + (1 - psi[[m]]) * !obs_det)
    })
    keep <- surveyed
    bad <- colSums(!is.finite(ll[keep, , drop = FALSE]))
    if (any(bad > 0)) message("  non-finite log-likelihood (cells per model): ", paste(names(bad), bad, collapse = ", "))
    stopifnot(all(is.finite(ll[keep, ])))
    cell_rows[[paste(sp, T_te)]] <- data.frame(species = sp, test_year = YEARS[T_te], cell = which(keep),
                                               blk = blk[keep], recorded = det_tr[keep], ll[keep, , drop = FALSE])
  }
}
cells <- do.call(rbind, cell_rows); rownames(cells) <- NULL

# contrasts per test year and pooled, block bootstrap over cells
rows <- list()
for (sp in species) for (grp in c("all", "recorded")) {
  d0 <- cells[cells$species == sp & (grp == "all" | cells$recorded), ]
  by_blk <- split(seq_len(nrow(d0)), d0$blk)
  boot_idx <- replicate(N_BOOT, unlist(by_blk[sample(length(by_blk), replace = TRUE)]), simplify = FALSE)
  for (ty in c(YEARS[TEST], NA)) {
    sel <- if (is.na(ty)) rep(TRUE, nrow(d0)) else d0$test_year == ty
    for (s in names(STATIC)) {
      dd <- (d0$dynamic - d0[[s]]) * sel
      b <- sapply(boot_idx, function(i) sum(dd[i]))
      rows[[paste(sp, grp, ty, s)]] <- data.frame(
        species = sp, cells = grp, test_year = if (is.na(ty)) "pooled" else as.character(ty),
        reference = s, reference_label = STATIC[[s]], n_cell_years = sum(sel),
        d_loglik = sum(dd), lo = quantile(b, 0.025), hi = quantile(b, 0.975),
        p_dynamic_better = mean(b > 0), row.names = NULL)
    }
  }
}
out <- do.call(rbind, rows); rownames(out) <- NULL
write.csv(out, file.path(out_dir, "rolling_validation.csv"), row.names = FALSE)
write.csv(do.call(rbind, score_rows), file.path(out_dir, "rolling_validation_scores.csv"), row.names = FALSE)
options(width = 200)
print(out[out$test_year == "pooled", ], digits = 3)

sessionInfo()
