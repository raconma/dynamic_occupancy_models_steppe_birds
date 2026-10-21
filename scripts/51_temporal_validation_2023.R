###############################################################################
# 51_temporal_validation_2023.R
#
# Does the dynamic model predict the next season better than a static
# distribution? Models are fitted to 2017-2022 and predict the independent
# 2023 detection histories (cells with at least one visit in 2023).
#
# Pre-specified (fixed before fitting):
#   dynamic        colext with the primary formulas (script 41) fitted to
#                  2017-2022; occupancy in 2023 = P(z_2022 = 1 | data) (1 - eps)
#                  + P(z_2022 = 0 | data) gamma, with gamma and eps of the
#                  2022-2023 transition
#   static_closed  (primary comparison) single-season occupancy model treating
#                  all 2017-2022 visits of a cell as replicates of one fixed
#                  state (psi and detection formulas of the primary model);
#                  occupancy in 2023 = P(z = 1 | data 2017-2022). The
#                  distribution is a snapshot that does not change.
#   static_sdm     (secondary) stacked single-season model (cell-years as
#                  sites); occupancy in 2023 = psi from habitat alone, as a
#                  correlative distribution model would predict
#   dynamic_uncond (secondary) dynamic model projected from initial occupancy
#                  with gamma and eps, without conditioning on the detections
# Each model predicts 2023 detection with its own detection submodel and the
# 2023 visit covariates. Scores over surveyed cells: held-out log-likelihood
# of the 2023 detection histories (sum), AUC and Brier score of the
# probability of at least one detection in 2023. Differences in
# log-likelihood from static_closed, and from static_sdm, get 95 % intervals by
# resampling 50 km blocks of cells (1,000 draws). Also reported by cells detected and not
# detected in 2017-2022.
# Caveat: the formulas were selected on 2017-2023 (script 41), so 2023 is not
# fully independent of model selection; the psi and detection formulas are
# shared by all models, the gamma and epsilon formulas only by the dynamic
# ones.
#
# Inputs:  results/{sp}_model_object.rds (script 41), cell_habitat.csv
# Outputs: results/robustness/temporal_validation_2023.csv
###############################################################################

library(here)
library(unmarked)

set.seed(42)

species <- c("otitar", "tettet", "pteori", "ptealc")
YEARS <- 2017:2023; T_years <- length(YEARS); T_tr <- T_years - 1
N_BOOT <- 1000; BLOCK_M <- 50000
hab <- read.csv(here("results", "core_spa", "cell_habitat.csv"))
out_dir <- here("results", "robustness"); dir.create(out_dir, showWarnings = FALSE)

rhs <- function(f) attr(terms(f), "term.labels")
auc <- function(p, o) {                                  # Mann-Whitney AUC
  r <- rank(p); n1 <- sum(o); n0 <- sum(!o)
  (sum(r[o]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
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
  sc <- siteCovs(u); ysc <- yearlySiteCovs(u); oc <- obsCovs(u)
  yr_row <- rep(seq_len(T_years), M)                         # yearly rows, site-major
  vis_yr <- rep(rep(seq_len(T_years), each = J), M)          # visit rows, site-major, year-major

  y_tr <- y[, seq_len(T_tr * J)]; y_te <- y[, T_tr * J + seq_len(J)]
  oc_tr <- oc[vis_yr <= T_tr, , drop = FALSE]; oc_te <- oc[vis_yr == T_years, , drop = FALSE]
  surveyed <- rowSums(!is.na(y_te)) > 0
  det_tr <- rowSums(y_tr == 1, na.rm = TRUE) > 0

  # dynamic model, 2017-2022
  umf_dyn <- unmarkedMultFrame(y = y_tr, siteCovs = sc, yearlySiteCovs = ysc[yr_row <= T_tr, , drop = FALSE],
                               obsCovs = oc_tr, numPrimary = T_tr)
  cf0 <- coef(mod)
  st <- c(starts_from(cf0, "psi", rhs(fl$psi)), starts_from(cf0, "col", rhs(fl$col)),
          starts_from(cf0, "ext", rhs(fl$ext)), starts_from(cf0, "p", rhs(fl$det)))
  f_dyn <- colext(psiformula = fl$psi, gammaformula = fl$col, epsilonformula = fl$ext,
                  pformula = fl$det, data = umf_dyn, starts = st, control = list(maxit = 3000))
  # static closed model: one state per cell over 2017-2022
  umf_cl <- unmarkedFrameOccu(y = y_tr, siteCovs = sc, obsCovs = oc_tr)
  f_cl <- occu(as.formula(paste(deparse1(fl$det), "~", deparse1(fl$psi[[2]]))), data = umf_cl,
               control = list(maxit = 3000))
  # static SDM: cell-years stacked as sites
  y_st <- do.call(rbind, lapply(seq_len(T_tr), function(t) y_tr[, (t - 1) * J + seq_len(J)]))
  oc_st <- do.call(rbind, lapply(seq_len(T_tr), function(t) oc[vis_yr == t, , drop = FALSE]))
  umf_st <- unmarkedFrameOccu(y = y_st, siteCovs = sc[rep(seq_len(M), T_tr), , drop = FALSE], obsCovs = oc_st)
  f_st <- occu(as.formula(paste(deparse1(fl$det), "~", deparse1(fl$psi[[2]]))), data = umf_st,
               control = list(maxit = 3000))

  # occupancy in 2023
  nd22 <- cbind(sc, ysc[yr_row == T_tr, , drop = FALSE])
  gam <- predict(f_dyn, type = "col", newdata = nd22)$Predicted
  eps <- predict(f_dyn, type = "ext", newdata = nd22)$Predicted
  # unconditional projection: psi in 2022 (psi_u22) and 2023 (psi_u)
  psi_u <- predict(f_dyn, type = "psi", newdata = sc)$Predicted
  for (t in seq_len(T_tr)) {
    if (t == T_tr) psi_u22 <- psi_u
    nd <- cbind(sc, ysc[yr_row == t, , drop = FALSE])
    psi_u <- psi_u * (1 - predict(f_dyn, type = "ext", newdata = nd)$Predicted) +
      (1 - psi_u) * predict(f_dyn, type = "col", newdata = nd)$Predicted
  }
  # posteriors: cells without any 2017-2022 visit are dropped by unmarked and
  # get the unconditional prediction (no data to condition on)
  has_tr <- rowSums(!is.na(y_tr)) > 0
  post_dyn <- ranef(f_dyn)@post[, 2, T_tr]; post_cl <- ranef(f_cl)@post[, 2, 1]
  stopifnot(length(post_dyn) %in% c(M, sum(has_tr)), length(post_cl) %in% c(M, sum(has_tr)))
  z22 <- psi_u22; z22[has_tr] <- if (length(post_dyn) == M) post_dyn[has_tr] else post_dyn
  z_cl <- predict(f_cl, type = "state", newdata = sc)$Predicted
  z_cl[has_tr] <- if (length(post_cl) == M) post_cl[has_tr] else post_cl
  psi_dyn <- z22 * (1 - eps) + (1 - z22) * gam
  message("  cells surveyed in 2023 without 2017-2022 visits: ", sum(surveyed & !has_tr))
  psi <- list(dynamic = psi_dyn, static_closed = z_cl,
              static_sdm = predict(f_st, type = "state", newdata = sc)$Predicted,
              dynamic_uncond = psi_u)
  fits <- list(dynamic = f_dyn, static_closed = f_cl, static_sdm = f_st, dynamic_uncond = f_dyn)

  # per-cell scores on 2023
  obs_det <- rowSums(y_te == 1, na.rm = TRUE) > 0
  score <- lapply(names(psi), function(m) {
    p <- matrix(predict(fits[[m]], type = "det", newdata = oc_te)$Predicted, M, J, byrow = TRUE)
    p[is.na(y_te)] <- NA
    lik_occ <- apply(p^y_te * (1 - p)^(1 - y_te), 1, prod, na.rm = TRUE)
    ll <- log(psi[[m]] * lik_occ + (1 - psi[[m]]) * !obs_det)
    pdet <- psi[[m]] * (1 - apply(1 - p, 1, prod, na.rm = TRUE))
    data.frame(ll = ll, pdet = pdet)
  }); names(score) <- names(psi)

  blk <- paste(floor(h$x / BLOCK_M), floor(h$y / BLOCK_M))
  for (grp in c("all", "detected_2017_2022", "not_detected_2017_2022")) {
    keep <- surveyed & switch(grp, all = TRUE, detected_2017_2022 = det_tr, not_detected_2017_2022 = !det_tr)
    if (sum(obs_det[keep]) < 2 || sum(!obs_det[keep]) < 2) next
    stopifnot(all(sapply(score, function(z) all(is.finite(z$ll[keep]) & is.finite(z$pdet[keep])))))
    by_blk <- split(which(keep), blk[keep])
    boot_idx <- replicate(N_BOOT, unlist(by_blk[sample(length(by_blk), replace = TRUE)]), simplify = FALSE)
    for (m in names(psi)) {
      s <- score[[m]]; s0 <- score$static_closed; s1 <- score$static_sdm
      d_ll <- s$ll - s0$ll; d_sdm <- s$ll - s1$ll
      boot <- sapply(boot_idx, function(i) sum(d_ll[i]))
      boot_sdm <- sapply(boot_idx, function(i) sum(d_sdm[i]))
      rows[[paste(sp, grp, m)]] <- data.frame(
        species = sp, cells = grp, model = m, n_surveyed = sum(keep), n_detected_2023 = sum(obs_det[keep]),
        loglik = sum(s$ll[keep]), d_loglik_vs_closed = sum(d_ll[keep]),
        d_lo = quantile(boot, 0.025), d_hi = quantile(boot, 0.975),
        d_loglik_vs_sdm = sum(d_sdm[keep]), d_sdm_lo = quantile(boot_sdm, 0.025), d_sdm_hi = quantile(boot_sdm, 0.975),
        auc = auc(s$pdet[keep], obs_det[keep]), brier = mean((s$pdet[keep] - obs_det[keep])^2),
        mean_psi_2023 = mean(psi[[m]][keep]),
        converged = fits[[m]]@opt$convergence == 0, row.names = NULL)
    }
  }
  print(do.call(rbind, rows[grepl(paste0("^", sp, " "), names(rows))]), digits = 3)
}
out <- do.call(rbind, rows); rownames(out) <- NULL
write.csv(out, file.path(out_dir, "temporal_validation_2023.csv"), row.names = FALSE)
print(out, digits = 3)

sessionInfo()
