###############################################################################
# 35_spa_turnover_contrast.R
#
# Static protection test: do SPA cells differ from matched unprotected cells
# in colonisation as well as in extinction? Script 23 compares extinction
# rates; here the colonisation rate per unoccupied cell-year is compared in
# the same way and with the same coefficient draws and cell resamples as the
# extinction rate, so both contrasts describe the same matched cells.
#
#   extinction rate   = sum(E) / sum(expected occupied cell-years)
#   colonisation rate = sum(C) / sum(expected unoccupied cell-years)
#
# on anchored recent transitions (2018->19 ... 2022->23, surveyed in both
# years) in cells with >= 1 detection. Under static protection SPA cells hold
# persistently occupied sites, with extinction and colonisation no higher
# (and colonisation possibly lower) than in matched cells; if SPAs acted as
# sources, colonisation of vacant SPA cells would be higher.
#
# Model variants: primary, eps_null, steprep_v3 (script 28). In both
# sandgrouse the alternative variants place colonisation at its boundary, so
# their colonisation contrasts are not informative (flagged by point_only).
# P. alchata is kept for colonisation only (extinction not identifiable).
#
# Inputs:  results/{sp}_model_object.rds and results/robustness/ (script 28),
#          colext inputs, results/core_spa/cell_habitat.csv (script 23)
# Output:  results/core_spa/spa_turnover_contrast.csv
###############################################################################

library(here)
library(unmarked)
library(MASS)

source(here("R", "colext_utils.R"))
source(here("R", "spa_matching.R"))
set.seed(42)

species <- c("otitar", "tettet", "pteori", "ptealc")
variants <- c("primary", "eps_null", "steprep_v3")
N_DRAWS <- 1000
hab <- read.csv(here("results", "core_spa", "cell_habitat.csv"))
q <- function(x, p) as.numeric(quantile(x, p, na.rm = TRUE))
summ <- function(x) { x <- x[is.finite(x)]; c(med = median(x), lo = q(x, .025), hi = q(x, .975)) }

rows <- list()
for (v in variants) for (sp in species) {
  message("=== ", v, " / ", sp)
  s <- spa_setup(sp, v, hab)
  m <- spa_matches(s)
  draws <- draw_coefs(s$mod, N_DRAWS)
  point_only <- attr(draws, "point_only")
  stats <- c("ext_spa", "ext_ctl", "col_spa", "col_ctl", "ext_diff", "col_diff",
             "ext_ratio", "col_ratio")
  out <- array(NA_real_, c(N_DRAWS, length(m), length(stats)),
               dimnames = list(NULL, names(m), stats))
  rec <- s$recent
  for (k in seq_len(N_DRAWS)) {
    fb <- forward_backward(colext_rates(draws[k, ], s$D), s$D)
    z <- fb$z1[, 1:(s$n_t - 1)]
    E <- rowSums((fb$ext * s$both)[, rec]); O <- rowSums((z * s$both)[, rec])
    C <- rowSums((fb$col * s$both)[, rec]); U <- rowSums(((1 - z) * s$both)[, rec])
    b <- sample(seq_along(s$tr), replace = TRUE)
    for (nm in names(m)) {
      t_i <- s$tr[b]; c_i <- m[[nm]][b]
      e_t <- sum(E[t_i]) / sum(O[t_i]); e_c <- sum(E[c_i]) / sum(O[c_i])
      c_t <- sum(C[t_i]) / sum(U[t_i]); c_c <- sum(C[c_i]) / sum(U[c_i])
      out[k, nm, ] <- c(e_t, e_c, c_t, c_c, e_t - e_c, c_t - c_c, e_t / e_c, c_t / c_c)
    }
  }
  # Naive rates on the same anchored transitions (no resampling)
  n0 <- s$naive[, rec]; n1 <- s$naive[, rec + 1]; a <- s$both[, rec]
  naive_rate <- function(i, from) {
    opp <- sum(n0[i, ] == from & a[i, ], na.rm = TRUE)
    ev <- sum(n0[i, ] == from & n1[i, ] == 1 - from & a[i, ], na.rm = TRUE)
    c(ev = ev, opp = opp)
  }
  for (nm in names(m)) {
    ext_t <- naive_rate(s$tr, 1); ext_c <- naive_rate(m[[nm]], 1)
    col_t <- naive_rate(s$tr, 0); col_c <- naive_rate(m[[nm]], 0)
    r <- data.frame(species = sp, variant = v, point_only = point_only, spec = nm,
                    n_spa = length(s$tr), n_unique_controls = length(unique(m[[nm]])),
                    naive_ext_spa = sprintf("%d/%d", ext_t[["ev"]], ext_t[["opp"]]),
                    naive_ext_ctl = sprintf("%d/%d", ext_c[["ev"]], ext_c[["opp"]]),
                    naive_col_spa = sprintf("%d/%d", col_t[["ev"]], col_t[["opp"]]),
                    naive_col_ctl = sprintf("%d/%d", col_c[["ev"]], col_c[["opp"]]))
    for (st in stats) {
      sm <- summ(out[, nm, st])
      r[[paste0(st, "_med")]] <- sm[["med"]]; r[[paste0(st, "_lo")]] <- sm[["lo"]]
      r[[paste0(st, "_hi")]] <- sm[["hi"]]
    }
    r$p_col_spa_lower <- mean(out[, nm, "col_diff"] < 0, na.rm = TRUE)
    r$n_col_undefined <- sum(!is.finite(out[, nm, "col_diff"]))
    rows[[paste(v, sp, nm)]] <- r
  }
}
res <- do.call(rbind, rows); rownames(res) <- NULL
write.csv(res, here("results", "core_spa", "spa_turnover_contrast.csv"), row.names = FALSE)
print(res[, c("species", "variant", "spec", "point_only", "naive_col_spa", "naive_col_ctl",
              "col_spa_med", "col_ctl_med", "col_diff_med", "col_diff_lo", "col_diff_hi",
              "ext_diff_med", "ext_diff_lo", "ext_diff_hi")], digits = 3)

sessionInfo()
