###############################################################################
# 47_stPGOcc_spa_trend.R
#
# Spatial check of the SPA results: spatio-temporal occupancy model
# (spOccupancy::stPGOcc; NNGP exponential spatial random effect, AR(1)
# temporal random effect) of the pooled detection histories, with
#   occupancy  ~ initial-occupancy covariates of the primary model +
#                pseudo-steppe + SPA cover + year + SPA cover x year
#   detection  ~ primary detection covariates + SPA cover (observers target
#                SPAs; script 44)
# The SPA x year coefficient compares the occupancy trend inside and outside
# SPAs once residual spatial autocorrelation is absorbed by the spatial
# random effect. The model does not separate colonisation from extinction;
# it complements the dynamic models of scripts 41-46.
#
# Usage:   Rscript scripts/47_stPGOcc_spa_trend.R [species]   (default otitar)
# Inputs:  results/{sp}_model_object.rds (script 41), the pooled input table,
#          cell_habitat.csv (script 23)
# Outputs: results/results_spatial/stPGOcc_spa_trend_{sp}.rds
#          results/results_spatial/stPGOcc_spa_trend_{sp}_summary.csv
###############################################################################

library(here)
library(unmarked)
library(spOccupancy)

set.seed(42)

args <- commandArgs(trailingOnly = TRUE)
sp <- if (length(args)) args[1] else "otitar"
YEARS <- 2017:2023; T_years <- length(YEARS)
N_NEIGHBORS <- 5
N_CHAINS <- 3; BATCH_LENGTH <- 25
N_BATCH <- as.integer(Sys.getenv("N_BATCH", "1600"))      # 40,000 iterations per chain
N_BURN <- as.integer(Sys.getenv("N_BURN", "20000"))       # env overrides for a quick test
N_THIN <- 20
# (a first run of 10,000 iterations left the spatial variance at Rhat 2.2)
out_dir <- here("results", "results_spatial"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

rhs <- function(f) attr(terms(f), "term.labels")
sc01 <- function(v) (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE)

mod <- readRDS(here("results", paste0(sp, "_model_object.rds")))
fl <- mod@formlist; u <- mod@data
tab <- read.csv(here("data", "processed_2023", sp, paste0(sp, "_occ_wide_dynamic.csv")))
hab <- read.csv(here("results", "core_spa", "cell_habitat.csv"))
h <- hab[match(tab$cells, hab$cells), ]
y2 <- getY(u); M <- nrow(y2); J <- ncol(y2) / T_years

# detection histories and detection covariates as site x year x visit arrays
to_arr <- function(m) aperm(array(m, c(M, J, T_years)), c(1, 3, 2))
y <- to_arr(y2)
oc <- obsCovs(u)
det_vars <- rhs(fl$det)
base_det <- unique(sub("^I\\((.*)\\^2\\)$", "\\1", det_vars))
det.covs <- lapply(setNames(base_det, base_det), function(v) to_arr(matrix(oc[[v]], M, J * T_years, byrow = TRUE)))
det.covs$spa_det <- to_arr(matrix(rep(h$spa_frac, each = J * T_years), M, J * T_years, byrow = TRUE))
det_formula <- reformulate(c(det_vars, "spa_det"))

# occupancy covariates: site-level vectors and the year matrix
sc <- siteCovs(u)
occ.covs <- as.list(sc[, rhs(fl$psi), drop = FALSE])
occ.covs$pseudo_steppe <- sc01(h$steppe_strict)
occ.covs$spa <- h$spa_frac
occ.covs$years <- matrix(YEARS - 2020, M, T_years, byrow = TRUE)
occ_formula <- reformulate(c(rhs(fl$psi), "pseudo_steppe", "spa", "years", "spa:years"))
coords <- cbind(h$x, h$y) / 1000                       # km, UTM 30N
stopifnot(!anyDuplicated(coords))

message(sp, ": ", M, " sites, ", T_years, " years, ", J, " visits; occupancy ",
        deparse(occ_formula), "; detection ", deparse(det_formula))
dd <- dist(coords[sample(M, min(M, 1000)), ])
d_min <- 3 / max(dd); d_max <- 3 / quantile(dd, 0.05)   # effective range 5 % to 100 % of max distance
t0 <- Sys.time()
fit <- stPGOcc(occ.formula = occ_formula, det.formula = det_formula,
               data = list(y = y, occ.covs = occ.covs, det.covs = det.covs, coords = coords),
               priors = list(phi.unif = c(d_min, d_max)),
               cov.model = "exponential", NNGP = TRUE, n.neighbors = N_NEIGHBORS,
               n.batch = N_BATCH, batch.length = BATCH_LENGTH, ar1 = TRUE,
               n.burn = N_BURN, n.thin = N_THIN, n.chains = N_CHAINS,
               n.omp.threads = 4, verbose = TRUE, n.report = 100)
message("fitted in ", round(as.numeric(Sys.time() - t0, units = "hours"), 2), " h")
saveRDS(fit, file.path(out_dir, paste0("stPGOcc_spa_trend_", sp, ".rds")))

q <- function(x) quantile(x, c(0.5, 0.025, 0.975))
beta <- fit$beta.samples
res <- data.frame(species = sp, parameter = colnames(beta),
                  t(apply(beta, 2, q)), rhat = fit$rhat$beta, ess = fit$ESS$beta, row.names = NULL)
names(res)[3:5] <- c("median", "lo", "hi")
theta <- fit$theta.samples
res <- rbind(res, data.frame(species = sp, parameter = colnames(theta), t(apply(theta, 2, q)),
                             rhat = fit$rhat$theta, ess = fit$ESS$theta, row.names = NULL,
                             check.names = FALSE) |> setNames(names(res)))
write.csv(res, file.path(out_dir, paste0("stPGOcc_spa_trend_", sp, "_summary.csv")), row.names = FALSE)
print(res, digits = 3)

sessionInfo()
