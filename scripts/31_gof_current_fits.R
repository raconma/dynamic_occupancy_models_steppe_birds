###############################################################################
# 31_gof_current_fits.R
#
# Goodness of fit of the primary colext models used in the Biological
# Conservation draft. Parametric bootstrap of the sum of squared errors
# (unmarked::parboot, default statistic SSE, 200 simulations), as in
# script 15.
#
# Script 15 stored bootstraps in results/{sp}_gof_parboot.rds, but two
# models were refitted afterwards. A stored bootstrap is reused only if its
# observed statistic equals the SSE of the current model object; otherwise
# the bootstrap is rerun and saved to results/{sp}_gof_parboot_current.rds.
#
# Output: results/manuscript_tables/tableS_gof.csv
###############################################################################

library(here)
library(unmarked)

set.seed(42)

species <- c("otitar", "ptealc", "pteori", "tettet")
N_SIM <- 200   # 500 in the earlier fits; 200 is ample for a P value and keeps refits with up to 30 visits per cell-year tractable
N_CORES <- max(1, min(4, parallel::detectCores() - 2))
tab_dir <- here("results", "manuscript_tables")
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

sse <- function(fm) {
  sum((getY(fm@data) - fitted(fm))^2, na.rm = TRUE)
}

rows <- list()
for (sp in species) {
  message("\n=== ", sp, " ===")
  mod <- readRDS(here("results", paste0(sp, "_model_object.rds")))
  sse_now <- sse(mod)
  stored_path <- here("results", paste0(sp, "_gof_parboot.rds"))
  stored <- if (file.exists(stored_path)) readRDS(stored_path) else NULL
  reuse <- !is.null(stored) && isTRUE(all.equal(unname(stored@t0[1]), sse_now,
                                                tolerance = 1e-6))
  if (reuse) {
    message("  Stored bootstrap matches the current fit; reused")
    gof <- stored
    source_file <- basename(stored_path)
  } else {
    message("  Stored bootstrap does not match the current fit; rerunning (",
            N_SIM, " simulations, ", N_CORES, " cores)")
    gof <- parboot(mod, nsim = N_SIM, parallel = N_CORES > 1, ncores = N_CORES)
    out <- here("results", paste0(sp, "_gof_parboot_current.rds"))
    saveRDS(gof, out)
    source_file <- basename(out)
  }
  t_star <- gof@t.star[, 1]
  rows[[sp]] <- data.frame(
    species = sp, sse_observed = sse_now, sse_sim_mean = mean(t_star),
    sse_sim_sd = sd(t_star), n_sim = length(t_star),
    p_value = mean(t_star >= gof@t0[1]), c_hat = gof@t0[1] / mean(t_star),
    reused_stored = reuse, source = source_file)
}

gof_df <- do.call(rbind, rows); rownames(gof_df) <- NULL
write.csv(gof_df, file.path(tab_dir, "tableS_gof.csv"), row.names = FALSE)
print(gof_df, digits = 3)
sessionInfo()
