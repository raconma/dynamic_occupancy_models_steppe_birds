###############################################################################
# 15_parboot_publication.R
#
# Runs publication-grade parametric bootstrap GOF (nsim = 500) and
# MacKenzie-Bailey GOF (nsim = 500) for each species' colext model.
#
# Loads fitted model objects from results/{sp}_model_object.rds.
# Saves results to results/{sp}_gof_parboot.rds and
# results/parboot_summary.csv.
#
# Estimated runtime: ~12 hours on 8-core / 16 GB machine (all 4 species).
#   ~22 seconds per parboot iteration per species.
#   500 iterations × 4 species = ~44,000 seconds.
#
# Safe to run overnight — memory usage is low (~200 MB per species).
###############################################################################

library(here)
library(unmarked)

NSIM_PARBOOT <- 500
NSIM_MB_GOF  <- 500

species_codes <- c("otitar", "ptealc", "pteori", "tettet")

results <- data.frame(
  species     = character(),
  chisq_obs   = numeric(),
  chisq_mean  = numeric(),
  chisq_pval  = numeric(),
  mb_chisq    = numeric(),
  mb_pval     = numeric(),
  c_hat       = numeric(),
  stringsAsFactors = FALSE
)

for (sp in species_codes) {
  message("\n===== ", toupper(sp), " =====")
  message("  Loading model from results/", sp, "_model_object.rds ...")

  mod_path <- here("results", paste0(sp, "_model_object.rds"))
  if (!file.exists(mod_path)) {
    message("  WARNING: model file not found — skipping ", sp)
    next
  }

  mod <- readRDS(mod_path)
  message("  Model loaded. AIC = ", round(mod@AIC, 2))

  # ---- Parametric bootstrap GOF ----
  message("  Running parboot GOF (nsim = ", NSIM_PARBOOT, ") ...")
  message("  Start: ", Sys.time())

  t0 <- proc.time()
  gof <- tryCatch(
    parboot(mod, nsim = NSIM_PARBOOT),
    error = function(e) {
      message("  ERROR in parboot: ", e$message)
      NULL
    }
  )
  elapsed <- (proc.time() - t0)[3]
  message("  Elapsed: ", round(elapsed / 60, 1), " min")

  if (!is.null(gof)) {
    saveRDS(gof, here("results", paste0(sp, "_gof_parboot.rds")))
    message("  Saved: results/", sp, "_gof_parboot.rds")

    # Extract p-value
    obs_stat <- gof@t0
    sim_stats <- gof@t.star
    p_val <- mean(sim_stats >= obs_stat)
    message("  Observed chi-sq: ", round(obs_stat, 2),
            "  |  Mean simulated: ", round(mean(sim_stats), 2),
            "  |  P-value: ", round(p_val, 4))
  }

  # ---- MacKenzie-Bailey GOF ----
  message("  Running MacKenzie-Bailey GOF (nsim = ", NSIM_MB_GOF, ") ...")
  t1 <- proc.time()
  mb_gof <- tryCatch(
    AICcmodavg::mb.gof.test(mod, nsim = NSIM_MB_GOF, plot.hist = FALSE),
    error = function(e) {
      message("  ERROR in mb.gof.test: ", e$message)
      NULL
    }
  )
  elapsed2 <- (proc.time() - t1)[3]
  message("  Elapsed: ", round(elapsed2 / 60, 1), " min")

  # Collect results
  row <- data.frame(
    species     = sp,
    chisq_obs   = if (!is.null(gof)) gof@t0 else NA,
    chisq_mean  = if (!is.null(gof)) mean(gof@t.star) else NA,
    chisq_pval  = if (!is.null(gof)) p_val else NA,
    mb_chisq    = if (!is.null(mb_gof)) mb_gof$chisq.table$chi.square[1] else NA,
    mb_pval     = if (!is.null(mb_gof)) mb_gof$p.value else NA,
    c_hat       = if (!is.null(mb_gof)) mb_gof$c.hat.est else NA,
    stringsAsFactors = FALSE
  )
  results <- rbind(results, row)

  if (!is.null(mb_gof)) {
    message("  MB chi-sq: ", round(mb_gof$chisq.table$chi.square[1], 2),
            "  |  c-hat: ", round(mb_gof$c.hat.est, 3),
            "  |  P-value: ", round(mb_gof$p.value, 4))
  }

  message("  Done with ", sp, " at ", Sys.time())
}

# Save summary table
write.csv(results, here("results", "parboot_summary.csv"), row.names = FALSE)
message("\n===== ALL DONE =====")
message("  Results saved to results/parboot_summary.csv")
message("  Individual parboot objects: results/{sp}_gof_parboot.rds")
message("  Finished at: ", Sys.time())
