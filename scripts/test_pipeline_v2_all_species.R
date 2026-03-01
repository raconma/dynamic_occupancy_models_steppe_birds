###############################################################################
# test_pipeline_v2_all_species.R
#
# Pipeline alternativo para las 4 especies de aves esteparias.
# Procesa cada especie y genera una tabla comparativa final.
#
# Cambios respecto al pipeline original:
#   1. site_vars = c("locality_id")  en vez de c("locality_id","observer_id")
#   2. min_obs = 2  en vez de 3
#   3. Filtro temporal: >= 3 anios con datos, early (2017-2019) + late (2020-2022)
#   4. Modelos progresivos: null -> p-only -> psi+p -> psi completo
#
# No usa variables dinamicas (GEE) ya que el cambio de sitios invalida el
# export original. Usa solo covariables estaticas + de observacion.
###############################################################################

library(here)
library(auk)
library(unmarked)
library(dplyr)
library(tidyr)
library(raster)
library(terra)
library(sf)

select <- dplyr::select
filter <- dplyr::filter
set.seed(42)

YEARS   <- 2017:2022
T_years <- length(YEARS)
J_reps  <- 10

# Species info
species_info <- list(
  otitar = list(code = "otitar", name = "Otis tarda (Great Bustard)"),
  ptealc = list(code = "ptealc", name = "Pterocles alchata (Pin-tailed Sandgrouse)"),
  pteori = list(code = "pteori", name = "Pterocles orientalis (Black-bellied Sandgrouse)"),
  tettet = list(code = "tettet", name = "Tetrax tetrax (Little Bustard)")
)

# Variables
vars_to_scale <- c("bio1", "tree_cover", "bio2", "grass_cover",
                    "topo_aspect", "topo_elev")
obs_covs <- c("time_observations_started", "duration_minutes",
              "effort_distance_km", "number_observers")
site_covs_wide <- c("locality_id", "latitude", "longitude", vars_to_scale)

base_data <- "/Users/gfandos/Documents/GitHub/dynamic_occupancy_models_steppe_birds/data-raw/data"
base_repo <- "/Users/gfandos/Documents/GitHub/dynamic_occupancy_models_steppe_birds"

# Load environmental rasters (shared)
cat("Cargando rasters ambientales...\n")
variables_environmental <- raster::stack(
  file.path(base_data, "environmental_data/environmental_data_occ/variables_spain.grd")
)
variables_aspect <- raster::stack(
  file.path(base_data, "topology_data/topo_aspect.asc")
)
variables_elev <- raster::stack(
  file.path(base_data, "topology_data/topo_elev.asc")
)

###############################################################################
# Helper function: run full pipeline for one species
###############################################################################
run_pipeline_v2 <- function(sp_code, sp_name) {

  cat("\n\n")
  cat(strrep("#", 70), "\n")
  cat("  ", sp_name, "\n")
  cat(strrep("#", 70), "\n\n")

  result <- list(species = sp_code, name = sp_name)

  # --- PASO 1: Cargar datos ---
  raw_path <- file.path(base_data, sp_code,
                         paste0("ebd_", sp_code, "_breeding_spain_zf.csv"))
  if (!file.exists(raw_path)) {
    cat("  ERROR: Archivo no encontrado: ", raw_path, "\n")
    result$error <- "File not found"
    return(result)
  }

  occ_raw <- read.csv(raw_path)
  occ_raw <- occ_raw %>% filter(year >= 2017, year <= 2022)
  result$n_raw <- nrow(occ_raw)
  cat("  Registros crudos: ", result$n_raw, "\n")

  # --- PASO 2: filter_repeat_visits ---
  occ <- filter_repeat_visits(
    occ_raw, min_obs = 2, max_obs = 10,
    annual_closure = TRUE, date_var = "observation_date",
    site_vars = c("locality_id")
  )
  result$n_filtered <- nrow(occ)
  result$n_sites_year <- n_distinct(occ$site)
  result$n_localities <- n_distinct(occ$locality_id)
  cat("  Registros filtrados: ", result$n_filtered, "\n")
  cat("  Sitios (loc x anio): ", result$n_sites_year, "\n")
  cat("  Localidades: ", result$n_localities, "\n")

  # --- PASO 3: Covariables ambientales ---
  coords <- occ[, c("longitude", "latitude")]

  env_ext <- as.data.frame(raster::extract(variables_environmental, coords, cellnumbers = TRUE))

  aspect_raw <- as.data.frame(terra::extract(variables_aspect, coords))
  elev_raw   <- as.data.frame(terra::extract(variables_elev, coords))
  if ("ID" %in% names(aspect_raw)) aspect_raw <- aspect_raw[, -which(names(aspect_raw) == "ID"), drop = FALSE]
  if ("ID" %in% names(elev_raw))   elev_raw   <- elev_raw[, -which(names(elev_raw) == "ID"), drop = FALSE]
  names(aspect_raw) <- "topo_aspect"
  names(elev_raw)   <- "topo_elev"

  occ_var <- occ %>%
    cbind(env_ext[, setdiff(names(env_ext), "ID")]) %>%
    cbind(aspect_raw) %>%
    cbind(elev_raw)

  if ("cell" %in% names(occ_var)) names(occ_var)[names(occ_var) == "cell"] <- "cells"
  occ_var <- occ_var %>% drop_na(bio12, tree_cover, topo_elev)
  result$n_with_covs <- nrow(occ_var)
  cat("  Con covariables: ", result$n_with_covs, "\n")

  # --- PASO 4: Estandarizar + formato wide ---
  scaling_params <- list()
  occ_var_std <- occ_var
  for (v in vars_to_scale) {
    vals <- occ_var_std[[v]]
    scaling_params[[v]] <- list(center = mean(vals, na.rm = TRUE),
                                 scale  = sd(vals, na.rm = TRUE))
    occ_var_std[[v]] <- (vals - scaling_params[[v]]$center) / scaling_params[[v]]$scale
  }

  occ_wide_list <- list()
  for (yr in YEARS) {
    occ_yr <- occ_var_std %>% filter(year == yr)
    if (nrow(occ_yr) == 0) next
    occ_wide_yr <- format_unmarked_occu(
      occ_yr, site_id = "site", response = "species_observed",
      site_covs = site_covs_wide, obs_covs = obs_covs
    ) %>%
      rename_with(~ paste0(., ".", yr), -locality_id)
    occ_wide_list[[as.character(yr)]] <- occ_wide_yr
  }

  occ_wide <- occ_wide_list[[1]]
  for (k in 2:length(occ_wide_list)) {
    occ_wide <- occ_wide %>%
      full_join(occ_wide_list[[k]], by = "locality_id", multiple = "all")
  }
  occ_wide <- occ_wide %>% distinct(locality_id, .keep_all = TRUE)

  cols_logical <- sapply(occ_wide, is.logical)
  occ_wide[, cols_logical] <- lapply(occ_wide[, cols_logical], as.numeric)

  first_year <- YEARS[1]
  for (v in c(vars_to_scale, "latitude", "longitude")) {
    base_col <- paste0(v, ".", first_year)
    if (base_col %in% names(occ_wide)) {
      occ_wide <- occ_wide %>% rename(!!v := !!base_col)
      year_cols <- paste0(v, ".", YEARS[-1])
      existing_cols <- year_cols[year_cols %in% names(occ_wide)]
      if (length(existing_cols) > 0) {
        occ_wide <- occ_wide %>%
          mutate(!!v := coalesce(!!sym(v), !!!syms(existing_cols)))
      }
    }
  }
  occ_wide_clean <- occ_wide %>% distinct(locality_id, .keep_all = TRUE)
  result$n_wide <- nrow(occ_wide_clean)

  # --- PASO 5: Matrices de deteccion ---
  y_list <- list()
  duration_list <- list()
  effort_list <- list()
  observers_list <- list()
  time_list <- list()

  for (yr in YEARS) {
    det_cols <- grep(paste0("^y\\.\\d+\\.", yr, "$"), names(occ_wide_clean), value = TRUE)
    rep_nums <- as.integer(sub(paste0("y\\.(\\d+)\\.", yr), "\\1", det_cols))
    det_cols <- det_cols[order(rep_nums)]

    n_have <- length(det_cols)
    if (n_have > 0) {
      y_yr <- as.matrix(occ_wide_clean[, det_cols])
    } else {
      y_yr <- matrix(NA, nrow = nrow(occ_wide_clean), ncol = 0)
    }
    if (n_have < J_reps) {
      y_yr <- cbind(y_yr, matrix(NA, nrow = nrow(occ_wide_clean), ncol = J_reps - n_have))
    } else if (n_have > J_reps) {
      y_yr <- y_yr[, 1:J_reps]
    }
    y_list[[as.character(yr)]] <- y_yr

    for (cov_info in list(
      list(prefix = "duration_minutes", lst = "duration_list"),
      list(prefix = "effort_distance_km", lst = "effort_list"),
      list(prefix = "number_observers", lst = "observers_list"),
      list(prefix = "time_observations_started", lst = "time_list")
    )) {
      cov_cols <- grep(paste0(cov_info$prefix, "\\.\\d+\\.", yr),
                       names(occ_wide_clean), value = TRUE)
      rep_n <- as.integer(sub(paste0(cov_info$prefix, "\\.(\\d+)\\.", yr), "\\1", cov_cols))
      cov_cols <- cov_cols[order(rep_n)]
      n_c <- length(cov_cols)
      if (n_c > 0) {
        cov_mat <- as.matrix(occ_wide_clean[, cov_cols])
      } else {
        cov_mat <- matrix(NA, nrow = nrow(occ_wide_clean), ncol = 0)
      }
      if (n_c < J_reps) {
        cov_mat <- cbind(cov_mat, matrix(NA, nrow = nrow(occ_wide_clean), ncol = J_reps - n_c))
      } else if (n_c > J_reps) {
        cov_mat <- cov_mat[, 1:J_reps]
      }
      if (cov_info$lst == "duration_list")  duration_list  <- c(duration_list, list(cov_mat))
      if (cov_info$lst == "effort_list")    effort_list    <- c(effort_list, list(cov_mat))
      if (cov_info$lst == "observers_list") observers_list <- c(observers_list, list(cov_mat))
      if (cov_info$lst == "time_list")      time_list      <- c(time_list, list(cov_mat))
    }
  }

  y.cross   <- do.call(cbind, y_list)
  duration  <- do.call(cbind, duration_list)
  effort    <- do.call(cbind, effort_list)
  observers <- do.call(cbind, observers_list)
  time_obs  <- do.call(cbind, time_list)

  n_sites <- nrow(y.cross)
  result$pct_na_pre <- round(100 * mean(is.na(y.cross)), 1)
  result$det_pre <- sum(y.cross == 1, na.rm = TRUE)

  # --- PASO 6: Filtro temporal ---
  has_data <- matrix(FALSE, n_sites, T_years)
  has_det  <- matrix(FALSE, n_sites, T_years)
  for (t in 1:T_years) {
    cols_t <- ((t-1) * J_reps + 1):(t * J_reps)
    has_data[, t] <- apply(y.cross[, cols_t], 1, function(x) any(!is.na(x)))
    has_det[, t]  <- apply(y.cross[, cols_t], 1, function(x) any(x == 1, na.rm = TRUE))
  }
  n_years_data <- rowSums(has_data)
  has_early <- rowSums(has_data[, 1:3]) > 0
  has_late  <- rowSums(has_data[, 4:6]) > 0
  ever_det  <- rowSums(has_det) > 0

  result$sites_pre <- n_sites
  result$detected_pre <- sum(ever_det)
  result$pct_1year <- round(100 * mean(n_years_data == 1), 1)
  result$sites_eligible <- sum(n_years_data >= 3 & has_early & has_late)

  keep <- n_years_data >= 3 & has_early & has_late

  if (sum(keep) < 10) {
    cat("  WARNING: Solo ", sum(keep), " sitios pasan el filtro. Relajando a >= 2 anios...\n")
    keep <- n_years_data >= 2 & has_early & has_late
    result$filter_relaxed <- TRUE
    if (sum(keep) < 10) {
      keep <- n_years_data >= 2
      result$filter_relaxed <- "very_relaxed"
    }
  } else {
    result$filter_relaxed <- FALSE
  }

  y_filt <- y.cross[keep, ]
  duration_filt  <- duration[keep, ]
  effort_filt    <- effort[keep, ]
  observers_filt <- observers[keep, ]
  time_filt      <- time_obs[keep, ]
  siteCovs_filt  <- occ_wide_clean[keep, vars_to_scale]

  n_filt <- nrow(y_filt)
  result$sites_post <- n_filt
  result$pct_na_post <- round(100 * mean(is.na(y_filt)), 1)
  result$det_post <- sum(y_filt == 1, na.rm = TRUE)
  result$detected_post <- sum(apply(y_filt, 1, function(x) any(x == 1, na.rm = TRUE)))

  cat("  Filtro temporal: ", result$sites_pre, " -> ", n_filt, " sitios\n")
  cat("  Detecciones: ", result$det_post, " en ", result$detected_post, " sitios\n")

  # Naive por anio
  naive_vec <- c()
  for (t in 1:T_years) {
    cols_t <- ((t-1) * J_reps + 1):(t * J_reps)
    naive <- mean(apply(y_filt[, cols_t], 1, function(x) any(x == 1, na.rm = TRUE)))
    naive_vec <- c(naive_vec, naive)
  }
  result$naive_by_year <- naive_vec
  names(result$naive_by_year) <- YEARS
  cat("  Naive: ", paste(sprintf("%d:%.1f%%", YEARS, naive_vec * 100), collapse = " "), "\n")

  # --- PASO 7: Ajustar modelos ---
  duration_sc  <- scale(duration_filt)
  effort_sc    <- scale(effort_filt)
  observers_sc <- scale(observers_filt)
  time_sc      <- scale(time_filt)

  siteCovs_df <- as.data.frame(siteCovs_filt)

  occ_umf <- tryCatch({
    unmarkedMultFrame(
      y = y_filt, siteCovs = siteCovs_df,
      obsCovs = list(duration = duration_sc, effort = effort_sc,
                     observers = observers_sc, time = time_sc),
      numPrimary = T_years
    )
  }, error = function(e) {
    cat("  ERROR construyendo UMF: ", e$message, "\n")
    return(NULL)
  })

  if (is.null(occ_umf)) {
    result$error <- "UMF construction failed"
    return(result)
  }

  cat("  UMF: ", n_filt, " sitios x ", T_years, " periodos x ", J_reps, " replicas\n\n")

  # Fit models
  fit_model <- function(formula_list, name) {
    cat("  Ajustando: ", name, "...\n")
    tryCatch({
      do.call(colext, c(formula_list, list(data = occ_umf)))
    }, error = function(e) {
      cat("    Error: ", e$message, "\n")
      NULL
    })
  }

  models <- list(
    "1: Null" = fit_model(
      list(psiformula = ~1, gammaformula = ~1, epsilonformula = ~1, pformula = ~1),
      "Null (~1,~1,~1,~1)"),
    "2: p=effort+obs" = fit_model(
      list(psiformula = ~1, gammaformula = ~1, epsilonformula = ~1,
           pformula = ~effort + observers),
      "p=effort+obs"),
    "3: psi+p" = fit_model(
      list(psiformula = ~grass_cover + tree_cover, gammaformula = ~1,
           epsilonformula = ~1, pformula = ~effort + observers),
      "psi+p"),
    "4: psi+p completo" = fit_model(
      list(psiformula = ~grass_cover + tree_cover, gammaformula = ~1,
           epsilonformula = ~1,
           pformula = ~effort + observers + duration + time),
      "psi+p completo"),
    "5: psi full" = fit_model(
      list(psiformula = ~bio1 + bio2 + tree_cover + grass_cover + topo_elev,
           gammaformula = ~1, epsilonformula = ~1,
           pformula = ~effort + observers + duration + time),
      "psi full")
  )

  # Extract results
  result$models <- list()
  for (mname in names(models)) {
    mod <- models[[mname]]
    if (is.null(mod)) {
      result$models[[mname]] <- list(error = TRUE)
      next
    }
    coefs <- coef(mod)
    s <- tryCatch(summary(mod), error = function(e) NULL)
    result$models[[mname]] <- list(
      psi1    = plogis(coefs["psi(Int)"]),
      gamma   = plogis(coefs["col(Int)"]),
      epsilon = plogis(coefs["ext(Int)"]),
      p       = plogis(coefs["p(Int)"]),
      aic     = mod@AIC,
      conv    = mod@opt$convergence,
      has_se  = !is.null(s),
      coefs   = coefs,
      summary = s
    )
  }

  # Load original model
  orig_path <- file.path(base_repo, "data/processed",
                          paste0(sp_code, "_model_object.rds"))
  if (file.exists(orig_path)) {
    mod_orig <- readRDS(orig_path)
    coefs_orig <- coef(mod_orig)
    result$original <- list(
      psi1    = plogis(coefs_orig["psi(Int)"]),
      gamma   = plogis(coefs_orig["col(Int)"]),
      epsilon = plogis(coefs_orig["ext(Int)"]),
      p       = plogis(coefs_orig["p(Int)"]),
      aic     = mod_orig@AIC,
      conv    = mod_orig@opt$convergence,
      n_sites = nrow(mod_orig@data@y),
      pct_na  = round(100 * mean(is.na(mod_orig@data@y)), 1),
      ext_logit = coefs_orig["ext(Int)"]
    )
    cat("\n  ORIGINAL: psi1=", round(result$original$psi1 * 100, 2), "%",
        " gamma=", round(result$original$gamma * 100, 2), "%",
        " epsilon=", round(result$original$epsilon * 100, 1), "%",
        " p=", round(result$original$p * 100, 1), "%\n")
  }

  return(result)
}

###############################################################################
# RUN ALL SPECIES
###############################################################################

all_results <- list()

for (sp_code in names(species_info)) {
  sp <- species_info[[sp_code]]
  res <- tryCatch(
    run_pipeline_v2(sp$code, sp$name),
    error = function(e) {
      cat("  CRITICAL ERROR for ", sp$code, ": ", e$message, "\n")
      list(species = sp$code, name = sp$name, error = e$message)
    }
  )
  all_results[[sp_code]] <- res
}

###############################################################################
# TABLA RESUMEN COMPARATIVA
###############################################################################

cat("\n\n")
cat(strrep("=", 90), "\n")
cat("  TABLA RESUMEN: PIPELINE v2 vs ORIGINAL (4 especies)\n")
cat(strrep("=", 90), "\n\n")

cat(sprintf("%-10s %-20s | %7s %7s %7s %7s | %6s | %5s\n",
            "Especie", "Pipeline", "psi1", "gamma", "epsil", "p", "Sites", "Conv"))
cat(strrep("-", 90), "\n")

for (sp_code in names(species_info)) {
  res <- all_results[[sp_code]]

  if (!is.null(res$error) && is.character(res$error)) {
    cat(sprintf("%-10s %-20s | %s\n", sp_code, "ERROR", res$error))
    next
  }

  # Original
  if (!is.null(res$original)) {
    cat(sprintf("%-10s %-20s | %6.1f%% %6.2f%% %6.1f%% %6.1f%% | %6d | %5s\n",
                sp_code, "ORIGINAL",
                res$original$psi1 * 100, res$original$gamma * 100,
                res$original$epsilon * 100, res$original$p * 100,
                res$original$n_sites,
                ifelse(res$original$epsilon > 0.99, "BNDRY", as.character(res$original$conv))))
  }

  # Best v2 model (lowest AIC with convergence = 0)
  best_name <- NULL
  best_aic <- Inf
  for (mname in names(res$models)) {
    m <- res$models[[mname]]
    if (is.null(m$error) && !isTRUE(m$error) && m$conv == 0 && m$aic < best_aic) {
      best_aic <- m$aic
      best_name <- mname
    }
  }

  if (!is.null(best_name)) {
    m <- res$models[[best_name]]
    cat(sprintf("%-10s %-20s | %6.1f%% %6.2f%% %6.1f%% %6.1f%% | %6d | %5d\n",
                "", paste0("v2: ", best_name),
                m$psi1 * 100, m$gamma * 100, m$epsilon * 100, m$p * 100,
                res$sites_post, m$conv))
  }

  # Also show null model v2
  m_null <- res$models[["1: Null"]]
  if (!is.null(m_null) && !isTRUE(m_null$error)) {
    cat(sprintf("%-10s %-20s | %6.1f%% %6.2f%% %6.1f%% %6.1f%% | %6d | %5d\n",
                "", "v2: Null",
                m_null$psi1 * 100, m_null$gamma * 100, m_null$epsilon * 100, m_null$p * 100,
                res$sites_post, m_null$conv))
  }
  cat(strrep("-", 90), "\n")
}

cat("\nReferencia biologica plausible: psi1=5-30%, gamma=1-15%, epsilon=5-25%, p=10-50%\n\n")

# Detailed model comparison per species
for (sp_code in names(species_info)) {
  res <- all_results[[sp_code]]
  if (!is.null(res$error) && is.character(res$error)) next

  cat("\n")
  cat(strrep("=", 70), "\n")
  cat("  DETALLE: ", species_info[[sp_code]]$name, "\n")
  cat(strrep("=", 70), "\n\n")

  cat("  Datos:\n")
  cat("    Registros crudos: ", res$n_raw, "\n")
  cat("    Tras filter_repeat_visits: ", res$n_filtered, "\n")
  cat("    Localidades: ", res$n_localities, "\n")
  cat("    Sitios pre-filtro: ", res$sites_pre, "\n")
  cat("    Sitios con 1 solo anio: ", res$pct_1year, "%\n")
  cat("    Sitios post-filtro: ", res$sites_post, "\n")
  cat("    Sitios con deteccion: ", res$detected_post, "\n")
  cat("    % NAs pre/post: ", res$pct_na_pre, "% / ", res$pct_na_post, "%\n")
  if (!is.null(res$filter_relaxed) && res$filter_relaxed != FALSE) {
    cat("    NOTA: Filtro relajado a ", res$filter_relaxed, "\n")
  }
  cat("    Naive por anio: ", paste(sprintf("%d:%.1f%%", YEARS, res$naive_by_year * 100), collapse = " "), "\n\n")

  cat(sprintf("  %-25s %7s %7s %7s %7s %9s %5s %5s\n",
              "Modelo", "psi1", "gamma", "epsil", "p", "AIC", "Conv", "SEs"))
  cat("  ", strrep("-", 75), "\n")

  for (mname in names(res$models)) {
    m <- res$models[[mname]]
    if (isTRUE(m$error)) {
      cat(sprintf("  %-25s  --- ERROR ---\n", mname))
      next
    }
    cat(sprintf("  %-25s %6.1f%% %6.2f%% %6.1f%% %6.1f%% %9.1f %5d %5s\n",
                mname,
                m$psi1 * 100, m$gamma * 100, m$epsilon * 100, m$p * 100,
                m$aic, m$conv, ifelse(m$has_se, "OK", "NO")))
  }

  # Print summary of best model
  best_name <- NULL
  best_aic <- Inf
  for (mname in names(res$models)) {
    m <- res$models[[mname]]
    if (!isTRUE(m$error) && m$conv == 0 && m$aic < best_aic) {
      best_aic <- m$aic
      best_name <- mname
    }
  }
  if (!is.null(best_name) && !is.null(res$models[[best_name]]$summary)) {
    cat("\n  Summary del mejor modelo (", best_name, "):\n\n")
    s <- res$models[[best_name]]$summary
    cat("  psi:\n"); print(s$psi)
    cat("  col:\n"); print(s$col)
    cat("  ext:\n"); print(s$ext)
    cat("  det:\n"); print(s$det)
  }
}

###############################################################################
# VEREDICTO GLOBAL
###############################################################################

cat("\n\n")
cat(strrep("=", 90), "\n")
cat("  VEREDICTO GLOBAL\n")
cat(strrep("=", 90), "\n\n")

for (sp_code in names(species_info)) {
  res <- all_results[[sp_code]]
  if (!is.null(res$error) && is.character(res$error)) {
    cat(sprintf("  %-10s: ERROR - %s\n", sp_code, res$error))
    next
  }

  # Get best model epsilon
  best_eps <- 1.0
  best_name <- "none"
  for (mname in names(res$models)) {
    m <- res$models[[mname]]
    if (!isTRUE(m$error) && m$conv == 0) {
      if (m$aic < (best_eps * 1e6)) {  # use aic for selection
        best_eps <- m$epsilon
        best_name <- mname
      }
    }
  }
  # Actually find best by AIC
  best_aic <- Inf
  for (mname in names(res$models)) {
    m <- res$models[[mname]]
    if (!isTRUE(m$error) && m$conv == 0 && m$aic < best_aic) {
      best_aic <- m$aic
      best_eps <- m$epsilon
      best_name <- mname
    }
  }

  orig_eps <- if (!is.null(res$original)) res$original$epsilon else NA

  if (is.na(orig_eps)) {
    cat(sprintf("  %-10s: epsilon v2 = %.1f%% (no original para comparar)\n",
                sp_code, best_eps * 100))
  } else if (best_eps < 0.30 && orig_eps > 0.90) {
    cat(sprintf("  %-10s: EXITO! epsilon: %.0f%% -> %.1f%% (%s)\n",
                sp_code, orig_eps * 100, best_eps * 100, best_name))
  } else if (best_eps < orig_eps) {
    cat(sprintf("  %-10s: MEJORA. epsilon: %.0f%% -> %.1f%% (%s)\n",
                sp_code, orig_eps * 100, best_eps * 100, best_name))
  } else {
    cat(sprintf("  %-10s: SIN MEJORA. epsilon: %.0f%% -> %.1f%%\n",
                sp_code, orig_eps * 100, best_eps * 100))
  }
}

cat("\n\nFIN del pipeline v2 multi-especie.\n")
