###############################################################################
# test_pipeline_v2_otitar.R
#
# Pipeline alternativo para Otis tarda.
# Cambios respecto al pipeline original:
#   1. site_vars = c("locality_id")  en vez de c("locality_id","observer_id")
#   2. min_obs = 2  en vez de 3
#   3. Filtro post-join: >= 3 anos con datos, early (2017-2019) + late (2020-2022)
#   4. Modelos progresivos: null -> p-only -> psi+p -> psi+p+epsilon simple
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

cat("\n")
cat(strrep("=", 70), "\n")
cat("  PIPELINE v2: Otis tarda\n")
cat(strrep("=", 70), "\n\n")

###############################################################################
# PASO 1: Cargar datos crudos (output de step 1)
###############################################################################

cat("--- PASO 1: Cargar datos crudos ---\n\n")

raw_path <- file.path(
  "/Users/gfandos/Documents/GitHub/dynamic_occupancy_models_steppe_birds",
  "data-raw/data/otitar/ebd_otitar_breeding_spain_zf.csv"
)
occ_raw <- read.csv(raw_path)
occ_raw <- occ_raw %>% filter(year >= 2017, year <= 2022)
cat("  Registros crudos: ", nrow(occ_raw), "\n\n")

###############################################################################
# PASO 2: filter_repeat_visits() con nueva configuracion
###############################################################################

cat("--- PASO 2: filter_repeat_visits (site_vars=locality_id, min_obs=2) ---\n\n")

occ <- filter_repeat_visits(
  occ_raw,
  min_obs  = 2,
  max_obs  = 10,
  annual_closure = TRUE,
  date_var = "observation_date",
  site_vars = c("locality_id")
)

cat("  Registros tras filtro: ", nrow(occ), "\n")
cat("  Sitios (locality x anio): ", n_distinct(occ$site), "\n")
cat("  Localidades unicas: ", n_distinct(occ$locality_id), "\n\n")

###############################################################################
# PASO 3: Extraer covariables ambientales
###############################################################################

cat("--- PASO 3: Extraer covariables ambientales ---\n\n")

base_data <- "/Users/gfandos/Documents/GitHub/dynamic_occupancy_models_steppe_birds/data-raw/data"

variables_environmental <- raster::stack(
  file.path(base_data, "environmental_data/environmental_data_occ/variables_spain.grd")
)
variables_aspect <- raster::stack(
  file.path(base_data, "topology_data/topo_aspect.asc")
)
variables_elev <- raster::stack(
  file.path(base_data, "topology_data/topo_elev.asc")
)

# Extraer valores en coordenadas de los sitios
coords <- occ[, c("longitude", "latitude")]

# Use raster::extract for RasterStack objects (more compatible with cells=TRUE)
env_ext <- as.data.frame(raster::extract(variables_environmental, coords, cellnumbers = TRUE))

# Extract topo variables - remove ID column and rename to standard names
aspect_raw <- as.data.frame(terra::extract(variables_aspect, coords))
elev_raw   <- as.data.frame(terra::extract(variables_elev, coords))

# Remove ID column if present, keep only value column(s)
if ("ID" %in% names(aspect_raw)) aspect_raw <- aspect_raw[, -which(names(aspect_raw) == "ID"), drop = FALSE]
if ("ID" %in% names(elev_raw))   elev_raw   <- elev_raw[, -which(names(elev_raw) == "ID"), drop = FALSE]

# Rename to standard names
names(aspect_raw) <- "topo_aspect"
names(elev_raw)   <- "topo_elev"

occ_var <- occ %>%
  cbind(env_ext[, setdiff(names(env_ext), "ID")]) %>%
  cbind(aspect_raw) %>%
  cbind(elev_raw)

# Renombrar la columna de celda del raster
if ("cell" %in% names(occ_var)) names(occ_var)[names(occ_var) == "cell"] <- "cells"

# Eliminar registros con NAs en variables clave
occ_var <- occ_var %>% drop_na(bio12, tree_cover, topo_elev)
cat("  Registros con covariables validas: ", nrow(occ_var), "\n")
cat("  Localidades con covariables: ", n_distinct(occ_var$locality_id), "\n\n")

###############################################################################
# PASO 4: Formato wide por anio + merge
###############################################################################

cat("--- PASO 4: Convertir a formato wide y merge por localidad ---\n\n")

vars_to_scale <- c("bio1", "tree_cover", "bio2", "grass_cover",
                    "topo_aspect", "topo_elev")
obs_covs <- c("time_observations_started", "duration_minutes",
              "effort_distance_km", "number_observers")
site_covs_wide <- c("locality_id",
                     "latitude", "longitude", vars_to_scale)

# Estandarizar covariables de sitio
scaling_params <- list()
occ_var_std <- occ_var
for (v in vars_to_scale) {
  vals <- occ_var_std[[v]]
  scaling_params[[v]] <- list(center = mean(vals, na.rm = TRUE),
                               scale  = sd(vals, na.rm = TRUE))
  occ_var_std[[v]] <- (vals - scaling_params[[v]]$center) / scaling_params[[v]]$scale
}

# Wide format por anio
occ_wide_list <- list()
for (yr in YEARS) {
  occ_yr <- occ_var_std %>% filter(year == yr)
  if (nrow(occ_yr) == 0) next

  occ_wide_yr <- format_unmarked_occu(
    occ_yr,
    site_id = "site",
    response = "species_observed",
    site_covs = site_covs_wide,
    obs_covs = obs_covs
  ) %>%
    rename_with(~ paste0(., ".", yr), -locality_id)
  occ_wide_list[[as.character(yr)]] <- occ_wide_yr
}

# Merge todos los anios por locality_id
occ_wide <- occ_wide_list[[1]]
for (k in 2:length(occ_wide_list)) {
  occ_wide <- occ_wide %>%
    full_join(occ_wide_list[[k]], by = "locality_id", multiple = "all")
}
occ_wide <- occ_wide %>% distinct(locality_id, .keep_all = TRUE)

# Convertir logicos a numerico
cols_logical <- sapply(occ_wide, is.logical)
occ_wide[, cols_logical] <- lapply(occ_wide[, cols_logical], as.numeric)

# Consolidar covariables de sitio (tomar primer valor no-NA entre anios)
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

# Deduplicar por localidad (ya que ya no tenemos la columna cells en wide)
occ_wide_clean <- occ_wide %>%
  distinct(locality_id, .keep_all = TRUE)

cat("  Sitios totales (celdas unicas): ", nrow(occ_wide_clean), "\n\n")

###############################################################################
# PASO 5: Construir matrices de deteccion y covariables
###############################################################################

cat("--- PASO 5: Construir matrices ---\n\n")

J_reps <- 10  # max replicas por periodo primario

# Columnas de deteccion: y.{1..10}.{year}
# (format_unmarked_occu names detection columns y.1, y.2, etc.)
det_pattern <- "^y\\."
det_cols_all <- grep(det_pattern, names(occ_wide_clean), value = TRUE)
cat("  Detection columns found: ", length(det_cols_all), "\n")
cat("  Sample: ", paste(head(det_cols_all, 5), collapse=", "), "\n")

# Organizar las columnas por anio
y_list <- list()
duration_list <- list()
effort_list <- list()
observers_list <- list()
time_list <- list()

for (yr in YEARS) {
  yr_suffix <- paste0("\\.", yr, "$")

  # Detection columns for this year: y.{rep}.{year}
  det_cols <- grep(paste0("^y\\.\\d+\\.", yr, "$"), names(occ_wide_clean), value = TRUE)
  # Sort by replicate number
  rep_nums <- as.integer(sub(paste0("y\\.(\\d+)\\.", yr), "\\1", det_cols))
  det_cols <- det_cols[order(rep_nums)]

  # Pad to J_reps columns
  n_have <- length(det_cols)
  if (n_have > 0) {
    y_yr <- as.matrix(occ_wide_clean[, det_cols])
  } else {
    y_yr <- matrix(NA, nrow = nrow(occ_wide_clean), ncol = 0)
  }
  if (n_have < J_reps) {
    pad <- matrix(NA, nrow = nrow(occ_wide_clean), ncol = J_reps - n_have)
    y_yr <- cbind(y_yr, pad)
  } else if (n_have > J_reps) {
    y_yr <- y_yr[, 1:J_reps]
  }
  y_list[[as.character(yr)]] <- y_yr

  # Observation covariates
  for (cov_info in list(
    list(prefix = "duration_minutes", out = "duration_list"),
    list(prefix = "effort_distance_km", out = "effort_list"),
    list(prefix = "number_observers", out = "observers_list"),
    list(prefix = "time_observations_started", out = "time_list")
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
      cov_mat <- cbind(cov_mat, matrix(NA, nrow = nrow(occ_wide_clean),
                                        ncol = J_reps - n_c))
    } else if (n_c > J_reps) {
      cov_mat <- cov_mat[, 1:J_reps]
    }
    # Store in appropriate list
    assign(cov_info$out,
           c(get(cov_info$out), list(cov_mat)),
           envir = .GlobalEnv)
  }
}

# Concatenar matrices por anio
y.cross <- do.call(cbind, y_list)
duration  <- do.call(cbind, duration_list)
effort    <- do.call(cbind, effort_list)
observers <- do.call(cbind, observers_list)
time_obs  <- do.call(cbind, time_list)

n_sites <- nrow(y.cross)
cat("  Matriz y: ", n_sites, " x ", ncol(y.cross), "\n")
cat("  % NAs: ", round(100 * mean(is.na(y.cross)), 1), "%\n")
cat("  Detecciones: ", sum(y.cross == 1, na.rm = TRUE), "\n")

###############################################################################
# PASO 6: Filtro temporal
###############################################################################

cat("\n--- PASO 6: Filtro temporal ---\n\n")

has_data <- matrix(FALSE, n_sites, T_years)
has_det <- matrix(FALSE, n_sites, T_years)
for (t in 1:T_years) {
  cols_t <- ((t-1) * J_reps + 1):(t * J_reps)
  has_data[, t] <- apply(y.cross[, cols_t], 1, function(x) any(!is.na(x)))
  has_det[, t]  <- apply(y.cross[, cols_t], 1, function(x) any(x == 1, na.rm = TRUE))
}
n_years_data <- rowSums(has_data)
has_early <- rowSums(has_data[, 1:3]) > 0
has_late  <- rowSums(has_data[, 4:6]) > 0
ever_det  <- rowSums(has_det) > 0

cat("  Antes del filtro:\n")
cat("    Sitios: ", n_sites, "\n")
cat("    Sitios detectados: ", sum(ever_det), "\n")
cat("    Sitios con 1 anio: ", sum(n_years_data == 1),
    " (", round(100 * mean(n_years_data == 1), 1), "%)\n")
cat("    Sitios con >= 3 anios + E+L: ",
    sum(n_years_data >= 3 & has_early & has_late), "\n")

# Aplicar filtro
keep <- n_years_data >= 3 & has_early & has_late

y_filt <- y.cross[keep, ]
duration_filt <- duration[keep, ]
effort_filt <- effort[keep, ]
observers_filt <- observers[keep, ]
time_filt <- time_obs[keep, ]
siteCovs_filt <- occ_wide_clean[keep, vars_to_scale]

n_filt <- nrow(y_filt)
cat("\n  Despues del filtro (>= 3 anios, E+L):\n")
cat("    Sitios: ", n_filt, "\n")
cat("    % NAs: ", round(100 * mean(is.na(y_filt)), 1), "%\n")
cat("    Detecciones: ", sum(y_filt == 1, na.rm = TRUE), "\n")
cat("    Sitios detectados: ",
    sum(apply(y_filt, 1, function(x) any(x == 1, na.rm = TRUE))), "\n")

# Naive por anio
cat("    Naive por anio: ")
for (t in 1:T_years) {
  cols_t <- ((t-1) * J_reps + 1):(t * J_reps)
  naive <- mean(apply(y_filt[, cols_t], 1, function(x) any(x == 1, na.rm = TRUE)))
  cat(sprintf("%d:%.1f%% ", YEARS[t], naive * 100))
}
cat("\n")

###############################################################################
# PASO 7: Construir unmarkedMultFrame y ajustar modelos
###############################################################################

cat("\n--- PASO 7: Ajustar modelos ---\n\n")

# Estandarizar covariables de observacion
duration_sc  <- scale(duration_filt)
effort_sc    <- scale(effort_filt)
observers_sc <- scale(observers_filt)
time_sc      <- scale(time_filt)

# Site covariates (ya estandarizadas en paso 4)
siteCovs_df <- as.data.frame(siteCovs_filt)

# Build UMF
occ_umf <- unmarkedMultFrame(
  y = y_filt,
  siteCovs = siteCovs_df,
  obsCovs = list(
    duration  = duration_sc,
    effort    = effort_sc,
    observers = observers_sc,
    time      = time_sc
  ),
  numPrimary = T_years
)

cat("  UMF construido: ", nrow(y_filt), " sitios, ",
    T_years, " periodos primarios, ", J_reps, " replicas\n\n")

# --- Modelo 1: Null ---
cat("  Ajustando Modelo 1: NULL (~1, ~1, ~1, ~1)...\n")
mod1 <- tryCatch({
  colext(psiformula = ~1, gammaformula = ~1,
         epsilonformula = ~1, pformula = ~1,
         data = occ_umf)
}, error = function(e) { cat("  Error: ", e$message, "\n"); NULL })

# --- Modelo 2: Solo deteccion ---
cat("  Ajustando Modelo 2: p = ~effort+observers (~1, ~1, ~1, ~effort+observers)...\n")
mod2 <- tryCatch({
  colext(psiformula = ~1, gammaformula = ~1,
         epsilonformula = ~1, pformula = ~effort + observers,
         data = occ_umf)
}, error = function(e) { cat("  Error: ", e$message, "\n"); NULL })

# --- Modelo 3: psi + deteccion ---
cat("  Ajustando Modelo 3: psi + p (~grass_cover+tree_cover, ~1, ~1, ~effort+observers)...\n")
mod3 <- tryCatch({
  colext(psiformula = ~grass_cover + tree_cover,
         gammaformula = ~1,
         epsilonformula = ~1,
         pformula = ~effort + observers,
         data = occ_umf)
}, error = function(e) { cat("  Error: ", e$message, "\n"); NULL })

# --- Modelo 4: psi + p + duracion deteccion ---
cat("  Ajustando Modelo 4: p completo (~grass_cover+tree_cover, ~1, ~1, ~effort+observers+duration+time)...\n")
mod4 <- tryCatch({
  colext(psiformula = ~grass_cover + tree_cover,
         gammaformula = ~1,
         epsilonformula = ~1,
         pformula = ~effort + observers + duration + time,
         data = occ_umf)
}, error = function(e) { cat("  Error: ", e$message, "\n"); NULL })

# --- Modelo 5: Modelo original de psi (sin dinamicas) ---
cat("  Ajustando Modelo 5: psi completo (~bio1+bio2+tree_cover+grass_cover+topo_elev, ~1, ~1, ~effort+observers+duration+time)...\n")
mod5 <- tryCatch({
  colext(psiformula = ~bio1 + bio2 + tree_cover + grass_cover + topo_elev,
         gammaformula = ~1,
         epsilonformula = ~1,
         pformula = ~effort + observers + duration + time,
         data = occ_umf)
}, error = function(e) { cat("  Error: ", e$message, "\n"); NULL })

###############################################################################
# PASO 8: Resultados
###############################################################################

cat("\n\n")
cat(strrep("=", 70), "\n")
cat("  RESULTADOS\n")
cat(strrep("=", 70), "\n\n")

models <- list(
  "1: Null" = mod1,
  "2: p=effort+obs" = mod2,
  "3: psi+p simple" = mod3,
  "4: psi+p completo" = mod4,
  "5: psi completo" = mod5
)

# Tabla resumen
cat(sprintf("%-25s %7s %7s %7s %7s %7s %5s\n",
            "Modelo", "psi1", "gamma", "epsil", "p(int)", "AIC", "Conv"))
cat(strrep("-", 75), "\n")

for (mname in names(models)) {
  mod <- models[[mname]]
  if (is.null(mod)) {
    cat(sprintf("%-25s  --- ERROR ---\n", mname))
    next
  }

  coefs <- coef(mod)
  psi_int <- plogis(coefs["psi(Int)"])
  col_int <- plogis(coefs["col(Int)"])
  ext_int <- plogis(coefs["ext(Int)"])
  det_int <- plogis(coefs["p(Int)"])

  cat(sprintf("%-25s %6.1f%% %6.2f%% %6.1f%% %6.1f%% %7.1f %5d\n",
              mname,
              psi_int * 100, col_int * 100, ext_int * 100, det_int * 100,
              mod@AIC, mod@opt$convergence))
}

cat("\nReferencia plausible: psi1=5-30%, gamma=1-15%, epsilon=5-25%, p=10-50%\n")

# Detalle de cada modelo
for (mname in names(models)) {
  mod <- models[[mname]]
  if (is.null(mod)) next

  cat("\n\n--- Detalle: ", mname, " ---\n\n")

  s <- tryCatch(summary(mod), error = function(e) NULL)
  if (!is.null(s)) {
    cat("psi:\n"); print(s$psi); cat("\n")
    cat("col:\n"); print(s$col); cat("\n")
    cat("ext:\n"); print(s$ext); cat("\n")
    cat("det:\n"); print(s$det); cat("\n")
    cat("AIC: ", mod@AIC, " | Convergence: ", mod@opt$convergence, "\n")
  } else {
    cat("  (No se pudo obtener summary - Hessiana singular)\n")
    coefs <- coef(mod)
    for (cn in names(coefs)) {
      cat(sprintf("  %s: %.4f -> prob=%.4f\n", cn, coefs[cn], plogis(coefs[cn])))
    }
  }
}

###############################################################################
# COMPARACION CON MODELO ORIGINAL
###############################################################################

cat("\n\n")
cat(strrep("=", 70), "\n")
cat("  COMPARACION CON MODELO ORIGINAL\n")
cat(strrep("=", 70), "\n\n")

mod_orig <- readRDS(
  "/Users/gfandos/Documents/GitHub/dynamic_occupancy_models_steppe_birds/data/processed/otitar_model_object.rds"
)
coefs_orig <- coef(mod_orig)

cat("MODELO ORIGINAL (pipeline actual):\n")
cat(sprintf("  psi1 (Int): logit=%.3f -> prob=%.4f (%.2f%%)\n",
            coefs_orig["psi(Int)"],
            plogis(coefs_orig["psi(Int)"]),
            plogis(coefs_orig["psi(Int)"]) * 100))
cat(sprintf("  gamma (Int): logit=%.3f -> prob=%.4f (%.2f%%)\n",
            coefs_orig["col(Int)"],
            plogis(coefs_orig["col(Int)"]),
            plogis(coefs_orig["col(Int)"]) * 100))
cat(sprintf("  epsilon (Int): logit=%.3f -> prob=%.4f (%.2f%%)\n",
            coefs_orig["ext(Int)"],
            plogis(coefs_orig["ext(Int)"]),
            plogis(coefs_orig["ext(Int)"]) * 100))
cat(sprintf("  p (Int): logit=%.3f -> prob=%.4f (%.2f%%)\n",
            coefs_orig["p(Int)"],
            plogis(coefs_orig["p(Int)"]),
            plogis(coefs_orig["p(Int)"]) * 100))
cat("  AIC: ", mod_orig@AIC, "\n")
cat("  Sitios: ", nrow(mod_orig@data@y), "\n")
cat("  Epsilon en boundary (logit=", round(coefs_orig["ext(Int)"], 1), ")\n\n")

# Best v2 model
best_v2 <- NULL
best_aic <- Inf
for (mname in names(models)) {
  mod <- models[[mname]]
  if (!is.null(mod) && mod@opt$convergence == 0 && mod@AIC < best_aic) {
    best_aic <- mod@AIC
    best_v2 <- mname
  }
}

if (!is.null(best_v2)) {
  mod_v2 <- models[[best_v2]]
  coefs_v2 <- coef(mod_v2)

  cat("MEJOR MODELO v2 (", best_v2, "):\n")
  cat(sprintf("  psi1 (Int): logit=%.3f -> prob=%.4f (%.2f%%)\n",
              coefs_v2["psi(Int)"],
              plogis(coefs_v2["psi(Int)"]),
              plogis(coefs_v2["psi(Int)"]) * 100))
  cat(sprintf("  gamma (Int): logit=%.3f -> prob=%.4f (%.2f%%)\n",
              coefs_v2["col(Int)"],
              plogis(coefs_v2["col(Int)"]),
              plogis(coefs_v2["col(Int)"]) * 100))
  cat(sprintf("  epsilon (Int): logit=%.3f -> prob=%.4f (%.2f%%)\n",
              coefs_v2["ext(Int)"],
              plogis(coefs_v2["ext(Int)"]),
              plogis(coefs_v2["ext(Int)"]) * 100))
  cat(sprintf("  p (Int): logit=%.3f -> prob=%.4f (%.2f%%)\n",
              coefs_v2["p(Int)"],
              plogis(coefs_v2["p(Int)"]),
              plogis(coefs_v2["p(Int)"]) * 100))
  cat("  AIC: ", mod_v2@AIC, "\n")
  cat("  Sitios: ", n_filt, "\n\n")

  # Resumen
  cat(strrep("=", 70), "\n")
  cat("  VEREDICTO\n")
  cat(strrep("=", 70), "\n\n")

  ext_orig <- plogis(coefs_orig["ext(Int)"])
  ext_v2 <- plogis(coefs_v2["ext(Int)"])

  if (ext_v2 < 0.30) {
    cat("  EXITO: epsilon bajo de ", round(ext_orig * 100, 0), "% a ",
        round(ext_v2 * 100, 1), "%\n")
    cat("  Los cambios en el filtrado FUNCIONAN.\n")
  } else {
    cat("  PARCIAL: epsilon = ", round(ext_v2 * 100, 1), "%\n")
    cat("  Mejora respecto al ", round(ext_orig * 100, 0),
        "% original, pero sigue alto.\n")
    cat("  Se recomienda tambien aumentar la escala espacial a 10 km.\n")
  }
}

cat("\nFIN del pipeline v2 para otitar.\n")
