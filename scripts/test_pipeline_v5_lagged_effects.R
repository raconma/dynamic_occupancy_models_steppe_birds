###############################################################################
# test_pipeline_v5_lagged_effects.R
#
# Pipeline v5 para Otis tarda: efectos retardados (lag-1) en gamma/epsilon.
#
# Extiende el pipeline v3 incorporando variables anuales con retardo temporal
# (lag-1): la covariable del ano t-1 predice la transicion del ano t al t+1.
#
# Hipotesis ecologica:
#   Las condiciones ambientales de un ano tienen efectos retardados sobre la
#   dinamica de ocupacion del siguiente. Ejemplo: una sequia en t-1 reduce
#   la productividad agricola y de pastizal, lo que afecta la disponibilidad
#   de habitat para la avutarda en t, influyendo en la colonizacion y
#   extincion medidas como transiciones t -> t+1.
#
# Variables lag-1 construidas:
#   Clima:
#   - NDVI_lag:     NDVI del ano anterior (productividad primaria)
#   - pr_lag:       precipitacion del ano anterior
#   - tmmx_lag:     temperatura maxima del ano anterior
#   Land cover:
#   - lc_grass_lag: % pastizal del ano anterior (IGBP clase 10)
#   - lc_crop_lag:  % cultivo del ano anterior (IGBP clase 12)
#
# Construccion del lag en colext:
#   En colext, la yearlySiteCov del ano t predice la transicion t -> t+1.
#   Para un efecto lag-1, queremos que el valor del ano t-1 prediga esa
#   misma transicion. Implementacion:
#     lag_mat[, 1] <- mat[, 1]  # 2017: sin lag real, usar mismo ano
#     lag_mat[, t] <- mat[, t-1]  # t >= 2: valor del ano anterior
#
# Estructura del pipeline:
#   Pasos 1-6: identicos al pipeline v3/v4/v4b (datos, filtro, UMF base)
#   Paso 7: Carga y preparacion de yearlySiteCovs con lag-1
#   Paso 8: Construccion del unmarkedMultFrame con yearlySiteCovs
#   Paso 9: Bateria de 12 modelos (baseline + lagged effects)
#   Paso 10: Seleccion de modelos, diagnosticos e interpretacion
#   Paso 11: Comparacion ocupacion predicha vs naive por anio
#   Paso 12: Exportar resultados a results/
#   Paso 13: Goodness-of-fit (parboot)
#
# Input:
#   - data-raw/data/otitar/ebd_otitar_breeding_spain_zf.csv (eBird zero-filled)
#   - data/processed/otitar_dynamic_variables.csv (GEE export anual)
#   - data-raw/data/environmental_data/ (covariables estaticas)
#   - data-raw/data/topology_data/ (topografia)
#
# Output:
#   - results/otitar_v5_model_selection.csv (tabla AIC 12 modelos)
#   - results/otitar_v5_best_model_coefficients.csv (coefs mejor modelo)
#   - results/otitar_v5_dynamic_scaling_params.csv (params escalado)
#
# Cambios heredados de v2:
#   1. site_vars = c("locality_id") en vez de c("locality_id","observer_id")
#   2. min_obs = 2 en vez de 3
#   3. Filtro temporal: >= 3 anos con datos, cubriendo early + late
###############################################################################

library(here)
library(auk)
library(unmarked)
library(dplyr)
library(tidyr)

select <- dplyr::select
filter <- dplyr::filter
set.seed(42)

YEARS   <- 2017:2022
T_years <- length(YEARS)

cat("\n")
cat(strrep("=", 70), "\n")
cat("  PIPELINE v5: Otis tarda - Efectos retardados (lag-1)\n")
cat("  Modelos con yearlySiteCovs lag-1 en gamma/epsilon\n")
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
# PASO 3: Extraer covariables ambientales (estaticas)
###############################################################################

cat("--- PASO 3: Extraer covariables ambientales estaticas ---\n\n")

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
cat("  Parametros de escalado guardados para ", length(scaling_params), " variables\n")

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

# Consolidar covariables de sitio
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

# Deduplicar por localidad
occ_wide_clean <- occ_wide %>%
  distinct(locality_id, .keep_all = TRUE)

cat("  Sitios totales (localidades unicas): ", nrow(occ_wide_clean), "\n\n")

###############################################################################
# PASO 5: Construir matrices de deteccion y covariables de observacion
###############################################################################

cat("--- PASO 5: Construir matrices ---\n\n")

J_reps <- 10

det_pattern <- "^y\\."
det_cols_all <- grep(det_pattern, names(occ_wide_clean), value = TRUE)
cat("  Detection columns found: ", length(det_cols_all), "\n")

y_list <- list()
duration_list <- list()
effort_list <- list()
observers_list <- list()
time_list <- list()

for (yr in YEARS) {
  yr_suffix <- paste0("\\.", yr, "$")

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
    pad <- matrix(NA, nrow = nrow(occ_wide_clean), ncol = J_reps - n_have)
    y_yr <- cbind(y_yr, pad)
  } else if (n_have > J_reps) {
    y_yr <- y_yr[, 1:J_reps]
  }
  y_list[[as.character(yr)]] <- y_yr

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
    assign(cov_info$out,
           c(get(cov_info$out), list(cov_mat)),
           envir = .GlobalEnv)
  }
}

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

keep <- n_years_data >= 3 & has_early & has_late

y_filt <- y.cross[keep, ]
duration_filt <- duration[keep, ]
effort_filt <- effort[keep, ]
observers_filt <- observers[keep, ]
time_filt <- time_obs[keep, ]
siteCovs_filt <- occ_wide_clean[keep, vars_to_scale]
locality_ids_filt <- occ_wide_clean$locality_id[keep]

n_filt <- nrow(y_filt)
cat("\n  Despues del filtro (>= 3 anios, E+L):\n")
cat("    Sitios: ", n_filt, "\n")
cat("    % NAs: ", round(100 * mean(is.na(y_filt)), 1), "%\n")
cat("    Detecciones: ", sum(y_filt == 1, na.rm = TRUE), "\n")
cat("    Sitios detectados: ",
    sum(apply(y_filt, 1, function(x) any(x == 1, na.rm = TRUE))), "\n")

cat("    Naive por anio: ")
for (t in 1:T_years) {
  cols_t <- ((t-1) * J_reps + 1):(t * J_reps)
  naive <- mean(apply(y_filt[, cols_t], 1, function(x) any(x == 1, na.rm = TRUE)))
  cat(sprintf("%d:%.1f%% ", YEARS[t], naive * 100))
}
cat("\n")

###############################################################################
# PASO 7: Cargar y preparar variables lag-1 (yearlySiteCovs)
###############################################################################

cat("\n--- PASO 7: Cargar variables y construir lag-1 ---\n\n")

dyn_path <- file.path(
  "/Users/gfandos/Documents/GitHub/dynamic_occupancy_models_steppe_birds",
  "data/processed/otitar_dynamic_variables.csv"
)
dyn_raw <- read.csv(dyn_path)
cat("  Archivo: otitar_dynamic_variables.csv\n")
cat("  Filas: ", nrow(dyn_raw), "\n")
cat("  Columnas: ", ncol(dyn_raw), "\n\n")

# --- 7.1: Hacer match entre locality_id del pipeline y del CSV dinamico ---
match_idx <- match(locality_ids_filt, dyn_raw$locality_id)
n_matched <- sum(!is.na(match_idx))
n_missing <- sum(is.na(match_idx))

cat(sprintf("  Match locality_id:\n"))
cat(sprintf("    Sitios filtrados en pipeline: %d\n", n_filt))
cat(sprintf("    Sitios con datos dinamicos: %d (%.1f%%)\n",
            n_matched, 100 * n_matched / n_filt))
cat(sprintf("    Sitios sin datos dinamicos: %d (%.1f%%)\n",
            n_missing, 100 * n_missing / n_filt))

dyn_matched <- dyn_raw[match_idx, ]

# --- 7.2: Funcion para construir lag-1 ---
# En colext, yearlySiteCov del ano t predice transicion t -> t+1
# Lag-1 = valor del ano t-1 predice la misma transicion
# Implementacion: desplazar la matriz 1 columna a la derecha
#   lag_mat[, 1] <- mat[, 1]  (2017: sin lag real, usar mismo ano)
#   lag_mat[, 2] <- mat[, 1]  (2018: usar valor de 2017)
#   lag_mat[, 3] <- mat[, 2]  (2019: usar valor de 2018)
#   ...etc.

make_lag1 <- function(mat) {
  lag_mat <- matrix(NA, nrow = nrow(mat), ncol = ncol(mat))
  lag_mat[, 1] <- mat[, 1]  # Ano 1: sin lag disponible, usar mismo valor
  for (t in 2:ncol(mat)) {
    lag_mat[, t] <- mat[, t - 1]
  }
  return(lag_mat)
}

# --- 7.3: Extraer variables climaticas ---
cat("  --- Variables climaticas ---\n")
climate_vars <- c("NDVI", "pr", "tmmx")

dyn_matrices_raw <- list()  # Matrices originales (sin lag)

for (dv in climate_vars) {
  mat <- matrix(NA, nrow = n_filt, ncol = T_years)
  for (t in seq_along(YEARS)) {
    col_name <- paste0(dv, "_", YEARS[t])
    if (col_name %in% names(dyn_matched)) {
      mat[, t] <- dyn_matched[[col_name]]
    }
  }

  # --- LIMPIEZA DE NAs (3 pasos, heredado de v4) ---
  n_na_raw <- sum(is.na(mat))

  # Paso A: Imputar NAs por media del sitio
  for (i in 1:nrow(mat)) {
    na_cols <- is.na(mat[i, ])
    if (any(na_cols) && !all(na_cols)) {
      mat[i, na_cols] <- mean(mat[i, !na_cols])
    }
  }
  # Paso B: Imputar NAs restantes con media global
  global_mean <- mean(mat, na.rm = TRUE)
  n_na_remaining <- sum(is.na(mat))
  if (n_na_remaining > 0) {
    mat[is.na(mat)] <- global_mean
    cat(sprintf("    %s: %d NAs imputados con media global (%.1f)\n",
                dv, n_na_remaining, global_mean))
  }
  # Paso C: Clamp outliers (percentil 1%)
  q01 <- quantile(mat, 0.01)
  n_clamped <- sum(mat < q01)
  if (n_clamped > 0) {
    mat[mat < q01] <- q01
    cat(sprintf("    %s: %d valores extremos clamped al P1%% (%.1f)\n",
                dv, n_clamped, q01))
  }

  dyn_matrices_raw[[dv]] <- mat

  cat(sprintf("    %s: %d x %d, NAs orig: %d, rango: [%.1f, %.1f]\n",
              dv, nrow(mat), ncol(mat), n_na_raw,
              min(mat), max(mat)))
}

# --- 7.4: Extraer variables de land cover ---
cat("\n  --- Variables de land cover ---\n")

lc_classes <- list(
  lc_grass = list(code = 10, name = "Grassland"),
  lc_crop  = list(code = 12, name = "Cropland")
)

for (lc_name in names(lc_classes)) {
  lc_code <- lc_classes[[lc_name]]$code
  mat <- matrix(NA, nrow = n_filt, ncol = T_years)
  for (t in seq_along(YEARS)) {
    col_name <- paste0("Land_Cover_Type_1_Percent_Class_", lc_code, "_", YEARS[t])
    if (col_name %in% names(dyn_matched)) {
      mat[, t] <- dyn_matched[[col_name]]
    }
  }
  # Land cover: 0 NAs confirmado
  dyn_matrices_raw[[lc_name]] <- mat

  cat(sprintf("    %s (IGBP %d): %d x %d, NAs: %d, rango: [%.1f, %.1f]\n",
              lc_name, lc_code, nrow(mat), ncol(mat), sum(is.na(mat)),
              min(mat, na.rm = TRUE), max(mat, na.rm = TRUE)))
}

# --- 7.5: Construir versiones lag-1 ---
cat("\n  --- Construyendo matrices lag-1 ---\n")
cat("    Logica: lag_mat[,t] = mat[,t-1] para t>=2; lag_mat[,1] = mat[,1]\n\n")

all_lag_vars <- c("NDVI", "pr", "tmmx", "lc_grass", "lc_crop")
lag_names <- paste0(all_lag_vars, "_lag")

dyn_matrices <- list()   # Matrices escaladas para el modelo
dyn_scaling <- list()     # Parametros de escalado

for (i in seq_along(all_lag_vars)) {
  dv <- all_lag_vars[i]
  lag_name <- lag_names[i]
  mat_raw <- dyn_matrices_raw[[dv]]

  # Construir lag-1
  mat_lag <- make_lag1(mat_raw)

  # Escalar la version lag
  mu <- mean(mat_lag, na.rm = TRUE)
  sd_val <- sd(mat_lag, na.rm = TRUE)
  if (sd_val > 0) {
    mat_lag_sc <- (mat_lag - mu) / sd_val
  } else {
    mat_lag_sc <- mat_lag * 0
    cat(sprintf("    AVISO: %s tiene varianza 0\n", lag_name))
  }
  dyn_scaling[[lag_name]] <- list(center = mu, scale = sd_val)
  dyn_matrices[[lag_name]] <- mat_lag_sc

  cat(sprintf("    %s -> %s: scaled range [%.2f, %.2f]\n",
              dv, lag_name, min(mat_lag_sc, na.rm = TRUE), max(mat_lag_sc, na.rm = TRUE)))
}

# Verificacion: la primera columna del lag debe ser igual a la original
cat("\n  Verificacion lag-1:\n")
for (i in seq_along(all_lag_vars)) {
  dv <- all_lag_vars[i]
  lag_name <- lag_names[i]
  # Comparar ano 2 del lag con ano 1 del original (ambos deben ser valor 2017)
  raw_yr1 <- dyn_matrices_raw[[dv]][1:3, 1]  # primeros 3 sitios, ano 1
  lag_yr2 <- make_lag1(dyn_matrices_raw[[dv]])[1:3, 2]  # primeros 3 sitios, lag ano 2
  match_ok <- all(raw_yr1 == lag_yr2, na.rm = TRUE)
  cat(sprintf("    %s: lag[,2] == orig[,1]? %s\n", lag_name,
              ifelse(match_ok, "OK", "ERROR")))
}

# --- 7.6: Convertir a formato yearlySiteCovs ---
yearly_df <- data.frame(row.names = 1:(n_filt * T_years))
for (lag_name in lag_names) {
  yearly_df[[lag_name]] <- as.vector(dyn_matrices[[lag_name]])
}

cat(sprintf("\n  yearlySiteCovs construido: %d filas (%d sitios x %d anios), %d columnas\n",
            nrow(yearly_df), n_filt, T_years, ncol(yearly_df)))
cat(sprintf("  Variables: %s\n", paste(names(yearly_df), collapse = ", ")))

###############################################################################
# PASO 8: Construir unmarkedMultFrame con yearlySiteCovs
###############################################################################

cat("\n--- PASO 8: Construir unmarkedMultFrame con yearlySiteCovs ---\n\n")

duration_sc  <- scale(duration_filt)
effort_sc    <- scale(effort_filt)
observers_sc <- scale(observers_filt)
time_sc      <- scale(time_filt)

siteCovs_df <- as.data.frame(siteCovs_filt)

occ_umf <- unmarkedMultFrame(
  y = y_filt,
  siteCovs = siteCovs_df,
  yearlySiteCovs = yearly_df,
  obsCovs = list(
    duration  = duration_sc,
    effort    = effort_sc,
    observers = observers_sc,
    time      = time_sc
  ),
  numPrimary = T_years
)

cat("  UMF construido: ", n_filt, " sitios, ",
    T_years, " periodos primarios, ", J_reps, " replicas\n")
cat("  yearlySiteCovs: ", paste(names(yearly_df), collapse = ", "), "\n")
cat("  siteCovs: ", paste(names(siteCovs_df), collapse = ", "), "\n\n")

###############################################################################
# PASO 9: Ajustar modelos (baseline + lag-1)
###############################################################################

cat("--- PASO 9: Ajustar 12 modelos ---\n\n")

check_model <- function(mod, model_name) {
  if (is.null(mod)) return(list(ok = FALSE, reason = "Error en ajuste"))
  if (mod@opt$convergence != 0) return(list(ok = FALSE, reason = "No converge"))
  s <- tryCatch(summary(mod), error = function(e) NULL)
  if (is.null(s)) return(list(ok = FALSE, reason = "Hessiana singular"))
  coefs <- coef(mod)
  gep_coefs <- coefs[grepl("^(col|ext|p)\\(", names(coefs))]
  boundary <- names(gep_coefs)[abs(gep_coefs) > 10]
  if (length(boundary) > 0)
    return(list(ok = FALSE, reason = paste("Boundary:", paste(boundary, collapse = ", "))))
  psi_coefs <- coefs[grepl("^psi\\(", names(coefs))]
  psi_boundary <- names(psi_coefs)[abs(psi_coefs) > 10]
  if (length(psi_boundary) > 0)
    return(list(ok = TRUE, reason = paste("OK (psi boundary:", paste(psi_boundary, collapse = ", "), ")")))
  return(list(ok = TRUE, reason = "OK"))
}

fit_model <- function(psi_f, gam_f, eps_f, p_f, data, label) {
  cat(sprintf("  %s\n", label))
  mod <- tryCatch({
    colext(psiformula = psi_f, gammaformula = gam_f,
           epsilonformula = eps_f, pformula = p_f, data = data)
  }, error = function(e) {
    cat(sprintf("    Error: %s\n", e$message))
    NULL
  })
  if (!is.null(mod)) {
    chk <- check_model(mod, label)
    cat(sprintf("    -> AIC=%.1f, %s\n", mod@AIC, chk$reason))
  }
  return(mod)
}

psi_full <- ~bio1 + bio2 + tree_cover + grass_cover + topo_elev
p_full   <- ~effort + observers + duration + time

# === MODELO 0: Baseline estatico (mejor del v3) ===
cat("  [m0] BASELINE: eps ~ grass_cover + tree_cover (mejor v3)\n")
m0 <- fit_model(psi_full, ~1, ~grass_cover + tree_cover, p_full, occ_umf,
                "[m0] Baseline estatico")

# === MODELOS CON EPSILON LAGGED ===
cat("\n  --- Epsilon con lag-1 ---\n")

m1 <- fit_model(psi_full, ~1, ~tree_cover + NDVI_lag, p_full, occ_umf,
                "[m1] eps ~ tree_cover + NDVI_lag")

m2 <- fit_model(psi_full, ~1, ~tree_cover + pr_lag, p_full, occ_umf,
                "[m2] eps ~ tree_cover + pr_lag")

m3 <- fit_model(psi_full, ~1, ~tree_cover + tmmx_lag, p_full, occ_umf,
                "[m3] eps ~ tree_cover + tmmx_lag")

m4 <- fit_model(psi_full, ~1, ~tree_cover + lc_grass_lag, p_full, occ_umf,
                "[m4] eps ~ tree_cover + lc_grass_lag")

# === MODELOS CON GAMMA LAGGED ===
cat("\n  --- Gamma con lag-1 ---\n")

m5 <- fit_model(psi_full, ~NDVI_lag, ~grass_cover + tree_cover, p_full, occ_umf,
                "[m5] gam ~ NDVI_lag")

m6 <- fit_model(psi_full, ~pr_lag, ~grass_cover + tree_cover, p_full, occ_umf,
                "[m6] gam ~ pr_lag")

m7 <- fit_model(psi_full, ~lc_grass_lag, ~grass_cover + tree_cover, p_full, occ_umf,
                "[m7] gam ~ lc_grass_lag")

# === MODELOS CON AMBOS LAGGED ===
cat("\n  --- Gamma + Epsilon con lag-1 ---\n")

m8 <- fit_model(psi_full, ~pr_lag, ~tree_cover + pr_lag, p_full, occ_umf,
                "[m8] gam ~ pr_lag, eps ~ tree_cover + pr_lag")

m9 <- fit_model(psi_full, ~NDVI_lag, ~tree_cover + NDVI_lag, p_full, occ_umf,
                "[m9] gam ~ NDVI_lag, eps ~ tree_cover + NDVI_lag")

m10 <- fit_model(psi_full, ~lc_grass_lag, ~tree_cover + lc_grass_lag, p_full, occ_umf,
                 "[m10] gam ~ lc_grass_lag, eps ~ tree_cover + lc_grass_lag")

# === MODELO COMBINADO ===
cat("\n  --- Modelo combinado lag ---\n")

m11 <- fit_model(psi_full, ~pr_lag + lc_grass_lag,
                 ~tree_cover + NDVI_lag + lc_grass_lag, p_full, occ_umf,
                 "[m11] gam ~ pr_lag+lc_grass_lag, eps ~ tree+NDVI_lag+lc_grass_lag")

###############################################################################
# PASO 10: Seleccion de modelos y diagnosticos
###############################################################################

cat("\n--- PASO 10: Seleccion de modelos ---\n\n")

models <- list(
  "m0: Baseline estatico"         = m0,
  "m1: eps~tree+NDVI_lag"         = m1,
  "m2: eps~tree+pr_lag"           = m2,
  "m3: eps~tree+tmmx_lag"         = m3,
  "m4: eps~tree+lc_grass_lag"     = m4,
  "m5: gam~NDVI_lag"              = m5,
  "m6: gam~pr_lag"                = m6,
  "m7: gam~lc_grass_lag"          = m7,
  "m8: ambos pr_lag"              = m8,
  "m9: ambos NDVI_lag"            = m9,
  "m10: ambos lc_grass_lag"       = m10,
  "m11: combinado lag"            = m11
)

model_formulas <- list(
  "m0: Baseline estatico"         = "psi(full) gam(~1) eps(~grass+tree) p(full)",
  "m1: eps~tree+NDVI_lag"         = "psi(full) gam(~1) eps(~tree+NDVI_lag) p(full)",
  "m2: eps~tree+pr_lag"           = "psi(full) gam(~1) eps(~tree+pr_lag) p(full)",
  "m3: eps~tree+tmmx_lag"         = "psi(full) gam(~1) eps(~tree+tmmx_lag) p(full)",
  "m4: eps~tree+lc_grass_lag"     = "psi(full) gam(~1) eps(~tree+lc_grass_lag) p(full)",
  "m5: gam~NDVI_lag"              = "psi(full) gam(~NDVI_lag) eps(~grass+tree) p(full)",
  "m6: gam~pr_lag"                = "psi(full) gam(~pr_lag) eps(~grass+tree) p(full)",
  "m7: gam~lc_grass_lag"          = "psi(full) gam(~lc_grass_lag) eps(~grass+tree) p(full)",
  "m8: ambos pr_lag"              = "psi(full) gam(~pr_lag) eps(~tree+pr_lag) p(full)",
  "m9: ambos NDVI_lag"            = "psi(full) gam(~NDVI_lag) eps(~tree+NDVI_lag) p(full)",
  "m10: ambos lc_grass_lag"       = "psi(full) gam(~lc_grass_lag) eps(~tree+lc_grass_lag) p(full)",
  "m11: combinado lag"            = "psi(full) gam(~pr_lag+lc_grass_lag) eps(~tree+NDVI_lag+lc_grass_lag) p(full)"
)

# 10.1: Tabla AIC
cat("  === TABLA AIC ===\n\n")

aic_table <- data.frame(
  Model = character(), Formula = character(), nPars = integer(),
  AIC = numeric(), Conv = integer(), Status = character(),
  stringsAsFactors = FALSE
)

for (mname in names(models)) {
  mod <- models[[mname]]
  chk <- check_model(mod, mname)
  if (!is.null(mod)) {
    aic_table <- rbind(aic_table, data.frame(
      Model = mname, Formula = model_formulas[[mname]],
      nPars = length(coef(mod)), AIC = round(mod@AIC, 1),
      Conv = mod@opt$convergence, Status = chk$reason,
      stringsAsFactors = FALSE))
  } else {
    aic_table <- rbind(aic_table, data.frame(
      Model = mname, Formula = model_formulas[[mname]],
      nPars = NA, AIC = NA, Conv = NA, Status = "Error en ajuste",
      stringsAsFactors = FALSE))
  }
}

valid_rows <- !is.na(aic_table$AIC) & grepl("^OK", aic_table$Status)
if (any(valid_rows)) {
  min_aic <- min(aic_table$AIC[valid_rows])
  aic_table$deltaAIC <- ifelse(!is.na(aic_table$AIC), round(aic_table$AIC - min_aic, 1), NA)
  delta_valid <- aic_table$deltaAIC[valid_rows]
  weights_raw <- exp(-0.5 * delta_valid)
  aic_table$Weight <- NA
  aic_table$Weight[valid_rows] <- round(weights_raw / sum(weights_raw), 4)
}

aic_table <- aic_table[order(aic_table$AIC, na.last = TRUE), ]

cat(sprintf("  %-35s %5s %7s %7s %6s  %-30s\n",
            "Modelo", "nPar", "AIC", "dAIC", "Peso", "Estado"))
cat(strrep("-", 100), "\n")
for (i in 1:nrow(aic_table)) {
  cat(sprintf("  %-35s %5s %7s %7s %6s  %-30s\n",
              aic_table$Model[i],
              ifelse(is.na(aic_table$nPars[i]), " -", as.character(aic_table$nPars[i])),
              ifelse(is.na(aic_table$AIC[i]), "   -", sprintf("%.1f", aic_table$AIC[i])),
              ifelse(is.na(aic_table$deltaAIC[i]), "   -", sprintf("%.1f", aic_table$deltaAIC[i])),
              ifelse(is.na(aic_table$Weight[i]), "  -", sprintf("%.4f", aic_table$Weight[i])),
              aic_table$Status[i]))
}

# 10.2: Seleccionar mejor modelo
best_name <- NULL; best_mod <- NULL; best_aic_val <- Inf
for (mname in names(models)) {
  mod <- models[[mname]]; chk <- check_model(mod, mname)
  if (chk$ok && !is.null(mod) && mod@AIC < best_aic_val) {
    best_aic_val <- mod@AIC; best_name <- mname; best_mod <- mod
  }
}
if (is.null(best_mod)) {
  for (mname in names(models)) {
    mod <- models[[mname]]
    if (!is.null(mod) && mod@opt$convergence == 0 && mod@AIC < best_aic_val) {
      best_aic_val <- mod@AIC; best_name <- mname; best_mod <- mod
    }
  }
}

cat(sprintf("\n  >>> MEJOR MODELO: %s (AIC=%.1f) <<<\n", best_name, best_aic_val))

baseline_aic <- ifelse(!is.null(m0), m0@AIC, NA)
if (!is.na(baseline_aic)) {
  delta_vs_baseline <- best_aic_val - baseline_aic
  cat(sprintf("  >>> dAIC vs baseline estatico: %.1f", delta_vs_baseline))
  if (delta_vs_baseline < -2) cat(" (MEJORA SUSTANCIAL)\n")
  else if (delta_vs_baseline < 2) cat(" (SIMILAR al baseline)\n")
  else cat(" (PEOR que baseline)\n")
}

# 10.3: Parametros del mejor modelo
if (!is.null(best_mod)) {
  cat("\n  === PARAMETROS DEL MEJOR MODELO ===\n\n")
  coefs_best <- coef(best_mod)
  cat("  Interceptos (escala probabilidad):\n")
  cat(sprintf("    psi1:    logit=%.3f -> %.2f%%\n", coefs_best["psi(Int)"], plogis(coefs_best["psi(Int)"]) * 100))
  cat(sprintf("    gamma:   logit=%.3f -> %.2f%%\n", coefs_best["col(Int)"], plogis(coefs_best["col(Int)"]) * 100))
  cat(sprintf("    epsilon: logit=%.3f -> %.2f%%\n", coefs_best["ext(Int)"], plogis(coefs_best["ext(Int)"]) * 100))
  cat(sprintf("    p:       logit=%.3f -> %.2f%%\n", coefs_best["p(Int)"], plogis(coefs_best["p(Int)"]) * 100))

  s_best <- tryCatch(summary(best_mod), error = function(e) NULL)
  if (!is.null(s_best)) {
    cat("\n  psi (ocupacion inicial):\n"); print(s_best$psi)
    cat("\n  col (colonizacion):\n"); print(s_best$col)
    cat("\n  ext (extincion):\n"); print(s_best$ext)
    cat("\n  det (deteccion):\n"); print(s_best$det)
  }

  # 10.5: Boundary detection
  cat("\n  === BOUNDARY DETECTION ===\n\n")
  for (mname in names(models)) {
    mod <- models[[mname]]
    if (is.null(mod)) next
    coefs <- coef(mod)
    boundary <- names(coefs)[abs(coefs) > 10]
    if (length(boundary) > 0)
      cat(sprintf("  %s: BOUNDARY en %s\n", mname,
                  paste(boundary, "=", round(coefs[boundary], 1), collapse = ", ")))
  }

  # 10.6: Numero de condicion
  cat("\n  === NUMERO DE CONDICION (mejor modelo) ===\n\n")
  tryCatch({
    vcov_mat <- vcov(best_mod)
    eig <- eigen(vcov_mat)$values
    cond_num <- max(abs(eig)) / min(abs(eig))
    cat(sprintf("  Numero de condicion: %.1f\n", cond_num))
    if (cond_num > 1e6) cat("  AVISO: Muy alto\n")
    else if (cond_num > 1e4) cat("  NOTA: Moderado\n")
    else cat("  OK: Aceptable\n")
  }, error = function(e) cat("  No se pudo calcular\n"))
}

###############################################################################
# PASO 11: Comparacion ocupacion predicha vs naive por anio
###############################################################################

cat("\n\n--- PASO 11: Ocupacion predicha vs naive por anio ---\n\n")

if (!is.null(best_mod)) {
  tryCatch({
    proj <- projected(best_mod)
    cat(sprintf("  %-6s %10s %10s %10s\n", "Anio", "Naive(%)", "Pred(%)", "Diff"))
    cat(strrep("-", 45), "\n")
    for (t in 1:T_years) {
      cols_t <- ((t-1) * J_reps + 1):(t * J_reps)
      naive <- mean(apply(y_filt[, cols_t], 1, function(x) any(x == 1, na.rm = TRUE))) * 100
      pred_occ <- if (is.list(proj)) mean(proj[[t]][, 2], na.rm = TRUE) * 100 else NA
      cat(sprintf("  %d   %8.1f%%  %8.1f%%  %+7.1f%%\n",
                  YEARS[t], naive, ifelse(is.na(pred_occ), NA, pred_occ),
                  ifelse(is.na(pred_occ), NA, pred_occ - naive)))
    }
  }, error = function(e) cat("  No se pudo calcular projected(): ", e$message, "\n"))
}

###############################################################################
# PASO 12: Exportar resultados
###############################################################################

cat("\n--- PASO 12: Exportar resultados ---\n\n")

results_dir <- file.path(
  "/Users/gfandos/Documents/GitHub/dynamic_occupancy_models_steppe_birds", "results"
)
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

model_sel_path <- file.path(results_dir, "otitar_v5_model_selection.csv")
write.csv(aic_table, model_sel_path, row.names = FALSE)
cat("  Tabla AIC guardada: ", model_sel_path, "\n")

if (!is.null(best_mod) && !is.null(s_best)) {
  coef_list <- list()
  for (sm_name in c("psi", "col", "ext", "det")) {
    sm <- s_best[[sm_name]]
    if (is.null(sm)) next
    sm_df <- as.data.frame(sm)
    sm_df$Parameter <- rownames(sm)
    sm_df$Submodel <- sm_name
    sm_df$Prob <- plogis(sm_df$Estimate)
    coef_list[[sm_name]] <- sm_df
  }
  coef_df <- do.call(rbind, coef_list)
  rownames(coef_df) <- NULL
  coef_path <- file.path(results_dir, "otitar_v5_best_model_coefficients.csv")
  write.csv(coef_df, coef_path, row.names = FALSE)
  cat("  Coeficientes guardados: ", coef_path, "\n")
}

scaling_df <- data.frame(
  variable = names(dyn_scaling),
  center = sapply(dyn_scaling, function(x) x$center),
  scale = sapply(dyn_scaling, function(x) x$scale),
  row.names = NULL
)
scaling_path <- file.path(results_dir, "otitar_v5_dynamic_scaling_params.csv")
write.csv(scaling_df, scaling_path, row.names = FALSE)
cat("  Parametros escalado: ", scaling_path, "\n")

###############################################################################
# PASO 13: Goodness-of-fit (parboot) del mejor modelo
###############################################################################

cat("\n--- PASO 13: Goodness-of-fit (parboot, nsim=50) ---\n\n")
cat("  NOTA: nsim=50 es bajo, solo para verificar. Para publicacion usar nsim>=500\n\n")

if (!is.null(best_mod)) {
  tryCatch({
    chisq_fn <- function(fm) {
      obs <- getY(fm@data); ex <- fitted(fm)
      ex[ex == 0] <- 1e-6
      return(c(chisq = sum((obs - ex)^2 / ex, na.rm = TRUE)))
    }
    pb <- parboot(best_mod, statistic = chisq_fn, nsim = 50, report = 10)
    cat("  Resultado parboot:\n")
    cat(sprintf("    Chi-cuadrado observado: %.1f\n", pb@t0))
    cat(sprintf("    Chi-cuadrado simulado (media): %.1f\n", mean(pb@t.star[,1])))
    p_gof <- mean(pb@t.star[,1] >= pb@t0)
    cat(sprintf("    p-value: %.3f\n", p_gof))
    if (p_gof > 0.05) cat("    => Modelo OK\n")
    else cat("    => AVISO: Posible falta de ajuste\n")
  }, error = function(e) cat("  Error en parboot: ", e$message, "\n"))
}

###############################################################################
# RESUMEN FINAL
###############################################################################

cat("\n\n")
cat(strrep("=", 70), "\n")
cat("  PIPELINE V5: EFECTOS RETARDADOS (LAG-1) — RESUMEN\n")
cat(strrep("=", 70), "\n\n")

if (!is.null(best_mod)) {
  cat(sprintf("  Mejor modelo: %s\n", best_name))
  cat(sprintf("  Formula: %s\n", model_formulas[[best_name]]))
  cat(sprintf("  AIC: %.1f\n", best_aic_val))
  if (!is.na(baseline_aic)) {
    cat(sprintf("  dAIC vs baseline estatico (v3): %.1f\n", best_aic_val - baseline_aic))
    if (best_aic_val - baseline_aic < -2) cat("  => Efectos retardados MEJORAN el modelo\n")
    else cat("  => Diferencia no sustancial con el baseline estatico\n")
  }
  coefs_best <- coef(best_mod)
  cat(sprintf("\n  Interceptos:\n"))
  cat(sprintf("    psi1: %.1f%%\n", plogis(coefs_best["psi(Int)"]) * 100))
  cat(sprintf("    gamma: %.2f%%\n", plogis(coefs_best["col(Int)"]) * 100))
  cat(sprintf("    epsilon: %.1f%%\n", plogis(coefs_best["ext(Int)"]) * 100))
  cat(sprintf("    p: %.1f%%\n", plogis(coefs_best["p(Int)"]) * 100))
}

cat("\n  Archivos generados:\n")
cat("    - results/otitar_v5_model_selection.csv\n")
cat("    - results/otitar_v5_best_model_coefficients.csv\n")
cat("    - results/otitar_v5_dynamic_scaling_params.csv\n")

cat("\n")
cat(strrep("=", 70), "\n")
cat("  FIN PIPELINE V5\n")
cat(strrep("=", 70), "\n\n")
