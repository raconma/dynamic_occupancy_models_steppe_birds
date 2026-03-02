###############################################################################
# test_pipeline_v4_dynamic_covs.R
#
# Pipeline v4 para Otis tarda: covariables dinamicas en gamma/epsilon.
#
# Extiende el pipeline v3 incorporando variables anuales (yearlySiteCovs)
# exportadas desde Google Earth Engine:
#   - NDVI anual: proxy de productividad y calidad del habitat
#   - EVI anual: similar al NDVI pero mas robusto en zonas de alta biomasa
#   - Precipitacion (pr): acumulada anual, driver hidrico clave para esteparias
#   - Temperatura maxima (tmmx): estres termico, afecta fenologia y cria
#   - Temperatura minima (tmmn): severidad invernal, puede afectar supervivencia
#
# Hipotesis ecologicas:
#   - Epsilon (extincion): anos secos (bajo NDVI/pr) -> mayor extincion local
#     porque el habitat estepario pierde calidad temporal. El tree_cover
#     (estatico) captura la idoneidad estructural del habitat.
#   - Gamma (colonizacion): anos productivos (alto NDVI/pr) -> mas colonizacion
#     porque la abundancia de recursos facilita la expansion.
#
# Estructura del pipeline:
#   Pasos 1-6: identicos al pipeline v3 (datos, filtro, UMF base)
#   Paso 7: Carga y preparacion de yearlySiteCovs desde GEE
#   Paso 8: Construccion del unmarkedMultFrame con yearlySiteCovs
#   Paso 9: Bateria de 11 modelos (baseline estatico + 10 dinamicos)
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
#   - results/otitar_v4_model_selection.csv (tabla AIC 11 modelos)
#   - results/otitar_v4_best_model_coefficients.csv (coefs mejor modelo)
#   - results/otitar_v4_dynamic_scaling_params.csv (params escalado)
#
# Limpieza de datos NDVI/EVI:
#   Los datos NDVI/EVI de GEE contienen valores negativos (hasta -1594)
#   y NAs (39 de 11928, 0.3%). El optimizador de colext falla con NAs
#   en yearlySiteCovs ("valor inicial en vmmin no es finito"). Solucion:
#     1. Imputar NAs por media del sitio (otros anios)
#     2. Imputar NAs residuales con media global (24 sitios sin datos)
#     3. Clamp extremos < percentil 1% (NDVI negativo -> P1%)
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
cat("  PIPELINE v4: Otis tarda - Covariables dinamicas\n")
cat("  Modelos con yearlySiteCovs (NDVI, EVI, precip, temp) en gamma/epsilon\n")
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

# Estandarizar covariables de sitio (guardar parametros para referencia)
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

# Deduplicar por localidad
occ_wide_clean <- occ_wide %>%
  distinct(locality_id, .keep_all = TRUE)

cat("  Sitios totales (localidades unicas): ", nrow(occ_wide_clean), "\n\n")

###############################################################################
# PASO 5: Construir matrices de deteccion y covariables de observacion
###############################################################################

cat("--- PASO 5: Construir matrices ---\n\n")

J_reps <- 10  # max replicas por periodo primario

det_pattern <- "^y\\."
det_cols_all <- grep(det_pattern, names(occ_wide_clean), value = TRUE)
cat("  Detection columns found: ", length(det_cols_all), "\n")

# Organizar las columnas por anio
y_list <- list()
duration_list <- list()
effort_list <- list()
observers_list <- list()
time_list <- list()

for (yr in YEARS) {
  yr_suffix <- paste0("\\.", yr, "$")

  # Detection columns: y.{rep}.{year}
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
locality_ids_filt <- occ_wide_clean$locality_id[keep]

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
# PASO 7: Cargar y preparar variables dinamicas (yearlySiteCovs)
###############################################################################

cat("\n--- PASO 7: Cargar variables dinamicas desde GEE ---\n\n")

dyn_path <- file.path(
  "/Users/gfandos/Documents/GitHub/dynamic_occupancy_models_steppe_birds",
  "data/processed/otitar_dynamic_variables.csv"
)
dyn_raw <- read.csv(dyn_path)
cat("  Archivo: otitar_dynamic_variables.csv\n")
cat("  Filas: ", nrow(dyn_raw), "\n")
cat("  Columnas: ", ncol(dyn_raw), "\n\n")

# Variables dinamicas disponibles por anio (columnas en formato VAR_YYYY)
dyn_vars <- c("NDVI", "EVI", "pr", "tmmn", "tmmx")
cat("  Variables dinamicas disponibles:\n")
for (dv in dyn_vars) {
  cols_dv <- grep(paste0("^", dv, "_\\d{4}$"), names(dyn_raw), value = TRUE)
  if (length(cols_dv) > 0) {
    vals <- unlist(dyn_raw[, cols_dv])
    cat(sprintf("    %s: %d columnas, rango [%.1f, %.1f], NAs: %d (%.1f%%)\n",
                dv, length(cols_dv),
                min(vals, na.rm = TRUE), max(vals, na.rm = TRUE),
                sum(is.na(vals)), 100 * mean(is.na(vals))))
  } else {
    cat(sprintf("    %s: NO ENCONTRADA\n", dv))
  }
}

# --- 7.1: Hacer match entre locality_id del pipeline y del CSV dinamico ---
# El CSV de GEE tiene locality_id como identificador

# Comprobar match
match_idx <- match(locality_ids_filt, dyn_raw$locality_id)
n_matched <- sum(!is.na(match_idx))
n_missing <- sum(is.na(match_idx))

cat(sprintf("\n  Match locality_id:\n"))
cat(sprintf("    Sitios filtrados en pipeline: %d\n", n_filt))
cat(sprintf("    Sitios con datos dinamicos: %d (%.1f%%)\n",
            n_matched, 100 * n_matched / n_filt))
cat(sprintf("    Sitios sin datos dinamicos: %d (%.1f%%)\n",
            n_missing, 100 * n_missing / n_filt))

# Si hay sitios sin match, los conservamos con NAs en las dinamicas
# (unmarked maneja NAs en yearlySiteCovs)

# --- 7.2: Construir matrices [n_filt x T_years] para cada variable ---
# Estas matrices tienen la misma dimension que las transiciones: n_sites filas,
# T_years columnas (una por anio). En colext, las yearlySiteCovs se aplican
# a las transiciones ENTRE anios, pero la convencion de unmarked es:
# - yearlySiteCovs tiene M * T filas (M sitios x T anios)
# - Orden: sitio_1_anio_1, sitio_2_anio_1, ..., sitio_M_anio_1,
#           sitio_1_anio_2, ..., sitio_M_anio_T

cat("\n  Construyendo matrices de yearlySiteCovs...\n")

# Extraer datos dinamicos para los sitios filtrados (manteniendo orden)
dyn_matched <- dyn_raw[match_idx, ]  # NAs donde no hay match

# Construir una lista de matrices [n_filt x T_years] por variable
dyn_matrices <- list()
dyn_scaling <- list()  # Guardar parametros de escalado

for (dv in dyn_vars) {
  mat <- matrix(NA, nrow = n_filt, ncol = T_years)
  for (t in seq_along(YEARS)) {
    col_name <- paste0(dv, "_", YEARS[t])
    if (col_name %in% names(dyn_matched)) {
      mat[, t] <- dyn_matched[[col_name]]
    }
  }

  # Diagnostico de NAs antes de limpieza
  n_na_raw <- sum(is.na(mat))
  pct_na_raw <- 100 * n_na_raw / length(mat)

  # --- LIMPIEZA DE NAs ---
  # Paso A: Imputar NAs por media del sitio (otros anios del mismo sitio)
  for (i in 1:nrow(mat)) {
    na_cols <- is.na(mat[i, ])
    if (any(na_cols) && !all(na_cols)) {
      mat[i, na_cols] <- mean(mat[i, !na_cols])
    }
  }
  # Paso B: Imputar NAs restantes (sitios sin ningun dato) con media global
  # Estos NAs causan "valor inicial en vmmin no es finito" en colext
  global_mean <- mean(mat, na.rm = TRUE)
  n_na_remaining <- sum(is.na(mat))
  if (n_na_remaining > 0) {
    mat[is.na(mat)] <- global_mean
    cat(sprintf("    %s: %d NAs imputados con media global (%.1f)\n",
                dv, n_na_remaining, global_mean))
  }

  # --- LIMPIEZA DE EXTREMOS ---
  # Clamp valores por debajo del percentil 1% (outliers como NDVI < 0)
  q01 <- quantile(mat, 0.01)
  n_clamped <- sum(mat < q01)
  if (n_clamped > 0) {
    mat[mat < q01] <- q01
    cat(sprintf("    %s: %d valores extremos clamped al P1%% (%.1f)\n",
                dv, n_clamped, q01))
  }

  # Escalar (centrar y dividir por sd)
  mu <- mean(mat)
  sd_val <- sd(mat)
  if (sd_val > 0) {
    mat_sc <- (mat - mu) / sd_val
  } else {
    mat_sc <- mat * 0
    cat(sprintf("    AVISO: %s tiene varianza 0, excluida\n", dv))
  }
  dyn_scaling[[dv]] <- list(center = mu, scale = sd_val)

  dyn_matrices[[dv]] <- mat_sc

  cat(sprintf("    %s: %d x %d, NAs orig: %d (%.1f%%), scaled range: [%.2f, %.2f]\n",
              dv, nrow(mat), ncol(mat), n_na_raw, pct_na_raw,
              min(mat_sc), max(mat_sc)))
}

# --- 7.3: Convertir matrices a formato yearlySiteCovs para unmarked ---
# unmarkedMultFrame espera un data.frame con M*T filas
# Orden: todos los sitios del anio 1, luego del anio 2, etc.

yearly_df <- data.frame(row.names = 1:(n_filt * T_years))
for (dv in dyn_vars) {
  # Apilar columnas: columna 1 (anio 1) de todos los sitios, luego columna 2, etc.
  yearly_df[[dv]] <- as.vector(dyn_matrices[[dv]])
}

cat(sprintf("\n  yearlySiteCovs construido: %d filas (%d sitios x %d anios), %d columnas\n",
            nrow(yearly_df), n_filt, T_years, ncol(yearly_df)))

# Verificar que no hay variables con demasiados NAs
for (dv in dyn_vars) {
  pct_na <- 100 * mean(is.na(yearly_df[[dv]]))
  if (pct_na > 50) {
    cat(sprintf("  AVISO: %s tiene %.1f%% NAs - considerar excluir del modelo\n", dv, pct_na))
  }
}

###############################################################################
# PASO 8: Construir unmarkedMultFrame con yearlySiteCovs
###############################################################################

cat("\n--- PASO 8: Construir unmarkedMultFrame con yearlySiteCovs ---\n\n")

# Estandarizar covariables de observacion
duration_sc  <- scale(duration_filt)
effort_sc    <- scale(effort_filt)
observers_sc <- scale(observers_filt)
time_sc      <- scale(time_filt)

# Site covariates (ya estandarizadas en paso 4)
siteCovs_df <- as.data.frame(siteCovs_filt)

# Build UMF con yearlySiteCovs
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
# PASO 9: Ajustar modelos (baseline + dinamicos)
###############################################################################

cat("--- PASO 9: Ajustar 11 modelos ---\n\n")

# Funcion para verificar modelo (heredada de v3)
check_model <- function(mod, model_name) {
  if (is.null(mod)) {
    return(list(ok = FALSE, reason = "Error en ajuste"))
  }
  if (mod@opt$convergence != 0) {
    return(list(ok = FALSE, reason = "No converge"))
  }
  s <- tryCatch(summary(mod), error = function(e) NULL)
  if (is.null(s)) {
    return(list(ok = FALSE, reason = "Hessiana singular"))
  }
  # Boundary en gamma, epsilon y p (no psi — es normal en especies raras)
  coefs <- coef(mod)
  gep_coefs <- coefs[grepl("^(col|ext|p)\\(", names(coefs))]
  boundary <- names(gep_coefs)[abs(gep_coefs) > 10]
  if (length(boundary) > 0) {
    return(list(ok = FALSE,
                reason = paste("Boundary:", paste(boundary, collapse = ", ")),
                boundary_params = boundary))
  }
  psi_coefs <- coefs[grepl("^psi\\(", names(coefs))]
  psi_boundary <- names(psi_coefs)[abs(psi_coefs) > 10]
  if (length(psi_boundary) > 0) {
    return(list(ok = TRUE,
                reason = paste("OK (psi boundary:", paste(psi_boundary, collapse = ", "), ")")))
  }
  return(list(ok = TRUE, reason = "OK"))
}

# Funcion auxiliar para ajustar modelo con tryCatch
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

# --- Formulas comunes ---
psi_full <- ~bio1 + bio2 + tree_cover + grass_cover + topo_elev
p_full   <- ~effort + observers + duration + time

# === MODELO 0: Baseline estatico (mejor del v3) ===
cat("  [m0] BASELINE: eps ~ grass_cover + tree_cover (mejor v3)\n")
m0 <- fit_model(psi_full, ~1, ~grass_cover + tree_cover, p_full, occ_umf,
                "[m0] Baseline estatico")

# === MODELOS CON EPSILON DINAMICO ===
cat("\n  --- Epsilon dinamico ---\n")

m1 <- fit_model(psi_full, ~1, ~tree_cover + NDVI, p_full, occ_umf,
                "[m1] eps ~ tree_cover + NDVI")

m2 <- fit_model(psi_full, ~1, ~tree_cover + pr, p_full, occ_umf,
                "[m2] eps ~ tree_cover + pr")

m3 <- fit_model(psi_full, ~1, ~tree_cover + tmmx, p_full, occ_umf,
                "[m3] eps ~ tree_cover + tmmx")

m4 <- fit_model(psi_full, ~1, ~tree_cover + NDVI + pr, p_full, occ_umf,
                "[m4] eps ~ tree_cover + NDVI + pr")

# === MODELOS CON GAMMA DINAMICO ===
cat("\n  --- Gamma dinamico ---\n")

m5 <- fit_model(psi_full, ~NDVI, ~grass_cover + tree_cover, p_full, occ_umf,
                "[m5] gam ~ NDVI")

m6 <- fit_model(psi_full, ~pr, ~grass_cover + tree_cover, p_full, occ_umf,
                "[m6] gam ~ pr")

m7 <- fit_model(psi_full, ~NDVI + pr, ~grass_cover + tree_cover, p_full, occ_umf,
                "[m7] gam ~ NDVI + pr")

# === MODELOS CON AMBOS DINAMICOS ===
cat("\n  --- Gamma + Epsilon dinamicos ---\n")

m8 <- fit_model(psi_full, ~NDVI, ~tree_cover + NDVI, p_full, occ_umf,
                "[m8] gam ~ NDVI, eps ~ tree_cover + NDVI")

m9 <- fit_model(psi_full, ~pr, ~tree_cover + pr, p_full, occ_umf,
                "[m9] gam ~ pr, eps ~ tree_cover + pr")

m10 <- fit_model(psi_full, ~NDVI + pr, ~tree_cover + NDVI + pr, p_full, occ_umf,
                 "[m10] gam ~ NDVI+pr, eps ~ tree_cover + NDVI + pr")

###############################################################################
# PASO 10: Seleccion de modelos y diagnosticos
###############################################################################

cat("\n--- PASO 10: Seleccion de modelos ---\n\n")

models <- list(
  "m0: Baseline estatico"       = m0,
  "m1: eps~tree+NDVI"           = m1,
  "m2: eps~tree+pr"             = m2,
  "m3: eps~tree+tmmx"           = m3,
  "m4: eps~tree+NDVI+pr"        = m4,
  "m5: gam~NDVI"                = m5,
  "m6: gam~pr"                  = m6,
  "m7: gam~NDVI+pr"             = m7,
  "m8: gam~NDVI,eps~tree+NDVI"  = m8,
  "m9: gam~pr,eps~tree+pr"      = m9,
  "m10: gam+eps NDVI+pr"        = m10
)

model_formulas <- list(
  "m0: Baseline estatico"       = "psi(full) gam(~1) eps(~grass+tree) p(full)",
  "m1: eps~tree+NDVI"           = "psi(full) gam(~1) eps(~tree+NDVI) p(full)",
  "m2: eps~tree+pr"             = "psi(full) gam(~1) eps(~tree+pr) p(full)",
  "m3: eps~tree+tmmx"           = "psi(full) gam(~1) eps(~tree+tmmx) p(full)",
  "m4: eps~tree+NDVI+pr"        = "psi(full) gam(~1) eps(~tree+NDVI+pr) p(full)",
  "m5: gam~NDVI"                = "psi(full) gam(~NDVI) eps(~grass+tree) p(full)",
  "m6: gam~pr"                  = "psi(full) gam(~pr) eps(~grass+tree) p(full)",
  "m7: gam~NDVI+pr"             = "psi(full) gam(~NDVI+pr) eps(~grass+tree) p(full)",
  "m8: gam~NDVI,eps~tree+NDVI"  = "psi(full) gam(~NDVI) eps(~tree+NDVI) p(full)",
  "m9: gam~pr,eps~tree+pr"      = "psi(full) gam(~pr) eps(~tree+pr) p(full)",
  "m10: gam+eps NDVI+pr"        = "psi(full) gam(~NDVI+pr) eps(~tree+NDVI+pr) p(full)"
)

# 10.1: Tabla AIC
cat("  === TABLA AIC ===\n\n")

aic_table <- data.frame(
  Model = character(),
  Formula = character(),
  nPars = integer(),
  AIC = numeric(),
  Conv = integer(),
  Status = character(),
  stringsAsFactors = FALSE
)

for (mname in names(models)) {
  mod <- models[[mname]]
  chk <- check_model(mod, mname)

  if (!is.null(mod)) {
    aic_table <- rbind(aic_table, data.frame(
      Model = mname,
      Formula = model_formulas[[mname]],
      nPars = length(coef(mod)),
      AIC = round(mod@AIC, 1),
      Conv = mod@opt$convergence,
      Status = chk$reason,
      stringsAsFactors = FALSE
    ))
  } else {
    aic_table <- rbind(aic_table, data.frame(
      Model = mname,
      Formula = model_formulas[[mname]],
      nPars = NA,
      AIC = NA,
      Conv = NA,
      Status = "Error en ajuste",
      stringsAsFactors = FALSE
    ))
  }
}

# deltaAIC y pesos (solo modelos validos y OK)
valid_rows <- !is.na(aic_table$AIC) & grepl("^OK", aic_table$Status)
if (any(valid_rows)) {
  min_aic <- min(aic_table$AIC[valid_rows])
  aic_table$deltaAIC <- ifelse(!is.na(aic_table$AIC),
                                round(aic_table$AIC - min_aic, 1), NA)
  delta_valid <- aic_table$deltaAIC[valid_rows]
  weights_raw <- exp(-0.5 * delta_valid)
  aic_table$Weight <- NA
  aic_table$Weight[valid_rows] <- round(weights_raw / sum(weights_raw), 4)
} else {
  # Si ningun modelo pasa check, calcular deltaAIC sobre todos los convergentes
  valid_rows <- !is.na(aic_table$AIC)
  if (any(valid_rows)) {
    min_aic <- min(aic_table$AIC[valid_rows])
    aic_table$deltaAIC <- ifelse(valid_rows, round(aic_table$AIC - min_aic, 1), NA)
    delta_valid <- aic_table$deltaAIC[valid_rows]
    weights_raw <- exp(-0.5 * delta_valid)
    aic_table$Weight <- NA
    aic_table$Weight[valid_rows] <- round(weights_raw / sum(weights_raw), 4)
  }
}

# Ordenar por AIC
aic_table <- aic_table[order(aic_table$AIC, na.last = TRUE), ]

# Imprimir tabla
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

# 10.2: Seleccionar mejor modelo (convergente, sin boundary en gam/eps/p)
best_name <- NULL
best_mod <- NULL
best_aic_val <- Inf

for (mname in names(models)) {
  mod <- models[[mname]]
  chk <- check_model(mod, mname)
  if (chk$ok && !is.null(mod) && mod@AIC < best_aic_val) {
    best_aic_val <- mod@AIC
    best_name <- mname
    best_mod <- mod
  }
}

# Fallback: si ningun modelo pasa check, tomar el de menor AIC convergente
if (is.null(best_mod)) {
  cat("\n  AVISO: Ningun modelo paso controles. Seleccionando menor AIC convergente...\n")
  for (mname in names(models)) {
    mod <- models[[mname]]
    if (!is.null(mod) && mod@opt$convergence == 0 && mod@AIC < best_aic_val) {
      best_aic_val <- mod@AIC
      best_name <- mname
      best_mod <- mod
    }
  }
}

cat(sprintf("\n  >>> MEJOR MODELO: %s (AIC=%.1f) <<<\n", best_name, best_aic_val))

# AIC del baseline para comparacion
baseline_aic <- ifelse(!is.null(m0), m0@AIC, NA)
if (!is.na(baseline_aic)) {
  delta_vs_baseline <- best_aic_val - baseline_aic
  cat(sprintf("  >>> ΔAIC vs baseline estatico: %.1f", delta_vs_baseline))
  if (delta_vs_baseline < -2) {
    cat(" (MEJORA SUSTANCIAL)\n")
  } else if (delta_vs_baseline < 2) {
    cat(" (SIMILAR al baseline)\n")
  } else {
    cat(" (PEOR que baseline)\n")
  }
}

# 10.3: Parametros del mejor modelo
if (!is.null(best_mod)) {
  cat("\n  === PARAMETROS DEL MEJOR MODELO ===\n\n")

  coefs_best <- coef(best_mod)
  cat("  Interceptos (escala probabilidad):\n")
  cat(sprintf("    psi1:    logit=%.3f -> %.2f%%\n",
              coefs_best["psi(Int)"], plogis(coefs_best["psi(Int)"]) * 100))
  cat(sprintf("    gamma:   logit=%.3f -> %.2f%%\n",
              coefs_best["col(Int)"], plogis(coefs_best["col(Int)"]) * 100))
  cat(sprintf("    epsilon: logit=%.3f -> %.2f%%\n",
              coefs_best["ext(Int)"], plogis(coefs_best["ext(Int)"]) * 100))
  cat(sprintf("    p:       logit=%.3f -> %.2f%%\n",
              coefs_best["p(Int)"], plogis(coefs_best["p(Int)"]) * 100))

  # Summary completo
  s_best <- tryCatch(summary(best_mod), error = function(e) NULL)
  if (!is.null(s_best)) {
    cat("\n  psi (ocupacion inicial):\n"); print(s_best$psi)
    cat("\n  col (colonizacion):\n"); print(s_best$col)
    cat("\n  ext (extincion):\n"); print(s_best$ext)
    cat("\n  det (deteccion):\n"); print(s_best$det)
  }

  # 10.4: Interpretacion ecologica de las covariables dinamicas
  cat("\n  === INTERPRETACION ECOLOGICA ===\n\n")

  if (!is.null(s_best)) {
    # Extraer coeficientes con SE, z y p-value de cada submodelo
    submodels <- list(col = "Colonizacion (gamma)", ext = "Extincion (epsilon)")
    for (sm_name in names(submodels)) {
      sm <- s_best[[sm_name]]
      if (is.null(sm)) next

      cat(sprintf("  %s:\n", submodels[[sm_name]]))
      for (i in 1:nrow(sm)) {
        param_name <- rownames(sm)[i]
        est <- sm[i, "Estimate"]
        se  <- sm[i, "SE"]
        z   <- sm[i, "z"]
        p   <- sm[i, "P(>|z|)"]

        # Identificar si es variable dinamica
        is_dynamic <- any(sapply(dyn_vars, function(dv) grepl(dv, param_name)))

        if (param_name == "(Intercept)") {
          cat(sprintf("    Intercepto: logit=%.3f (prob=%.1f%%), SE=%.3f\n",
                      est, plogis(est) * 100, se))
        } else {
          # Direccion del efecto
          if (sm_name == "ext") {
            direction <- ifelse(est > 0,
                                "aumenta extincion (negativo para la especie)",
                                "reduce extincion (positivo para la especie)")
          } else {
            direction <- ifelse(est > 0,
                                "aumenta colonizacion (positivo para la especie)",
                                "reduce colonizacion (negativo para la especie)")
          }

          signif_label <- ifelse(p < 0.001, "***",
                          ifelse(p < 0.01, "**",
                          ifelse(p < 0.05, "*",
                          ifelse(p < 0.1, ".", "ns"))))

          cat(sprintf("    %s: coef=%.3f (SE=%.3f), z=%.2f, p=%.4f %s\n",
                      param_name, est, se, z, p, signif_label))
          cat(sprintf("      -> %s\n", direction))

          # Interpretacion especifica para variables dinamicas
          if (is_dynamic) {
            if (grepl("NDVI|EVI", param_name)) {
              if (sm_name == "ext" && est < 0) {
                cat("      Ecol: Anos con mayor productividad vegetal reducen la extincion local\n")
              } else if (sm_name == "ext" && est > 0) {
                cat("      Ecol: Anos con mayor productividad vegetal aumentan la extincion (inesperado)\n")
              } else if (sm_name == "col" && est > 0) {
                cat("      Ecol: Anos productivos facilitan la colonizacion de nuevos sitios\n")
              } else if (sm_name == "col" && est < 0) {
                cat("      Ecol: Anos productivos reducen la colonizacion (inesperado)\n")
              }
            }
            if (grepl("^pr$|^pr\\)", param_name) || param_name == "pr") {
              if (sm_name == "ext" && est < 0) {
                cat("      Ecol: Anos lluviosos reducen la extincion (mas recursos hidricos)\n")
              } else if (sm_name == "col" && est > 0) {
                cat("      Ecol: Anos lluviosos facilitan la colonizacion\n")
              }
            }
            if (grepl("tmmx", param_name)) {
              if (sm_name == "ext" && est > 0) {
                cat("      Ecol: Veranos mas calurosos aumentan la extincion (estres termico)\n")
              } else if (sm_name == "ext" && est < 0) {
                cat("      Ecol: Temperaturas altas reducen la extincion (inesperado para esteparias)\n")
              }
            }
          }
        }
      }
      cat("\n")
    }
  }

  # 10.5: Boundary detection en todos los modelos
  cat("  === BOUNDARY DETECTION ===\n\n")
  for (mname in names(models)) {
    mod <- models[[mname]]
    if (is.null(mod)) next
    coefs <- coef(mod)
    boundary <- names(coefs)[abs(coefs) > 10]
    if (length(boundary) > 0) {
      cat(sprintf("  %s: BOUNDARY en %s\n", mname,
                  paste(boundary, "=", round(coefs[boundary], 1), collapse = ", ")))
    }
  }

  # 10.6: Numero de condicion
  cat("\n  === NUMERO DE CONDICION (mejor modelo) ===\n\n")
  tryCatch({
    vcov_mat <- vcov(best_mod)
    eig <- eigen(vcov_mat)$values
    cond_num <- max(abs(eig)) / min(abs(eig))
    cat(sprintf("  Numero de condicion: %.1f\n", cond_num))
    if (cond_num > 1e6) {
      cat("  AVISO: Muy alto - posible colinealidad severa\n")
    } else if (cond_num > 1e4) {
      cat("  NOTA: Moderado - revisar correlaciones\n")
    } else {
      cat("  OK: Aceptable\n")
    }
  }, error = function(e) {
    cat("  No se pudo calcular (vcov no disponible)\n")
  })

  # 10.7: Resumen intercepts de todos los modelos
  cat("\n  === RESUMEN TODOS LOS MODELOS ===\n\n")
  cat(sprintf("  %-35s %7s %7s %7s %7s %7s %5s\n",
              "Modelo", "psi1", "gamma", "epsil", "p(int)", "AIC", "Conv"))
  cat(strrep("-", 105), "\n")

  for (mname in names(models)) {
    mod <- models[[mname]]
    if (is.null(mod)) {
      cat(sprintf("  %-35s  --- ERROR ---\n", mname))
      next
    }
    coefs <- coef(mod)
    psi_int <- plogis(coefs["psi(Int)"])
    col_int <- plogis(coefs["col(Int)"])
    ext_int <- plogis(coefs["ext(Int)"])
    det_int <- plogis(coefs["p(Int)"])
    cat(sprintf("  %-35s %6.1f%% %6.2f%% %6.1f%% %6.1f%% %7.1f %5d\n",
                mname,
                psi_int * 100, col_int * 100, ext_int * 100, det_int * 100,
                mod@AIC, mod@opt$convergence))
  }
}

###############################################################################
# PASO 11: Comparacion ocupacion predicha vs naive por anio
###############################################################################

cat("\n\n--- PASO 11: Ocupacion predicha vs naive por anio ---\n\n")

if (!is.null(best_mod)) {
  tryCatch({
    # Prediccion de ocupacion suavizada (projected)
    proj <- projected(best_mod)

    # projected() devuelve una lista de matrices [n_sites x 2] por anio
    # Columna 1 = P(no ocupado), Columna 2 = P(ocupado)
    cat(sprintf("  %-6s %10s %10s %10s\n", "Anio", "Naive(%)", "Pred(%)", "Diff"))
    cat(strrep("-", 45), "\n")

    for (t in 1:T_years) {
      cols_t <- ((t-1) * J_reps + 1):(t * J_reps)
      naive <- mean(apply(y_filt[, cols_t], 1,
                          function(x) any(x == 1, na.rm = TRUE))) * 100

      # Probabilidad media de ocupacion predicha
      if (is.list(proj)) {
        pred_occ <- mean(proj[[t]][, 2], na.rm = TRUE) * 100
      } else {
        pred_occ <- NA
      }

      cat(sprintf("  %d   %8.1f%%  %8.1f%%  %+7.1f%%\n",
                  YEARS[t], naive,
                  ifelse(is.na(pred_occ), NA, pred_occ),
                  ifelse(is.na(pred_occ), NA, pred_occ - naive)))
    }
  }, error = function(e) {
    cat("  No se pudo calcular projected(): ", e$message, "\n")
  })
}

###############################################################################
# PASO 12: Exportar resultados
###############################################################################

cat("\n--- PASO 12: Exportar resultados ---\n\n")

results_dir <- file.path(
  "/Users/gfandos/Documents/GitHub/dynamic_occupancy_models_steppe_birds",
  "results"
)
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

# 12.1: Tabla de seleccion de modelos
model_sel_path <- file.path(results_dir, "otitar_v4_model_selection.csv")
write.csv(aic_table, model_sel_path, row.names = FALSE)
cat("  Tabla AIC guardada: ", model_sel_path, "\n")

# 12.2: Coeficientes del mejor modelo
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

  coef_path <- file.path(results_dir, "otitar_v4_best_model_coefficients.csv")
  write.csv(coef_df, coef_path, row.names = FALSE)
  cat("  Coeficientes guardados: ", coef_path, "\n")
}

# 12.3: Parametros de escalado de las variables dinamicas
scaling_df <- data.frame(
  variable = names(dyn_scaling),
  center = sapply(dyn_scaling, function(x) x$center),
  scale = sapply(dyn_scaling, function(x) x$scale),
  row.names = NULL
)
scaling_path <- file.path(results_dir, "otitar_v4_dynamic_scaling_params.csv")
write.csv(scaling_df, scaling_path, row.names = FALSE)
cat("  Parametros escalado: ", scaling_path, "\n")

###############################################################################
# PASO 13: Goodness-of-fit (parboot) del mejor modelo
###############################################################################

cat("\n--- PASO 13: Goodness-of-fit (parboot, nsim=50) ---\n\n")
cat("  NOTA: nsim=50 es bajo, solo para verificar. Para publicacion usar nsim>=500\n\n")

if (!is.null(best_mod)) {
  tryCatch({
    # Funcion de estadistico chi-cuadrado para parboot
    chisq_fn <- function(fm) {
      obs <- getY(fm@data)
      ex  <- fitted(fm)
      # Reemplazar 0s en expected para evitar division por 0
      ex[ex == 0] <- 1e-6
      ts <- sum((obs - ex)^2 / ex, na.rm = TRUE)
      return(c(chisq = ts))
    }

    pb <- parboot(best_mod, statistic = chisq_fn, nsim = 50, report = 10)

    cat("  Resultado parboot:\n")
    cat(sprintf("    Chi-cuadrado observado: %.1f\n", pb@t0))
    cat(sprintf("    Chi-cuadrado simulado (media): %.1f\n", mean(pb@t.star[,1])))
    cat(sprintf("    p-value: %.3f\n", mean(pb@t.star[,1] >= pb@t0)))

    p_gof <- mean(pb@t.star[,1] >= pb@t0)
    if (p_gof > 0.05) {
      cat("    => Modelo OK (no se rechaza H0 de buen ajuste)\n")
    } else {
      cat("    => AVISO: Posible falta de ajuste (p < 0.05)\n")
    }
  }, error = function(e) {
    cat("  Error en parboot: ", e$message, "\n")
    cat("  (Esto puede ocurrir con modelos complejos / datos con muchos NAs)\n")
  })
}

###############################################################################
# RESUMEN FINAL
###############################################################################

cat("\n\n")
cat(strrep("=", 70), "\n")
cat("  PIPELINE V4: COVARIABLES DINAMICAS — RESUMEN\n")
cat(strrep("=", 70), "\n\n")

if (!is.null(best_mod)) {
  # Extraer formulas del mejor modelo
  best_formula <- model_formulas[[best_name]]

  cat(sprintf("  Mejor modelo: %s\n", best_name))
  cat(sprintf("  Formula: %s\n", best_formula))
  cat(sprintf("  AIC: %.1f\n", best_aic_val))

  if (!is.na(baseline_aic)) {
    cat(sprintf("  ΔAIC vs baseline estatico (v3): %.1f\n", best_aic_val - baseline_aic))
    if (best_aic_val - baseline_aic < -2) {
      cat("  => Las covariables dinamicas MEJORAN el modelo\n")
    } else if (best_aic_val - baseline_aic > 2) {
      cat("  => Las covariables dinamicas NO mejoran el modelo (baseline estatico mejor)\n")
    } else {
      cat("  => Diferencia no sustancial con el baseline estatico\n")
    }
  }

  cat(sprintf("\n  Interceptos:\n"))
  coefs_best <- coef(best_mod)
  cat(sprintf("    psi1 (ocupacion inicial):  %.1f%%\n", plogis(coefs_best["psi(Int)"]) * 100))
  cat(sprintf("    gamma (colonizacion):      %.2f%%\n", plogis(coefs_best["col(Int)"]) * 100))
  cat(sprintf("    epsilon (extincion):       %.1f%%\n", plogis(coefs_best["ext(Int)"]) * 100))
  cat(sprintf("    p (deteccion):             %.1f%%\n", plogis(coefs_best["p(Int)"]) * 100))

  # Efectos clave de las variables dinamicas
  cat("\n  Efectos clave:\n")
  if (!is.null(s_best)) {
    for (sm_name in c("col", "ext")) {
      sm <- s_best[[sm_name]]
      if (is.null(sm)) next
      for (i in 1:nrow(sm)) {
        param <- rownames(sm)[i]
        is_dynamic <- any(sapply(dyn_vars, function(dv) grepl(dv, param)))
        if (is_dynamic) {
          est <- sm[i, "Estimate"]
          p_val <- sm[i, "P(>|z|)"]
          signif_label <- ifelse(p_val < 0.05, "SIGNIFICATIVO", "no significativo")
          submodel_label <- ifelse(sm_name == "col", "gamma", "epsilon")
          cat(sprintf("    %s en %s: coef=%.3f (p=%.4f) — %s\n",
                      param, submodel_label, est, p_val, signif_label))
        }
      }
    }
  }

  cat(sprintf("\n  Conclusion: "))
  if (!is.na(baseline_aic) && best_aic_val - baseline_aic < -2) {
    cat("Las covariables dinamicas mejoran sustancialmente el modelo.\n")
    cat("  La variabilidad interanual del clima/vegetacion es importante\n")
    cat("  para explicar la dinamica de ocupacion de Otis tarda.\n")
  } else if (!is.na(baseline_aic) && best_name == "m0: Baseline estatico") {
    cat("El baseline estatico sigue siendo el mejor modelo.\n")
    cat("  Las covariables dinamicas no aportan poder explicativo adicional.\n")
    cat("  La dinamica de Otis tarda parece determinada por factores\n")
    cat("  estructurales del habitat mas que por la variabilidad interanual.\n")
  } else {
    cat("Resultados mixtos. Revisar interpretacion ecologica.\n")
  }
}

cat("\n  Archivos generados:\n")
cat("    - results/otitar_v4_model_selection.csv\n")
cat("    - results/otitar_v4_best_model_coefficients.csv\n")
cat("    - results/otitar_v4_dynamic_scaling_params.csv\n")

cat("\n")
cat(strrep("=", 70), "\n")
cat("  FIN PIPELINE V4\n")
cat(strrep("=", 70), "\n\n")
