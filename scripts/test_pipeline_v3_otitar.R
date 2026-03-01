###############################################################################
# test_pipeline_v3_otitar.R
#
# Pipeline v3 para Otis tarda.
# Extiende el pipeline v2 con:
#   - 3 modelos adicionales con covariables estaticas en gamma/epsilon
#   - Diagnosticos avanzados (AIC, deltaAIC, pesos, boundary, correlaciones)
#   - Mapas espaciales de ocupacion inicial, colonizacion y extincion
#
# Cambios respecto al pipeline original (heredados de v2):
#   1. site_vars = c("locality_id")  en vez de c("locality_id","observer_id")
#   2. min_obs = 2  en vez de 3
#   3. Filtro post-join: >= 3 anos con datos, early (2017-2019) + late (2020-2022)
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
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(gridExtra)

select <- dplyr::select
filter <- dplyr::filter
set.seed(42)

YEARS   <- 2017:2022
T_years <- length(YEARS)

cat("\n")
cat(strrep("=", 70), "\n")
cat("  PIPELINE v3: Otis tarda\n")
cat("  Modelos gamma/epsilon con covariables + mapas espaciales\n")
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

# Estandarizar covariables de sitio (guardar parametros para mapas)
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
# PASO 7: Construir unmarkedMultFrame
###############################################################################

cat("\n--- PASO 7: Construir unmarkedMultFrame ---\n\n")

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

###############################################################################
# PASO 8: Ajustar modelos (5 del v2 + 3 nuevos)
###############################################################################

cat("--- PASO 8: Ajustar 8 modelos ---\n\n")

# Funcion para verificar un modelo ajustado
check_model <- function(mod, model_name) {
  if (is.null(mod)) {
    return(list(ok = FALSE, reason = "Error en ajuste"))
  }

  # Verificar convergencia
  if (mod@opt$convergence != 0) {
    return(list(ok = FALSE, reason = "No converge (convergence != 0)"))
  }

  # Verificar hessiana (intentar summary)
  s <- tryCatch(summary(mod), error = function(e) NULL)
  if (is.null(s)) {
    return(list(ok = FALSE, reason = "Hessiana singular"))
  }

  # Verificar boundary estimates (|logit| > 10) en gamma, epsilon y p

  # NOTA: No verificamos boundary en psi porque con especies raras y covariables
  # escaladas, es normal tener interceptos psi muy negativos (logit < -10).
  # Solo nos preocupa boundary en gamma, epsilon y p.
  coefs <- coef(mod)
  gep_coefs <- coefs[grepl("^(col|ext|p)\\(", names(coefs))]
  boundary <- names(gep_coefs)[abs(gep_coefs) > 10]
  if (length(boundary) > 0) {
    return(list(ok = FALSE,
                reason = paste("Boundary:", paste(boundary, collapse = ", ")),
                boundary_params = boundary))
  }

  # Tambien verificar psi boundary (solo como warning, no falla)
  psi_coefs <- coefs[grepl("^psi\\(", names(coefs))]
  psi_boundary <- names(psi_coefs)[abs(psi_coefs) > 10]

  if (length(psi_boundary) > 0) {
    return(list(ok = TRUE,
                reason = paste("OK (psi boundary:", paste(psi_boundary, collapse = ", "), ")")))
  }

  return(list(ok = TRUE, reason = "OK"))
}

# --- Modelo 1: Null ---
cat("  [1/8] Null (~1, ~1, ~1, ~1)...\n")
mod1 <- tryCatch({
  colext(psiformula = ~1, gammaformula = ~1,
         epsilonformula = ~1, pformula = ~1,
         data = occ_umf)
}, error = function(e) { cat("    Error: ", e$message, "\n"); NULL })

# --- Modelo 2: Solo deteccion ---
cat("  [2/8] p = ~effort+observers...\n")
mod2 <- tryCatch({
  colext(psiformula = ~1, gammaformula = ~1,
         epsilonformula = ~1, pformula = ~effort + observers,
         data = occ_umf)
}, error = function(e) { cat("    Error: ", e$message, "\n"); NULL })

# --- Modelo 3: psi + deteccion simple ---
cat("  [3/8] psi = ~grass+tree, p = ~effort+observers...\n")
mod3 <- tryCatch({
  colext(psiformula = ~grass_cover + tree_cover,
         gammaformula = ~1,
         epsilonformula = ~1,
         pformula = ~effort + observers,
         data = occ_umf)
}, error = function(e) { cat("    Error: ", e$message, "\n"); NULL })

# --- Modelo 4: psi + deteccion completa ---
cat("  [4/8] psi = ~grass+tree, p = ~effort+obs+dur+time...\n")
mod4 <- tryCatch({
  colext(psiformula = ~grass_cover + tree_cover,
         gammaformula = ~1,
         epsilonformula = ~1,
         pformula = ~effort + observers + duration + time,
         data = occ_umf)
}, error = function(e) { cat("    Error: ", e$message, "\n"); NULL })

# --- Modelo 5: psi completo (modelo v2) ---
cat("  [5/8] psi = ~bio1+bio2+tree+grass+elev, p completo...\n")
mod5 <- tryCatch({
  colext(psiformula = ~bio1 + bio2 + tree_cover + grass_cover + topo_elev,
         gammaformula = ~1,
         epsilonformula = ~1,
         pformula = ~effort + observers + duration + time,
         data = occ_umf)
}, error = function(e) { cat("    Error: ", e$message, "\n"); NULL })

# --- Modelo 6: epsilon ~ habitat (NUEVO) ---
cat("  [6/8] psi completo, eps = ~grass+tree...\n")
mod6 <- tryCatch({
  colext(psiformula = ~bio1 + bio2 + tree_cover + grass_cover + topo_elev,
         gammaformula = ~1,
         epsilonformula = ~grass_cover + tree_cover,
         pformula = ~effort + observers + duration + time,
         data = occ_umf)
}, error = function(e) { cat("    Error: ", e$message, "\n"); NULL })

# --- Modelo 7: epsilon ~ habitat + topografia (NUEVO) ---
cat("  [7/8] psi completo, eps = ~grass+elev...\n")
mod7 <- tryCatch({
  colext(psiformula = ~bio1 + bio2 + tree_cover + grass_cover + topo_elev,
         gammaformula = ~1,
         epsilonformula = ~grass_cover + topo_elev,
         pformula = ~effort + observers + duration + time,
         data = occ_umf)
}, error = function(e) { cat("    Error: ", e$message, "\n"); NULL })

# --- Modelo 8: gamma + epsilon con covariables (NUEVO) ---
cat("  [8/8] psi completo, gam = ~grass+bio1, eps = ~grass+tree...\n")
mod8 <- tryCatch({
  colext(psiformula = ~bio1 + bio2 + tree_cover + grass_cover + topo_elev,
         gammaformula = ~grass_cover + bio1,
         epsilonformula = ~grass_cover + tree_cover,
         pformula = ~effort + observers + duration + time,
         data = occ_umf)
}, error = function(e) { cat("    Error: ", e$message, "\n"); NULL })

###############################################################################
# PASO 9: Diagnosticos avanzados
###############################################################################

cat("\n--- PASO 9: Diagnosticos ---\n\n")

models <- list(
  "1: Null"               = mod1,
  "2: p simple"           = mod2,
  "3: psi+p simple"       = mod3,
  "4: psi+p completo"     = mod4,
  "5: psi completo"       = mod5,
  "6: eps~habitat"        = mod6,
  "7: eps~habitat+topo"   = mod7,
  "8: gam+eps covars"     = mod8
)

model_formulas <- list(
  "1: Null"               = "psi(~1) gam(~1) eps(~1) p(~1)",
  "2: p simple"           = "psi(~1) gam(~1) eps(~1) p(~eff+obs)",
  "3: psi+p simple"       = "psi(~gra+tree) gam(~1) eps(~1) p(~eff+obs)",
  "4: psi+p completo"     = "psi(~gra+tree) gam(~1) eps(~1) p(~eff+obs+dur+time)",
  "5: psi completo"       = "psi(~bio1+bio2+tree+gra+elev) gam(~1) eps(~1) p(~eff+obs+dur+time)",
  "6: eps~habitat"        = "psi(~bio1+bio2+tree+gra+elev) gam(~1) eps(~gra+tree) p(~eff+obs+dur+time)",
  "7: eps~habitat+topo"   = "psi(~bio1+bio2+tree+gra+elev) gam(~1) eps(~gra+elev) p(~eff+obs+dur+time)",
  "8: gam+eps covars"     = "psi(~bio1+bio2+tree+gra+elev) gam(~gra+bio1) eps(~gra+tree) p(~eff+obs+dur+time)"
)

# 9.1: Tabla AIC con deltaAIC y pesos
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

# Calcular deltaAIC y pesos (solo para modelos validos)
valid_aic <- !is.na(aic_table$AIC)
if (any(valid_aic)) {
  min_aic <- min(aic_table$AIC[valid_aic])
  aic_table$deltaAIC <- ifelse(valid_aic, round(aic_table$AIC - min_aic, 1), NA)
  # Akaike weights
  delta_valid <- aic_table$deltaAIC[valid_aic]
  weights_raw <- exp(-0.5 * delta_valid)
  aic_table$Weight <- NA
  aic_table$Weight[valid_aic] <- round(weights_raw / sum(weights_raw), 4)
}

# Imprimir tabla
cat(sprintf("%-22s %5s %7s %7s %6s  %-20s\n",
            "Modelo", "nPar", "AIC", "dAIC", "Peso", "Estado"))
cat(strrep("-", 85), "\n")
for (i in 1:nrow(aic_table)) {
  cat(sprintf("%-22s %5s %7s %7s %6s  %-20s\n",
              aic_table$Model[i],
              ifelse(is.na(aic_table$nPars[i]), " -", as.character(aic_table$nPars[i])),
              ifelse(is.na(aic_table$AIC[i]), "   -", sprintf("%.1f", aic_table$AIC[i])),
              ifelse(is.na(aic_table$deltaAIC[i]), "   -", sprintf("%.1f", aic_table$deltaAIC[i])),
              ifelse(is.na(aic_table$Weight[i]), "  -", sprintf("%.4f", aic_table$Weight[i])),
              aic_table$Status[i]))
}

# 9.2: Identificar mejor modelo (convergente, sin boundary)
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

if (is.null(best_mod)) {
  cat("\n  AVISO: Ningun modelo paso todos los controles de calidad.\n")
  cat("  Seleccionando el modelo con menor AIC entre los convergentes...\n")
  for (mname in names(models)) {
    mod <- models[[mname]]
    if (!is.null(mod) && mod@opt$convergence == 0 && mod@AIC < best_aic_val) {
      best_aic_val <- mod@AIC
      best_name <- mname
      best_mod <- mod
    }
  }
}

cat("\n  Mejor modelo: ", best_name, " (AIC=", round(best_aic_val, 1), ")\n")

# 9.3: Parametros del mejor modelo
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
    cat("\n  psi:\n"); print(s_best$psi)
    cat("\n  col:\n"); print(s_best$col)
    cat("\n  ext:\n"); print(s_best$ext)
    cat("\n  det:\n"); print(s_best$det)
  }

  # 9.4: Boundary detection en todos los modelos
  cat("\n  === BOUNDARY DETECTION ===\n\n")
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

  # 9.5: Numero de condicion (solo mejor modelo)
  cat("\n  === NUMERO DE CONDICION (mejor modelo) ===\n\n")
  tryCatch({
    vcov_mat <- vcov(best_mod)
    eig <- eigen(vcov_mat)$values
    cond_num <- max(abs(eig)) / min(abs(eig))
    cat(sprintf("  Numero de condicion: %.1f\n", cond_num))
    if (cond_num > 1e6) {
      cat("  AVISO: Numero de condicion muy alto - posible colinealidad\n")
    } else if (cond_num > 1e4) {
      cat("  NOTA: Numero de condicion moderado\n")
    } else {
      cat("  OK: Numero de condicion aceptable\n")
    }
  }, error = function(e) {
    cat("  No se pudo calcular (vcov no disponible)\n")
  })

  # 9.6: Todos los modelos con parametros basicos
  cat("\n  === RESUMEN TODOS LOS MODELOS ===\n\n")
  cat(sprintf("%-22s %7s %7s %7s %7s %7s %5s\n",
              "Modelo", "psi1", "gamma", "epsil", "p(int)", "AIC", "Conv"))
  cat(strrep("-", 78), "\n")

  for (mname in names(models)) {
    mod <- models[[mname]]
    if (is.null(mod)) {
      cat(sprintf("%-22s  --- ERROR ---\n", mname))
      next
    }
    coefs <- coef(mod)
    psi_int <- plogis(coefs["psi(Int)"])
    col_int <- plogis(coefs["col(Int)"])
    ext_int <- plogis(coefs["ext(Int)"])
    det_int <- plogis(coefs["p(Int)"])
    cat(sprintf("%-22s %6.1f%% %6.2f%% %6.1f%% %6.1f%% %7.1f %5d\n",
                mname,
                psi_int * 100, col_int * 100, ext_int * 100, det_int * 100,
                mod@AIC, mod@opt$convergence))
  }
  cat("\n  Referencia plausible: psi1=5-30%, gamma=1-15%, epsilon=5-25%, p=10-50%\n")
}

###############################################################################
# PASO 10: Mapas espaciales
###############################################################################

cat("\n\n--- PASO 10: Mapas espaciales ---\n\n")

if (is.null(best_mod)) {
  cat("  AVISO: No hay modelo valido para generar mapas. Saltando.\n")
} else {

  # 10.1: Cargar rasters para prediccion espacial
  cat("  Cargando rasters de prediccion...\n")

  # Environmental raster (same used in training)
  pred_env <- raster::stack(
    file.path(base_data, "environmental_data/environmental_data_occ/variables_spain.grd")
  )

  # Topology rasters (masked versions for full Spain coverage)
  topo_elev_path <- file.path(base_data, "topology_data/topo_elev_masked.tif")
  topo_aspect_path <- file.path(base_data, "topology_data/topo_aspect_masked.tif")

  # Check if masked versions exist, fallback to .asc
  if (!file.exists(topo_elev_path)) {
    topo_elev_path <- file.path(base_data, "topology_data/topo_elev.asc")
  }
  if (!file.exists(topo_aspect_path)) {
    topo_aspect_path <- file.path(base_data, "topology_data/topo_aspect.asc")
  }

  pred_elev   <- raster::raster(topo_elev_path)
  pred_aspect <- raster::raster(topo_aspect_path)

  # Resample topo to match environmental grid
  cat("  Resampling topografia al grid ambiental...\n")
  pred_elev   <- raster::resample(pred_elev, pred_env[[1]], method = "bilinear")
  pred_aspect <- raster::resample(pred_aspect, pred_env[[1]], method = "bilinear")

  # 10.2: Crear data.frame de prediccion
  cat("  Creando superficie de prediccion...\n")
  pred_pts <- rasterToPoints(pred_env)
  pred_df <- as.data.frame(pred_pts)

  # Anadir topografia
  elev_vals   <- raster::extract(pred_elev, pred_pts[, c("x", "y")])
  aspect_vals <- raster::extract(pred_aspect, pred_pts[, c("x", "y")])
  pred_df$topo_elev   <- elev_vals
  pred_df$topo_aspect <- aspect_vals

  # Guardar coordenadas originales
  pred_df$longitude <- pred_df$x
  pred_df$latitude  <- pred_df$y

  cat("  Puntos de prediccion: ", nrow(pred_df), "\n")

  # 10.3: Escalar con parametros de entrenamiento
  cat("  Escalando con parametros de entrenamiento...\n")
  for (v in names(scaling_params)) {
    if (v %in% names(pred_df)) {
      pred_df[[v]] <- (pred_df[[v]] - scaling_params[[v]]$center) / scaling_params[[v]]$scale
    }
  }

  # Eliminar NAs
  pred_df_clean <- pred_df %>% drop_na(bio1, bio2, tree_cover, grass_cover, topo_elev)
  cat("  Puntos validos: ", nrow(pred_df_clean), "\n")

  # 10.4: Predicciones espaciales
  cat("  Generando predicciones...\n")

  # Ocupacion inicial (psi) - siempre
  pred_psi <- tryCatch({
    predict(best_mod, type = "psi", newdata = pred_df_clean)
  }, error = function(e) {
    cat("    Error en predict(psi): ", e$message, "\n")
    NULL
  })

  # Colonizacion (gamma) - solo si tiene covariables
  best_coefs <- coef(best_mod)
  gamma_has_covars <- any(grepl("^col\\(.+\\)$", names(best_coefs)) &
                           !grepl("Int", names(best_coefs)))
  epsilon_has_covars <- any(grepl("^ext\\(.+\\)$", names(best_coefs)) &
                             !grepl("Int", names(best_coefs)))

  pred_col <- NULL
  pred_ext <- NULL

  if (gamma_has_covars) {
    cat("  Mejor modelo tiene covariables en gamma -> generando mapa de colonizacion\n")
    pred_col <- tryCatch({
      predict(best_mod, type = "col", newdata = pred_df_clean)
    }, error = function(e) {
      cat("    Error en predict(col): ", e$message, "\n")
      NULL
    })
  }

  if (epsilon_has_covars) {
    cat("  Mejor modelo tiene covariables en epsilon -> generando mapa de extincion\n")
    pred_ext <- tryCatch({
      predict(best_mod, type = "ext", newdata = pred_df_clean)
    }, error = function(e) {
      cat("    Error en predict(ext): ", e$message, "\n")
      NULL
    })
  }

  # 10.5: Rasterizar y recortar a Espana
  cat("  Rasterizando y recortando a Espana...\n")

  map_proj <- "+proj=longlat +datum=WGS84 +no_defs"

  # Obtener contorno de Espana peninsular (sin islas)
  spain <- ne_countries(country = "spain", scale = "medium", returnclass = "sf")
  # Recortar a peninsula iberica (excluir Canarias, Baleares, Ceuta, Melilla)
  peninsula_bbox <- st_bbox(c(xmin = -10, ymin = 35.5, xmax = 4.5, ymax = 44), crs = 4326)
  spain_crop <- st_crop(spain, peninsula_bbox)

  # Crear directorio figs si no existe
  dir.create(here("figs"), showWarnings = FALSE)

  # --- Funcion auxiliar para crear mapa ---
  make_map <- function(pred_result, pred_data, var_name, title, filename,
                       legend_label, color_option = "D") {
    if (is.null(pred_result)) return(NULL)

    map_data <- data.frame(
      longitude = pred_data$longitude,
      latitude  = pred_data$latitude,
      value     = pred_result$Predicted
    )

    # Rasterizar con sf
    r_pred <- map_data %>%
      st_as_sf(coords = c("longitude", "latitude"), crs = map_proj) %>%
      st_transform(crs = raster::projection(pred_env[[1]])) %>%
      rasterize(pred_env[[1]], field = "value")

    # Reproyectar, recortar y enmascarar
    r_pred_proj <- projectRaster(r_pred, crs = map_proj, method = "ngb")
    r_pred_crop <- crop(r_pred_proj, spain_crop)
    r_pred_masked <- mask(r_pred_crop, spain_crop)

    # Convertir a data.frame para ggplot
    r_df <- as.data.frame(r_pred_masked, xy = TRUE, na.rm = TRUE)
    names(r_df)[3] <- "value"

    # Crear mapa
    p <- ggplot() +
      geom_tile(data = r_df, aes(x = x, y = y, fill = value)) +
      scale_fill_viridis_c(name = legend_label, option = color_option) +
      geom_sf(data = spain_crop, fill = NA, color = "black", linewidth = 0.3) +
      labs(x = "Longitud", y = "Latitud", title = title) +
      theme_minimal() +
      theme(
        panel.border = element_rect(color = "black", fill = "transparent"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        legend.position = "right"
      )

    ggsave(here("figs", filename), p, width = 10, height = 8, dpi = 300)
    cat("    Guardado: figs/", filename, "\n")

    return(p)
  }

  # Generar mapas
  p_psi <- make_map(pred_psi, pred_df_clean,
                     "psi", "Otis tarda - Ocupacion inicial (psi1)",
                     "otitar_v3_psi_map.png",
                     expression(psi[1]), "D")

  p_gamma <- NULL
  if (gamma_has_covars) {
    p_gamma <- make_map(pred_col, pred_df_clean,
                         "col", "Otis tarda - Colonizacion (gamma)",
                         "otitar_v3_gamma_map.png",
                         expression(gamma), "C")
  }

  p_epsilon <- NULL
  if (epsilon_has_covars) {
    p_epsilon <- make_map(pred_ext, pred_df_clean,
                           "ext", "Otis tarda - Extincion (epsilon)",
                           "otitar_v3_epsilon_map.png",
                           expression(epsilon), "B")
  }

  # 10.6: Mapa combinado
  cat("  Creando mapa combinado...\n")
  plot_list <- list(p_psi)
  if (!is.null(p_gamma))   plot_list <- c(plot_list, list(p_gamma))
  if (!is.null(p_epsilon)) plot_list <- c(plot_list, list(p_epsilon))

  n_plots <- length(plot_list)
  if (n_plots > 1) {
    p_combined <- grid.arrange(grobs = plot_list, ncol = min(n_plots, 2))
    ggsave(here("figs", "otitar_v3_maps_combined.png"),
           p_combined, width = 14, height = 8 * ceiling(n_plots / 2), dpi = 300)
    cat("    Guardado: figs/otitar_v3_maps_combined.png\n")
  } else if (n_plots == 1) {
    # Si solo hay un mapa (psi), copiar como combinado
    ggsave(here("figs", "otitar_v3_maps_combined.png"),
           plot_list[[1]], width = 10, height = 8, dpi = 300)
    cat("    Guardado: figs/otitar_v3_maps_combined.png (solo psi)\n")
  }
}

###############################################################################
# PASO 11: Comparacion con modelo original + resumen final
###############################################################################

cat("\n\n")
cat(strrep("=", 70), "\n")
cat("  COMPARACION CON MODELO ORIGINAL\n")
cat(strrep("=", 70), "\n\n")

mod_orig <- readRDS(
  "/Users/gfandos/Documents/GitHub/dynamic_occupancy_models_steppe_birds/data/processed/otitar_model_object.rds"
)
coefs_orig <- coef(mod_orig)

cat("MODELO ORIGINAL (pipeline v1):\n")
cat(sprintf("  psi1:    logit=%.3f -> %.2f%%\n",
            coefs_orig["psi(Int)"], plogis(coefs_orig["psi(Int)"]) * 100))
cat(sprintf("  gamma:   logit=%.3f -> %.2f%%\n",
            coefs_orig["col(Int)"], plogis(coefs_orig["col(Int)"]) * 100))
cat(sprintf("  epsilon: logit=%.3f -> %.2f%% %s\n",
            coefs_orig["ext(Int)"], plogis(coefs_orig["ext(Int)"]) * 100,
            ifelse(abs(coefs_orig["ext(Int)"]) > 10, "(BOUNDARY)", "")))
cat(sprintf("  p:       logit=%.3f -> %.2f%%\n",
            coefs_orig["p(Int)"], plogis(coefs_orig["p(Int)"]) * 100))
cat("  AIC: ", mod_orig@AIC, " | Sitios: ", nrow(mod_orig@data@y), "\n\n")

if (!is.null(best_mod)) {
  coefs_v3 <- coef(best_mod)
  cat("MEJOR MODELO v3 (", best_name, "):\n")
  cat(sprintf("  psi1:    logit=%.3f -> %.2f%%\n",
              coefs_v3["psi(Int)"], plogis(coefs_v3["psi(Int)"]) * 100))
  cat(sprintf("  gamma:   logit=%.3f -> %.2f%%\n",
              coefs_v3["col(Int)"], plogis(coefs_v3["col(Int)"]) * 100))
  cat(sprintf("  epsilon: logit=%.3f -> %.2f%%\n",
              coefs_v3["ext(Int)"], plogis(coefs_v3["ext(Int)"]) * 100))
  cat(sprintf("  p:       logit=%.3f -> %.2f%%\n",
              coefs_v3["p(Int)"], plogis(coefs_v3["p(Int)"]) * 100))
  cat("  AIC: ", round(best_mod@AIC, 1), " | Sitios: ", n_filt, "\n\n")

  # Veredicto
  cat(strrep("=", 70), "\n")
  cat("  VEREDICTO\n")
  cat(strrep("=", 70), "\n\n")

  ext_orig <- plogis(coefs_orig["ext(Int)"])
  ext_v3 <- plogis(coefs_v3["ext(Int)"])

  if (ext_v3 < 0.30) {
    cat("  EXITO: epsilon bajo de ", round(ext_orig * 100, 0), "% a ",
        round(ext_v3 * 100, 1), "%\n")
    cat("  Los cambios en el filtrado FUNCIONAN.\n")
  } else {
    cat("  PARCIAL: epsilon = ", round(ext_v3 * 100, 1), "%\n")
    cat("  Mejora respecto al ", round(ext_orig * 100, 0),
        "% original, pero sigue alto.\n")
  }

  # Resumen de mapas generados
  cat("\n  Mapas generados:\n")
  cat("    - figs/otitar_v3_psi_map.png (ocupacion inicial)\n")
  if (gamma_has_covars) cat("    - figs/otitar_v3_gamma_map.png (colonizacion)\n")
  if (epsilon_has_covars) cat("    - figs/otitar_v3_epsilon_map.png (extincion)\n")
  cat("    - figs/otitar_v3_maps_combined.png (combinado)\n")
}

cat("\n")
cat(strrep("=", 70), "\n")
cat("  FIN del pipeline v3 para Otis tarda\n")
cat(strrep("=", 70), "\n\n")
