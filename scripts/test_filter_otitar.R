###############################################################################
# test_filter_otitar.R
#
# Prueba rapida: comparar filter_repeat_visits() con diferentes configuraciones
# para Otis tarda. Evaluar si el cambio de site_vars mejora la estructura
# de datos lo suficiente para obtener parametros plausibles en un modelo null.
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
J_reps  <- 10

cat("\n")
cat(strrep("=", 70), "\n")
cat("  TEST: filter_repeat_visits() para Otis tarda\n")
cat(strrep("=", 70), "\n\n")

###############################################################################
# 1. CARGAR DATOS CRUDOS
###############################################################################

raw_path <- "/Users/gfandos/Documents/GitHub/dynamic_occupancy_models_steppe_birds/data-raw/data/otitar/ebd_otitar_breeding_spain_zf.csv"
occ_raw <- read.csv(raw_path)
occ_raw <- occ_raw %>% filter(year >= 2017, year <= 2022)

cat("Datos crudos: ", nrow(occ_raw), " registros\n")
cat("  Localidades unicas: ", n_distinct(occ_raw$locality_id), "\n")
cat("  Observadores unicos: ", n_distinct(occ_raw$observer_id), "\n")
cat("  Combos localidad+observador: ",
    n_distinct(paste(occ_raw$locality_id, occ_raw$observer_id)), "\n")
cat("  Detecciones (species_observed==TRUE): ",
    sum(occ_raw$species_observed == TRUE), "\n\n")

###############################################################################
# 2. COMPARAR CONFIGURACIONES DE filter_repeat_visits()
###############################################################################

configs <- list(
  # Config A: ACTUAL (locality_id + observer_id, min_obs=3)
  A_actual = list(
    label = "ACTUAL: locality+observer, min_obs=3",
    site_vars = c("locality_id", "observer_id"),
    min_obs = 3, max_obs = 10
  ),
  # Config B: Solo locality_id, min_obs=3
  B_locality = list(
    label = "locality_id solamente, min_obs=3",
    site_vars = c("locality_id"),
    min_obs = 3, max_obs = 10
  ),
  # Config C: Solo locality_id, min_obs=2 (mas permisivo)
  C_loc_min2 = list(
    label = "locality_id solamente, min_obs=2",
    site_vars = c("locality_id"),
    min_obs = 2, max_obs = 10
  ),
  # Config D: Solo locality_id, min_obs=2, max_obs=30 (aun mas permisivo)
  D_loc_loose = list(
    label = "locality_id solamente, min_obs=2, max_obs=30",
    site_vars = c("locality_id"),
    min_obs = 2, max_obs = 30
  )
)

results <- list()

for (cfg_name in names(configs)) {
  cfg <- configs[[cfg_name]]

  cat(strrep("-", 60), "\n")
  cat("CONFIG ", cfg_name, ": ", cfg$label, "\n")
  cat(strrep("-", 60), "\n")

  occ <- tryCatch({
    filter_repeat_visits(
      occ_raw,
      min_obs = cfg$min_obs,
      max_obs = cfg$max_obs,
      annual_closure = TRUE,
      date_var = "observation_date",
      site_vars = cfg$site_vars
    )
  }, error = function(e) {
    cat("  ERROR: ", e$message, "\n")
    return(NULL)
  })

  if (is.null(occ)) next

  n_sites <- n_distinct(occ$site)
  cat("  Sitios totales: ", n_sites, "\n")
  cat("  Registros totales: ", nrow(occ), "\n")

  # Detecciones
  n_det <- sum(occ$species_observed == TRUE)
  cat("  Detecciones: ", n_det, "\n")
  cat("  Tasa deteccion bruta: ",
      round(100 * n_det / nrow(occ), 2), "%\n")

  # Por anio
  cat("\n  Sitios con datos por anio:\n")
  year_summary <- occ %>%
    group_by(year) %>%
    summarise(
      n_lists = n(),
      n_sites = n_distinct(site),
      n_det = sum(species_observed == TRUE),
      .groups = "drop"
    )
  for (i in 1:nrow(year_summary)) {
    cat(sprintf("    %d: %4d sitios, %5d listas, %3d detecciones\n",
                year_summary$year[i],
                year_summary$n_sites[i],
                year_summary$n_lists[i],
                year_summary$n_det[i]))
  }

  # Cobertura temporal por sitio
  years_per_site <- occ %>%
    group_by(site) %>%
    summarise(
      n_years = n_distinct(year),
      has_early = any(year <= 2019),
      has_late = any(year >= 2020),
      has_det = any(species_observed == TRUE),
      .groups = "drop"
    )

  cat("\n  Cobertura temporal:\n")
  for (ny in 1:6) {
    n <- sum(years_per_site$n_years == ny)
    cat(sprintf("    %d anos: %4d sitios (%5.1f%%)\n",
                ny, n, 100 * n / n_sites))
  }

  n_3plus <- sum(years_per_site$n_years >= 3)
  n_early_late <- sum(years_per_site$has_early & years_per_site$has_late)
  n_det_sites <- sum(years_per_site$has_det)

  cat(sprintf("\n  Sitios con >= 3 anos: %d (%.1f%%)\n",
              n_3plus, 100 * n_3plus / n_sites))
  cat(sprintf("  Sitios con datos temprano Y tardio: %d (%.1f%%)\n",
              n_early_late, 100 * n_early_late / n_sites))
  cat(sprintf("  Sitios con al menos 1 deteccion: %d (%.1f%%)\n",
              n_det_sites, 100 * n_det_sites / n_sites))

  # Filtro: >= 3 anos con datos tempranos y tardios
  sites_good <- years_per_site %>%
    filter(n_years >= 3, has_early, has_late)
  cat(sprintf("  Sitios que pasan filtro temporal (>=3y, early+late): %d (%.1f%%)\n",
              nrow(sites_good), 100 * nrow(sites_good) / n_sites))

  n_det_good <- sum(sites_good$has_det)
  cat(sprintf("  De ellos con detecciones: %d (%.1f%%)\n",
              n_det_good, 100 * n_det_good / nrow(sites_good)))

  results[[cfg_name]] <- list(
    n_sites = n_sites,
    n_det = n_det,
    n_3plus = n_3plus,
    n_early_late = n_early_late,
    n_good = nrow(sites_good),
    n_det_good = n_det_good,
    occ = occ,
    sites_good = sites_good
  )

  cat("\n")
}

###############################################################################
# 3. TABLA COMPARATIVA
###############################################################################

cat("\n")
cat(strrep("=", 70), "\n")
cat("  TABLA COMPARATIVA\n")
cat(strrep("=", 70), "\n\n")

cat(sprintf("%-42s %7s %6s %8s %7s %6s\n",
            "Configuracion", "Sitios", "Dets", ">=3anos", "E+L", "E+L+det"))
cat(strrep("-", 70), "\n")
for (cfg_name in names(results)) {
  r <- results[[cfg_name]]
  cat(sprintf("%-42s %7d %6d %8d %7d %6d\n",
              configs[[cfg_name]]$label,
              r$n_sites, r$n_det, r$n_3plus,
              r$n_early_late, r$n_det_good))
}

###############################################################################
# 4. ELEGIR MEJOR CONFIG Y AJUSTAR MODELO NULL
###############################################################################

cat("\n\n")
cat(strrep("=", 70), "\n")
cat("  AJUSTE DE MODELO NULL CON DATOS FILTRADOS\n")
cat(strrep("=", 70), "\n\n")

# Usar la config mas prometedora: la que maximice sitios con buena cobertura
# Probaremos con Config C (locality, min_obs=2) que deberia dar mas sitios
best_cfg <- "C_loc_min2"
cat("Usando configuracion: ", configs[[best_cfg]]$label, "\n\n")

occ_best <- results[[best_cfg]]$occ
good_sites <- results[[best_cfg]]$sites_good

# --- 4a. SIN filtro temporal (como el modelo actual pero con site_vars corregido)
cat("--- 4a. Modelo NULL sin filtro temporal (site_vars corregido) ---\n\n")

occ_nofilt <- occ_best

# Construir detection history en formato wide
build_umf <- function(occ_data, label) {

  cat("  Construyendo UMF para: ", label, "\n")

  # Formatear para unmarked usando format_unmarked_occu por anio
  occ_wide_list <- list()
  for (yr in YEARS) {
    occ_yr <- occ_data %>% filter(year == yr)
    if (nrow(occ_yr) == 0) next
    occ_wide_yr <- tryCatch({
      format_unmarked_occu(
        occ_yr,
        site_id = "site",
        response = "species_observed",
        site_covs = c("latitude", "longitude"),
        obs_covs = c("time_observations_started", "duration_minutes",
                      "effort_distance_km", "number_observers")
      )
    }, error = function(e) {
      cat("    Error en format_unmarked_occu para ", yr, ": ", e$message, "\n")
      return(NULL)
    })
    if (!is.null(occ_wide_yr)) {
      occ_wide_yr <- occ_wide_yr %>%
        rename_with(~ paste0(., ".", yr), -c(1))  # rename all except site id
      occ_wide_list[[as.character(yr)]] <- occ_wide_yr
    }
  }

  # Merge
  if (length(occ_wide_list) == 0) {
    cat("    No hay datos suficientes\n")
    return(NULL)
  }

  # Get the site ID column name
  id_col <- names(occ_wide_list[[1]])[1]

  occ_wide <- occ_wide_list[[1]]
  for (k in 2:length(occ_wide_list)) {
    occ_wide <- occ_wide %>%
      full_join(occ_wide_list[[k]], by = id_col, multiple = "all")
  }
  occ_wide <- occ_wide %>% distinct(!!sym(id_col), .keep_all = TRUE)

  # Convert logical to numeric
  cols_logical <- sapply(occ_wide, is.logical)
  occ_wide[, cols_logical] <- lapply(occ_wide[, cols_logical], as.numeric)

  n_sites <- nrow(occ_wide)
  cat("  Sitios en wide: ", n_sites, "\n")

  # Build detection matrix: find species_observed columns
  det_cols <- grep("^species_observed\\.", names(occ_wide), value = TRUE)
  y <- as.matrix(occ_wide[, det_cols])

  # Ensure exactly J_reps * T_years columns
  # We may have fewer than 10 reps per year depending on max_obs
  reps_per_year <- table(sub(".*\\.(\\d{4})$", "\\1", det_cols))
  actual_J <- as.integer(max(reps_per_year))

  # Check dimensions
  expected_cols <- actual_J * T_years
  if (ncol(y) != expected_cols) {
    cat("  Warning: y tiene ", ncol(y), " cols, esperadas ", expected_cols, "\n")
    # Pad with NAs if needed
    if (ncol(y) < expected_cols) {
      pad <- matrix(NA, nrow = n_sites, ncol = expected_cols - ncol(y))
      y <- cbind(y, pad)
    } else {
      y <- y[, 1:expected_cols]
    }
  }

  cat("  Matriz y: ", nrow(y), " x ", ncol(y), "\n")
  cat("  % NAs: ", round(100 * mean(is.na(y)), 1), "%\n")
  cat("  Detecciones: ", sum(y == 1, na.rm = TRUE), "\n")
  cat("  Sitios detectados: ", sum(apply(y, 1, function(x) any(x == 1, na.rm = TRUE))), "\n")

  # Naive occupancy per year
  cat("  Naive por anio: ")
  for (t in 1:T_years) {
    cols_t <- ((t-1) * actual_J + 1):(t * actual_J)
    if (max(cols_t) <= ncol(y)) {
      naive <- mean(apply(y[, cols_t, drop=FALSE], 1,
                          function(x) any(x == 1, na.rm = TRUE)))
      cat(sprintf("%d:%.1f%% ", YEARS[t], naive * 100))
    }
  }
  cat("\n")

  # Temporal coverage
  has_data_yr <- matrix(FALSE, n_sites, T_years)
  for (t in 1:T_years) {
    cols_t <- ((t-1) * actual_J + 1):min(t * actual_J, ncol(y))
    has_data_yr[, t] <- apply(y[, cols_t, drop=FALSE], 1,
                               function(x) any(!is.na(x)))
  }
  n_years_data <- rowSums(has_data_yr)
  cat("  Sitios con >=3 anos datos: ",
      sum(n_years_data >= 3), " (", round(100 * mean(n_years_data >= 3), 1), "%)\n")
  cat("  Sitios con 1 ano datos: ",
      sum(n_years_data == 1), " (", round(100 * mean(n_years_data == 1), 1), "%)\n")

  # Build UMF (null model: no covariates needed)
  umf <- unmarkedMultFrame(
    y = y,
    numPrimary = T_years
  )

  return(list(umf = umf, y = y, actual_J = actual_J,
              n_sites = n_sites, has_data_yr = has_data_yr))
}

# Build for unfiltered
umf_nofilt <- build_umf(occ_nofilt, "sin filtro temporal")

if (!is.null(umf_nofilt)) {
  cat("\n  Ajustando modelo null...\n")
  mod_null_nofilt <- tryCatch({
    colext(
      psiformula = ~ 1,
      gammaformula = ~ 1,
      epsilonformula = ~ 1,
      pformula = ~ 1,
      data = umf_nofilt$umf
    )
  }, error = function(e) {
    cat("  Error: ", e$message, "\n")
    NULL
  })

  if (!is.null(mod_null_nofilt)) {
    cat("\n  RESULTADOS MODELO NULL (sin filtro temporal):\n")
    coefs <- coef(mod_null_nofilt)
    cat(sprintf("    psi1: logit=%.3f -> prob=%.4f (%.1f%%)\n",
                coefs["psi(Int)"], plogis(coefs["psi(Int)"]),
                plogis(coefs["psi(Int)"]) * 100))
    cat(sprintf("    gamma: logit=%.3f -> prob=%.4f (%.1f%%)\n",
                coefs["col(Int)"], plogis(coefs["col(Int)"]),
                plogis(coefs["col(Int)"]) * 100))
    cat(sprintf("    epsilon: logit=%.3f -> prob=%.4f (%.1f%%)\n",
                coefs["ext(Int)"], plogis(coefs["ext(Int)"]),
                plogis(coefs["ext(Int)"]) * 100))
    cat(sprintf("    p: logit=%.3f -> prob=%.4f (%.1f%%)\n",
                coefs["p(Int)"], plogis(coefs["p(Int)"]),
                plogis(coefs["p(Int)"]) * 100))
    cat(sprintf("    AIC: %.1f\n", mod_null_nofilt@AIC))
    cat(sprintf("    Convergence: %d\n", mod_null_nofilt@opt$convergence))
  }
}

# --- 4b. CON filtro temporal (>= 3 anos, early + late)
cat("\n\n--- 4b. Modelo NULL con filtro temporal (>=3 anos, early+late) ---\n\n")

occ_filt <- occ_best %>%
  filter(site %in% good_sites$site)
cat("  Registros tras filtro: ", nrow(occ_filt), "\n")
cat("  Sitios tras filtro: ", n_distinct(occ_filt$site), "\n\n")

umf_filt <- build_umf(occ_filt, "con filtro temporal")

if (!is.null(umf_filt)) {
  cat("\n  Ajustando modelo null...\n")
  mod_null_filt <- tryCatch({
    colext(
      psiformula = ~ 1,
      gammaformula = ~ 1,
      epsilonformula = ~ 1,
      pformula = ~ 1,
      data = umf_filt$umf
    )
  }, error = function(e) {
    cat("  Error: ", e$message, "\n")
    NULL
  })

  if (!is.null(mod_null_filt)) {
    cat("\n  RESULTADOS MODELO NULL (con filtro temporal):\n")
    coefs <- coef(mod_null_filt)
    cat(sprintf("    psi1: logit=%.3f -> prob=%.4f (%.1f%%)\n",
                coefs["psi(Int)"], plogis(coefs["psi(Int)"]),
                plogis(coefs["psi(Int)"]) * 100))
    cat(sprintf("    gamma: logit=%.3f -> prob=%.4f (%.1f%%)\n",
                coefs["col(Int)"], plogis(coefs["col(Int)"]),
                plogis(coefs["col(Int)"]) * 100))
    cat(sprintf("    epsilon: logit=%.3f -> prob=%.4f (%.1f%%)\n",
                coefs["ext(Int)"], plogis(coefs["ext(Int)"]),
                plogis(coefs["ext(Int)"]) * 100))
    cat(sprintf("    p: logit=%.3f -> prob=%.4f (%.1f%%)\n",
                coefs["p(Int)"], plogis(coefs["p(Int)"]),
                plogis(coefs["p(Int)"]) * 100))
    cat(sprintf("    AIC: %.1f\n", mod_null_filt@AIC))
    cat(sprintf("    Convergence: %d\n", mod_null_filt@opt$convergence))

    # SEs
    s <- summary(mod_null_filt)
    cat("\n  Coeficientes con SE:\n")
    cat("    psi:\n"); print(s$psi)
    cat("    col:\n"); print(s$col)
    cat("    ext:\n"); print(s$ext)
    cat("    det:\n"); print(s$det)
  }
}

# --- 4c. TAMBIEN probar Config A actual con filtro temporal, para comparar
cat("\n\n--- 4c. Modelo NULL: Config ACTUAL (loc+obs) + filtro temporal ---\n\n")

occ_actual <- results[["A_actual"]]$occ
good_sites_actual <- results[["A_actual"]]$sites_good

occ_actual_filt <- occ_actual %>%
  filter(site %in% good_sites_actual$site)
cat("  Registros tras filtro: ", nrow(occ_actual_filt), "\n")
cat("  Sitios tras filtro: ", n_distinct(occ_actual_filt$site), "\n\n")

umf_actual_filt <- build_umf(occ_actual_filt, "actual + filtro temporal")

if (!is.null(umf_actual_filt)) {
  cat("\n  Ajustando modelo null...\n")
  mod_null_actual_filt <- tryCatch({
    colext(
      psiformula = ~ 1,
      gammaformula = ~ 1,
      epsilonformula = ~ 1,
      pformula = ~ 1,
      data = umf_actual_filt$umf
    )
  }, error = function(e) {
    cat("  Error: ", e$message, "\n")
    NULL
  })

  if (!is.null(mod_null_actual_filt)) {
    cat("\n  RESULTADOS MODELO NULL (actual + filtro temporal):\n")
    coefs <- coef(mod_null_actual_filt)
    cat(sprintf("    psi1: logit=%.3f -> prob=%.4f (%.1f%%)\n",
                coefs["psi(Int)"], plogis(coefs["psi(Int)"]),
                plogis(coefs["psi(Int)"]) * 100))
    cat(sprintf("    gamma: logit=%.3f -> prob=%.4f (%.1f%%)\n",
                coefs["col(Int)"], plogis(coefs["col(Int)"]),
                plogis(coefs["col(Int)"]) * 100))
    cat(sprintf("    epsilon: logit=%.3f -> prob=%.4f (%.1f%%)\n",
                coefs["ext(Int)"], plogis(coefs["ext(Int)"]),
                plogis(coefs["ext(Int)"]) * 100))
    cat(sprintf("    p: logit=%.3f -> prob=%.4f (%.1f%%)\n",
                coefs["p(Int)"], plogis(coefs["p(Int)"]),
                plogis(coefs["p(Int)"]) * 100))
    cat(sprintf("    AIC: %.1f\n", mod_null_actual_filt@AIC))
    cat(sprintf("    Convergence: %d\n", mod_null_actual_filt@opt$convergence))
  }
}

###############################################################################
# 5. TAMBIEN: refiltrar directamente desde el modelo original
###############################################################################

cat("\n\n")
cat(strrep("=", 70), "\n")
cat("  BONUS: Refiltrar datos del modelo ORIGINAL por cobertura temporal\n")
cat(strrep("=", 70), "\n\n")

mod_orig <- readRDS("/Users/gfandos/Documents/GitHub/dynamic_occupancy_models_steppe_birds/data/processed/otitar_model_object.rds")
y_orig <- mod_orig@data@y
n_orig <- nrow(y_orig)

# Identificar sitios con >= 3 anos de datos
has_data_orig <- matrix(FALSE, n_orig, T_years)
for (t in 1:T_years) {
  cols_t <- ((t-1) * 10 + 1):(t * 10)
  has_data_orig[, t] <- apply(y_orig[, cols_t], 1,
                               function(x) any(!is.na(x)))
}
n_years_orig <- rowSums(has_data_orig)

# Filtro: >= 3 anos, con al menos 1 temprano y 1 tardio
has_early_orig <- rowSums(has_data_orig[, 1:3]) > 0
has_late_orig <- rowSums(has_data_orig[, 4:6]) > 0
keep <- n_years_orig >= 3 & has_early_orig & has_late_orig

cat("Modelo original: ", n_orig, " sitios\n")
cat("Sitios que pasan filtro: ", sum(keep), " (", round(100 * mean(keep), 1), "%)\n")

# Filtrar y reconstruir UMF
y_filt2 <- y_orig[keep, ]
cat("Sitios filtrados: ", nrow(y_filt2), "\n")
cat("% NAs: ", round(100 * mean(is.na(y_filt2)), 1), "%\n")
cat("Detecciones: ", sum(y_filt2 == 1, na.rm = TRUE), "\n")
cat("Sitios detectados: ", sum(apply(y_filt2, 1, function(x)
  any(x == 1, na.rm = TRUE))), "\n")

# Naive per year
cat("Naive por anio: ")
for (t in 1:T_years) {
  cols_t <- ((t-1) * 10 + 1):(t * 10)
  naive <- mean(apply(y_filt2[, cols_t], 1,
                      function(x) any(x == 1, na.rm = TRUE)))
  cat(sprintf("%d:%.1f%% ", YEARS[t], naive * 100))
}
cat("\n\n")

# Fit null model
umf_filt2 <- unmarkedMultFrame(y = y_filt2, numPrimary = T_years)

mod_null_filt2 <- tryCatch({
  colext(
    psiformula = ~ 1, gammaformula = ~ 1,
    epsilonformula = ~ 1, pformula = ~ 1,
    data = umf_filt2
  )
}, error = function(e) {
  cat("Error: ", e$message, "\n")
  NULL
})

if (!is.null(mod_null_filt2)) {
  cat("RESULTADOS MODELO NULL (datos originales refiltrados):\n")
  coefs <- coef(mod_null_filt2)
  cat(sprintf("  psi1: logit=%.3f -> prob=%.4f (%.1f%%)\n",
              coefs["psi(Int)"], plogis(coefs["psi(Int)"]),
              plogis(coefs["psi(Int)"]) * 100))
  cat(sprintf("  gamma: logit=%.3f -> prob=%.4f (%.1f%%)\n",
              coefs["col(Int)"], plogis(coefs["col(Int)"]),
              plogis(coefs["col(Int)"]) * 100))
  cat(sprintf("  epsilon: logit=%.3f -> prob=%.4f (%.1f%%)\n",
              coefs["ext(Int)"], plogis(coefs["ext(Int)"]),
              plogis(coefs["ext(Int)"]) * 100))
  cat(sprintf("  p: logit=%.3f -> prob=%.4f (%.1f%%)\n",
              coefs["p(Int)"], plogis(coefs["p(Int)"]),
              plogis(coefs["p(Int)"]) * 100))
  cat(sprintf("  AIC: %.1f\n", mod_null_filt2@AIC))
  cat(sprintf("  Convergence: %d\n", mod_null_filt2@opt$convergence))

  s2 <- summary(mod_null_filt2)
  cat("\n  Coeficientes con SE:\n")
  cat("  psi:\n"); print(s2$psi)
  cat("  col:\n"); print(s2$col)
  cat("  ext:\n"); print(s2$ext)
  cat("  det:\n"); print(s2$det)
}

###############################################################################
# 6. RESUMEN FINAL
###############################################################################

cat("\n\n")
cat(strrep("=", 70), "\n")
cat("  RESUMEN: COMPARACION MODELOS NULL\n")
cat(strrep("=", 70), "\n\n")

cat("Interpretacion de resultados plausibles:\n")
cat("  psi1: 5-30% (estas especies son escasas pero presentes)\n")
cat("  gamma: 1-15% (colonizacion baja, son sedentarias)\n")
cat("  epsilon: 5-25% (extincion local moderada)\n")
cat("  p: 10-50% (deteccion variable)\n\n")

cat("Si los parametros del modelo null con datos filtrados siguen siendo\n")
cat("implausibles, el problema es mas profundo que el filtrado y requiere\n")
cat("cambio de escala espacial o enfoque analitico diferente.\n")
