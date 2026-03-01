###############################################################################
# test_filter_otitar_v2.R
#
# Prueba rapida para Otis tarda. Dos enfoques:
#   PARTE A: Analizar la cobertura temporal a nivel de LOCALIDAD en los datos
#            crudos, evaluando el efecto de quitar observer_id de site_vars.
#   PARTE B: Tomar el modelo original ya ajustado, filtrar los sitios por
#            cobertura temporal (>= 3 anos), y re-ajustar un modelo null.
#            Esto nos dice inmediatamente si el filtro temporal basta.
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
cat("  TEST v2: Otis tarda - diagnostico de cobertura temporal\n")
cat(strrep("=", 70), "\n\n")

###############################################################################
# PARTE A: Analisis de cobertura a nivel de LOCALIDAD
###############################################################################

cat(strrep("=", 70), "\n")
cat("  PARTE A: Cobertura temporal por LOCALIDAD en datos crudos\n")
cat(strrep("=", 70), "\n\n")

raw_path <- "/Users/gfandos/Documents/GitHub/dynamic_occupancy_models_steppe_birds/data-raw/data/otitar/ebd_otitar_breeding_spain_zf.csv"
occ_raw <- read.csv(raw_path)
occ_raw <- occ_raw %>% filter(year >= 2017, year <= 2022)

cat("Datos crudos totales: ", nrow(occ_raw), " checklists\n")
cat("Localidades unicas: ", n_distinct(occ_raw$locality_id), "\n")
cat("Observadores unicos: ", n_distinct(occ_raw$observer_id), "\n\n")

# --- A1: Cobertura temporal por LOCALIDAD (sin considerar observer_id) ---
cat("--- A1: Cobertura por LOCALIDAD (ignorando observer_id) ---\n\n")

loc_summary <- occ_raw %>%
  group_by(locality_id) %>%
  summarise(
    n_checklists = n(),
    n_years = n_distinct(year),
    n_observers = n_distinct(observer_id),
    has_early = any(year <= 2019),
    has_late = any(year >= 2020),
    n_det = sum(species_observed == TRUE),
    years_list = paste(sort(unique(year)), collapse = ","),
    .groups = "drop"
  )

n_locs <- nrow(loc_summary)
cat("Total localidades: ", n_locs, "\n\n")

cat("Distribucion de checklists por localidad:\n")
cat("  Media: ", round(mean(loc_summary$n_checklists), 1), "\n")
cat("  Mediana: ", median(loc_summary$n_checklists), "\n")
cat("  Max: ", max(loc_summary$n_checklists), "\n\n")

cat("Cobertura temporal por localidad:\n")
for (ny in 1:6) {
  n <- sum(loc_summary$n_years == ny)
  n_det <- sum(loc_summary$n_years == ny & loc_summary$n_det > 0)
  cat(sprintf("  %d anos: %5d localidades (%5.1f%%), %4d con detecciones\n",
              ny, n, 100 * n / n_locs, n_det))
}

n_3plus <- sum(loc_summary$n_years >= 3)
n_el <- sum(loc_summary$has_early & loc_summary$has_late)
n_good <- sum(loc_summary$n_years >= 3 & loc_summary$has_early & loc_summary$has_late)
n_good_det <- sum(loc_summary$n_years >= 3 & loc_summary$has_early &
                    loc_summary$has_late & loc_summary$n_det > 0)

cat(sprintf("\nLocalidades con >= 3 anos: %d (%.1f%%)\n", n_3plus, 100 * n_3plus / n_locs))
cat(sprintf("Localidades con datos early+late: %d (%.1f%%)\n", n_el, 100 * n_el / n_locs))
cat(sprintf("Localidades 'buenas' (>=3y, E+L): %d (%.1f%%)\n", n_good, 100 * n_good / n_locs))
cat(sprintf("De ellas con al menos 1 deteccion: %d (%.1f%%)\n",
            n_good_det, 100 * n_good_det / max(n_good, 1)))

# --- A2: Comparar replicas: locality_id vs locality_id+observer_id ---
cat("\n--- A2: Efecto de quitar observer_id ---\n\n")

# Con observer_id (actual): replicas = checklists del MISMO observador en la misma localidad-ano
loc_obs_yr <- occ_raw %>%
  group_by(locality_id, observer_id, year) %>%
  summarise(n_reps = n(), .groups = "drop")

cat("CON observer_id (actual):\n")
cat("  Replicas por sitio-ano (media): ", round(mean(loc_obs_yr$n_reps), 1), "\n")
cat("  Replicas por sitio-ano (mediana): ", median(loc_obs_yr$n_reps), "\n")
cat("  Sitios-ano con >= 3 replicas: ",
    sum(loc_obs_yr$n_reps >= 3),
    " (", round(100 * mean(loc_obs_yr$n_reps >= 3), 1), "%)\n")

# Sin observer_id (propuesto): replicas = TODAS las checklists en la misma localidad-ano
loc_yr <- occ_raw %>%
  group_by(locality_id, year) %>%
  summarise(
    n_reps = n(),
    n_observers = n_distinct(observer_id),
    .groups = "drop"
  )

cat("\nSIN observer_id (propuesto):\n")
cat("  Replicas por localidad-ano (media): ", round(mean(loc_yr$n_reps), 1), "\n")
cat("  Replicas por localidad-ano (mediana): ", median(loc_yr$n_reps), "\n")
cat("  Localidades-ano con >= 3 replicas: ",
    sum(loc_yr$n_reps >= 3),
    " (", round(100 * mean(loc_yr$n_reps >= 3), 1), "%)\n")
cat("  Localidades-ano con >= 2 replicas: ",
    sum(loc_yr$n_reps >= 2),
    " (", round(100 * mean(loc_yr$n_reps >= 2), 1), "%)\n")
cat("  Media de observadores por localidad-ano: ",
    round(mean(loc_yr$n_observers), 1), "\n")

# A3: Localidades con >= 2 replicas por ano, y cuantas tienen buena cobertura temporal
cat("\n--- A3: Localidades con >= 2 replicas Y buena cobertura temporal ---\n\n")

# Para cada localidad, cuantos anos tienen >= 2 replicas?
loc_good_reps <- loc_yr %>%
  filter(n_reps >= 2) %>%
  group_by(locality_id) %>%
  summarise(
    n_years_good = n(),
    has_early = any(year <= 2019),
    has_late = any(year >= 2020),
    total_reps = sum(n_reps),
    .groups = "drop"
  )

cat("Localidades con >= 2 replicas en al menos 1 ano: ", nrow(loc_good_reps), "\n")
for (ny in 1:6) {
  n <- sum(loc_good_reps$n_years_good == ny)
  cat(sprintf("  %d anos con >=2 reps: %5d localidades\n", ny, n))
}
n_temporal_good <- sum(loc_good_reps$n_years_good >= 3 &
                         loc_good_reps$has_early & loc_good_reps$has_late)
cat(sprintf("\nLocalidades con >=2 reps en >=3 anos (E+L): %d\n", n_temporal_good))

# Chequear detecciones en estas localidades
good_loc_ids <- loc_good_reps$locality_id[
  loc_good_reps$n_years_good >= 3 & loc_good_reps$has_early & loc_good_reps$has_late
]
det_in_good <- occ_raw %>%
  filter(locality_id %in% good_loc_ids, species_observed == TRUE) %>%
  n_distinct("locality_id")

cat("De ellas con al menos 1 deteccion de O. tarda: ",
    sum(occ_raw$locality_id %in% good_loc_ids & occ_raw$species_observed == TRUE), " detecciones\n")
cat("En ", n_distinct(occ_raw$locality_id[occ_raw$locality_id %in% good_loc_ids &
                                           occ_raw$species_observed == TRUE]),
    " localidades distintas\n")


###############################################################################
# PARTE B: Refiltrar el modelo ORIGINAL por cobertura temporal
###############################################################################

cat("\n\n")
cat(strrep("=", 70), "\n")
cat("  PARTE B: Re-ajustar modelo NULL filtrando por cobertura temporal\n")
cat(strrep("=", 70), "\n\n")

mod_orig <- readRDS("/Users/gfandos/Documents/GitHub/dynamic_occupancy_models_steppe_birds/data/processed/otitar_model_object.rds")
y_orig <- mod_orig@data@y
n_orig <- nrow(y_orig)
J <- 10  # replicas por periodo primario

cat("Modelo original: ", n_orig, " sitios x ", ncol(y_orig), " observaciones\n\n")

# Calcular cobertura temporal por sitio
has_data <- matrix(FALSE, n_orig, T_years)
has_det <- matrix(FALSE, n_orig, T_years)
n_reps <- matrix(0L, n_orig, T_years)
for (t in 1:T_years) {
  cols_t <- ((t-1) * J + 1):(t * J)
  has_data[, t] <- apply(y_orig[, cols_t], 1, function(x) any(!is.na(x)))
  has_det[, t] <- apply(y_orig[, cols_t], 1, function(x) any(x == 1, na.rm = TRUE))
  n_reps[, t] <- apply(y_orig[, cols_t], 1, function(x) sum(!is.na(x)))
}
n_years <- rowSums(has_data)
has_early <- rowSums(has_data[, 1:3]) > 0
has_late <- rowSums(has_data[, 4:6]) > 0
ever_det <- rowSums(has_det) > 0

cat("Distribucion original de cobertura temporal:\n")
for (ny in 0:6) {
  n <- sum(n_years == ny)
  n_d <- sum(n_years == ny & ever_det)
  cat(sprintf("  %d anos: %5d sitios (%5.1f%%), %3d con detecciones\n",
              ny, n, 100 * n / n_orig, n_d))
}

###############################################################################
# B1: Filtros progresivos
###############################################################################

filters <- list(
  "Original (sin filtro)" = rep(TRUE, n_orig),
  ">= 2 anos" = n_years >= 2,
  ">= 3 anos" = n_years >= 3,
  ">= 3 anos + E+L" = n_years >= 3 & has_early & has_late,
  ">= 4 anos" = n_years >= 4,
  ">= 4 anos + E+L" = n_years >= 4 & has_early & has_late
)

cat("\n")
cat(sprintf("%-25s %6s %5s %6s %6s %6s %6s %6s\n",
            "Filtro", "Sites", "Dets", "NA%", "Naiv1", "NaivT", "S.det", "%det"))
cat(strrep("-", 85), "\n")

null_results <- list()

for (fname in names(filters)) {
  keep <- filters[[fname]]

  if (sum(keep) < 10) {
    cat(sprintf("%-25s %6d  -- demasiado pocos sitios --\n", fname, sum(keep)))
    next
  }

  y_f <- y_orig[keep, ]
  n_f <- nrow(y_f)
  na_pct <- round(100 * mean(is.na(y_f)), 1)
  total_det <- sum(y_f == 1, na.rm = TRUE)
  sites_det <- sum(apply(y_f, 1, function(x) any(x == 1, na.rm = TRUE)))

  # Naive year 1
  naive1 <- mean(apply(y_f[, 1:J], 1, function(x) any(x == 1, na.rm = TRUE)))
  # Naive total (any year)
  naiveT <- mean(apply(y_f, 1, function(x) any(x == 1, na.rm = TRUE)))

  cat(sprintf("%-25s %6d %5d %5.1f%% %5.1f%% %5.1f%% %5d %5.1f%%\n",
              fname, n_f, total_det, na_pct,
              naive1 * 100, naiveT * 100,
              sites_det, 100 * sites_det / n_f))

  # Ajustar modelo null
  umf_f <- unmarkedMultFrame(y = y_f, numPrimary = T_years)
  mod_f <- tryCatch({
    colext(psiformula = ~1, gammaformula = ~1,
           epsilonformula = ~1, pformula = ~1,
           data = umf_f)
  }, error = function(e) NULL)

  if (!is.null(mod_f)) {
    null_results[[fname]] <- mod_f
  }
}

###############################################################################
# B2: Mostrar resultados de todos los modelos null
###############################################################################

cat("\n\n")
cat(strrep("=", 70), "\n")
cat("  RESULTADOS: Modelos null con diferentes filtros temporales\n")
cat(strrep("=", 70), "\n\n")

cat(sprintf("%-25s %7s %7s %7s %7s %5s %5s\n",
            "Filtro", "psi1", "gamma", "epsil", "p", "AIC", "Conv"))
cat(strrep("-", 70), "\n")

for (fname in names(null_results)) {
  mod_f <- null_results[[fname]]
  coefs <- coef(mod_f)

  psi_p <- plogis(coefs["psi(Int)"])
  col_p <- plogis(coefs["col(Int)"])
  ext_p <- plogis(coefs["ext(Int)"])
  det_p <- plogis(coefs["p(Int)"])

  cat(sprintf("%-25s %6.1f%% %6.2f%% %6.1f%% %6.1f%% %5.0f %5d\n",
              fname,
              psi_p * 100, col_p * 100, ext_p * 100, det_p * 100,
              mod_f@AIC, mod_f@opt$convergence))
}

cat("\nReferencia biologica plausible:\n")
cat("  psi1: 5-30% | gamma: 1-15% | epsilon: 5-25% | p: 10-50%\n\n")

###############################################################################
# B3: Detalle del mejor modelo
###############################################################################

# Buscar el filtro con resultados mas plausibles
cat("\n")
cat(strrep("=", 70), "\n")
cat("  DETALLE: Modelos null individuales\n")
cat(strrep("=", 70), "\n")

for (fname in names(null_results)) {
  mod_f <- null_results[[fname]]
  cat("\n--- ", fname, " ---\n")
  s <- summary(mod_f)
  cat("  psi (ocupacion inicial):\n"); print(s$psi)
  cat("  col (colonizacion):\n");     print(s$col)
  cat("  ext (extincion):\n");        print(s$ext)
  cat("  det (deteccion):\n");        print(s$det)
  cat("  AIC: ", mod_f@AIC, "\n")
  cat("  Convergence: ", mod_f@opt$convergence, "\n")

  # Chequear vcov para correlaciones altas
  vc <- tryCatch(vcov(mod_f), error = function(e) NULL)
  if (!is.null(vc)) {
    se_v <- sqrt(diag(vc))
    cor_m <- vc / (se_v %o% se_v)
    high <- which(abs(cor_m) > 0.8 & upper.tri(cor_m), arr.ind = TRUE)
    if (nrow(high) > 0) {
      cat("  Correlaciones > 0.8:\n")
      for (r in 1:nrow(high)) {
        cat(sprintf("    %s <-> %s: r = %.3f\n",
                    rownames(cor_m)[high[r, 1]],
                    colnames(cor_m)[high[r, 2]],
                    cor_m[high[r, 1], high[r, 2]]))
      }
    } else {
      cat("  Sin correlaciones > 0.8 entre parametros\n")
    }
  }
}

###############################################################################
# CONCLUSIONES
###############################################################################

cat("\n\n")
cat(strrep("=", 70), "\n")
cat("  CONCLUSIONES\n")
cat(strrep("=", 70), "\n\n")

cat("1. EFECTO de filter_repeat_visits():\n")
cat("   - annual_closure=TRUE crea sitios POR ANO, asi que cada sitio\n")
cat("     tiene datos en exactamente 1 ano por definicion.\n")
cat("   - La cobertura temporal multianual viene del full_join posterior\n")
cat("     por locality_id en step 2.\n")
cat("   - Quitar observer_id de site_vars AUMENTA las replicas por\n")
cat("     localidad-ano (pooling observadores), pero no cambia la\n")
cat("     cobertura temporal directamente.\n\n")

cat("2. EFECTO del filtro temporal en el modelo null:\n")
cat("   - Comparar las lineas de arriba para ver si filtrar por >= 3 anos\n")
cat("     produce parametros plausibles.\n")
cat("   - Si epsilon sigue > 50% incluso con el filtro, el problema\n")
cat("     requiere cambio de escala espacial.\n\n")

cat("3. SIGUIENTE PASO:\n")
cat("   - Si el filtro temporal mejora epsilon significativamente,\n")
cat("     implementar el filtro en el pipeline completo.\n")
cat("   - Si no, probar escala espacial de 10 km.\n")
