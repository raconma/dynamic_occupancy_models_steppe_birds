###############################################################################
# test_sensitivity_tree_cover.R
#
# Analisis de sensibilidad: colocacion de tree_cover en psi vs epsilon
#
# Pregunta: tree_cover tiene un efecto enorme en psi (estimacion ~ -12,
# prob ~ 0%) y tambien aparece como predictor significativo de epsilon.
# Pero si los sitios con alto tree_cover nunca estan ocupados (psi ~ 0),
# epsilon alli es irrelevante. ¿Es el efecto en epsilon redundante?
#
# Variantes:
#   s0: ACTUAL (v3 best) - tree_cover en psi Y epsilon
#   s1: tree_cover solo en psi
#   s2: tree_cover solo en epsilon
#   s3: tree_cover en ninguno
#   s4: tree_cover en psi, solo tree (sin grass) en epsilon
#   s5: tree_cover solo en epsilon (sin grass), no en psi
#
# Todos con gam(~1) y p(~effort+observers+duration+time)
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
cat("  SENSIBILIDAD: tree_cover en psi vs epsilon\n")
cat("  ¿Es redundante tener tree_cover en ambos submodelos?\n")
cat(strrep("=", 70), "\n\n")

###############################################################################
# PASO 1: Cargar datos crudos
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
# PASO 2: filter_repeat_visits()
###############################################################################

cat("--- PASO 2: filter_repeat_visits ---\n\n")

occ <- filter_repeat_visits(
  occ_raw,
  min_obs  = 2,
  max_obs  = 10,
  annual_closure = TRUE,
  date_var = "observation_date",
  site_vars = c("locality_id")
)

cat("  Registros tras filtro: ", nrow(occ), "\n")
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
cat("  Registros con covariables validas: ", nrow(occ_var), "\n\n")

###############################################################################
# PASO 4: Formato wide por anio + merge
###############################################################################

cat("--- PASO 4: Convertir a formato wide ---\n\n")

vars_to_scale <- c("bio1", "tree_cover", "bio2", "grass_cover",
                    "topo_aspect", "topo_elev")
obs_covs <- c("time_observations_started", "duration_minutes",
              "effort_distance_km", "number_observers")
site_covs_wide <- c("locality_id",
                     "latitude", "longitude", vars_to_scale)

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
    occ_yr,
    site_id = "site",
    response = "species_observed",
    site_covs = site_covs_wide,
    obs_covs = obs_covs
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

occ_wide_clean <- occ_wide %>%
  distinct(locality_id, .keep_all = TRUE)

cat("  Sitios totales: ", nrow(occ_wide_clean), "\n\n")

###############################################################################
# PASO 5: Construir matrices de deteccion
###############################################################################

cat("--- PASO 5: Construir matrices ---\n\n")

J_reps <- 10

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

###############################################################################
# PASO 6: Filtro temporal
###############################################################################

cat("--- PASO 6: Filtro temporal ---\n\n")

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

keep <- n_years_data >= 3 & has_early & has_late

y_filt <- y.cross[keep, ]
duration_filt <- duration[keep, ]
effort_filt <- effort[keep, ]
observers_filt <- observers[keep, ]
time_filt <- time_obs[keep, ]
siteCovs_filt <- occ_wide_clean[keep, vars_to_scale]

n_filt <- nrow(y_filt)
cat("  Sitios tras filtro: ", n_filt, "\n\n")

###############################################################################
# PASO 7: Construir unmarkedMultFrame (sin yearlySiteCovs)
###############################################################################

cat("--- PASO 7: Construir UMF estatico ---\n\n")

duration_sc  <- scale(duration_filt)
effort_sc    <- scale(effort_filt)
observers_sc <- scale(observers_filt)
time_sc      <- scale(time_filt)

siteCovs_df <- as.data.frame(siteCovs_filt)

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

cat("  UMF: ", n_filt, " sitios, ", T_years, " periodos, ", J_reps, " replicas\n\n")

###############################################################################
# PASO 8: Ajustar 6 variantes de tree_cover
###############################################################################

cat("--- PASO 8: Ajustar 6 variantes de tree_cover ---\n\n")

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

p_full <- ~effort + observers + duration + time

# s0: ACTUAL (v3 best) - tree_cover en psi Y epsilon
s0 <- fit_model(
  ~bio1 + bio2 + tree_cover + grass_cover + topo_elev, ~1,
  ~grass_cover + tree_cover, p_full, occ_umf,
  "[s0] ACTUAL: psi(tree+) eps(grass+tree)")

# s1: tree_cover solo en psi
s1 <- fit_model(
  ~bio1 + bio2 + tree_cover + grass_cover + topo_elev, ~1,
  ~grass_cover, p_full, occ_umf,
  "[s1] tree solo en psi: psi(tree+) eps(grass)")

# s2: tree_cover solo en epsilon
s2 <- fit_model(
  ~bio1 + bio2 + grass_cover + topo_elev, ~1,
  ~grass_cover + tree_cover, p_full, occ_umf,
  "[s2] tree solo en eps: psi(sin tree) eps(grass+tree)")

# s3: tree_cover en ninguno
s3 <- fit_model(
  ~bio1 + bio2 + grass_cover + topo_elev, ~1,
  ~grass_cover, p_full, occ_umf,
  "[s3] sin tree: psi(sin tree) eps(grass)")

# s4: tree en psi y epsilon, pero eps solo tree (sin grass)
s4 <- fit_model(
  ~bio1 + bio2 + tree_cover + grass_cover + topo_elev, ~1,
  ~tree_cover, p_full, occ_umf,
  "[s4] psi(tree+) eps(solo tree)")

# s5: tree solo en epsilon, sin grass en eps
s5 <- fit_model(
  ~bio1 + bio2 + grass_cover + topo_elev, ~1,
  ~tree_cover, p_full, occ_umf,
  "[s5] psi(sin tree) eps(solo tree)")

###############################################################################
# PASO 9: Tabla AIC comparativa
###############################################################################

cat("\n--- PASO 9: Comparacion AIC ---\n\n")

models <- list(
  "s0: psi(tree+) eps(grass+tree)"  = s0,
  "s1: psi(tree+) eps(grass)"       = s1,
  "s2: psi(-) eps(grass+tree)"      = s2,
  "s3: psi(-) eps(grass)"           = s3,
  "s4: psi(tree+) eps(tree)"        = s4,
  "s5: psi(-) eps(tree)"            = s5
)

model_desc <- list(
  "s0: psi(tree+) eps(grass+tree)"  = "tree en AMBOS (actual v3)",
  "s1: psi(tree+) eps(grass)"       = "tree solo en PSI",
  "s2: psi(-) eps(grass+tree)"      = "tree solo en EPS",
  "s3: psi(-) eps(grass)"           = "tree en NINGUNO",
  "s4: psi(tree+) eps(tree)"        = "tree en ambos, eps sin grass",
  "s5: psi(-) eps(tree)"            = "tree solo en eps, sin grass"
)

aic_table <- data.frame(
  Model = character(), Description = character(), nPars = integer(),
  AIC = numeric(), Conv = integer(), Status = character(),
  stringsAsFactors = FALSE
)

for (mname in names(models)) {
  mod <- models[[mname]]
  chk <- check_model(mod, mname)
  if (!is.null(mod)) {
    aic_table <- rbind(aic_table, data.frame(
      Model = mname, Description = model_desc[[mname]],
      nPars = length(coef(mod)), AIC = round(mod@AIC, 1),
      Conv = mod@opt$convergence, Status = chk$reason,
      stringsAsFactors = FALSE))
  } else {
    aic_table <- rbind(aic_table, data.frame(
      Model = mname, Description = model_desc[[mname]],
      nPars = NA, AIC = NA, Conv = NA, Status = "Error",
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

cat(sprintf("  %-40s %5s %7s %7s %6s  %-30s\n",
            "Modelo", "nPar", "AIC", "dAIC", "Peso", "Estado"))
cat(strrep("-", 105), "\n")
for (i in 1:nrow(aic_table)) {
  cat(sprintf("  %-40s %5s %7s %7s %6s  %-30s\n",
              aic_table$Model[i],
              ifelse(is.na(aic_table$nPars[i]), " -", as.character(aic_table$nPars[i])),
              ifelse(is.na(aic_table$AIC[i]), "   -", sprintf("%.1f", aic_table$AIC[i])),
              ifelse(is.na(aic_table$deltaAIC[i]), "   -", sprintf("%.1f", aic_table$deltaAIC[i])),
              ifelse(is.na(aic_table$Weight[i]), "  -", sprintf("%.4f", aic_table$Weight[i])),
              aic_table$Status[i]))
}

###############################################################################
# PASO 10: Coeficientes de tree_cover en cada modelo
###############################################################################

cat("\n\n--- PASO 10: Coeficientes de tree_cover en cada modelo ---\n\n")

for (mname in names(models)) {
  mod <- models[[mname]]
  if (is.null(mod)) { cat(sprintf("  %s: NO AJUSTADO\n", mname)); next }

  s <- tryCatch(summary(mod), error = function(e) NULL)
  if (is.null(s)) { cat(sprintf("  %s: HESSIANA SINGULAR\n", mname)); next }

  cat(sprintf("  === %s ===\n", mname))

  # tree_cover en psi?
  psi_df <- as.data.frame(s$psi)
  psi_df$Parameter <- rownames(psi_df)
  tree_psi <- psi_df[psi_df$Parameter == "tree_cover", ]
  if (nrow(tree_psi) > 0) {
    cat(sprintf("    psi(tree_cover): est=%.3f, SE=%.3f, z=%.2f, p=%.4f\n",
                tree_psi$Estimate, tree_psi$SE, tree_psi$z, tree_psi$`P(>|z|)`))
  } else {
    cat("    psi(tree_cover): NO incluido\n")
  }

  # tree_cover en eps?
  ext_df <- as.data.frame(s$ext)
  ext_df$Parameter <- rownames(ext_df)
  tree_ext <- ext_df[ext_df$Parameter == "tree_cover", ]
  if (nrow(tree_ext) > 0) {
    cat(sprintf("    eps(tree_cover): est=%.3f, SE=%.3f, z=%.2f, p=%.4f\n",
                tree_ext$Estimate, tree_ext$SE, tree_ext$z, tree_ext$`P(>|z|)`))
  } else {
    cat("    eps(tree_cover): NO incluido\n")
  }

  # grass_cover en eps?
  grass_ext <- ext_df[ext_df$Parameter == "grass_cover", ]
  if (nrow(grass_ext) > 0) {
    cat(sprintf("    eps(grass_cover): est=%.3f, SE=%.3f, z=%.2f, p=%.4f\n",
                grass_ext$Estimate, grass_ext$SE, grass_ext$z, grass_ext$`P(>|z|)`))
  }

  cat("\n")
}

###############################################################################
# PASO 11: Comparaciones clave
###############################################################################

cat("\n--- PASO 11: Interpretacion ---\n\n")

if (!is.null(s0) && !is.null(s1)) {
  cat(sprintf("  1. ¿tree_cover en eps aporta algo (cuando ya esta en psi)?\n"))
  cat(sprintf("     s0 (ambos) vs s1 (solo psi): dAIC = %.1f\n",
              s0@AIC - s1@AIC))
  if (s0@AIC < s1@AIC - 2) cat("     => SI, tree en eps MEJORA\n\n")
  else if (s0@AIC > s1@AIC + 2) cat("     => tree en eps EMPEORA\n\n")
  else cat("     => Diferencia no sustancial\n\n")
}

if (!is.null(s2) && !is.null(s3)) {
  cat(sprintf("  2. ¿tree_cover en eps aporta (cuando NO esta en psi)?\n"))
  cat(sprintf("     s2 (solo eps) vs s3 (ninguno): dAIC = %.1f\n",
              s2@AIC - s3@AIC))
  if (s2@AIC < s3@AIC - 2) cat("     => SI, tree en eps MEJORA cuando no esta en psi\n\n")
  else if (s2@AIC > s3@AIC + 2) cat("     => tree en eps EMPEORA\n\n")
  else cat("     => Diferencia no sustancial\n\n")
}

if (!is.null(s0) && !is.null(s2)) {
  cat(sprintf("  3. ¿tree_cover en psi aporta algo?\n"))
  cat(sprintf("     s0 (ambos) vs s2 (solo eps): dAIC = %.1f\n",
              s0@AIC - s2@AIC))
  if (s0@AIC < s2@AIC - 2) cat("     => SI, tree en psi MEJORA\n\n")
  else if (s0@AIC > s2@AIC + 2) cat("     => tree en psi EMPEORA\n\n")
  else cat("     => Diferencia no sustancial\n\n")
}

if (!is.null(s0) && !is.null(s3)) {
  cat(sprintf("  4. ¿tree_cover es util en general?\n"))
  cat(sprintf("     s0 (ambos) vs s3 (ninguno): dAIC = %.1f\n",
              s0@AIC - s3@AIC))
  if (s0@AIC < s3@AIC - 2) cat("     => SI, tree_cover es UTIL en algun submodelo\n\n")
  else cat("     => tree_cover no aporta nada\n\n")
}

###############################################################################
# PASO 12: Exportar resultados
###############################################################################

cat("\n--- PASO 12: Exportar ---\n\n")

results_dir <- file.path(
  "/Users/gfandos/Documents/GitHub/dynamic_occupancy_models_steppe_birds", "results"
)
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

write.csv(aic_table, file.path(results_dir, "otitar_sensitivity_tree_cover.csv"),
          row.names = FALSE)
cat("  Guardado: results/otitar_sensitivity_tree_cover.csv\n")

###############################################################################
# RESUMEN
###############################################################################

cat("\n\n")
cat(strrep("=", 70), "\n")
cat("  SENSIBILIDAD tree_cover — RESUMEN\n")
cat(strrep("=", 70), "\n\n")

best_idx <- which.min(aic_table$AIC)
cat(sprintf("  Mejor modelo: %s (AIC=%.1f)\n",
            aic_table$Model[best_idx], aic_table$AIC[best_idx]))
cat(sprintf("  Descripcion: %s\n", aic_table$Description[best_idx]))

cat("\n")
cat(strrep("=", 70), "\n")
cat("  FIN ANALISIS DE SENSIBILIDAD\n")
cat(strrep("=", 70), "\n\n")
