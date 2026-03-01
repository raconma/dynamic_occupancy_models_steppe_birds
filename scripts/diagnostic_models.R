###############################################################################
# diagnostic_models.R
#
# Purpose: Comprehensive diagnostic of 4 colext dynamic occupancy models.
#          Evaluates detection, raw data patterns, identifiability, and
#          potential causes of biologically implausible results.
###############################################################################

library(unmarked)

cat("\n##########################################################################\n")
cat("# INFORME DIAGNOSTICO - MODELOS DE OCUPACION DINAMICA (colext)\n")
cat("# Fecha:", as.character(Sys.time()), "\n")
cat("##########################################################################\n\n")

# --- Path to model objects ---
base_path <- "/Users/gfandos/Documents/GitHub/dynamic_occupancy_models_steppe_birds/data/processed"
species_codes <- c("otitar", "ptealc", "pteori", "tettet")
species_names <- c("Otis tarda", "Pterocles alchata", "Pterocles orientalis", "Tetrax tetrax")

YEARS <- 2017:2022
T_years <- length(YEARS)
J_reps <- 10

###############################################################################
# MAIN DIAGNOSTIC LOOP
###############################################################################

for (idx in seq_along(species_codes)) {
  sp <- species_codes[idx]
  sp_name <- species_names[idx]

  cat("\n\n")
  cat(strrep("=", 80), "\n")
  cat("DIAGNOSTICO:", sp, "(", sp_name, ")\n")
  cat(strrep("=", 80), "\n")

  # Load model
  mod_path <- file.path(base_path, paste0(sp, "_model_object.rds"))
  if (!file.exists(mod_path)) {
    cat("  ERROR: No se encontro el archivo:", mod_path, "\n")
    next
  }
  mod <- readRDS(mod_path)
  umf <- mod@data

  ############################################################################
  # 1. COEFICIENTES DEL MODELO
  ############################################################################
  cat("\n--- 1. COEFICIENTES DEL MODELO ---\n\n")

  # Extract summary
  s <- summary(mod)

  cat("  SUBMODELO: psi (ocupacion inicial)\n")
  psi_coefs <- s$psi
  print(psi_coefs)
  cat("\n")

  cat("  SUBMODELO: gamma (colonizacion)\n")
  col_coefs <- s$col
  print(col_coefs)
  cat("\n")

  cat("  SUBMODELO: epsilon (extincion)\n")
  ext_coefs <- s$ext
  print(ext_coefs)
  cat("\n")

  cat("  SUBMODELO: p (deteccion)\n")
  det_coefs <- s$det
  print(det_coefs)
  cat("\n")

  cat("  AIC:", mod@AIC, "\n")

  # Check for convergence warnings
  opt <- mod@opt
  cat("  Convergence code:", opt$convergence, "\n")
  if (opt$convergence != 0) {
    cat("  *** WARNING: MODELO NO CONVERGIO ***\n")
    cat("  Mensaje:", opt$message, "\n")
  }

  ############################################################################
  # 2. EVALUACION DE LA DETECCION (p)
  ############################################################################
  cat("\n--- 2. EVALUACION DE LA DETECCION ---\n\n")

  # Predicted detection across all observations
  p_pred <- tryCatch({
    predict(mod, type = "det")
  }, error = function(e) {
    cat("  Error al predecir deteccion:", e$message, "\n")
    NULL
  })

  if (!is.null(p_pred)) {
    p_vals <- p_pred$Predicted
    cat("  p media:", round(mean(p_vals, na.rm = TRUE), 4), "\n")
    cat("  p mediana:", round(median(p_vals, na.rm = TRUE), 4), "\n")
    cat("  p rango: [", round(min(p_vals, na.rm = TRUE), 4), ",",
        round(max(p_vals, na.rm = TRUE), 4), "]\n")
    cat("  p SD:", round(sd(p_vals, na.rm = TRUE), 4), "\n")

    # Quantiles
    p_q <- quantile(p_vals, probs = c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE)
    cat("  p cuantiles (5%, 25%, 50%, 75%, 95%):",
        paste(round(p_q, 4), collapse = ", "), "\n")

    # Interpretation
    p_mean <- mean(p_vals, na.rm = TRUE)
    if (p_mean < 0.05) {
      cat("  *** INTERPRETACION: MUY BAJA - el modelo NO puede separar ",
          "no-deteccion de ausencia real ***\n")
    } else if (p_mean < 0.1) {
      cat("  *** INTERPRETACION: BAJA - resultados poco fiables, ",
          "confusion psi/p probable ***\n")
    } else if (p_mean < 0.2) {
      cat("  INTERPRETACION: MODERADA-BAJA - resultados inciertos\n")
    } else {
      cat("  INTERPRETACION: ADECUADA\n")
    }

    # Detection covariate effects
    cat("\n  Efecto de covariables de deteccion:\n")
    det_df <- det_coefs
    for (r in seq_len(nrow(det_df))) {
      var_name <- rownames(det_df)[r]
      est <- det_df[r, "Estimate"]
      se <- det_df[r, "SE"]
      z <- det_df[r, "z"]
      p_val <- det_df[r, "P(>|z|)"]
      sig <- ifelse(p_val < 0.001, "***", ifelse(p_val < 0.01, "**",
             ifelse(p_val < 0.05, "*", "ns")))
      cat(sprintf("    %-25s Estimate=%7.4f  SE=%7.4f  z=%7.3f  p=%s %s\n",
                  var_name, est, se, z, format.pval(p_val, digits = 3), sig))
    }
  }

  ############################################################################
  # 3. DATOS CRUDOS
  ############################################################################
  cat("\n--- 3. DATOS CRUDOS ---\n\n")

  y <- umf@y
  n_sites <- nrow(y)
  n_total_cols <- ncol(y)

  cat("  Dimension de la matriz y:", n_sites, "sitios x", n_total_cols,
      "observaciones\n")
  cat("  Periodos primarios:", T_years, "| Replicas por periodo:", J_reps, "\n")

  # Total NAs in y matrix
  total_cells <- prod(dim(y))
  na_cells <- sum(is.na(y))
  cat("  Total celdas en y:", total_cells, "\n")
  cat("  Celdas NA:", na_cells, "(", round(100 * na_cells / total_cells, 1), "%)\n")

  # Per-year analysis
  cat("\n  Ocupacion naive y detecciones por anio:\n")
  cat(sprintf("  %-6s %8s %12s %12s %10s %8s\n",
              "Anio", "Naive%", "Detecciones", "Sitios_det", "Rep_validas", "NA%"))

  naive_by_year <- numeric(T_years)
  det_by_year <- numeric(T_years)
  sites_det_by_year <- numeric(T_years)
  valid_reps_by_year <- numeric(T_years)
  na_pct_by_year <- numeric(T_years)

  for (t in 1:T_years) {
    cols <- ((t - 1) * J_reps + 1):(t * J_reps)
    y_year <- y[, cols]

    # Naive occupancy: at least one detection
    any_det <- apply(y_year, 1, function(x) any(x == 1, na.rm = TRUE))
    naive_by_year[t] <- mean(any_det, na.rm = TRUE)

    # Total detections
    det_by_year[t] <- sum(y_year == 1, na.rm = TRUE)

    # Sites with at least one detection
    sites_det_by_year[t] <- sum(any_det, na.rm = TRUE)

    # Mean valid (non-NA) replicates per site
    valid_per_site <- apply(y_year, 1, function(x) sum(!is.na(x)))
    valid_reps_by_year[t] <- mean(valid_per_site)

    # NA percentage this year
    na_pct_by_year[t] <- mean(is.na(y_year)) * 100

    cat(sprintf("  %-6d %7.1f%% %12d %12d %10.1f %7.1f%%\n",
                YEARS[t],
                naive_by_year[t] * 100,
                det_by_year[t],
                sites_det_by_year[t],
                valid_reps_by_year[t],
                na_pct_by_year[t]))
  }

  cat("\n  Total sitios:", n_sites, "\n")
  cat("  Detecciones totales:", sum(y == 1, na.rm = TRUE), "\n")
  cat("  Ceros totales:", sum(y == 0, na.rm = TRUE), "\n")
  cat("  Ratio detecciones/no-detecciones:",
      round(sum(y == 1, na.rm = TRUE) / sum(y == 0, na.rm = TRUE), 4), "\n")

  ############################################################################
  # 4. COMPARACION NAIVE vs MODELO
  ############################################################################
  cat("\n--- 4. COMPARACION NAIVE vs MODELO ---\n\n")

  # psi1 predicted
  psi_pred <- tryCatch({
    predict(mod, type = "psi")
  }, error = function(e) {
    cat("  Error al predecir psi:", e$message, "\n")
    NULL
  })

  if (!is.null(psi_pred)) {
    psi1_mean <- mean(psi_pred$Predicted, na.rm = TRUE)
    psi1_median <- median(psi_pred$Predicted, na.rm = TRUE)
    psi1_range <- range(psi_pred$Predicted, na.rm = TRUE)

    cat("  psi1 estimado (media):", round(psi1_mean, 4), "\n")
    cat("  psi1 estimado (mediana):", round(psi1_median, 4), "\n")
    cat("  psi1 estimado (rango): [", round(psi1_range[1], 4), ",",
        round(psi1_range[2], 4), "]\n")
    cat("  Naive anio 1:", round(naive_by_year[1], 4), "\n")

    discrepancy <- abs(psi1_mean - naive_by_year[1])
    cat("  Discrepancia |psi1 - naive|:", round(discrepancy, 4), "\n")

    if (psi1_mean < 0.05 && naive_by_year[1] > 0.05) {
      cat("  *** PROBLEMA: psi1 muy bajo pero hay detecciones en anio 1 ***\n")
      cat("  Posible confusion entre deteccion y ocupacion inicial\n")
    }
    if (psi1_mean > naive_by_year[1] * 3) {
      cat("  *** PROBLEMA: psi1 mucho mayor que naive - modelo sobreestima ***\n")
    }
  }

  # Colonization and extinction predicted
  col_pred <- tryCatch({
    predict(mod, type = "col")
  }, error = function(e) NULL)

  ext_pred <- tryCatch({
    predict(mod, type = "ext")
  }, error = function(e) NULL)

  if (!is.null(col_pred)) {
    cat("\n  Colonizacion estimada (media):", round(mean(col_pred$Predicted, na.rm = TRUE), 4), "\n")
    cat("  Colonizacion estimada (rango): [",
        round(min(col_pred$Predicted, na.rm = TRUE), 4), ",",
        round(max(col_pred$Predicted, na.rm = TRUE), 4), "]\n")
  }

  if (!is.null(ext_pred)) {
    cat("  Extincion estimada (media):", round(mean(ext_pred$Predicted, na.rm = TRUE), 4), "\n")
    cat("  Extincion estimada (rango): [",
        round(min(ext_pred$Predicted, na.rm = TRUE), 4), ",",
        round(max(ext_pred$Predicted, na.rm = TRUE), 4), "]\n")
  }

  ############################################################################
  # 5. PATRONES SOSPECHOSOS EN LOS DATOS
  ############################################################################
  cat("\n--- 5. PATRONES SOSPECHOSOS ---\n\n")

  # 5a. Sites detected only in late years
  first_det_year <- apply(y, 1, function(row) {
    for (t in 1:T_years) {
      cols <- ((t - 1) * J_reps + 1):(t * J_reps)
      if (any(row[cols] == 1, na.rm = TRUE)) return(t)
    }
    return(NA)
  })

  last_det_year <- apply(y, 1, function(row) {
    for (t in T_years:1) {
      cols <- ((t - 1) * J_reps + 1):(t * J_reps)
      if (any(row[cols] == 1, na.rm = TRUE)) return(t)
    }
    return(NA)
  })

  ever_detected <- !is.na(first_det_year)
  n_ever_detected <- sum(ever_detected)

  cat("  Sitios con al menos una deteccion:", n_ever_detected,
      "(", round(100 * n_ever_detected / n_sites, 1), "%)\n")
  cat("  Sitios NUNCA detectados:", sum(!ever_detected),
      "(", round(100 * sum(!ever_detected) / n_sites, 1), "%)\n")

  # Late-only sites: first detection in years 4-6 but NOT in years 1-3
  late_only <- first_det_year >= 4 & !is.na(first_det_year)
  n_late_only <- sum(late_only, na.rm = TRUE)
  cat("  Sitios detectados SOLO en anios tardios (4-6):", n_late_only,
      "(", round(100 * n_late_only / n_ever_detected, 1), "% de sitios detectados)\n")

  # Sites detected only in year 1
  early_only <- last_det_year == 1 & first_det_year == 1 & !is.na(first_det_year)
  n_early_only <- sum(early_only, na.rm = TRUE)
  cat("  Sitios detectados SOLO en anio 1:", n_early_only,
      "(", round(100 * n_early_only / n_ever_detected, 1), "%)\n")

  # First detection year distribution
  cat("\n  Distribucion de primera deteccion por anio:\n")
  for (t in 1:T_years) {
    n_first <- sum(first_det_year == t, na.rm = TRUE)
    cat(sprintf("    Anio %d (%d): %d sitios (%.1f%%)\n",
                t, YEARS[t], n_first,
                100 * n_first / n_ever_detected))
  }

  # 5b. Singletons: sites with exactly one detection across all years
  total_det_per_site <- apply(y, 1, function(x) sum(x == 1, na.rm = TRUE))
  singletons <- sum(total_det_per_site == 1)
  cat("\n  Singletons (1 sola deteccion en todo el periodo):", singletons,
      "(", round(100 * singletons / n_sites, 1), "% de todos los sitios)\n")

  if (n_ever_detected > 0) {
    cat("  Singletons como % de sitios detectados:",
        round(100 * singletons / n_ever_detected, 1), "%\n")
  }

  # Sites with 2 detections
  n_2det <- sum(total_det_per_site == 2)
  cat("  Sitios con exactamente 2 detecciones:", n_2det, "\n")

  # 5c. Concentration of detections
  det_sorted <- sort(total_det_per_site[total_det_per_site > 0], decreasing = TRUE)
  if (length(det_sorted) > 0) {
    total_det <- sum(det_sorted)
    # Top 10% of sites
    top10_n <- max(1, round(0.1 * length(det_sorted)))
    top10_det <- sum(det_sorted[1:top10_n])
    cat("  Concentracion: top 10% de sitios ocupados tienen",
        round(100 * top10_det / total_det, 1), "% de las detecciones\n")

    # Top 25%
    top25_n <- max(1, round(0.25 * length(det_sorted)))
    top25_det <- sum(det_sorted[1:top25_n])
    cat("  Concentracion: top 25% de sitios ocupados tienen",
        round(100 * top25_det / total_det, 1), "% de las detecciones\n")

    cat("  Detecciones maximas por sitio:", max(det_sorted), "\n")
    cat("  Detecciones medianas por sitio (solo sitios detectados):",
        round(median(det_sorted), 1), "\n")
  }

  # 5d. Temporal coverage: do sites have more surveys in later years?
  cat("\n  Cobertura de muestreo (replicas validas) por anio:\n")
  for (t in 1:T_years) {
    cols <- ((t - 1) * J_reps + 1):(t * J_reps)
    valid <- apply(y[, cols], 1, function(x) sum(!is.na(x)))
    sites_with_any <- sum(valid > 0)
    cat(sprintf("    %d: media=%.1f replicas, %d sitios con datos (%.0f%%), max=%d\n",
                YEARS[t], mean(valid), sites_with_any,
                100 * sites_with_any / n_sites, max(valid)))
  }

  # 5e. Sites with data only in some years
  has_data_year <- matrix(FALSE, n_sites, T_years)
  for (t in 1:T_years) {
    cols <- ((t - 1) * J_reps + 1):(t * J_reps)
    has_data_year[, t] <- apply(y[, cols], 1, function(x) any(!is.na(x)))
  }
  n_years_with_data <- rowSums(has_data_year)
  cat("\n  Sitios por numero de anios con datos:\n")
  for (ny in 0:T_years) {
    n_ny <- sum(n_years_with_data == ny)
    cat(sprintf("    %d anios: %d sitios (%.1f%%)\n",
                ny, n_ny, 100 * n_ny / n_sites))
  }

  ############################################################################
  # 6. IDENTIFICABILIDAD DEL MODELO
  ############################################################################
  cat("\n--- 6. IDENTIFICABILIDAD DEL MODELO ---\n\n")

  # Variance-covariance matrix
  vc <- tryCatch({
    vcov(mod)
  }, error = function(e) {
    cat("  Error al calcular vcov:", e$message, "\n")
    NULL
  })

  if (!is.null(vc)) {
    # Correlation matrix
    se_vec <- sqrt(diag(vc))
    cor_mat <- vc / (se_vec %o% se_vec)

    # Find high correlations (>0.7)
    n_params <- nrow(cor_mat)
    high_cors <- list()
    for (i in 1:(n_params - 1)) {
      for (j in (i + 1):n_params) {
        if (!is.na(cor_mat[i, j]) && abs(cor_mat[i, j]) > 0.7) {
          high_cors[[length(high_cors) + 1]] <- list(
            par1 = rownames(cor_mat)[i],
            par2 = colnames(cor_mat)[j],
            cor = cor_mat[i, j]
          )
        }
      }
    }

    if (length(high_cors) > 0) {
      cat("  *** CORRELACIONES ALTAS ENTRE PARAMETROS (|r| > 0.7): ***\n")
      for (hc in high_cors) {
        cat(sprintf("    %s <-> %s : r = %.4f %s\n",
                    hc$par1, hc$par2, hc$cor,
                    ifelse(abs(hc$cor) > 0.9, "*** CRITICA ***", "")))
      }
    } else {
      cat("  No se encontraron correlaciones > 0.7 entre parametros\n")
    }

    # Check for extreme SE values
    cat("\n  Parametros con SE extremos:\n")
    all_coefs <- coef(mod)
    all_se <- se_vec
    names(all_se) <- names(all_coefs)

    for (k in seq_along(all_coefs)) {
      if (all_se[k] > 10) {
        cat(sprintf("    %s: coef=%.4f, SE=%.4f *** SE MUY GRANDE ***\n",
                    names(all_coefs)[k], all_coefs[k], all_se[k]))
      }
    }

    # Check for boundary parameters
    cat("\n  Parametros en boundary (probabilidad ~0 o ~1):\n")
    # For logit-scale parameters, check if back-transformed they're near 0 or 1
    for (k in seq_along(all_coefs)) {
      p_val <- plogis(all_coefs[k])
      if (p_val < 0.001 || p_val > 0.999) {
        cat(sprintf("    %s: logit=%.4f -> prob=%.6f %s\n",
                    names(all_coefs)[k], all_coefs[k], p_val,
                    ifelse(p_val < 0.001, "(~0)", "(~1)")))
      }
    }
  }

  # Check the intercepts specifically
  cat("\n  Intercepts en escala de probabilidad:\n")
  all_coefs <- coef(mod)

  # psi intercept
  psi_int <- all_coefs[grep("^psi\\(Int", names(all_coefs))]
  if (length(psi_int) > 0) {
    cat(sprintf("    psi(Intercept): logit=%.4f -> prob=%.4f\n",
                psi_int, plogis(psi_int)))
  }

  # col intercept
  col_int <- all_coefs[grep("^col\\(Int", names(all_coefs))]
  if (length(col_int) > 0) {
    cat(sprintf("    col(Intercept): logit=%.4f -> prob=%.4f\n",
                col_int, plogis(col_int)))
  }

  # ext intercept
  ext_int <- all_coefs[grep("^ext\\(Int", names(all_coefs))]
  if (length(ext_int) > 0) {
    cat(sprintf("    ext(Intercept): logit=%.4f -> prob=%.4f\n",
                ext_int, plogis(ext_int)))
  }

  # det intercept
  det_int <- all_coefs[grep("^p\\(Int", names(all_coefs))]
  if (length(det_int) > 0) {
    cat(sprintf("    p(Intercept): logit=%.4f -> prob=%.4f\n",
                det_int, plogis(det_int)))
  }

  ############################################################################
  # 7. HIPOTESIS DIAGNOSTICAS
  ############################################################################
  cat("\n--- 7. HIPOTESIS DIAGNOSTICAS ---\n\n")

  # Calculate key metrics
  if (!is.null(p_pred)) p_mean_val <- mean(p_pred$Predicted, na.rm = TRUE) else p_mean_val <- NA
  if (!is.null(psi_pred)) psi1_val <- mean(psi_pred$Predicted, na.rm = TRUE) else psi1_val <- NA
  if (!is.null(ext_pred)) ext_val <- mean(ext_pred$Predicted, na.rm = TRUE) else ext_val <- NA
  if (!is.null(col_pred)) col_val <- mean(col_pred$Predicted, na.rm = TRUE) else col_val <- NA

  # H1: Detection confounded with initial occupancy
  cat("  H1 - Deteccion confundida con psi1:\n")
  if (!is.na(p_mean_val) && p_mean_val < 0.15) {
    cat("    *** PROBABLE *** p=", round(p_mean_val, 3),
        " es baja. Con deteccion baja y variable,\n")
    cat("    el modelo atribuye toda la ocupacion a colonizacion tardia\n")
    cat("    en lugar de presencia inicial no detectada.\n")
  } else {
    cat("    p=", round(p_mean_val, 3), " - deteccion moderada/adecuada,\n")
    cat("    menos probable esta hipotesis.\n")
  }

  # H2: Temporal imbalance in sampling
  cat("\n  H2 - Desbalance temporal de muestreo:\n")
  late_surveys <- sum(!is.na(y[, ((4 - 1) * J_reps + 1):(T_years * J_reps)]))
  early_surveys <- sum(!is.na(y[, 1:(3 * J_reps)]))
  ratio_late_early <- late_surveys / max(early_surveys, 1)
  cat("    Observaciones validas anios 1-3:", early_surveys, "\n")
  cat("    Observaciones validas anios 4-6:", late_surveys, "\n")
  cat("    Ratio tardio/temprano:", round(ratio_late_early, 2), "\n")
  if (ratio_late_early > 1.5) {
    cat("    *** PROBABLE *** Hay", round(ratio_late_early, 1),
        "x mas datos en anios tardios.\n")
    cat("    Sitios 'aparecen' como colonizaciones cuando en realidad\n")
    cat("    simplemente no se muestrearon antes.\n")
  } else {
    cat("    Ratio equilibrado - menos probable esta hipotesis.\n")
  }

  # H3: Spatial scale
  cat("\n  H3 - Escala espacial inadecuada:\n")
  cat("    Analisis requiere datos espaciales no disponibles en el modelo objeto.\n")
  cat("    Celdas de 2.5 km pueden ser demasiado pequenas para especies con\n")
  cat("    home ranges grandes (avutarda: ~10-50 km2).\n")

  # H4: Site filtering bias
  cat("\n  H4 - Sesgo de filtrado de sitios:\n")
  sites_all_years <- sum(n_years_with_data == T_years)
  sites_few_years <- sum(n_years_with_data <= 2)
  cat("    Sitios con datos todos los anios:", sites_all_years,
      "(", round(100 * sites_all_years / n_sites, 1), "%)\n")
  cat("    Sitios con datos en 1-2 anios:", sites_few_years,
      "(", round(100 * sites_few_years / n_sites, 1), "%)\n")
  cat("    Singletons:", singletons, "(", round(100 * singletons / n_sites, 1), "%)\n")
  if (singletons / n_sites > 0.3) {
    cat("    *** PROBABLE *** >30% de singletons sugiere sitios con\n")
    cat("    detecciones esporadicas retenidos por el filtrado.\n")
  }

  ############################################################################
  # 8. DIAGNOSTICO PROBABLE
  ############################################################################
  cat("\n--- 8. DIAGNOSTICO PROBABLE ---\n\n")

  # Score each hypothesis
  problems <- character(0)

  if (!is.na(p_mean_val) && p_mean_val < 0.15) {
    problems <- c(problems, "Deteccion muy baja (p<0.15)")
  }
  if (!is.na(ext_val) && ext_val > 0.4) {
    problems <- c(problems, paste0("Extincion implausiblemente alta (",
                                    round(ext_val * 100, 0), "%)"))
  }
  if (!is.na(psi1_val) && naive_by_year[1] > 0.05 && psi1_val < 0.05) {
    problems <- c(problems, "psi1 ~ 0 pero naive > 5%")
  }
  if (ratio_late_early > 1.5) {
    problems <- c(problems, paste0("Desbalance temporal (ratio=",
                                    round(ratio_late_early, 1), ")"))
  }
  if (singletons / n_sites > 0.2) {
    problems <- c(problems, paste0("Alta proporcion de singletons (",
                                    round(100 * singletons / n_sites, 0), "%)"))
  }
  if (n_late_only / max(n_ever_detected, 1) > 0.3) {
    problems <- c(problems, paste0("Muchos sitios aparecen solo en anios tardios (",
                                    round(100 * n_late_only / n_ever_detected, 0), "%)"))
  }

  if (length(problems) > 0) {
    cat("  Problemas identificados:\n")
    for (prob in problems) {
      cat("    - ", prob, "\n")
    }
  } else {
    cat("  No se identificaron problemas graves.\n")
  }

  ############################################################################
  # 9. RECOMENDACION
  ############################################################################
  cat("\n--- 9. RECOMENDACION ---\n\n")

  if (!is.na(p_mean_val) && p_mean_val < 0.15) {
    cat("  1. REDUCIR complejidad del modelo de deteccion o usar ",
        "informative priors.\n")
    cat("  2. Considerar un modelo null para p para evaluar si las ",
        "covariables empeoran la identificabilidad.\n")
  }
  if (ratio_late_early > 1.5) {
    cat("  3. BALANCEAR temporalmente: filtrar sitios sin datos en anios ",
        "tempranos o\n")
    cat("     requerir un minimo de datos por anio.\n")
  }
  if (singletons / n_sites > 0.2) {
    cat("  4. CAMBIAR filtrado de sitios: requerir minimo 2 detecciones y\n")
    cat("     presencia en multiples anios.\n")
  }
  if (!is.na(ext_val) && ext_val > 0.4) {
    cat("  5. EVALUAR escala espacial: probar celdas de 5 o 10 km.\n")
    cat("  6. CONSIDERAR que extincion alta puede reflejar movimiento entre\n")
    cat("     celdas, no extincion real.\n")
  }

  cat("\n")
}

###############################################################################
# CONCLUSION GENERAL
###############################################################################
cat("\n\n")
cat(strrep("#", 80), "\n")
cat("# CONCLUSION GENERAL\n")
cat(strrep("#", 80), "\n\n")

cat("Se analizaron los 4 modelos de ocupacion dinamica (colext) para aves\n")
cat("esteparias en Espana. Los resultados de cada especie se presentaron arriba.\n\n")

cat("PATRON COMUN ESPERADO:\n")
cat("- Tasas de extincion implausiblemente altas (>30-60%)\n")
cat("- Colonizacion muy baja (<1%)\n")
cat("- Resultado neto: poblaciones en colapso, lo cual contradice la realidad\n\n")

cat("CAUSAS MAS PROBABLES (en orden de prioridad):\n\n")
cat("1. CONFUSION PSI/P: Con deteccion baja (p<0.15), el modelo no puede\n")
cat("   distinguir entre 'no detectado pero presente' y 'ausente'. Esto\n")
cat("   subestima psi1 y sobreestima gamma (colonizacion aparente).\n\n")
cat("2. DESBALANCE TEMPORAL: Si hay mas listas eBird en 2020-2022 que en\n")
cat("   2017-2018, los sitios 'aparecen' en anios tardios. El modelo\n")
cat("   interpreta esto como colonizacion + extincion alta (turnover).\n\n")
cat("3. ESCALA ESPACIAL: Celdas de 2.5 km son demasiado pequenas para\n")
cat("   especies con home ranges de 10-50 km2. Los animales se mueven\n")
cat("   entre celdas -> el modelo ve 'extincion' donde hay movimiento.\n\n")
cat("4. FILTRADO DE SITIOS: El filtro '3-10 avistamientos' puede retener\n")
cat("   sitios con detecciones esporadicas y eliminar sitios con ocupacion\n")
cat("   estable, sesgando hacia alta rotacion.\n\n")

cat("RECOMENDACION PRIORITARIA:\n")
cat("1. Aumentar la escala espacial a 5 o 10 km\n")
cat("2. Eliminar el filtro '3-10 avistamientos' por sitio\n")
cat("3. Requerir que cada sitio tenga datos en al menos 3 de los 6 anios\n")
cat("4. Evaluar un modelo con p constante (~ 1) como baseline\n")
cat("5. Si p < 0.1, considerar modelos Bayesianos con priors informativos\n")
cat("   para p basados en protocolos de campo conocidos.\n\n")

cat("FIN DEL INFORME DIAGNOSTICO\n")
