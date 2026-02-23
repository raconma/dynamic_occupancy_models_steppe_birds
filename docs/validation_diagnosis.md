# Diagnóstico de la validación — 5_validation_otitar.R

**Fecha:** 2026-02-23
**Problema reportado por el alumno:**
1. Spearman ρ y R² "casi perfectos" → sospecha de sobreajuste
2. Moran's I significativo → autocorrelación espacial en residuos
3. Duda: ¿comparar predicciones de hábitat potencial vs datos de censo del atlas es válido?

---

## DIAGNÓSTICO

He identificado **4 problemas**, de los cuales **2 son errores metodológicos** y **2 son limitaciones esperables** que se deben documentar pero no invalidan el estudio.

---

### PROBLEMA 1 (ERROR) — Spearman y R² calculados sobre deciles, no sobre datos crudos

**Síntoma:** Spearman ρ ≈ 1.0, R² ≈ 0.98+

**Causa:** El calibration plot (líneas 109-125 del script actual) hace esto:

```r
atlas_deciles <- atlas_valid %>%
  mutate(prob_bin = ntile(pred_prob_mean, 10)) %>%    # ← agrupa en 10 bins
  group_by(prob_bin) %>%
  summarise(mean_prob = mean(pred_prob_mean),          # ← media de predicciones
            obs_prev  = mean(as.numeric(OTITAR)),      # ← prevalencia observada
            n = n())

spearman_corr <- cor.test(atlas_deciles$mean_prob,
                          atlas_deciles$obs_prev, method="spearman")
r2 <- summary(lm(obs_prev ~ mean_prob, data=atlas_deciles))$r.squared
```

**El problema:** Está calculando la correlación sobre **10 puntos agregados** (los deciles), no sobre los ~2000+ polígonos individuales del atlas. Al agrupar en deciles:
- Se eliminan el 99.5% de la varianza individual
- Cualquier modelo mínimamente razonable producirá una correlación casi perfecta entre medias de bins
- Es un artefacto estadístico, **no** sobreajuste real

Esto es como decir "la media de mis notas altas es mayor que la media de mis notas bajas" — siempre será verdad, con cualquier modelo.

**Gravedad:** El calibration plot con deciles está bien como **visualización**, pero las métricas Spearman y R² deben calcularse sobre los datos **individuales** (polígono a polígono).

**Fix:**
```r
# CORRECTO: correlación sobre datos individuales
spearman_raw <- cor.test(atlas_valid$pred_prob_mean,
                         as.numeric(atlas_valid$OTITAR),
                         method = "spearman")

# El calibration plot con deciles está bien SOLO como gráfico
# Pero reportar rho y R² de los datos individuales, no de los deciles
```

**Valores esperados con la corrección:** Spearman ρ ≈ 0.3-0.6 (que sería un resultado bueno y honesto para un modelo de ocupación con datos de citizen science).

---

### PROBLEMA 2 (ERROR) — El escalado de las predicciones es incorrecto

**Síntoma:** Potencialmente afecta a todas las métricas de validación.

**Causa (ya diagnosticada en Phase 1 del audit):** En `4_occupancy_otitar.R`, línea 497:

```r
pred_surface_std <- pred_surface %>%
  mutate_at(c("bio1", "bio2", "tree_cover", "grass_cover", "topo_elev"),
            ~(scale(.) %>% as.vector))
```

Esto re-escala la superficie de predicción usando **su propia media y sd**, no las del dataset de entrenamiento. Los coeficientes del modelo se estimaron con una escala diferente.

**Consecuencia para la validación:** El mapa de ocupación que se valida (`occ_otitar_prediction.csv`) tiene probabilidades sesgadas. La validación evalúa un mapa incorrecto.

**Gravedad:** Alta. Este bug ya fue corregido en la nueva versión del script 4 (Phase 2), pero las predicciones actuales que el alumno está validando siguen teniendo este error.

**Fix:** Regenerar las predicciones con el script 4 corregido (que usa los scaling_params guardados) y después re-ejecutar la validación.

---

### PROBLEMA 3 (ESPERADO/DOCUMENTAR) — Moran's I significativo en residuos

**Síntoma:** Autocorrelación espacial significativa en los residuos del modelo.

**Causa:** Esto es **esperable y no necesariamente un error**. Hay dos razones combinadas:

1. **Las aves esteparias tienen distribuciones geográficamente agrupadas.** Las mesetas castellanas concentran la mayoría de poblaciones. Un modelo de ocupación basado solo en variables ambientales no captura procesos espaciales como dispersión limitada, barreras geográficas, o historia biogeográfica. Esto deja autocorrelación residual.

2. **El atlas del MITECO y las predicciones del modelo miden cosas diferentes:**
   - **Predicciones del modelo:** Probabilidad de ocupación en función de hábitat (potencial)
   - **Atlas:** Presencia/ausencia observada en censos (realizada)

   Puede haber hábitat adecuado sin la especie presente (por ejemplo, si hay barreras de dispersión, extinción local histórica, o simplemente no se ha colonizado aún). Esto genera residuos positivos agrupados espacialmente (zonas de hábitat potencial vacías).

**¿Invalida los resultados?** No, pero hay que:
- **Reportarlo honestamente** en el paper (sección de limitaciones)
- **No usar métricas que asuman independencia** (el Spearman crudo ya viola este supuesto)
- Considerar si el modelo se beneficiaría de un **spatial random effect** (pero esto añade complejidad excesiva para un TFM)

**Redacción sugerida para el paper:**
> "Moran's I tests revealed significant positive spatial autocorrelation in model residuals (I = X.XX, p < 0.001), indicating that the model underpredicts occupancy in core range areas (e.g., central Castilla y León) and overpredicts in peripheral areas. This is consistent with the expectation that purely environmental models cannot capture dispersal limitations and historical biogeographic processes that shape species' realized distributions. The use of spatial block cross-validation (blockCV) for threshold optimization partially accounts for this spatial dependence."

---

### PROBLEMA 4 (METODOLÓGICO, MENOR) — El spatial block CV no se usa para Spearman/R²

**Síntoma:** El blockCV (líneas 47-96) se usa correctamente para calcular AUC, TSS, RMSE en cada fold. Pero el calibration plot (Spearman, R²) se calcula sobre **todos los datos** (`atlas_valid`), sin usar los folds.

**Causa:** Inconsistencia en el script: los threshold metrics respetan el CV espacial, pero las métricas continuas no.

**Fix:** O bien:
- (a) Calcular Spearman/R² por fold y reportar la media ± sd (consistente con el resto)
- (b) Dejar el calibration plot sobre todos los datos pero decir explícitamente que es un análisis exploratorio, no una métrica de validación formal

---

## RESUMEN DE ACCIONES

| # | Acción | Prioridad | Esfuerzo |
|---|--------|-----------|----------|
| 1 | **Calcular Spearman/R² sobre datos individuales**, no sobre deciles | ALTA | 5 min |
| 2 | **Regenerar predicciones** con script 4 corregido (scaling fix) | ALTA | 4-8 h (modelo) |
| 3 | **Documentar Moran's I** como limitación esperada en el paper | MEDIA | 30 min redacción |
| 4 | **Calcular Spearman por fold** o documentar que es exploratorio | BAJA | 15 min |

---

## RESPUESTA A LA PREGUNTA DEL ALUMNO

> "¿El hecho de comparar predicciones de hábitat potencial vs datos de censo del atlas puede causar la autocorrelación?"

**Sí, parcialmente.** Pero no es la única causa ni es un error per se:

1. **Diferencia conceptual:** El modelo predice **idoneidad de hábitat** (¿es buen hábitat?), mientras que el atlas registra **presencia real** (¿está la especie ahí?). Puede haber hábitat bueno sin la especie (por dispersión limitada, extinción local) y presencia en hábitat marginal (por inercia poblacional). Esto genera residuos espacialmente correlacionados.

2. **Temporal mismatch:** Las predicciones son para 2017-2022, el atlas puede tener datos de periodos diferentes.

3. **Resolución espacial:** Los polígonos del atlas (cuadrículas 10×10 km) promedian la heterogeneidad interna.

4. **El Moran's I significativo no invalida el modelo**, simplemente indica que quedan procesos espaciales no capturados. Esto es normal y esperable en modelos de distribución basados solo en ambiente. Para un paper de `colext()` centrado en dinámicas de colonización/extinción, es una limitación honesta, no un fallo fatal.
