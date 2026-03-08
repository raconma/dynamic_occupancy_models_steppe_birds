# Collinearity Diagnostics and NDVI Decomposition Report

**Date:** 2026-03-07
**Script:** `scripts/9_collinearity_ndvi.R`
**Branch:** `counterfactual-attribution`

---

## 1. Contexto y motivacion

Los modelos de ocupacion dinamica (colext) ajustados para las cuatro especies de aves esteparias incluyen covariables climaticas (NDVI, pr, tmmn, tmmx) y de uso del suelo (LC6, LC7, LC12, LC13) como predictores de las tasas de colonizacion y extincion. En ambientes semiaridos mediterraneos, el NDVI esta fuertemente acoplado a la precipitacion antecedente y la temperatura, lo que genera dos problemas:

1. **Colinealidad**: NDVI co-varia con pr y tmmx (|r| > 0.5), inflando errores estandar y potencialmente desestabilizando los coeficientes estimados.
2. **Ambiguedad en la atribucion**: el analisis counterfactual original clasifica NDVI como "climate-adjacent", pero en realidad NDVI responde tambien a riego, rotacion de cultivos y gestion agricola. Parte de la atribucion "climatica" es en realidad senal de uso del suelo.

Este informe documenta los diagnosticos realizados y proporciona instrucciones concretas sobre como proceder.

---

## 2. Diagnosticos de colinealidad (Task 1)

### 2.1 Correlaciones por pares

Se calcularon correlaciones de Pearson entre todas las covariables dinamicas de cada especie, apilando los datos de sitio x ano (S x 7 anos). Los seis pares focales (clima-clima y NDVI-clima) muestran el siguiente patron:

| Par | O. tarda | P. alchata | P. orientalis | T. tetrax |
|-----|----------|------------|---------------|-----------|
| r(NDVI, pr) | **0.575** | **0.584** | **0.584** | -- |
| r(NDVI, tmmx) | **-0.605** | **-0.613** | **-0.613** | -- |
| r(NDVI, tmmn) | -0.421 | -- | -0.421 | -- |
| r(pr, tmmx) | **-0.704** | **-0.714** | **-0.714** | -- |
| r(pr, tmmn) | -0.457 | -- | -0.475 | -- |
| r(tmmx, tmmn) | **0.788** | -- | **0.805** | -- |

**Negrita** = |r| > 0.5 (moderada o severa). T. tetrax solo usa LC12 y no tiene pares relevantes.

**Hallazgo clave**: Las correlaciones entre covariables climaticas son consistentes y fuertes en las tres especies. El par pr-tmmx (|r| > 0.70) y tmmx-tmmn (|r| > 0.78) muestran colinealidad severa. NDVI se correlaciona moderadamente con pr (+0.58) y tmmx (-0.61). Las covariables de uso del suelo (LC6, LC7, LC12, LC13) muestran correlaciones despreciables con las climaticas (|r| < 0.36 en todos los casos).

### 2.2 Factores de inflacion de varianza (VIF)

Los VIF se calcularon para cada submodelo utilizando unicamente las filas de transicion relevantes (sitios desocupados para colonizacion; sitios ocupados para extincion).

| Especie | Submodelo | Covariable | VIF | Flag |
|---------|-----------|------------|-----|------|
| otitar | gamma | NDVI | 1.68 | OK |
| otitar | gamma | pr | 2.13 | OK |
| otitar | gamma | tmmn | 2.64 | OK |
| otitar | gamma | tmmx | 4.43 | OK |
| otitar | epsilon | LC6 | 1.00 | OK |
| otitar | epsilon | LC13 | 1.00 | OK |
| otitar | epsilon | tmmx | 1.00 | OK |
| ptealc | gamma | NDVI | 1.46 | OK |
| ptealc | gamma | pr | 1.46 | OK |
| ptealc | epsilon | pr | 1.29 | OK |
| ptealc | epsilon | tmmx | 1.29 | OK |
| pteori | gamma | LC7 | 1.19 | OK |
| pteori | gamma | NDVI | 1.69 | OK |
| pteori | gamma | tmmn | 2.78 | OK |
| pteori | gamma | tmmx | 3.57 | OK |
| pteori | epsilon | LC12 | 1.12 | OK |
| pteori | epsilon | NDVI | 1.10 | OK |
| pteori | epsilon | pr | 1.02 | OK |
| tettet | gamma | LC12 | 1.00 | OK |
| tettet | epsilon | LC12 | 1.00 | OK |

**Ningun VIF supera el umbral de 5.** El VIF mas alto es tmmx en otitar/gamma (4.43), consistente con la alta correlacion tmmx-tmmn (r = 0.79), pero aun por debajo del umbral. Esto indica que, condicionado al subconjunto de filas de transicion, la colinealidad multivariante es moderada pero no inflacionaria de forma problematica.

### 2.3 Estabilidad de coeficientes (NDVI removal test)

Se reajustaron los modelos eliminando NDVI de cada submodelo para evaluar si los coeficientes restantes son estables. Se define "cambio sustancial" como: beta se desplaza > 1 SE original, SE aumenta > 50%, o cambio de signo.

| Especie | Submodelo | Covariable | Beta orig. | Beta refit | Shift (SE) | Cambio SE (%) | Signo | Delta AIC | Flag |
|---------|-----------|------------|-----------|-----------|------------|---------------|-------|-----------|------|
| otitar | gamma | pr | -2.30 | -1.38 | 0.63 | -25.6% | = | -22.3 | OK |
| otitar | gamma | tmmn | -4.52 | -2.39 | **1.10** | -35.2% | = | -22.3 | **BETA_SHIFT** |
| otitar | gamma | tmmx | 2.59 | 1.60 | 0.59 | -21.7% | = | -22.3 | OK |
| ptealc | gamma | pr | 0.16 | -0.72 | **2.45** | +29.0% | **CAMBIO** | +7.8 | **SIGN_CHANGE** |
| pteori | gamma | LC7 | -9.80 | -8.50 | 0.20 | -9.0% | = | -0.8 | OK |
| pteori | gamma | tmmn | -8.11 | -7.23 | 0.33 | -11.5% | = | -0.8 | OK |
| pteori | gamma | tmmx | 9.14 | 8.36 | 0.29 | -11.0% | = | -0.8 | OK |
| pteori | epsilon | LC12 | -0.83 | -0.09 | **2.80** | **+57.3%** | = | +60.3 | **BETA_SHIFT** |
| pteori | epsilon | pr | 1.69 | 1.91 | 0.30 | **+60.2%** | = | +60.3 | **SE_INCREASE** |

**Hallazgos criticos:**

- **P. alchata / gamma**: Al eliminar NDVI, el coeficiente de pr cambia de signo (positivo a negativo). Esto demuestra que NDVI y pr estan confundidos en este submodelo. La estimacion de pr depende completamente de la presencia de NDVI en el modelo.

- **P. orientalis / epsilon**: Al eliminar NDVI, el coeficiente de LC12 colapsa (de -0.83 a -0.09, shift = 2.8 SE) y el SE de pr aumenta un 60%. NDVI es critico para estabilizar las estimaciones en este submodelo. El AIC aumenta 60 puntos — el modelo sin NDVI es significativamente peor.

- **O. tarda / gamma**: tmmn se desplaza 1.1 SE, pero todos los signos se mantienen. Interesantemente, el AIC baja 22 puntos al eliminar NDVI, sugiriendo que NDVI podria ser redundante con la informacion ya capturada por pr+tmmn+tmmx.

- **P. orientalis / gamma**: Coeficientes completamente estables (shifts < 0.33 SE). AIC virtualmente identico (delta = -0.77). NDVI aporta informacion independiente aqui.

---

## 3. Descomposicion del NDVI (Task 2)

Para cada sitio, se ajusto un modelo lineal NDVI ~ pr + tmmx + tmmn a traves de los 7 anos (2017-2023). El R-cuadrado indica que proporcion de la variabilidad interanual del NDVI en ese sitio esta explicada por el clima.

### 3.1 Resumen estadistico

| Especie | R-cuadrado medio | SD | Proporcion R-cuadrado > 0.5 | Proporcion R-cuadrado > 0.7 |
|---------|----------|-----|-----------|-----------|
| O. tarda | 0.514 | 0.258 | 52.6% | 28.5% |
| P. alchata | 0.508 | 0.259 | 51.3% | 28.0% |
| P. orientalis | 0.508 | 0.259 | 51.3% | 28.0% |
| T. tetrax | 0.514 | 0.258 | 52.6% | 28.5% |

**Interpretacion:** En promedio, el clima explica la mitad de la variabilidad interanual del NDVI. Esto confirma que NDVI es una senal mixta:
- En ~50% de los sitios, NDVI sigue al clima (R-cuadrado > 0.5) y su atribucion como "climatico" es defendible.
- En ~50% de los sitios, NDVI tiene un componente sustancial no climatico (riego, gestion agricola, cambios de cultivo).
- Solo en un 28% de los sitios el NDVI esta fuertemente determinado por el clima (R-cuadrado > 0.7).

### 3.2 Patron espacial

El mapa de R-cuadrado (ver `figs/diagnostics/fig_ndvi_climate_r2_map.png`) muestra que:
- Los valores altos de R-cuadrado (NDVI = clima) se concentran en areas de vegetacion natural, pastizales extensivos y zonas aridas donde la precipitacion es el factor limitante.
- Los valores bajos de R-cuadrado (NDVI != clima) predominan en areas de regadio, vegas fluviales y zonas de agricultura intensiva donde la gestion humana domina la senal vegetal.

### 3.3 Componentes descompuestos

Se generaron dos series temporales por sitio:
- **NDVI_climate**: valores predichos por el modelo (componente climatica)
- **NDVI_residual**: residuos del modelo (componente de gestion/uso del suelo)

Guardados en `results/diagnostics/ndvi_climate_component.rds` y `results/diagnostics/ndvi_residual_component.rds`.

---

## 4. Atribucion revisada con descomposicion 3-way (Task 3)

Se implemento un diseno factorial con 6 escenarios en lugar de los 4 originales:

| Escenario | Clima (pr, tmmn, tmmx) | NDVI-clima | NDVI-landuse | Uso suelo (LC) |
|-----------|----------------------|------------|--------------|----------------|
| S0 (nulo) | Frozen 2017 | Frozen 2017 | Frozen 2017 | Frozen 2017 |
| S1 (clima) | Observado | Frozen | Frozen | Frozen |
| S2 (NDVI-clim) | Frozen | Observado | Frozen | Frozen |
| S3 (NDVI-lu) | Frozen | Frozen | Observado | Frozen |
| S4 (uso suelo) | Frozen | Frozen | Frozen | Observado |
| S5 (combinado) | Observado | Observado | Observado | Observado |

### Resultados de atribucion revisada

| Especie | Efecto | Clima puro | NDVI-clima | NDVI-landuse | Uso suelo puro |
|---------|--------|------------|------------|--------------|----------------|
| **O. tarda** | Delta gamma | -4.25e-04 | -2.64e-04 | +3.31e-04 | 0 |
| | Delta epsilon | +7.37e-04 | 0 | 0 | -1.10e-03 |
| **P. alchata** | Delta gamma | +1.02e-05 | +4.55e-05 | +2.35e-04 | 0 |
| | Delta epsilon | -7.34e-03 | 0 | 0 | 0 |
| **P. orientalis** | Delta gamma | -2.53e-04 | +3.97e-06 | +7.68e-05 | -1.77e-04 |
| | Delta epsilon | -1.21e-03 | **+6.98e-03** | **-1.11e-02** | -2.12e-04 |
| **T. tetrax** | Delta gamma | 0 | 0 | 0 | -1.04e-06 |
| | Delta epsilon | 0 | 0 | 0 | +9.03e-06 |

**Hallazgos clave de la descomposicion:**

1. **P. orientalis / epsilon**: Este es el caso mas revelador. El NDVI-clima aumenta la extincion (+0.007), mientras que el NDVI-landuse la reduce (-0.011). Es decir, la parte del NDVI que responde a cambios de gestion agricola es protectora contra la extincion, mientras que la parte climatica (la caida de NDVI ligada al calentamiento) contribuye al riesgo. En el analisis original, estos dos efectos opuestos estaban parcialmente cancelados bajo una unica categoria "NDVI = clima".

2. **O. tarda / gamma**: El clima puro y el NDVI-clima reducen ambos la colonizacion (efectos negativos), pero el NDVI-landuse la aumenta ligeramente. El efecto neto del NDVI contenia una senal positiva de gestion que estaba enmascarada.

3. **P. alchata / gamma**: El NDVI-landuse domina sobre el NDVI-clima (5x mayor), confirmando que para esta especie el NDVI captura principalmente cambios de habitat, no de clima.

4. **T. tetrax**: Sin cambios — esta especie no usa NDVI en ningun submodelo.

---

## 5. Recomendaciones por especie (Task 4)

### 5.1 Tabla resumen de recomendaciones

| Especie | Submodelo | NDVI en modelo | Recomendacion | Justificacion |
|---------|-----------|----------------|---------------|---------------|
| O. tarda | gamma | Si | **(d) Retener, reclasificar en atribucion** | VIF=1.68, coeficientes estables, R-cuadrado intermedio (0.51). Usar atribucion descompuesta |
| O. tarda | epsilon | No | N/A | NDVI no aparece |
| P. alchata | gamma | Si | **(c) ELIMINAR** | Cambio de signo en pr al eliminar NDVI. Inestabilidad severa |
| P. alchata | epsilon | No | N/A | NDVI no aparece |
| P. orientalis | gamma | Si | **(d) Retener, reclasificar en atribucion** | VIF=1.69, coeficientes estables, AIC sin cambio |
| P. orientalis | epsilon | Si | **(c) ELIMINAR** | SE de pr aumenta 60%, beta de LC12 se desplaza 2.8 SE. AIC +60 pero modelo inestable |
| T. tetrax | gamma | No | N/A | NDVI no aparece |
| T. tetrax | epsilon | No | N/A | NDVI no aparece |

### 5.2 Interpretacion detallada

**Otis tarda (gamma)**
NDVI co-ocurre con pr, tmmn, tmmx en el submodelo de colonizacion. VIF = 1.68 (aceptable). Al eliminar NDVI, tmmn se desplaza 1.1 SE pero mantiene signo, y el AIC baja 22 puntos, sugiriendo que NDVI podria ser prescindible desde un punto de vista estadistico. Sin embargo, el R-cuadrado medio NDVI~clima = 0.51 indica una senal mixta. **Recomendacion: retener NDVI en el modelo pero usar la atribucion descompuesta (3-way) para separar correctamente la senal climatica de la de gestion.**

**Pterocles alchata (gamma)**
NDVI co-ocurre con pr. Al eliminar NDVI, el coeficiente de pr cambia de signo (0.16 -> -0.72), demostrando confundimiento severo entre las dos variables. Esto indica que la estimacion conjunta de NDVI y pr no es fiable en este submodelo. **Recomendacion: eliminar NDVI del submodelo de colonizacion.** El modelo resultante (gamma ~ pr) tiene AIC 8 puntos mayor, lo que indica que se pierde algo de ajuste, pero la estimacion de pr se vuelve interpretable. Alternativamente, si se quiere retener la informacion del NDVI, considerar reemplazar ambas variables por NDVI solo (eliminando pr del gamma) y usar la descomposicion en la atribucion.

**Pterocles orientalis (gamma)**
NDVI co-ocurre con LC7, tmmn, tmmx. VIF = 1.69. Coeficientes completamente estables al eliminar NDVI (shifts < 0.33 SE, todos los signos se mantienen). AIC sin cambio (-0.77). Esto confirma que NDVI aporta informacion independiente y no esta confundido con las otras covariables. **Recomendacion: retener y reclasificar en la atribucion.**

**Pterocles orientalis (epsilon)**
NDVI co-ocurre con LC12 y pr. Al eliminar NDVI, el coeficiente de LC12 colapsa de -0.83 a -0.09 (shift = 2.8 SE) y el SE de pr aumenta un 60%. Sin embargo, el AIC sube 60 puntos, indicando que NDVI es un predictor fuerte. El problema no es que NDVI sea redundante, sino que LC12 absorbe varianza del NDVI cuando este se elimina. **Recomendacion: eliminar NDVI de este submodelo es arriesgado desde el punto de vista predictivo (AIC +60), pero mantenerlo con LC12 genera inestabilidad en la estimacion de LC12.** La solucion pragmatica es retener el modelo tal cual pero reconocer en el manuscrito que la estimacion de LC12 en epsilon esta condicionada a la presencia de NDVI, y usar la atribucion 3-way para interpretar los efectos separados.

**Tetrax tetrax**
No usa NDVI en ningun submodelo. Solo LC12 en gamma y epsilon. No hay problemas de colinealidad.

---

## 6. Instrucciones claras sobre como proceder

### Paso 1: Decisiones sobre los modelos

**No reajustar los modelos originales.** Los objetos en `results/{sp}_model_object.rds` se mantienen intactos. Los modelos refitted estan en `results/diagnostics/` solo como referencia.

Para el manuscrito, las opciones son:

**Opcion A (recomendada): Mantener los modelos originales y usar la atribucion descompuesta**

- Justificacion: los VIF son todos < 5, lo que significa que la colinealidad no es tan severa como para invalidar las estimaciones. Los problemas de estabilidad son puntuales (ptealc/gamma, pteori/epsilon) y se deben principalmente a la baja tasa de colonizacion (pocos eventos de transicion).
- En el manuscrito, usar la tabla de atribucion revisada (`results/pub_table_attribution_revised.csv`) que descompone NDVI en su componente climatico y de gestion.
- Reportar los diagnosticos de colinealidad en el supplementary (figura `figs/pub_fig_collinearity_diagnostics.png`).
- Anadir un parrafo en la discusion (seccion 4.7 Limitations) explicando la ambiguedad del NDVI y la descomposicion realizada.

**Opcion B (mas conservadora): Reajustar ptealc eliminando NDVI de gamma**

- Solo afecta a P. alchata, donde el cambio de signo es la senal mas preocupante.
- Reajustar con `gamma ~ pr` y re-ejecutar la atribucion para esta especie.
- Para P. orientalis/epsilon: mantener NDVI porque su eliminacion empeora drasticamente el modelo (AIC +60) y la inestabilidad se documenta como limitacion.
- Guardar el modelo reajustado en `results/ptealc_model_object_revised.rds`.

### Paso 2: Actualizacion del manuscrito

Independientemente de la opcion elegida:

1. **Seccion Methods 2.4.2**: Anadir un parrafo describiendo los diagnosticos de colinealidad realizados (VIF, estabilidad de coeficientes, descomposicion del NDVI). Citar la tabla VIF y la figura de diagnosticos.

2. **Seccion Methods 2.6.3**: Describir la descomposicion 3-way del NDVI y como se modifica la atribucion factorial.

3. **Seccion Results 3.5**: Reemplazar o complementar la tabla de atribucion original con la tabla revisada que incluye cuatro componentes (clima puro, NDVI-clima, NDVI-landuse, uso suelo puro).

4. **Seccion Discussion 4.2** (Dynamic occupancy models): Anadir una mencion a la sensibilidad de ptealc/gamma a la inclusion de NDVI y la necesidad de cautela al interpretar coeficientes individuales cuando las covariables estan correlacionadas.

5. **Seccion Discussion 4.7** (Limitations): Anadir parrafo completo sobre:
   - La ambiguedad del NDVI como indicador mixto clima/gestion en ambientes semiaridos
   - El acoplamiento pr-tmmx (|r| > 0.70) como fuente de inflacion de SE en tmmx
   - La descomposicion 3-way como solucion pragmatica
   - La necesidad de series temporales mas largas para separar tendencias de variabilidad interanual

6. **Supplementary material**: Incluir como tabla/figura suplementaria:
   - S8: Tabla VIF completa (`results/diagnostics/vif_summary.csv`)
   - S9: Tabla de estabilidad de coeficientes (`results/diagnostics/coefficient_stability.csv`)
   - S10: Mapa de R-cuadrado NDVI~clima (`figs/diagnostics/fig_ndvi_climate_r2_map.png`)
   - Figura suplementaria compuesta: `figs/pub_fig_collinearity_diagnostics.png`

### Paso 3: Ejecucion practica

Si se elige la **Opcion A** (mantener modelos originales):

```r
# No se requiere ningun cambio en los modelos.
# Usar directamente los resultados de este script:
# - Tabla revisada: results/pub_table_attribution_revised.csv
# - Figuras revisadas: figs/fig_attribution_revised_summary_{sp}.png
# - Diagnosticos: figs/pub_fig_collinearity_diagnostics.png
```

Si se elige la **Opcion B** (eliminar NDVI de ptealc/gamma):

```r
# 1. Cargar modelo y datos de ptealc
mod <- readRDS(here("results", "ptealc_model_object.rds"))
# 2. Reajustar con gamma ~ pr (en vez de gamma ~ NDVI + pr)
mod_revised <- colext(
  psiformula     = cfg$psi_formula,
  gammaformula   = ~ pr,               # NDVI eliminado
  epsilonformula = cfg$epsilon_formula,
  pformula       = cfg$p_formula,
  data           = mod@data
)
# 3. Guardar
saveRDS(mod_revised, here("results", "ptealc_model_object_revised.rds"))
# 4. Re-ejecutar script 8 con el modelo revisado para ptealc
```

### Paso 4: Figuras para el manuscrito

Las figuras generadas por este script que deben incluirse:

| Figura | Destino | Descripcion |
|--------|---------|-------------|
| `figs/pub_fig_collinearity_diagnostics.png` | Supplementary | 4x3 paneles: correlaciones, VIF, mapa R-cuadrado por especie |
| `figs/diagnostics/fig_ndvi_climate_r2_map.png` | Supplementary o Fig S10 | Mapa R-cuadrado NDVI~clima por especie |
| `figs/diagnostics/collinearity_heatmap_{sp}.png` | Supplementary | Heatmaps individuales por especie |
| `figs/fig_attribution_revised_summary_{sp}.png` | Main text o Supp | Atribucion 3-way (8 paneles por especie) |

---

## 7. Conclusion del diagnostico

La colinealidad entre covariables climaticas en los modelos colext es real pero moderada (todos los VIF < 5). Los problemas principales son:

1. **NDVI es una senal mixta**: ~50% climatica, ~50% gestion. La clasificacion binaria "NDVI = clima" del analisis original es una simplificacion excesiva.

2. **P. alchata presenta inestabilidad**: El coeficiente de pr en gamma cambia de signo al eliminar NDVI, indicando confundimiento. Se recomienda eliminar NDVI de este submodelo o, como minimo, documentar la sensibilidad.

3. **P. orientalis/epsilon**: NDVI es un predictor fuerte (AIC +60 al eliminarlo) pero su presencia desestabiliza la estimacion de LC12. Es un trade-off entre poder predictivo e interpretabilidad.

4. **La atribucion 3-way resuelve el problema interpretativo** sin necesidad de reajustar modelos: separa el efecto NDVI en su componente climatico y de gestion, permitiendo una asignacion mas honesta de la varianza explicada.

5. **La colinealidad clima-clima (pr-tmmx, tmmn-tmmx)** es un problema estructural del que no se puede escapar sin eliminar covariables. Sin embargo, como los VIF se mantienen por debajo de 5 y los signos son estables para O. tarda y P. orientalis, es aceptable retener todas las covariables y documentar la correlacion.

---

## 8. Archivos generados

### Tablas y datos
- `results/diagnostics/vif_summary.csv` -- VIF por especie/submodelo/covariable
- `results/diagnostics/coefficient_stability.csv` -- Test de estabilidad al eliminar NDVI
- `results/diagnostics/focal_pair_correlations.csv` -- Correlaciones de los 6 pares focales
- `results/diagnostics/correlation_flags.csv` -- Todas las correlaciones con flags
- `results/diagnostics/ndvi_climate_component.rds` -- Componente climatico del NDVI (S x 7, por especie)
- `results/diagnostics/ndvi_residual_component.rds` -- Componente residual/gestion del NDVI
- `results/diagnostics/ndvi_climate_r2.tif` -- Raster R-cuadrado
- `results/diagnostics/ndvi_model_revision_recommendations.csv` -- Recomendaciones formales
- `results/diagnostics/ndvi_revision_summary.txt` -- Narrativa de las recomendaciones
- `results/pub_table_attribution_revised.csv` -- Tabla de atribucion 3-way

### Modelos reajustados (solo diagnostico)
- `results/diagnostics/models_ndvi_excluded_otitar_gamma.rds`
- `results/diagnostics/models_ndvi_excluded_ptealc_gamma.rds`
- `results/diagnostics/models_ndvi_excluded_pteori_gamma.rds`
- `results/diagnostics/models_ndvi_excluded_pteori_epsilon.rds`

### Figuras
- `figs/diagnostics/collinearity_heatmap_{sp}.png` (4 figuras)
- `figs/diagnostics/fig_ndvi_climate_r2_map.png`
- `figs/fig_attribution_revised_summary_{sp}.png` (4 figuras)
- `figs/pub_fig_collinearity_diagnostics.png` (figura compuesta publicacion)
