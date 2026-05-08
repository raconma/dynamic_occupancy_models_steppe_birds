# stepRep — entrega para refitting de los modelos

Documento de integración para que el operador del modelo (Raúl) pueda
añadir la covariable de representatividad esteparia al submodelo de
detección con la **mínima fricción posible** dentro de la pipeline ya
existente (`scripts/1..4` + `R/model_configs.R`).

## TL;DR — opción A: solo correr (recomendada)

```
Rscript scripts/4_occupancy_models_v2.R
```

Eso ya hace todo:

1. Si `data/derived/stepRep_cellyear_*.csv` no existen, regenerarlas con
   `Rscript scripts/build_stepRep.R` antes (~5 min).
2. `4_occupancy_models_v2.R` detecta automáticamente si las columnas
   `stepRep_*_<year>` no están en los CSVs `_occ_wide_dynamic.csv` y
   en ese caso `source()` a `scripts/3b_add_stepRep.R`.
3. Después fija `cfg$p_formula <- update(cfg$p_formula, ~ . + stepRep_obs)`
   inline (no toca `R/model_configs.R`) y ajusta colext con `stepRep`
   añadido al submodelo de detección.
4. Todos los outputs llevan tag `_v2`
   (`results/{sp}_v2_model_object.rds`, `figs/{sp}_v2_*.png`, etc.) para
   no machacar los fits sin stepRep.

Si quieres cambiar de variante (`stepRep_strict_1km`, `stepRep_broad_500m`,
`stepRep_broad_1km`), edita la constante `STEPREP_VARIANT` al principio
del script y re-corre.

## Opción B: aplicar los patches a mano sobre `scripts/4_occupancy_models.R`

Si prefiere modificar el script v1 (porque ya lo tiene tuneado o
quiere mantener una sola ruta de entrada), los patches mínimos están
descritos abajo. La opción A es equivalente.

## 1. Por qué este patrón

`scripts/4_occupancy_models.R` ya construye el `unmarkedMultFrame` a
partir de `{sp}_occ_wide_dynamic.csv` extrayendo cada covariable
dinámica con dos formas paralelas: una matriz `<var>` (n_sites × T_years)
en `yearlySiteCovs` y una matriz expandida `<var>_obs` (n_sites × T·J_reps)
en `obsCovs`. Ej.: `NDVI` en yearly + `NDVI_obs` en obs. `stepRep` encaja
en ese mismo molde:

- `stepRep` (yearly) → captura per-year scaling para `train_dyn_scale`,
  igual que NDVI/pr/tmmn/tmmx.
- `stepRep_obs` (expandida) → entra en `pformula` como
  `~ ... + stepRep_obs`, como `NDVI_obs`.

La integración no toca ni el formato de los datos ni la estructura del
`unmarkedMultFrame`. Y el helper de response curves del script
(`plot_response_curves`) acepta cualquier covariable que esté en
`obsCovs`, así que `stepRep_obs` aparece automáticamente en la figura
de detección.

Verificado: la columna `cells` que usa la pipeline es el índice de la
malla WorldClim 5 km, idéntica al `cells` que produce `build_stepRep.R`.
Para *Pterocles alchata*, el solapamiento `pipeline ∩ stepRep` es
**4 131 / 4 131 celdas** (100 %).

Nota: tras unir, ~38 % de cell-years quedan inicialmente como NA porque
ese cell tuvo visitas en algunos años pero no en otros. Se imputan
primero por la media temporal del propio cell (siempre disponible, ya
que cada site tiene ≥ 1 año con stepRep real); el fallback de mediana
global queda en 0 NAs en las 4 especies. Esto es suficiente para
unmarked y preserva la variación temporal en los cells con datos en
varios años (que es donde el modelo aprende el efecto).

## 2. Patches a `scripts/4_occupancy_models.R` (opción B)

Cinco adiciones dentro del loop por especie. Líneas de referencia
del fichero v1 (commit `4ec50ca`).

### Patch 1 — extraer la matriz yearly (después de `tmmx`, ~línea 207)

```r
  tmmx <- occ_wide_clean[, paste0("tmmx_", YEARS)]
  stepRep <- occ_wide_clean[, paste0("stepRep_strict_500m_", YEARS)]   # añadir
```

### Patch 2 — estandarizar (después de `tmmx <- scale(tmmx)`, ~línea 227)

```r
  tmmx <- scale(tmmx)
  stepRep <- scale(stepRep)                                            # añadir
```

### Patch 3 — capturar el scaling per-year (~línea 234–243)

Añadir `stepRep = stepRep` al final de la lista `.dyn_mats`:

```r
  .dyn_mats <- list(
    EVI = EVI, NDVI = NDVI, pr = pr, tmmn = tmmn, tmmx = tmmx,
    Land_Cover_Type_1_Percent_Class_0  = Land_Cover_Type_1_Percent_Class_0,
    # ... resto igual ...
    Land_Cover_Type_1_Percent_Class_14 = Land_Cover_Type_1_Percent_Class_14,
    stepRep = stepRep                                                  # añadir
  )
```

Esto guarda center/scale por año en `{sp}_train_dyn_scale.rds`,
necesario para los scripts de atribución (8, 10) y predicción.

### Patch 4 — expandir a observation-level (~línea 270, junto a `NDVI_obs`)

```r
  NDVI_obs        <- expand_matrix(NDVI, J_reps)
  pr_obs          <- expand_matrix(pr, J_reps)
  topo_aspect_obs <- expand_matrix(siteCovs[, "topo_aspect", drop = FALSE], J_reps * T_years)
  topo_elev_obs   <- expand_matrix(siteCovs[, "topo_elev", drop = FALSE], J_reps * T_years)
  stepRep_obs     <- expand_matrix(stepRep, J_reps)                    # añadir
```

### Patch 5 — añadir a `yearlySiteCovs` y `obsCovs` (~línea 276–294)

```r
  occ_umf <- unmarkedMultFrame(
    y = y.cross,
    siteCovs = data.frame(siteCovs),
    yearlySiteCovs = list(
      years  = years_df, EVI = EVI, ...,
      NDVI = NDVI, pr = pr, tmmn = tmmn, tmmx = tmmx,
      stepRep = stepRep                                                # añadir
    ),
    obsCovs = list(
      duration = duration, effort = effort, observers = observers, time = time,
      NDVI_obs = NDVI_obs, pr_obs = pr_obs,
      topo_aspect_obs = topo_aspect_obs, topo_elev_obs = topo_elev_obs,
      stepRep_obs = stepRep_obs                                        # añadir
    ),
    numPrimary = T_years
  )
```

## 3. Patch a `R/model_configs.R`

Editar el `p_formula` de cada especie. Para *Pterocles alchata*
(línea 86):

```r
      p_formula       = ~ effort + observers + time,
```

cambiar a:

```r
      p_formula       = ~ effort + observers + time + stepRep_obs,
```

Replicar análogamente para las otras tres especies (líneas 27, 48, 67)
cuando se quiera extender. La opción A (`4_occupancy_models_v2.R`)
hace este paso inline con `update(cfg$p_formula, ~ . + stepRep_obs)`
para que `R/model_configs.R` no se toque.

## 4. Cómo correr la sensibilidad

`scripts/3b_add_stepRep.R` añade 4 variantes (28 columnas en total):
strict 500 m, strict 1 km, broad 500 m, broad 1 km.

**Opción A (v2 script):** editar la constante al principio del script
y re-ejecutar:

```r
STEPREP_VARIANT <- "stepRep_strict_1km"   # o stepRep_broad_500m / _broad_1km
```

**Opción B (patches manuales):** cambiar la línea del Patch 1 a la
variante deseada:

```r
  stepRep <- occ_wide_clean[, paste0("stepRep_strict_1km_",  YEARS)]
```

y re-ejecutar el step 4 con cada variante. El nombre `stepRep` en
la fórmula no cambia.

## 5. Lectura del coeficiente

- **Dirección esperada**: positiva. Mayor stepRep en una celda-año
  significa que más checklists caen sobre hábitat compatible →
  mayor *p* de detección por checklist.
- **Magnitud orientativa**: para que el efecto sobre γ y ε sea visible,
  el coeficiente sobre `stepRep` (logit) debería ser de al menos
  ~0.1–0.3 por SD. Por debajo de eso la covariable está cerca de ser
  irrelevante a este nivel de agregación.
- **Test de robustez**: comparar γ y ε con/sin `stepRep`:
  - Si γ baja y ε sube al añadir `stepRep` → la firma temporal
    descrita en `stepRep_diagnostics.md` §5 estaba inflando γ y
    deflactando ε. Es lo esperable dado que stepRep sube monotónicamente
    en celdas focales 2017–2023.
  - Si γ y ε apenas se mueven → `stepRep` no está absorbiendo sesgo
    a este nivel de agregación. Reportarla como auxiliar igualmente.

## 6. Banderas rojas

1. La tendencia temporal en celdas focales **sube** (≈ +0.10–0.14 sobre
   la escala 0–1 a lo largo de 2017–2023, ver §5 del informe
   diagnóstico). Esto va al revés de la hipótesis "newcomers diluyen el
   muestreo". El sesgo es real y temporal, pero su signo es opuesto al
   asumido inicialmente: los checklists tardíos en celdas focales son
   **más informativos por visita** que los tempranos.
2. La correlación de Spearman entre `log(n_checklists)` y `stepRep` por
   celda-año es **prácticamente cero** (-0.005 a -0.05). Si los modelos
   ya controlan por `effort`/`observers`/`duration`, eso **no** captura
   el sesgo: el sesgo opera vía estructura espacial intra-celda, no vía
   esfuerzo agregado.
3. ~38 % de los cell-years en las tablas del modelo se imputan
   (row-mean). Para los cells con datos en varios años la imputación
   conserva la variación; para los cells con un sólo año observado
   queda un valor constante a lo largo de los 7 años. Si Raúl quiere
   ser conservador puede filtrar `n_checklists ≥ 3` por celda-año antes
   de imputar (cambio menor en `scripts/3b_add_stepRep.R`); reportar
   ambas versiones como sensibilidad.

## 7. Diagnósticos previos al fit

`reports/stepRep_diagnostics.md` (mismo commit) contiene:

- §1: provenencia y máscara CORINE.
- §3: frecuencia de clases CLC bajo checklists (sesgo birder
  cuantificado).
- §5: tendencia temporal por celdas focales vs peninsulares (figura).
- §6: mapa de stepRep medio por celda Iberia 5 km.
- §7: correlación con esfuerzo (figura).
- §10: borrador de Methods y Results en prosa narrativa, listo para
  insertar en el manuscrito GCB con marcas `[REF NEEDED]`.

## 8. Reproducibilidad

Pipeline completo (orden):

```
scripts/build_stepRep.R           # CORINE → focal rasters → cell-year table  (~5 min)
scripts/4_occupancy_models_v2.R   # auto-llama a 3b si hace falta y ajusta colext + stepRep
                                  # outputs en results/{sp}_v2_*.rds y figs/{sp}_v2_*.png
```

Si quieres ver el join intermedio:
```
scripts/3b_add_stepRep.R          # cell-year → wide model table              (~10 s)
```

Para stPGOcc el script v2 análogo:
```
scripts/18_stPGOcc_production_run_v2.R   # auto-3b + stepRep_obs en det.covs/det.formula
```

Outputs separados de v1: `results/results_spatial/results_production_v2/`
y `figs/{sp}_spatial_diagnostics_production_v2.png`. Compara
automáticamente con el fit v1 (`results_production/{sp}_stPGOcc.rds`)
si existe y reporta el delta de WAIC.

Commits que cierran el entregable:
- `4ec50ca` — Add steppe-representativeness covariate (CLC2018) for detection sub-model
- `0cf9f43` — Add stepRep handover doc with concrete colext / stPGOcc integration
- `a21c9ef` — Add 3b_add_stepRep.R integrator + minimal-friction recipe
- `a2ff874` — Add 4_occupancy_models_v2.R drop-in script with stepRep wired
- (este) — Add 18_stPGOcc_production_run_v2.R drop-in script with stepRep wired
