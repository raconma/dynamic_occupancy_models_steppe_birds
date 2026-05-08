# stepRep — entrega para refitting de los modelos

Documento de integración para que el operador del modelo (Raúl) pueda
añadir la covariable de representatividad esteparia al submodelo de
detección con la **mínima fricción posible** dentro de la pipeline ya
existente (`scripts/1..4` + `R/model_configs.R`).

## TL;DR — los 3 pasos

1. **Asegurarse de tener las tablas de la covariable.** Si no las tiene,
   regenerar (~5 min) con `Rscript scripts/build_stepRep.R`. Si Guille
   se las pasa por correo / Zenodo, copiarlas a
   `data/derived/stepRep_cellyear_{otitar,ptealc,pteori,tettet}.csv`.

2. **Correr el script integrador (una vez):**
   ```
   Rscript scripts/3b_add_stepRep.R
   ```
   Este paso añade 28 columnas yearly-site-cov a cada
   `data/processed_2023/{sp}/{sp}_occ_wide_dynamic.csv` (4 variantes ×
   7 años: strict/broad × 500 m/1 km × 2017–2023). Es idempotente y
   no destructivo: si los `stepRep_*_<year>` ya están, no hace nada;
   y en cualquier caso siempre se puede regenerar el dynamic CSV
   re-corriendo `scripts/3_prepare_dynamic_variables.R`.

3. **Aplicar los 3 patches mínimos a `scripts/4_occupancy_models.R` y
   `R/model_configs.R`** (detallados abajo) y re-ejecutar el step 4 como
   siempre. Los patches añaden 4 líneas y modifican 1.

A partir de ahí, todo el pipeline downstream
(`5_validation.R`, `15_parboot_publication.R`,
`18_stPGOcc_production_run.R`, etc.) sigue funcionando sin cambios.

## 1. Por qué este patrón

`scripts/4_occupancy_models.R` ya construye el `unmarkedMultFrame` a
partir de `{sp}_occ_wide_dynamic.csv` extrayendo cada covariable yearly
con `occ_wide_clean[, paste0("<varname>_", YEARS)]` (ver bloque líneas
~215–225, donde se cargan `EVI`, `NDVI`, `pr`, `tmmn`, `tmmx` y los
`Land_Cover_Type_1_Percent_Class_*`). `stepRep` encaja exactamente en
ese mismo molde porque es también una **yearly site covariate**: un
valor por (cell, año). La integración minimiza fricción porque:

- No cambia el sitio donde se cargan los datos (`occ_wide_dynamic.csv`).
- No cambia el formato de las columnas (sigue siendo `<varname>_<year>`).
- No cambia la estructura del `unmarkedMultFrame`.
- No requiere tocar nada en steps 1–3 ni en R/.

Verificado: la columna `cells` que usa la pipeline es el índice de la
malla WorldClim 5 km, idéntica al `cells` que produce `build_stepRep.R`.
Para *Pterocles alchata*, el solapamiento `pipeline ∩ stepRep` es
**4 131 / 4 131 celdas** (100 %).

Nota: tras unir, ~38 % de cell-years quedan inicialmente como NA porque
ese cell tuvo visitas en algunos años pero no en otros. Se imputan
primero por la media temporal del propio cell (siempre disponible, ya
que cada site tiene ≥ 1 año con stepRep real); el fallback de mediana
global queda en 0 NAs en los 4 species. Esto es suficiente para
unmarked y preserva la variación temporal en los cells con datos en
varios años (que es donde el modelo aprende el efecto).

## 2. Patches a `scripts/4_occupancy_models.R`

Tres adiciones, todas dentro del loop por especie. Las líneas de
referencia son del fichero actual (commit `4ec50ca`).

### Patch 1 — extraer la matriz (después de `tmmx`, ~línea 223)

Buscar el final del bloque que extrae las yearly site covs:

```r
  tmmx <- occ_wide_clean[, paste0("tmmx_", YEARS)]
```

Añadir inmediatamente debajo:

```r
  stepRep <- occ_wide_clean[, paste0("stepRep_strict_500m_", YEARS)]
```

### Patch 2 — estandarizar (después de `tmmx <- scale(tmmx)`, ~línea 235)

Buscar el final del bloque de estandarización:

```r
  tmmx <- scale(tmmx)
```

Añadir inmediatamente debajo:

```r
  stepRep <- scale(stepRep)
```

### Patch 3 — capturar parámetros de escalado (~línea 240–253)

Dentro del bloque que construye `train_dyn_scale` con la lista
`.dyn_mats`, añadir `stepRep = stepRep` a la lista:

```r
  .dyn_mats <- list(
    EVI = EVI, NDVI = NDVI, pr = pr, tmmn = tmmn, tmmx = tmmx,
    Land_Cover_Type_1_Percent_Class_0  = Land_Cover_Type_1_Percent_Class_0,
    # ... resto igual ...
    Land_Cover_Type_1_Percent_Class_14 = Land_Cover_Type_1_Percent_Class_14,
    stepRep = stepRep                       # <-- añadir
  )
```

Esto guarda center/scale por año en `{sp}_train_dyn_scale.rds`, igual
que el resto. Necesario para que los scripts de atribución (8, 10) y
de predicción puedan reescalar consistentemente.

### Patch 4 — añadir al `yearlySiteCovs` (~línea 276–287)

```r
  occ_umf <- unmarkedMultFrame(
    y = y.cross,
    siteCovs = data.frame(siteCovs),
    yearlySiteCovs = list(
      years  = years_df,
      EVI    = EVI,
      Land_Cover_Type_1_Percent_Class_0  = Land_Cover_Type_1_Percent_Class_0,
      # ... resto igual ...
      NDVI = NDVI, pr = pr, tmmn = tmmn, tmmx = tmmx,
      stepRep = stepRep                     # <-- añadir
    ),
    obsCovs = list( ... ),
    numPrimary = T_years
  )
```

## 3. Patch a `R/model_configs.R`

Editar el `p_formula` de cada especie. Para el piloto en *Pterocles
alchata* (línea 86 actual):

```r
      p_formula       = ~ effort + observers + time,
```

cambiar a:

```r
      p_formula       = ~ effort + observers + time + stepRep,
```

Replicar análogamente para las otras tres especies (líneas 27, 48, 67)
una vez confirmado el efecto en ptealc.

## 4. Cómo correr la sensibilidad

`scripts/3b_add_stepRep.R` añade 4 variantes (28 columnas en total):
strict 500 m, strict 1 km, broad 500 m, broad 1 km. Para reportar la
sensibilidad sólo hay que cambiar la línea del Patch 1:

```r
  stepRep <- occ_wide_clean[, paste0("stepRep_strict_1km_",  YEARS)]   # buffer 1 km
  stepRep <- occ_wide_clean[, paste0("stepRep_broad_500m_",  YEARS)]   # incluye dehesa
  stepRep <- occ_wide_clean[, paste0("stepRep_broad_1km_",   YEARS)]   # ambas combinadas
```

y re-ejecutar el step 4 con cada variante. El resto de patches son
idénticos. El nombre `stepRep` en la fórmula no cambia.

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

```
scripts/build_stepRep.R   # CORINE → focal rasters → cell-year table  (~5 min)
scripts/3b_add_stepRep.R  # cell-year → wide model table              (~10 s)
scripts/4_occupancy_models.R    # fit colext (con los patches arriba)
scripts/18_stPGOcc_production_run.R   # fit stPGOcc (análogo, det.covs)
```

Commits que cierran el entregable:
- `4ec50ca` — Add steppe-representativeness covariate (CLC2018) for detection sub-model
- `0cf9f43` — Add stepRep handover doc with concrete colext / stPGOcc integration
- (este) — Add 3b_add_stepRep.R integrator + minimal-friction recipe
