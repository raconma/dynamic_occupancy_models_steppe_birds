# stepRep — paquete de entrega para refitting de los modelos

Documento de integración pensado para que el operador del modelo (Raúl) pueda
añadir la covariable de representatividad esteparia al submodelo de detección
sin necesidad de re-ejecutar `scripts/build_stepRep.R`.

## 1. Qué se entrega

Cuatro tablas (una por especie), formato CSV o RDS, cada una en
`data/derived/`:

```
stepRep_cellyear_otitar.{csv,rds}    # Otis tarda
stepRep_cellyear_ptealc.{csv,rds}    # Pterocles alchata   (piloto recomendado)
stepRep_cellyear_pteori.{csv,rds}    # Pterocles orientalis
stepRep_cellyear_tettet.{csv,rds}    # Tetrax tetrax
```

Estructura (9 columnas, ~42–46 k filas por especie):

| columna                | tipo  | rango   | descripción |
|------------------------|-------|---------|-------------|
| `cells`                | int   | —       | índice de celda en la malla WorldClim 5 km (mismo `cells` que ya usa la pipeline) |
| `year`                 | int   | 2017–2023 | año |
| `n_checklists`         | int   | ≥ 1     | nº de checklists eBird de la especie en esa celda-año |
| `stepRep_strict_point` | num   | 0–1     | proporción de checklists cuyo píxel CAE en hábitat estricto (211/231/321/333) |
| `stepRep_broad_point`  | num   | 0–1     | idem con máscara amplia (+ 242/243/244) |
| `stepRep_strict_500m`  | num   | 0–1     | **principal** — proporción media de buffer 500 m en hábitat estricto |
| `stepRep_strict_1km`   | num   | 0–1     | sensibilidad — buffer 1 km, máscara estricta |
| `stepRep_broad_500m`   | num   | 0–1     | sensibilidad — buffer 500 m, máscara amplia |
| `stepRep_broad_1km`    | num   | 0–1     | sensibilidad — buffer 1 km, máscara amplia |

Variable principal sugerida para el ajuste piloto: **`stepRep_strict_500m`**.

Acompañar con: [reports/stepRep_diagnostics.md](stepRep_diagnostics.md)
(documentación + figuras), por si Raúl necesita el contexto antes de tocar
código.

## 2. Quick start (R)

```r
library(here)
library(data.table)

sp <- "ptealc"  # piloto
stepRep <- fread(here("data","derived",
                      paste0("stepRep_cellyear_", sp, ".csv")))
stepRep
#       cells  year  n_checklists  stepRep_strict_point ...  stepRep_strict_500m ...
# 1:    1234   2017             3                 0.000               0.045 ...
# 2:    1234   2018             5                 0.200               0.183 ...
```

Clave de unión con la tabla de modelos: **`cells`** (indices WorldClim) +
**`year`**.

## 3. Integración en `colext` (unmarked)

`stepRep_strict_500m` es una **yearly site covariate** (un valor por
(celda, año)) y entra en el `pformula`. Encaja directamente en el mismo
slot que tu pipeline ya usa para `Land_Cover_Type_1_Percent_Class_*`,
`NDVI`, `pr`, `tmmn`, `tmmx` (cf. `scripts/4_occupancy_models.R:276–286`).

```r
library(unmarked); library(data.table)

YEARS <- 2017:2023
sp    <- "ptealc"

# (1) Cargar la covariable y poner en wide: filas = cells, cols = years
stepRep_long <- fread(here("data","derived",
                           paste0("stepRep_cellyear_", sp, ".csv")))
stepRep_wide <- dcast(stepRep_long, cells ~ year,
                      value.var = "stepRep_strict_500m")

# (2) Alinear con el orden de filas del occ_wide_clean (sites del modelo)
stepRep_wide <- stepRep_wide[match(occ_wide_clean$cells, cells)]
stepRep_mat  <- as.matrix(stepRep_wide[, as.character(YEARS), with = FALSE])

# (3) Tratar NAs: cell-years sin checklists pueden quedar NA porque no
#     existen en stepRep_long. Decisión recomendada: imputar con la media
#     temporal del propio cell (intuitivo y conserva los efectos fijos);
#     si el cell no tiene NINGUNA observación en ningún año, asignar la
#     mediana global. unmarked descarta filas con covariables NA en p,
#     así que no imputar perderia esas visitas.
row_means <- rowMeans(stepRep_mat, na.rm = TRUE)
for (j in seq_len(ncol(stepRep_mat))) {
  na_mask <- is.na(stepRep_mat[, j])
  stepRep_mat[na_mask, j] <- row_means[na_mask]
}
stepRep_mat[is.na(stepRep_mat)] <- median(stepRep_long$stepRep_strict_500m,
                                          na.rm = TRUE)

# (4) Estandarizar (centrar y escalar). Guardar parametros para prediccion
#     consistente, igual que se hace ya con NDVI/pr/tmmn/tmmx en el script 4.
mu_step    <- mean(stepRep_mat)
sigma_step <- sd(stepRep_mat)
stepRep_mat_z <- (stepRep_mat - mu_step) / sigma_step
colnames(stepRep_mat_z) <- as.character(YEARS)

# (5) Anadir al unmarkedMultFrame existente:
occ_umf <- unmarkedMultFrame(
  y = y.cross,
  siteCovs = data.frame(siteCovs),
  yearlySiteCovs = list(
    years  = years_df,
    EVI    = EVI,
    Land_Cover_Type_1_Percent_Class_0  = Land_Cover_Type_1_Percent_Class_0,
    # ... resto igual ...
    NDVI = NDVI, pr = pr, tmmn = tmmn, tmmx = tmmx,
    stepRep = stepRep_mat_z              # <-- nueva yearly site cov
  ),
  obsCovs = list(
    duration = duration, effort = effort, observers = observers, time = time,
    NDVI_obs = NDVI_obs, pr_obs = pr_obs,
    topo_aspect_obs = topo_aspect_obs, topo_elev_obs = topo_elev_obs
  ),
  numPrimary = T_years
)

# (6) Anadir al pformula del modelo final. Para Pterocles alchata el
#     formula actual en R/model_configs.R es:
#       ~ effort + observers + time
#     Pasaria a:
#       ~ effort + observers + time + stepRep
Mod.final <- colext(
  psiformula     = cfg$psi_formula,
  gammaformula   = cfg$gamma_formula,
  epsilonformula = cfg$epsilon_formula,
  pformula       = ~ effort + observers + time + stepRep,
  data           = occ_umf
)

# (7) Guardar parametros de escalado igual que el resto de dynamic vars:
train_dyn_scale[["stepRep"]]$center <- mu_step
train_dyn_scale[["stepRep"]]$scale  <- sigma_step
saveRDS(train_dyn_scale,
        here("results", paste0(sp, "_train_dyn_scale.rds")))
```

## 4. Integración en `stPGOcc` (spOccupancy)

`stPGOcc` espera `det.covs` como lista; las covariables que varían por
(celda, año) se pasan como matrices `n_sites x n_primary_years` —misma
forma que `stepRep_mat_z` arriba. Anadir:

```r
det.covs <- list(
  effort    = effort_array,    # ya en uso
  duration  = duration_array,
  observers = observers_array,
  time      = time_array,
  stepRep   = stepRep_mat_z    # <-- aqui
)

fit <- stPGOcc(
  occ.formula = ~ ...,
  det.formula = ~ effort + duration + observers + time + stepRep,
  data        = list(y = y_array, occ.covs = occ.covs, det.covs = det.covs,
                     coords = coords, ...),
  n.batch = ..., batch.length = ..., n.neighbors = 10,
  ...
)
```

Equivalente a `colext`: una vez confirmado el efecto en `ptealc`, extender
a las otras tres especies.

## 5. Lectura del coeficiente

- **Direccion esperada**: positiva. Mayor stepRep en una celda-año implica
  que mas checklists caen en hábitat compatible -> mayor p de detección.
- **Magnitud orientativa**: el efecto sobre la probabilidad de detección
  por checklist deberia ser sustancial (al menos 0.1 logit por SD de
  stepRep) si la firma temporal observada en celdas focales es real.
- **Test de robustez**: comparar gamma y epsilon con/sin stepRep:
  - Si gamma baja y epsilon sube al anadir stepRep -> el sesgo
    diagnosticado en `reports/stepRep_diagnostics.md` §5 estaba inflando
    gamma y deflactando epsilon en los modelos previos. Ese es el
    resultado esperable dada la tendencia ascendente de stepRep en
    celdas focales 2017–2023.
  - Si los puntos estimados de gamma/epsilon apenas se mueven, stepRep
    no esta absorbiendo sesgo y la covariable puede dejarse fuera del
    modelo final, pero sigue siendo defendible reportarla como auxiliar.
- **Sensibilidades**: re-ajustar el modelo con `stepRep_strict_1km` y con
  `stepRep_broad_500m` y reportar (Tabla suplementaria). La firma de
  resultados debería mantener signo y orden de magnitud.

## 6. Banderas rojas

1. La tendencia temporal en celdas focales **sube** durante 2017–2023
   (~0.10–0.14 sobre la escala 0–1, ver §5 del informe diagnóstico). Esto
   es lo opuesto a la hipótesis inicial. La covariable sigue motivada,
   pero la dirección esperada del coeficiente y el sentido del sesgo
   detectado en gamma/epsilon va al revés de lo que pensábamos.
2. La correlación de Spearman entre `log(n_checklists)` y stepRep por
   celda-año es **prácticamente cero** (-0.005 a -0.05). Si tus modelos
   ya controlan por effort, **eso no captura este sesgo**: el sesgo
   opera vía estructura espacial intra-celda, no vía esfuerzo.
3. La covariable no es independiente del esfuerzo en cell-years muy
   pequeñas (n=1: stepRep_strict_500m solo puede tomar valores en
   {0, focal_strict_500m(checklist único)}). Si n_checklists < 2 te
   incomoda, podemos reportar dos versiones del modelo: con todas las
   celda-año y con n_checklists ≥ 3.
4. Imputar NAs por media-de-celda (pasos 3 del Quick start) afecta a
   <5% de celda-años. Probar también con la imputación a 0 (ningún
   checklist visitó hábitat estepario ese año) para ver sensibilidad.

## 7. Reproducibilidad

- Script generador: [scripts/build_stepRep.R](../scripts/build_stepRep.R).
  Es idempotente (skip-if-exists). Tiempo end-to-end ~5–6 min en M2 con
  CORINE recortada en disco.
- Inputs requeridos para regenerar: CORINE 2018 v2020_20u1 (100 m, EPSG:3035),
  los CSVs zero-fill de eBird por especie y la malla WorldClim. Documentado
  al inicio del script.
- Versionado: el commit
  `4ec50ca Add steppe-representativeness covariate (CLC2018) for detection sub-model`
  cierra el entregable. Los outputs binarios (rasters, RDS, CSVs derivadas)
  están en `data/derived/` y no van al repositorio (se regeneran).
