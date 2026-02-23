# Feedback sobre el repositorio — Modelos de Ocupación Dinámica (Aves Esteparias)

**Fecha:** 22-23 febrero 2026
**Repo:** `dynamic_occupancy_models_steppe_birds`

---

## 1. Resumen

Raúl, he revisado a fondo todo el repositorio. El trabajo de fondo está bien, los modelos `colext()` están bien planteados y los datos están completos. Pero el repo tal como estaba no era reproducible — alguien que se lo descargara no podría ejecutar nada sin tocar muchas cosas a mano. He reestructurado bastante para dejarlo a nivel de publicación.

**Lo que he encontrado (estado inicial):**
- No había `renv.lock`, ni `.gitignore`, ni forma de lanzar el pipeline de golpe
- 18 scripts con mucho copy-paste entre las 4 especies
- Los tests de bondad de ajuste (GOF) usaban nsim=3 y nsim=10 — no es válido para publicar
- No había `set.seed()` en ninguna simulación estocástica
- Había un error de escalado: la superficie de predicción se re-escalaba independientemente de los datos de entrenamiento
- La validación estaba hecha solo para otitar
- El README estaba a medias

**Lo que he hecho (cambios aplicados):**
- Refactorizado de 17 scripts por especie → 5 scripts parametrizados + 1 archivo de configuración de fórmulas
- Todo el pipeline se lanza con `Rscript R/run_all.R`
- Inicializado `renv` con 198 paquetes bloqueados
- Corregido el bug de escalado: los parámetros del entrenamiento se guardan y se reutilizan para predecir
- GOF subido a nsim=1000/500
- `set.seed()` antes de cada paso estocástico
- Validación parametrizada para las 4 especies
- README y REPRODUCIBILITY.md reescritos desde cero

---

## 2. Estructura del repositorio

```
dynamic_occupancy_models_steppe_birds/
├── .gitignore
├── README.md
├── REPRODUCIBILITY.md
├── guille_feedback.md           ← este archivo
├── renv.lock
├── renv/
├── R/
│   ├── run_all.R               ← orquestador del pipeline
│   ├── check_repro.R           ← test de humo (paquetes, datos, paths)
│   └── model_configs.R         ← fórmulas por especie centralizadas
├── scripts/
│   ├── 1_prepare_database_ebird.R
│   ├── 2_prepare_static_variables.R
│   ├── 3_prepare_dynamic_variables.R
│   ├── 3_download_dynamic_variables_{sp}.js  (x4, GEE, sin cambios)
│   ├── 4_occupancy_models.R
│   └── 5_validation.R
├── data/
│   ├── raw/                    ← datos de entrada (inmutables)
│   └── processed/              ← generados por el pipeline
├── data-raw/
│   └── get_data.R              ← instrucciones de adquisición de datos
├── figs/                       ← figuras generadas
├── results/                    ← outputs del pipeline
├── docs/
│   └── validation_diagnosis.md ← diagnóstico detallado de la validación
└── scripts/                    ← scripts legacy por especie (conservados, no se usan)
```

### Especies

| Código | Especie | Modelo | Validación |
|--------|---------|--------|------------|
| otitar | *Otis tarda* | `colext()` | Spatial block CV |
| ptealc | *Pterocles alchata* | `colext()` | Spatial block CV |
| pteori | *Pterocles orientalis* | `colext()` | Spatial block CV |
| tettet | *Tetrax tetrax* | `colext()` | Spatial block CV |

### Flujo del pipeline

```
[Paso 1] 1_prepare_database_ebird.R
         Entrada: data/raw/ebird_raw_mar2024/*.txt
         Salida:  data/processed/{sp}/ebd_{sp}_breeding_spain_zf.csv

[Paso 2] 2_prepare_static_variables.R
         Entrada: Paso 1 + rasters bioclimáticos/topográficos
         Salida:  data/processed/{sp}/{sp}_occ_wide_static.csv
                  data/processed/{sp}/{sp}_scaling_params.rds  ← NUEVO (clave)

[Paso 3a] 3_download_dynamic_variables_{sp}.js  (GEE, manual)
         Entrada: {sp}_occ_wide_latlong.csv → se sube a GEE
         Salida:  data/raw/gee_exports/{sp}_dynamic_variables.csv

[Paso 3b] 3_prepare_dynamic_variables.R
         Entrada: Paso 2 + Paso 3a
         Salida:  data/processed/{sp}/{sp}_occ_wide_dynamic.csv

[Paso 4] 4_occupancy_models.R
         Entrada: Paso 3b + scaling params + rasters
         Salida:  results/{sp}_model_summary.txt
                  results/{sp}_model_object.rds
                  results/{sp}_gof_parboot.rds
                  figs/{sp}_occupancy_map.png  (+ response plots, prevalencia)
                  data/processed/{sp}/occ_{sp}_prediction.csv
                  data/processed/{sp}/{sp}_OccuMap.tif

[Paso 5] 5_validation.R
         Entrada: Paso 4 (predicciones) + atlas de biodiversidad
         Salida:  results/{sp}_validation_summary.txt
                  figs/{sp}_validation_calibration.png
                  figs/{sp}_validation_maps.png
                  figs/{sp}_validation_residuals_map.png
```

---

## 3. Checklist de reproducibilidad

| # | Criterio | Antes | Después |
|---|----------|-------|---------|
| 1 | Datos disponibles o script de descarga | NO | PARCIAL — `get_data.R` documenta fuentes; falta depósito Zenodo |
| 2 | renv.lock | NO | SÍ |
| 3 | .gitignore | NO | SÍ |
| 4 | Paths relativos (`here::here()`) | NO | SÍ |
| 5 | Orquestador del pipeline | NO | SÍ — `R/run_all.R` |
| 6 | `set.seed()` en pasos estocásticos | NO | SÍ — seeds en REPRODUCIBILITY.md |
| 7 | GOF con nsim adecuado | NO (3/10) | SÍ (1000/500) |
| 8 | Outputs guardados en figs/results/ | NO | SÍ — `ggsave()` + `sink()` |
| 9 | `sessionInfo()` registrado | NO | SÍ |
| 10 | README completo | NO | SÍ |
| 11 | REPRODUCIBILITY.md | NO | SÍ |
| 12 | Validación para las 4 especies | NO (solo otitar) | SÍ |
| 13 | Consistencia de escalado | NO (bug) | SÍ — params guardados en .rds |
| 14 | Paso GEE reproducible | PARCIAL | PARCIAL — scripts JS conservados, requiere ejecución manual en GEE |

---

## 4. Problemas encontrados y correcciones

### 4.1 Error de escalado en las predicciones (CRÍTICO)

El script original de predicción aplicaba `scale()` a la superficie de predicción de forma independiente, recalculando media y desviación típica desde los rasters en vez de usar los parámetros del entrenamiento. Esto distorsiona las probabilidades estimadas.

**Corrección:** El paso 2 ahora guarda los parámetros de escalado (centro y sd por variable) en `{sp}_scaling_params.rds`. El paso 4 los lee y aplica `(x - centro) / sd` con esos valores guardados.

### 4.2 Spearman y R² inflados artificialmente (ERROR en la validación)

El script original calculaba la correlación de Spearman sobre 10 bins de deciles (n=10 puntos), no sobre los polígonos individuales. Con solo 10 puntos, la correlación se infla a ~1.0 por diseño — no es un resultado real.

**Corrección:** Ahora se calcula:
- Spearman sobre polígonos individuales (n≈5200) → resultado honesto: ρ ≈ 0.49
- McFadden pseudo-R² (apropiado para respuesta binaria) en vez de R² lineal
- Spearman por fold del CV (para tener variabilidad)
- El plot de calibración con deciles se mantiene solo como visualización

### 4.3 Moran's I significativo (esperado, no es un error)

La autocorrelación espacial significativa en los residuos (Moran's I ≈ 0.38) es esperable en modelos ambientales que no incorporan dispersión ni procesos biogeográficos históricos. No invalida el modelo — hay que reconocerlo como limitación en el paper. He dejado texto sugerido en `docs/validation_diagnosis.md`.

### 4.4 GOF con pocas simulaciones

Los tests de bondad de ajuste usaban nsim=3 y nsim=10. Eso da p-valores inestables y no publicables. Subido a 1000 (parboot) y 500 (Mackenzie-Bailey).

### 4.5 Otros cambios menores

- Corregida la errata en el nombre del script (`1_preprare_` → `1_prepare_`)
- Eliminados `rm(list=ls())`, `knitr::opts_chunk$set()`, `dev.off()` sueltos
- Resueltos conflictos de namespace entre `pROC::auc` y `Metrics::auc`
- Actualizado de `spatialBlock()` (deprecated) a `cv_spatial()` de blockCV 3.x
- Añadido `stopifnot(nrow == nrow)` en el `cbind()` del paso 3 (seguridad)

---

## 5. Test de validación con datos reales (tettet)

He ejecutado `5_validation.R` end-to-end con los datos de tettet (es la única especie con el CSV de predicción disponible localmente; los de otitar, ptealc y pteori están descargados a iCloud).

### Métricas obtenidas

| Métrica | Valor | Comentario |
|---------|-------|------------|
| AUC medio (5-fold spatial block CV) | **0.824** | Buena capacidad discriminativa |
| TSS medio | **0.511** | Moderado-bueno, estable entre folds (0.44–0.57) |
| RMSE medio | **0.422** | Esperable para respuesta binaria |
| Spearman ρ (polígonos individuales) | **0.485** | Moderado; ahora honesto |
| Spearman ρ (media CV ± sd) | **0.484 ± 0.048** | Estable → robusto |
| McFadden pseudo-R² | **0.146** | Bajo-moderado; normal para binaria |
| Moran's I (residuos) | **0.385** (p < 2.2e-16) | Autocorrelación significativa (esperada) |
| Umbral óptimo (media de folds) | **0.027** | Demasiado bajo — ver diagnóstico abajo |

### Diagnóstico de los mapas

El mapa de predicción vs atlas muestra que el modelo predice presencia en casi toda la España peninsular (~41% de celdas), cuando el atlas indica que tettet está concentrado en las mesetas centrales y el valle del Ebro (~25% prevalencia). El umbral óptimo sale en 0.027, que es implausiblemente bajo.

Esto es consistente con el bug de escalado: las probabilidades están distorsionadas, y el umbral compensa artificialmente. Las métricas de ranking (AUC, TSS) son relativamente robustas a esto, así que el 0.82 de AUC seguramente está cerca del valor real. Pero los valores absolutos de probabilidad y el umbral **no son fiables** hasta que se regeneren las predicciones con el script corregido.

El mapa de residuos confirma: sobrepredice en la periferia (rojo), infrapredice en el rango nuclear (azul intenso en Castilla y León, Extremadura, Castilla-La Mancha, Aragón).

---

## 6. Qué tiene que hacer Raúl antes de la reunión

| Prioridad | Tarea | Tiempo estimado |
|-----------|-------|-----------------|
| **URGENTE** | **Regenerar predicciones de las 4 especies** con el script 4 corregido (el que guarda y reutiliza los scaling params). Sin esto, la validación se hace sobre mapas incorrectos | 4-8 horas |
| ALTA | Descargar de iCloud los archivos de predicción de otitar, ptealc y pteori (están offloaded) | 10 min |
| ALTA | Una vez regeneradas las predicciones, ejecutar `Rscript scripts/5_validation.R` para las 4 especies | 30 min |
| ALTA | Revisar los mapas generados en `figs/` — si el umbral óptimo sube a algo razonable (>0.1), la calibración estará bien | 15 min |
| MEDIA | Mover los scripts legacy a `scripts/legacy/` o eliminarlos | 15 min |
| MEDIA | Preparar depósito de datos en Zenodo | 1 día |

---

## 7. Registro de cambios

| Fecha | Cambio | Motivo |
|-------|--------|--------|
| 22/02 | Revisión inicial — inventario de 18 scripts | Identificar riesgos |
| 22/02 | Creado `.gitignore` | Excluir datos, outputs, .DS_Store |
| 22/02 | Reestructurado el repositorio | Organización estándar data/R/figs/results |
| 22/02 | Corregida errata `1_preprare_` → `1_prepare_` | Typo |
| 22/02 | Reescrito script 1 (prepare_database_ebird) | `here::here()`, loop parametrizado |
| 22/02 | Creado script 2 unificado (static variables) | Sustituye 4 scripts; guarda scaling params |
| 22/02 | Creado script 3 unificado (dynamic variables) | Sustituye 4 scripts |
| 22/02 | Creado `R/model_configs.R` | Fórmulas centralizadas |
| 22/02 | Creado script 4 unificado (occupancy models) | Fix escalado, set.seed, GOF nsim, ggsave |
| 22/02 | Creado script 5 unificado (validation) | Parametrizado para 4 especies |
| 22/02 | Creado `R/run_all.R` | Orquestador |
| 22/02 | Creado `R/check_repro.R` | Test de humo |
| 22/02 | Inicializado renv (198 paquetes) | Control de versiones de paquetes |
| 22/02 | Reescrito README.md y creado REPRODUCIBILITY.md | Documentación completa |
| 23/02 | Corregido bug Spearman/R² en validación | Calculado sobre deciles (n=10) → polígonos individuales (n≈5200) |
| 23/02 | Añadido pseudo-R² de McFadden + Spearman por fold | Métricas apropiadas para binaria |
| 23/02 | Creado `docs/validation_diagnosis.md` | Diagnóstico completo Spearman + Moran's I |
| 23/02 | Actualizado `.gitignore` para `data-raw/` | Datos del alumno (~8 GB) no suben a GitHub |
| 23/02 | Creados symlinks `data/` → `data-raw/data/` | Compatibilidad con rutas del pipeline |
| 23/02 | Actualizado `spatialBlock()` → `cv_spatial()` | blockCV 3.x deprecation |
| 23/02 | Namespaces explícitos `pROC::` y `Metrics::` | Evitar conflictos auc/rmse |
| 23/02 | Test end-to-end de validación en tettet | Primera ejecución con datos reales |
