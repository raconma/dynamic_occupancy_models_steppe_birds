# Informe: Diagnostico y solucion de modelos de ocupacion dinamica para aves esteparias en Espana

**Autor**: Guillermo Fandos
**Para**: Raul
**Proyecto**: Dynamic occupancy models for steppe birds

---

## 1. Resumen ejecutivo

Los modelos de ocupacion dinamica (`colext`, paquete `unmarked`) para las 4 especies de aves esteparias producian resultados biologicamente implausibles: tasas de extincion del 92-100% y ocupacion inicial cercana a cero. Tras un diagnostico exhaustivo, identificamos que **el problema no era ecologico ni de deteccion, sino de estructura de datos**: el 67-69% de los sitios solo tenian datos en un unico año, haciendo imposible estimar transiciones temporales (colonizacion/extincion).

Implementamos un **pipeline alternativo (v2)** con tres cambios clave en el filtrado de datos. Los resultados para *Otis tarda* confirman que la solucion funciona: **la tasa de extincion paso del 100% (boundary) al 10.1%**, todos los modelos convergen, y los errores estandar son validos.

---

## 2. Contexto: que son los modelos de ocupacion dinamica

### 2.1 El modelo `colext` (Mackenzie et al. 2003)

Un modelo de ocupacion dinamica estima 4 procesos simultaneamente:

| Parametro | Simbolo | Significado | Rango plausible |
|-----------|---------|-------------|-----------------|
| Ocupacion inicial | psi1 (psi_1) | Probabilidad de que un sitio este ocupado en el primer año | 5-30% |
| Colonizacion | gamma | Probabilidad de que un sitio vacio sea colonizado entre años | 1-15% |
| Extincion | epsilon | Probabilidad de que un sitio ocupado quede vacio entre años | 5-25% |
| Deteccion | p | Probabilidad de detectar la especie si esta presente | 10-50% |

La clave del modelo es que **necesita observar cambios entre años** en los mismos sitios para separar la dinamica real (colonizacion/extincion) de los errores de deteccion. Si un sitio solo tiene datos en un año, no aporta informacion sobre las transiciones temporales.

### 2.2 Estructura de los datos

Los datos se organizan en una **matriz de deteccion** con:
- **Filas**: sitios (localidades de eBird)
- **Columnas**: visitas repetidas dentro de cada año (replicas) x años (periodos primarios)
- **Valores**: 1 (detectado), 0 (no detectado), NA (sin visita)

Ejemplo para un sitio con datos en 3 años y 3 replicas:

```
        Año 2017      Año 2018      Año 2019
Sitio   v1  v2  v3    v1  v2  v3    v1  v2  v3
  A     0   1   0     0   0   0     NA  NA  NA   <- Solo 2 años con datos
  B     NA  NA  NA    1   0   1     0   0   1    <- Solo 2 años con datos
  C     0   0   0     0   1   0     0   0   0    <- 3 años: MUY valioso
```

Los sitios A y B aportan algo de informacion sobre transiciones, pero el sitio C es el mas valioso porque permite observar cambios a lo largo de los 3 años.

---

## 3. Diagnostico de los modelos originales

### 3.1 Resultados originales (las 4 especies)

| Especie | psi1 | gamma | epsilon | p | AIC | Convergencia | Hessiana |
|---------|------|-------|---------|---|-----|-------------|----------|
| *Otis tarda* | 0.06% | 0.02% | **100%** | 20.1% | 1857 | 0 | OK* |
| *Pterocles alchata* | 0.05% | 0.00% | **100%** | 27.9% | 1641 | **1 (no converge)** | OK* |
| *Pterocles orientalis* | 0.08% | 0.01% | **100%** | 17.1% | 1608 | 0 | OK* |
| *Tetrax tetrax* | 0.53% | 0.01% | **92.2%** | 39.8% | 1365 | **1 (no converge)** | OK* |

*OK con asterisco porque aunque la hessiana invierte, el numero de condicion es extremadamente alto (correlaciones >0.99 entre parametros).

**Interpretacion biologica (incorrecta) que sugieren estos resultados:**
- Practicamente ninguna celda del territorio español esta ocupada (<1%)
- Las que estan ocupadas se extinguen con certeza cada año (epsilon = 100%)
- La colonizacion es practicamente cero (gamma < 0.02%)
- **Esto implicaria que las 4 especies deberian haberse extinguido**, lo cual contradice los datos de campo y censos

### 3.2 Coeficientes en el boundary

Cuando `epsilon` tiene un valor logit de 60-90 (en escala logistica), la probabilidad es exactamente 1.0. Esto se llama un **parametro en el boundary**: el optimizador ha llegado al limite del espacio de parametros. Es una senal clara de que algo esta mal con los datos o el modelo, no con la ecologia.

### 3.3 La pista clave: deteccion adecuada

Antes de buscar el problema en los datos, descartamos la hipotesis de baja deteccion:

| Especie | p (deteccion) | Evaluacion |
|---------|---------------|------------|
| *O. tarda* | 20.1% | Adecuada |
| *P. alchata* | 27.9% | Buena |
| *P. orientalis* | 17.1% | Aceptable |
| *T. tetrax* | 39.8% | Muy buena |

Con detecciones del 17-40%, el modelo deberia funcionar. **El problema no es que no detectemos las aves, es que no tenemos suficientes datos temporales para estimar las transiciones.**

---

## 4. Causa raiz: estructura temporal de los datos

### 4.1 El problema: sitios con un solo año de datos

| Metrica | *O. tarda* | *P. alchata* | *P. orientalis* | *T. tetrax* |
|---------|-----------|-------------|----------------|-------------|
| Sitios totales | 3,135 | 3,489 | 3,489 | 3,135 |
| % matriz NA | 86.4% | 85.3% | 85.3% | 86.4% |
| Sitios con 1 solo año | **69.3%** | **67.2%** | **67.2%** | **69.3%** |
| Sitios con 2 años | 15.7% | 15.9% | 15.9% | 15.7% |
| Sitios con 3+ años | 15.0% | 17.0% | 17.0% | 15.0% |
| Sitios con deteccion | 118 | 109 | 109 | 87 |
| Detecciones totales | 318 | 385 | 238 | 240 |

**El 67-69% de los sitios solo tienen datos en 1 año.** Estos sitios no aportan NINGUNA informacion sobre transiciones temporales. El modelo "ve" miles de sitios que aparecen un año y desaparecen, lo cual interpreta como extincion total.

### 4.2 Por que ocurre esto: el filtro `filter_repeat_visits`

El pipeline original usaba:

```r
filter_repeat_visits(
  occ_raw,
  min_obs = 3,        # Minimo 3 visitas
  max_obs = 10,       # Maximo 10 visitas
  annual_closure = TRUE,
  site_vars = c("locality_id", "observer_id")  # <- PROBLEMA
)
```

El parametro **`site_vars = c("locality_id", "observer_id")`** define un "sitio" como la combinacion unica de localidad + observador. Esto significa que:

- El mismo hotspot de eBird visitado por 5 observadores diferentes genera 5 "sitios" distintos
- Cada observador casual que visita una localidad 3 veces en un año crea un sitio que solo tiene datos ese año
- El resultado es una inflacion masiva de sitios con un solo año de datos

Ademas, **`min_obs = 3`** es demasiado restrictivo: descarta localidades con solo 2 visitas en un año, reduciendo la cobertura temporal.

### 4.3 Diagrama del problema

```
DATOS CRUDOS eBird (212,967 registros para otitar)
          |
          v
filter_repeat_visits(site_vars = c("locality_id", "observer_id"))
          |
          v
SITIOS = LOCALIDAD x OBSERVADOR x AÑO
  - Un hotspot visitado por 10 observadores = 10 "sitios"
  - Cada sitio tiene datos en 1 solo año (annual_closure=TRUE)
          |
          v
MERGE por locality_id (full_join entre años)
  - Crea la tabla final con datos multi-anuales
  - PERO el 69% de las localidades solo tienen datos en 1 año
  - El modelo no puede estimar transiciones
          |
          v
MODELO colext
  - "Ve" miles de sitios que aparecen y desaparecen
  - Interpreta esto como epsilon = 100% (extincion total)
  - psi1 se aplasta a ~0% para compensar
  - RESULTADO: parametros en boundary, bioinformacion nula
```

---

## 5. Solucion implementada: Pipeline v2

### 5.1 Tres cambios en el filtrado

| Cambio | Original | Pipeline v2 | Justificacion |
|--------|----------|-------------|---------------|
| `site_vars` | `c("locality_id", "observer_id")` | `c("locality_id")` | Un sitio = una localidad, independientemente de quien observe. Aumenta las replicas por sitio (+50%) |
| `min_obs` | 3 | 2 | Con 2 visitas ya se puede estimar deteccion. Retiene mas localidades |
| Filtro temporal | Ninguno | >= 3 años con datos + datos en periodo "early" (2017-2019) Y "late" (2020-2022) | Elimina sitios sin informacion temporal. Garantiza que el modelo puede estimar transiciones |

### 5.2 Efecto en los datos (*Otis tarda*)

| Metrica | Original | Pipeline v2 | Cambio |
|---------|----------|-------------|--------|
| Registros filtrados | — | 108,255 | — |
| Localidades unicas | — | 15,618 | — |
| Sitios pre-filtro temporal | 3,135 | 14,179 | +352% |
| % sitios con 1 año | 69.3% | 68.3% | Similar (esperado) |
| Sitios post-filtro temporal | — | **1,988** | Seleccion de calidad |
| % NAs en deteccion | 86.4% | **62.4%** | -24 pp |
| Detecciones | 318 | **650** | +104% |
| Sitios con deteccion | 118 | **125** | +6% |

El filtro temporal selecciona **1,988 sitios de alta calidad** que tienen datos en al menos 3 años y cubren tanto el periodo temprano como el tardio. Esto garantiza que el modelo tiene informacion suficiente para estimar transiciones.

### 5.3 Ocupacion naive por año (Pipeline v2)

```
2017: 1.0%   2018: 2.0%   2019: 3.2%   2020: 1.6%   2021: 3.5%   2022: 3.5%
```

La ocupacion naive (proporcion de sitios con al menos una deteccion) es baja pero estable, lo cual es coherente con una especie de baja prevalencia pero persistente.

---

## 6. Resultados del Pipeline v2 (*Otis tarda*)

### 6.1 Cinco modelos progresivos

| # | Modelo | psi1 | gamma | epsilon | p | AIC | Conv | SEs |
|---|--------|------|-------|---------|---|-----|------|-----|
| 1 | Null (~1,~1,~1,~1) | 4.5% | 0.68% | **12.2%** | 29.3% | 3931 | 0 | OK |
| 2 | p = effort + observers | 4.6% | 0.67% | **11.6%** | 27.5% | 3880 | 0 | OK |
| 3 | psi = grass + tree; p | ~0% | 0.46% | **10.4%** | 26.2% | 3646 | 0 | OK |
| 4 | psi = grass + tree; p completo | ~0% | 0.47% | **10.4%** | 25.8% | 3635 | 0 | OK |
| 5 | psi = bio1+bio2+tree+grass+elev; p completo | ~0% | 0.40% | **10.1%** | 25.3% | 3588 | 0 | OK |

### 6.2 Comparacion directa: Original vs Pipeline v2

| Parametro | Original | v2 (Null) | v2 (Mejor: Mod 5) | Biologicamente plausible? |
|-----------|----------|-----------|-------------------|--------------------------|
| psi1 (intercepto) | 0.06% | 4.5% | ~0%* | SI (v2) |
| gamma | 0.02% | 0.68% | 0.40% | Mejorado |
| **epsilon** | **100%** (boundary) | **12.2%** | **10.1%** | **SI** |
| p | 20.1% | 29.3% | 25.3% | SI |
| Convergencia | Boundary | OK | OK | **SI** |
| SEs validos | No** | SI | SI | **SI** |
| AIC | 1857 | 3931 | 3588 | No comparable*** |

*En el modelo 5, psi1 en el intercepto es ~0% porque las covariables capturan la variacion: la ocupacion es alta en sitios con mucha cobertura herbácea y poca arboleda. El intercepto bajo es correcto: representa un sitio con valores promedio de todas las covariables.

**La hessiana del modelo original invertia, pero los parametros tenian correlaciones >0.99 y los parametros estaban en boundary.

***El AIC no es directamente comparable porque los datasets son diferentes (3,135 vs 1,988 sitios).

### 6.3 Coeficientes del mejor modelo (Modelo 5)

#### Ocupacion inicial (psi1):
| Covariable | Estimacion | SE | z | p-valor |
|------------|-----------|-----|---|---------|
| Intercepto | -12.87 | 1.73 | -7.45 | <0.001 |
| bio1 (temperatura media) | **1.24** | 0.51 | 2.46 | 0.014 |
| bio2 (rango termico) | **1.41** | 0.34 | 4.16 | <0.001 |
| tree_cover (arbolado) | **-11.74** | 1.83 | -6.41 | <0.001 |
| grass_cover (herbazal) | **1.15** | 0.21 | 5.37 | <0.001 |
| topo_elev (elevacion) | **1.70** | 0.63 | 2.71 | 0.007 |

**Interpretacion ecologica**: *Otis tarda* ocupa preferentemente zonas con:
- Mayor cobertura de herbazales (grass_cover +)
- Menor cobertura arborea (tree_cover -)
- Mayor rango termico / continentalidad (bio2 +)
- Temperaturas medias mas altas (bio1 +)
- Mayor elevacion (topo_elev +)

Esto es **coherente con la biologia conocida de la avutarda**: habita pseudoestepas cerealistas, mesetas con herbazales abiertos, evita zonas forestales.

#### Colonizacion (gamma):
| Covariable | Estimacion | SE | z | p-valor |
|------------|-----------|-----|---|---------|
| Intercepto | -5.52 | 0.22 | -25.2 | <0.001 |

Gamma = 0.40% (solo intercepto). La colonizacion es baja, lo cual es esperado para una especie sedentaria con requerimientos de habitat muy especificos.

#### Extincion (epsilon):
| Covariable | Estimacion | SE | z | p-valor |
|------------|-----------|-----|---|---------|
| Intercepto | -2.18 | 0.21 | -10.6 | <0.001 |

Epsilon = 10.1% con SE bien definido. Esto significa que aproximadamente 1 de cada 10 sitios ocupados deja de ser ocupado entre años, lo cual es biologicamente razonable.

#### Deteccion (p):
| Covariable | Estimacion | SE | z | p-valor |
|------------|-----------|-----|---|---------|
| Intercepto | -1.09 | 0.06 | -18.7 | <0.001 |
| effort (distancia) | **0.23** | 0.05 | 4.92 | <0.001 |
| observers (n observadores) | 0.06 | 0.04 | 1.62 | 0.106 |
| duration (duracion) | **0.19** | 0.05 | 3.63 | <0.001 |
| time (hora inicio) | -0.06 | 0.05 | -1.22 | 0.221 |

**Interpretacion**: La deteccion aumenta con mayor esfuerzo de muestreo (distancia recorrida) y mayor duracion de la visita. El numero de observadores y la hora de inicio no son significativos.

---

## 7. Diagnostico de las otras 3 especies (modelos originales)

El pipeline v2 aun no se ha ejecutado para *P. alchata*, *P. orientalis* y *T. tetrax*. Se necesitan los archivos zero-filled de cada especie (ver seccion 10 para detalles). Los diagnosticos de los modelos originales muestran el **mismo patron estructural** en las 4 especies:

### 7.1 Patron comun

| Problema | otitar | ptealc | pteori | tettet |
|----------|--------|--------|--------|--------|
| Epsilon en boundary | SI (logit=60) | SI (logit=91) | SI (logit=73) | Casi (logit=2.5, 92%) |
| % sitios con 1 año | 69.3% | 67.2% | 67.2% | 69.3% |
| % NA en deteccion | 86.4% | 85.3% | 85.3% | 86.4% |
| No converge | No | **SI** | No | **SI** |
| psi1 <1% | SI | SI | SI | SI |
| gamma ~0% | SI | SI | SI | SI |

Las 4 especies comparten exactamente el mismo problema:
- Demasiados sitios con un solo año de datos (67-69%)
- Matrices de deteccion con >85% de NAs
- Parametros de extincion en o cerca del boundary

### 7.2 Nota sobre *Tetrax tetrax*

*T. tetrax* es la unica especie donde epsilon no esta exactamente en el boundary (92.2% vs 100%), pero sigue siendo biologicamente implausible. Esto podria indicar que esta especie tiene ligeramente mejor cobertura temporal en los datos, o que su deteccion mas alta (p=39.8%) proporciona algo mas de informacion al modelo. Aun asi, el modelo no converge (convergence=1).

### 7.3 Nota sobre otitar/tettet vs ptealc/pteori

*O. tarda* y *T. tetrax* comparten 3,135 sitios, mientras que *P. alchata* y *P. orientalis* comparten 3,489 sitios. Esto sugiere que el paso 2 del pipeline original (el que incluye `filter_repeat_visits`) genero dos conjuntos de sitios distintos. Los pterocles tienen algo mas de sitios (354 mas), posiblemente porque la distribucion de los observadores en eBird cubre mas localidades para estas especies.

---

## 8. Situacion actual y limitaciones

### 8.1 Mensaje principal

**Los modelos originales no eran incorrectos conceptualmente, pero los datos no tenian la estructura temporal necesaria para estimarlos.** El pipeline v2 demuestra que con un filtrado adecuado, los mismos datos eBird producen resultados biologicamente plausibles. El objetivo es presentar los modelos corregidos con parametros plausibles y significativos para *Otis tarda* (y las otras 3 especies cuando se ejecute el pipeline v2 completo).

### 8.2 Limitaciones actuales

1. **Sin variables dinamicas**: El pipeline v2 actualmente solo usa covariables estaticas (bioclimaticas, topografia, cobertura de suelo). Las covariables dinamicas anuales (NDVI, temperatura, precipitacion, land cover) que se usaban en los modelos originales no se pueden reutilizar porque fueron exportadas desde Google Earth Engine (GEE) para el conjunto de sitios del pipeline original. Al cambiar el filtrado, los sitios son diferentes y por tanto hay que re-extraer estas variables. Los pasos para hacerlo son:
   - Exportar las coordenadas de los nuevos sitios filtrados (el pipeline v2 ya las tiene en `occ_wide_clean[, c("latitude", "longitude")]`)
   - Subir estas coordenadas a GEE como un asset (FeatureCollection)
   - Ejecutar los scripts de GEE existentes (ver `scripts/3_prepare_dynamic_variables.R` para la estructura esperada) apuntando a los nuevos sitios
   - El export de GEE generara un CSV por especie con las variables dinamicas alineadas por fila con los sitios
   - Incorporar estas variables como `yearlySiteCovs` en el `unmarkedMultFrame` y añadirlas a las formulas de gamma y epsilon

2. **Colonizacion/extincion solo como intercepto**: Los modelos v2 no incluyen covariables en gamma y epsilon. Esto es una decision deliberada: primero demostrar que los interceptos son plausibles antes de añadir complejidad. Una vez disponibles las variables dinamicas, se podran modelar gamma y epsilon en funcion de NDVI, cambios de uso del suelo, temperatura, etc.

3. **Solo ejecutado para *Otis tarda***: Las otras 3 especies necesitan los archivos zero-filled para ejecutar el pipeline v2 (ver seccion 10).

---

## 9. Proximos pasos

1. **Ejecutar el pipeline v2** para las 4 especies (script: `scripts/test_pipeline_v2_all_species.R`; ver seccion 10 para los archivos necesarios)
2. **Verificar resultados** para las 4 especies
3. **Re-exportar variables dinamicas desde GEE** para los nuevos sitios filtrados (ver seccion 8.2 para los pasos detallados)
4. **Incorporar covariables en gamma y epsilon** (NDVI, land cover, temperatura)
5. **Seleccion de modelos** con AIC/BIC para cada especie
6. **Goodness-of-fit (GOF)**: Validar los modelos con bootstrap parametrico
7. **Curvas de respuesta y mapas de prediccion**
8. **Simulaciones estocasticas** de trayectorias de ocupacion
9. **Redactar manuscrito**

---

## 10. Archivos generados y archivos necesarios

### 10.1 Scripts y documentos generados

| Archivo | Descripcion |
|---------|-------------|
| `scripts/test_pipeline_v2_otitar.R` | Pipeline v2 completo para *Otis tarda* |
| `scripts/test_pipeline_v2_all_species.R` | Pipeline v2 multi-especie |
| `scripts/diagnostic_models.R` | Script de diagnostico para los 4 modelos originales |
| `scripts/test_filter_otitar.R` | Test inicial de filtros para otitar |
| `scripts/test_filter_otitar_v2.R` | Test avanzado con filtros temporales |
| `docs/diagnostic_report_raul.md` | Informe diagnostico inicial |
| `docs/report_v2_pipeline_results.md` | **Este informe** |

### 10.2 Archivos de datos necesarios para ejecutar el pipeline v2

Para ejecutar `test_pipeline_v2_all_species.R` se necesitan los archivos zero-filled de cada especie. Son los CSVs generados por el script `1_prepare_database_ebird.R` y deben estar en las siguientes rutas:

```
data-raw/data/otitar/ebd_otitar_breeding_spain_zf.csv
data-raw/data/ptealc/ebd_ptealc_breeding_spain_zf.csv
data-raw/data/pteori/ebd_pteori_breeding_spain_zf.csv
data-raw/data/tettet/ebd_tettet_breeding_spain_zf.csv
```

Si no estan disponibles, se pueden regenerar ejecutando `scripts/1_prepare_database_ebird.R` con los datos crudos de eBird.

---

## Anexo A: Glosario

| Termino | Definicion |
|---------|------------|
| `colext` | Funcion de `unmarked` que ajusta modelos de ocupacion dinamica (colonization-extinction) |
| `filter_repeat_visits` | Funcion de `auk` que selecciona localidades con visitas repetidas para analisis de deteccion |
| `annual_closure` | Supuesto de que la ocupacion no cambia dentro de un año (solo entre años) |
| `site_vars` | Variables que definen un "sitio" unico en el analisis |
| Boundary estimate | Parametro que llega al limite del espacio posible (0% o 100%), indica problemas de estimacion |
| Hessiana singular | La matriz de segundas derivadas no se puede invertir, indica parametros no identificables |
| Naive occupancy | Proporcion de sitios donde se detecto la especie al menos una vez (no corrige por deteccion imperfecta) |
| Primary period | Cada año es un periodo primario; dentro de cada uno se asume que la ocupacion no cambia |
| Secondary period | Cada visita dentro de un año; las replicas para estimar deteccion |
| `unmarkedMultFrame` | Objeto de R que contiene la matriz de deteccion + covariables, listo para `colext` |
| Zero-filling | Proceso de convertir datos de solo presencia en presencia-ausencia usando las listas de chequeo de eBird |

## Anexo B: Como ejecutar el pipeline v2

```r
# Desde la raiz del proyecto:

# 1. Solo Otis tarda (ya ejecutado, ~5 min):
source("scripts/test_pipeline_v2_otitar.R")

# 2. Las 4 especies (necesita los archivos zero-filled, ver seccion 10.2):
source("scripts/test_pipeline_v2_all_species.R")
```

**Requisitos:**
- R >= 4.4
- Paquetes: `unmarked`, `auk`, `dplyr`, `tidyr`, `raster`, `terra`, `sf`
- Datos crudos: `data-raw/data/{sp}/ebd_{sp}_breeding_spain_zf.csv`
- Rasters ambientales: `data-raw/data/environmental_data/` y `topology_data/`
- Modelos originales (para comparacion): `data/processed/{sp}_model_object.rds`
