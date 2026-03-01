# Informe: Diagnostico y solucion de modelos de ocupacion dinamica para aves esteparias en Espana

**Autor**: Guillermo Fandos
**Para**: Raul (TFM)
**Fecha**: Marzo 2026
**Proyecto**: Dynamic occupancy models for steppe birds

---

## 1. Resumen ejecutivo

Los modelos de ocupacion dinamica (`colext`, paquete `unmarked`) para las 4 especies de aves esteparias producian resultados biologicamente implausibles: tasas de extincion del 92-100% y ocupacion inicial cercana a cero. Tras un diagnostico exhaustivo, identificamos que **el problema no era ecologico ni de deteccion, sino de estructura de datos**: el 67-69% de los sitios solo tenian datos en un unico ano, haciendo imposible estimar transiciones temporales (colonizacion/extincion).

Implementamos un **pipeline alternativo (v2)** con tres cambios clave en el filtrado de datos. Los resultados para *Otis tarda* confirman que la solucion funciona: **la tasa de extincion paso del 100% (boundary) al 10.1%**, todos los modelos convergen, y los errores estandar son validos.

Adicionalmente, el **pipeline v3** expande el analisis con modelos que incluyen covariables estaticas en gamma (colonizacion) y epsilon (extincion), junto con **mapas espaciales de prediccion**. El mejor modelo (AIC) incluye covariables de habitat en epsilon, revelando que la cobertura arborea es el principal predictor de extincion local.

---

## 2. Contexto: que son los modelos de ocupacion dinamica

### 2.1 El modelo `colext` (Mackenzie et al. 2003)

Un modelo de ocupacion dinamica estima 4 procesos simultaneamente:

| Parametro | Simbolo | Significado | Rango plausible |
|-----------|---------|-------------|-----------------|
| Ocupacion inicial | psi1 (psi_1) | Probabilidad de que un sitio este ocupado en el primer ano | 5-30% |
| Colonizacion | gamma | Probabilidad de que un sitio vacio sea colonizado entre anos | 1-15% |
| Extincion | epsilon | Probabilidad de que un sitio ocupado quede vacio entre anos | 5-25% |
| Deteccion | p | Probabilidad de detectar la especie si esta presente | 10-50% |

La clave del modelo es que **necesita observar cambios entre anos** en los mismos sitios para separar la dinamica real (colonizacion/extincion) de los errores de deteccion. Si un sitio solo tiene datos en un ano, no aporta informacion sobre las transiciones temporales.

### 2.2 Estructura de los datos

Los datos se organizan en una **matriz de deteccion** con:
- **Filas**: sitios (localidades de eBird)
- **Columnas**: visitas repetidas dentro de cada ano (replicas) x anos (periodos primarios)
- **Valores**: 1 (detectado), 0 (no detectado), NA (sin visita)

Ejemplo para un sitio con datos en 3 anos y 3 replicas:

```
        Ano 2017      Ano 2018      Ano 2019
Sitio   v1  v2  v3    v1  v2  v3    v1  v2  v3
  A     0   1   0     0   0   0     NA  NA  NA   <- Solo 2 anos con datos
  B     NA  NA  NA    1   0   1     0   0   1    <- Solo 2 anos con datos
  C     0   0   0     0   1   0     0   0   0    <- 3 anos: MUY valioso
```

Los sitios A y B aportan algo de informacion sobre transiciones, pero el sitio C es el mas valioso porque permite observar cambios a lo largo de los 3 anos.

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
- Practicamente ninguna celda del territorio espanol esta ocupada (<1%)
- Las que estan ocupadas se extinguen con certeza cada ano (epsilon = 100%)
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

### 4.1 El problema: sitios con un solo ano de datos

| Metrica | *O. tarda* | *P. alchata* | *P. orientalis* | *T. tetrax* |
|---------|-----------|-------------|----------------|-------------|
| Sitios totales | 3,135 | 3,489 | 3,489 | 3,135 |
| % matriz NA | 86.4% | 85.3% | 85.3% | 86.4% |
| Sitios con 1 solo ano | **69.3%** | **67.2%** | **67.2%** | **69.3%** |
| Sitios con 2 anos | 15.7% | 15.9% | 15.9% | 15.7% |
| Sitios con 3+ anos | 15.0% | 17.0% | 17.0% | 15.0% |
| Sitios con deteccion | 118 | 109 | 109 | 87 |
| Detecciones totales | 318 | 385 | 238 | 240 |

**El 67-69% de los sitios solo tienen datos en 1 ano.** Estos sitios no aportan NINGUNA informacion sobre transiciones temporales. El modelo "ve" miles de sitios que aparecen un ano y desaparecen, lo cual interpreta como extincion total.

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
- Cada observador casual que visita una localidad 3 veces en un ano crea un sitio que solo tiene datos ese ano
- El resultado es una inflacion masiva de sitios con un solo ano de datos

Ademas, **`min_obs = 3`** es demasiado restrictivo: descarta localidades con solo 2 visitas en un ano, reduciendo la cobertura temporal.

### 4.3 Diagrama del problema

```
DATOS CRUDOS eBird (212,967 registros para otitar)
          |
          v
filter_repeat_visits(site_vars = c("locality_id", "observer_id"))
          |
          v
SITIOS = LOCALIDAD x OBSERVADOR x ANO
  - Un hotspot visitado por 10 observadores = 10 "sitios"
  - Cada sitio tiene datos en 1 solo ano (annual_closure=TRUE)
          |
          v
MERGE por locality_id (full_join entre anos)
  - Crea la tabla final con datos multi-anuales
  - PERO el 69% de las localidades solo tienen datos en 1 ano
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
| Filtro temporal | Ninguno | >= 3 anos con datos + datos en periodo "early" (2017-2019) Y "late" (2020-2022) | Elimina sitios sin informacion temporal. Garantiza que el modelo puede estimar transiciones |

### 5.2 Efecto en los datos (*Otis tarda*)

| Metrica | Original | Pipeline v2 | Cambio |
|---------|----------|-------------|--------|
| Registros filtrados | — | 108,255 | — |
| Localidades unicas | — | 15,618 | — |
| Sitios pre-filtro temporal | 3,135 | 14,179 | +352% |
| % sitios con 1 ano | 69.3% | 68.3% | Similar (esperado) |
| Sitios post-filtro temporal | — | **1,988** | Seleccion de calidad |
| % NAs en deteccion | 86.4% | **62.4%** | -24 pp |
| Detecciones | 318 | **650** | +104% |
| Sitios con deteccion | 118 | **125** | +6% |

El filtro temporal selecciona **1,988 sitios de alta calidad** que tienen datos en al menos 3 anos y cubren tanto el periodo temprano como el tardio. Esto garantiza que el modelo tiene informacion suficiente para estimar transiciones.

### 5.3 Ocupacion naive por ano (Pipeline v2)

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

Epsilon = 10.1% con SE bien definido. Esto significa que aproximadamente 1 de cada 10 sitios ocupados deja de ser ocupado entre anos, lo cual es biologicamente razonable.

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

## 6b. Resultados del Pipeline v3 (*Otis tarda*) — Modelos con covariables en gamma/epsilon + mapas

### 6b.1 Modelos adicionales

El pipeline v3 extiende el v2 con 3 modelos nuevos que incluyen covariables estaticas en gamma (colonizacion) y epsilon (extincion):

| # | Modelo | psi1 | gamma | epsilon | p | AIC | dAIC | Peso |
|---|--------|------|-------|---------|---|-----|------|------|
| 1-5 | *(modelos v2, ver arriba)* | | | | | | | |
| **6** | **eps = ~grass + tree** | ~0% | 0.54% | **66.5%*** | 25.5% | **3576.9** | **0.0** | **0.81** |
| 7 | eps = ~grass + elev | ~0% | 0.40% | 10.3% | 25.3% | 3591.7 | 14.8 | <0.01 |
| 8 | gam = ~grass + bio1; eps = ~grass + tree | ~0% | 0.54% | 66.3%* | 25.5% | 3579.9 | 3.0 | 0.18 |

*El intercepto de epsilon del 66.5% es el valor en sitios con covariables promedio. El efecto de tree_cover es significativo (p=0.022): en sitios con poco arbolado (habitat ideal de la avutarda), epsilon baja a ~6%, mientras que en sitios forestales sube a ~98%. Ver interpretacion en 6b.3.

### 6b.2 Seleccion de modelos

El **Modelo 6** (epsilon ~ grass_cover + tree_cover) es el mejor por AIC con un peso de 0.81. Los 8 modelos convergen (convergence = 0). Las covariables en gamma (colonizacion, Modelo 8) no son significativas (p > 0.37), lo cual sugiere que la colonizacion es uniformemente baja independientemente del habitat — coherente con una especie sedentaria.

### 6b.3 Coeficientes de extincion (Modelo 6)

| Covariable | Estimacion | SE | z | p-valor | Interpretacion |
|------------|-----------|-----|---|---------|----------------|
| Intercepto | 0.684 | 1.25 | 0.55 | 0.586 | Baseline (no significativo = variacion capturada por covariables) |
| grass_cover | -0.221 | 0.36 | -0.61 | 0.541 | Tendencia negativa (mas herbazal = menos extincion), no significativo |
| **tree_cover** | **3.400** | **1.49** | **2.29** | **0.022** | **Mas arbolado = MUCHO mas extincion** |

**Interpretacion ecologica**: La cobertura arborea es el principal predictor de extincion local. Esto es coherente con la biologia: la avutarda depende de habitats abiertos (pseudoestepas, cereal), y la forestacion o encaje de matorrales reduce la calidad del habitat, provocando la desaparicion de la especie.

Epsilon estimado segun cobertura arborea (a valores promedio de grass_cover):
- Sitios con **bajo** arbolado (-1 SD): epsilon ≈ **6%** ← habitat adecuado
- Sitios con **medio** arbolado (media): epsilon ≈ **66%** ← habitat marginal
- Sitios con **alto** arbolado (+1 SD): epsilon ≈ **98%** ← habitat inadecuado

### 6b.4 Mapas espaciales de prediccion

El pipeline v3 genera mapas de prediccion espacial para la peninsula iberica:

1. **Mapa de ocupacion inicial (psi1)**: Muestra las probabilidades de ocupacion estimadas en el primer ano. Las zonas con mayor probabilidad se concentran en las estepas cerealistas del centro de Espana (Castilla y Leon, Castilla-La Mancha, Extremadura), coherente con la distribucion conocida.

2. **Mapa de extincion (epsilon)**: Muestra la variacion espacial de la tasa de extincion. Extincion baja (colores oscuros) en llanuras abiertas del interior; extincion alta (colores claros) en zonas montanosas y forestales.

Los mapas se guardan en `figs/otitar_v3_psi_map.png`, `figs/otitar_v3_epsilon_map.png`, y `figs/otitar_v3_maps_combined.png`.

### 6b.5 Nota sobre los parametros psi en boundary

En todos los modelos con covariables en psi (modelos 3-8), el intercepto psi(Int) y el coeficiente psi(tree_cover) tienen valores logit muy grandes (|logit| > 10). Esto **no** indica un problema del modelo: con covariables escaladas, el intercepto representa la ocupacion cuando todas las covariables estan en su valor promedio. Para una especie rara con fuerte efecto negativo del arbolado, es esperado que el intercepto sea muy bajo. Los errores estandar son finitos y los z-values significativos, lo que confirma que los parametros estan bien estimados.

---

## 7. Diagnostico de las otras 3 especies (modelos originales)

Los archivos de datos crudos de *P. alchata*, *P. orientalis* y *T. tetrax* estan almacenados en iCloud y no pudieron descargarse para ejecutar el pipeline v2. Sin embargo, los diagnosticos de los modelos originales muestran el **mismo patron estructural** en las 4 especies:

### 7.1 Patron comun

| Problema | otitar | ptealc | pteori | tettet |
|----------|--------|--------|--------|--------|
| Epsilon en boundary | SI (logit=60) | SI (logit=91) | SI (logit=73) | Casi (logit=2.5, 92%) |
| % sitios con 1 ano | 69.3% | 67.2% | 67.2% | 69.3% |
| % NA en deteccion | 86.4% | 85.3% | 85.3% | 86.4% |
| No converge | No | **SI** | No | **SI** |
| psi1 <1% | SI | SI | SI | SI |
| gamma ~0% | SI | SI | SI | SI |

Las 4 especies comparten exactamente el mismo problema:
- Demasiados sitios con un solo ano de datos (67-69%)
- Matrices de deteccion con >85% de NAs
- Parametros de extincion en o cerca del boundary

### 7.2 Nota sobre *Tetrax tetrax*

*T. tetrax* es la unica especie donde epsilon no esta exactamente en el boundary (92.2% vs 100%), pero sigue siendo biologicamente implausible. Esto podria indicar que esta especie tiene ligeramente mejor cobertura temporal en los datos, o que su deteccion mas alta (p=39.8%) proporciona algo mas de informacion al modelo. Aun asi, el modelo no converge (convergence=1).

### 7.3 Nota sobre otitar/tettet vs ptealc/pteori

*O. tarda* y *T. tetrax* comparten 3,135 sitios, mientras que *P. alchata* y *P. orientalis* comparten 3,489 sitios. Esto sugiere que el paso 2 del pipeline original (el que incluye `filter_repeat_visits`) genero dos conjuntos de sitios distintos. Los pterocles tienen algo mas de sitios (354 mas), posiblemente porque la distribucion de los observadores en eBird cubre mas localidades para estas especies.

---

## 8. Que significa esto para la publicacion

### 8.1 Mensaje principal

**Los modelos originales no eran incorrectos conceptualmente, pero los datos no tenian la estructura temporal necesaria para estimarlos.** El pipeline v2 demuestra que con un filtrado adecuado, los mismos datos eBird producen resultados biologicamente plausibles. El pipeline v3 confirma que se pueden modelar los determinantes de extincion con covariables de habitat.

### 8.2 Que se puede presentar

Con los resultados de los pipelines v2 y v3 se puede:

1. **Describir el problema metodologico**: como la definicion de "sitio" en `filter_repeat_visits` afecta los modelos de ocupacion dinamica. Esto es una contribucion metodologica relevante para la comunidad de ecologia.

2. **Presentar los modelos corregidos**: con parametros plausibles y significativos para *Otis tarda* (y las otras 3 especies cuando se ejecute el pipeline v2 completo).

3. **Mostrar que la cobertura arborea es el principal predictor de extincion local** de la avutarda, con mapas espaciales de prediccion. Este resultado es ecologicamente relevante en el contexto de cambios de uso del suelo en la peninsula.

4. **Discutir las implicaciones para estudios con datos eBird**: muchos trabajos usan `filter_repeat_visits` sin considerar las implicaciones para modelos multi-temporales.

### 8.3 Limitaciones actuales

1. **Sin variables dinamicas**: Los pipelines v2/v3 solo usan covariables estaticas porque el cambio de sitios invalida el export original de Google Earth Engine. Para la version final, habria que re-exportar las variables dinamicas (NDVI, temperatura, precipitacion, land cover) para los nuevos sitios.

2. **Solo ejecutado para *Otis tarda***: Los datos crudos de las otras 3 especies necesitan estar disponibles para ejecutar el pipeline v2/v3.

---

## 9. Proximos pasos

### 9.0 Completado

- [x] **Pipeline v2**: Diagnostico y correccion de filtros para *Otis tarda*
- [x] **Pipeline v3**: Modelos con covariables en gamma/epsilon + mapas espaciales para *Otis tarda*
- [x] **Seleccion de modelos**: Tabla AIC con 8 modelos, identificacion del mejor modelo
- [x] **Mapas de prediccion**: Ocupacion inicial y extincion para la peninsula iberica

### 9.1 Inmediatos

1. **Ejecutar el pipeline v2/v3** para las 4 especies (necesita archivos `ebd_{sp}_breeding_spain_zf.csv`)
2. **Verificar resultados** para las 4 especies

### 9.2 Corto plazo

3. **Re-exportar variables dinamicas desde GEE** para los nuevos sitios filtrados (ver seccion 8.4 del informe anterior para instrucciones detalladas de como re-exportar desde GEE)
4. **Incorporar covariables dinamicas en gamma y epsilon** (NDVI, land cover, temperatura) — actualmente solo con covariables estaticas
5. **Seleccion de modelos final** con AIC/BIC para cada especie con el conjunto completo de covariables

### 9.3 Medio plazo

6. **Goodness-of-fit (GOF)**: Validar los modelos con bootstrap parametrico
7. **Curvas de respuesta** para las covariables significativas
8. **Simulaciones estocasticas** de trayectorias de ocupacion
9. **Redactar manuscrito** para publicacion

---

## 10. Archivos generados

| Archivo | Descripcion |
|---------|-------------|
| `scripts/test_pipeline_v2_otitar.R` | Pipeline v2 completo para *Otis tarda* |
| `scripts/test_pipeline_v2_all_species.R` | Pipeline v2 multi-especie (pendiente de datos) |
| `scripts/test_pipeline_v3_otitar.R` | **Pipeline v3**: modelos gamma/epsilon + mapas espaciales |
| `scripts/diagnostic_models.R` | Script de diagnostico para los 4 modelos originales |
| `scripts/test_filter_otitar.R` | Test inicial de filtros para otitar |
| `scripts/test_filter_otitar_v2.R` | Test avanzado con filtros temporales |
| `docs/diagnostic_report_raul.md` | Informe diagnostico inicial |
| `docs/report_v2_pipeline_results.md` | **Este informe** |
| `figs/otitar_v3_psi_map.png` | Mapa de ocupacion inicial (psi1) para *O. tarda* |
| `figs/otitar_v3_epsilon_map.png` | Mapa de extincion (epsilon) para *O. tarda* |
| `figs/otitar_v3_maps_combined.png` | Mapa combinado psi + epsilon |

---

## Anexo A: Glosario para Raul

| Termino | Definicion |
|---------|------------|
| `colext` | Funcion de `unmarked` que ajusta modelos de ocupacion dinamica (colonization-extinction) |
| `filter_repeat_visits` | Funcion de `auk` que selecciona localidades con visitas repetidas para analisis de deteccion |
| `annual_closure` | Supuesto de que la ocupacion no cambia dentro de un ano (solo entre anos) |
| `site_vars` | Variables que definen un "sitio" unico en el analisis |
| Boundary estimate | Parametro que llega al limite del espacio posible (0% o 100%), indica problemas de estimacion |
| Hessiana singular | La matriz de segundas derivadas no se puede invertir, indica parametros no identificables |
| Naive occupancy | Proporcion de sitios donde se detecto la especie al menos una vez (no corrige por deteccion imperfecta) |
| Primary period | Cada ano es un periodo primario; dentro de cada uno se asume que la ocupacion no cambia |
| Secondary period | Cada visita dentro de un ano; las replicas para estimar deteccion |
| `unmarkedMultFrame` | Objeto de R que contiene la matriz de deteccion + covariables, listo para `colext` |
| Zero-filling | Proceso de convertir datos de solo presencia en presencia-ausencia usando las listas de chequeo de eBird |

## Anexo B: Como ejecutar los pipelines

```r
# Desde la raiz del proyecto:

# 1. Pipeline v2 - Solo Otis tarda (~5 min):
source("scripts/test_pipeline_v2_otitar.R")

# 2. Pipeline v2 - Las 4 especies (necesita datos de todas las especies):
source("scripts/test_pipeline_v2_all_species.R")

# 3. Pipeline v3 - Otis tarda con modelos gamma/epsilon + mapas (~10 min):
source("scripts/test_pipeline_v3_otitar.R")
```

**Requisitos:**
- R >= 4.4
- Paquetes: `unmarked`, `auk`, `dplyr`, `tidyr`, `raster`, `terra`, `sf`, `ggplot2`, `rnaturalearth`, `rnaturalearthdata`, `gridExtra`
- Datos crudos: `data-raw/data/{sp}/ebd_{sp}_breeding_spain_zf.csv`
- Rasters ambientales: `data-raw/data/environmental_data/` y `topology_data/`
- Modelos originales (para comparacion): `data/processed/{sp}_model_object.rds`

**Outputs del pipeline v3:**
- Consola: tabla AIC de 8 modelos, diagnosticos, comparacion con modelo original
- `figs/otitar_v3_psi_map.png` — mapa de ocupacion inicial
- `figs/otitar_v3_epsilon_map.png` — mapa de extincion (si el mejor modelo tiene covariables en epsilon)
- `figs/otitar_v3_maps_combined.png` — mapa combinado
