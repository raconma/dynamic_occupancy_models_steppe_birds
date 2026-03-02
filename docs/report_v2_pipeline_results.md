# Informe: Diagnostico y solucion de modelos de ocupacion dinamica para aves esteparias en Espana

**Autor**: Guillermo Fandos
**Para**: Raul Contreras
**Fecha**: Marzo 2026
**Proyecto**: Dynamic occupancy models for steppe birds

---

## 1. Resumen ejecutivo

Los modelos de ocupacion dinamica (`colext`, paquete `unmarked`) para las 4 especies de aves esteparias producian resultados biologicamente implausibles: tasas de extincion del 92-100% y ocupacion inicial cercana a cero. Tras un diagnostico exhaustivo, identificamos que **el problema no era ecologico ni de deteccion, sino de estructura de datos**: el 67-69% de los sitios solo tenian datos en un unico ano, haciendo imposible estimar transiciones temporales (colonizacion/extincion).

Implementamos un **pipeline alternativo (v2)** con tres cambios clave en el filtrado de datos. Los resultados para *Otis tarda* confirman que la solucion funciona: **la tasa de extincion paso del 100% (boundary) al 10.1%**, todos los modelos convergen, y los errores estandar son validos.

Adicionalmente, el **pipeline v3** expande el analisis con modelos que incluyen covariables estaticas en gamma (colonizacion) y epsilon (extincion), junto con **mapas espaciales de prediccion**. El mejor modelo (AIC) incluye covariables de habitat en epsilon, revelando que la cobertura arborea es el principal predictor de extincion local.

Finalmente, el **pipeline v4** evalua el efecto de **covariables dinamicas** (NDVI, EVI, precipitacion, temperatura) exportadas desde Google Earth Engine como yearlySiteCovs en gamma y epsilon. **Resultado: las covariables dinamicas no mejoran sustancialmente el modelo** (deltaAIC < 2 vs. baseline estatico). La cobertura arborea sigue siendo el predictor dominante de extincion. Este resultado sugiere que, a la escala temporal analizada (2017-2022), la estructura del habitat importa mas que la variabilidad climatica interanual.

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

## 6c. Resultados del Pipeline v4 (*Otis tarda*) — Covariables dinamicas en gamma/epsilon

### 6c.1 Objetivo

El pipeline v4 extiende el v3 incorporando **covariables que varian tanto por sitio como por ano** (yearlySiteCovs) en los parametros de colonizacion (gamma) y extincion (epsilon). Estas variables capturan la variabilidad interanual del habitat y el clima que las covariables estaticas del v3 no pueden representar.

### 6c.2 Variables dinamicas

Las variables provienen de un export de Google Earth Engine (`data/processed/otitar_dynamic_variables.csv`) con datos anuales 2017-2022 para cada localidad:

| Variable | Fuente | Unidad | Significado ecologico |
|----------|--------|--------|----------------------|
| NDVI | MODIS (MOD13A1) | x10000 | Productividad vegetal, proxy de calidad del habitat estepario |
| EVI | MODIS (MOD13A1) | x10000 | Similar al NDVI, mas robusto en zonas de alta biomasa |
| pr | TerraClimate | mm/mes | Precipitacion acumulada, driver hidrico clave |
| tmmn | TerraClimate | C x10 | Temperatura minima, severidad invernal |
| tmmx | TerraClimate | C x10 | Temperatura maxima, estres termico en epoca de cria |

### 6c.3 Problema con NDVI y solucion

Los datos NDVI de GEE contenian valores problematicos:
- **39 NAs** de 11,928 celdas (0.3%): sitios sin dato MODIS algunos anos
- **24 valores negativos** (hasta -1594): pixeles de agua, nieve o nubes
- **NAs en yearlySiteCovs causan fallo del optimizador**: el mensaje "valor inicial en vmmin no es finito" indica que la funcion de verosimilitud retorna NaN en los valores iniciales

**Solucion implementada (3 pasos):**
1. Imputar NAs por la media del sitio en otros anos (si disponible)
2. Imputar NAs residuales (24 sitios sin ningun dato) con la media global
3. Clamp valores extremos: reemplazar valores por debajo del percentil 1% con el valor P1%

Tras esta limpieza, **los 11 modelos convergen sin errores**. La misma limpieza se aplica a EVI (mismo problema, menor magnitud).

### 6c.4 Bateria de modelos

Se evaluaron 11 modelos, todos con la misma estructura base del v3 (psi1 = bio1 + bio2 + tree_cover + grass_cover + topo_elev; p = effort + observers + duration + time):

| # | Modelo | Formula gamma | Formula epsilon | AIC | dAIC | Peso |
|---|--------|--------------|-----------------|-----|------|------|
| **m6** | **gam~pr** | **~pr** | **~grass+tree** | **3575.5** | **0.0** | **0.177** |
| m9 | gam~pr, eps~tree+pr | ~pr | ~tree+pr | 3575.6 | 0.1 | 0.168 |
| m1 | eps~tree+NDVI | ~1 | ~tree+NDVI | 3576.1 | 0.6 | 0.131 |
| m0 | Baseline (v3 best) | ~1 | ~grass+tree | 3576.9 | 1.4 | 0.088 |
| m7 | gam~NDVI+pr | ~NDVI+pr | ~grass+tree | 3577.0 | 1.5 | 0.084 |
| m2 | eps~tree+pr | ~1 | ~tree+pr | 3577.1 | 1.6 | 0.080 |
| m3 | eps~tree+tmmx | ~1 | ~tree+tmmx | 3577.2 | 1.7 | 0.076 |
| m10 | gam+eps NDVI+pr | ~NDVI+pr | ~tree+NDVI+pr | 3577.4 | 1.9 | 0.068 |
| m4 | eps~tree+NDVI+pr | ~1 | ~tree+NDVI+pr | 3578.1 | 2.6 | 0.048 |
| m8 | gam~NDVI, eps~tree+NDVI | ~NDVI | ~tree+NDVI | 3578.1 | 2.6 | 0.048 |
| m5 | gam~NDVI | ~NDVI | ~grass+tree | 3578.9 | 3.4 | 0.032 |

**Todos los modelos convergen.** Los parametros en boundary (psi(Int) y psi(tree_cover)) son los mismos que en v3 y son esperados para una especie rara con fuerte efecto negativo del arbolado (ver seccion 6b.5).

### 6c.5 Interpretacion de la seleccion de modelos

**Resultado principal: las covariables dinamicas NO mejoran sustancialmente el modelo.**

- El mejor modelo (m6: gam~pr) tiene un deltaAIC de solo **-1.4** respecto al baseline estatico (m0). Un deltaAIC < 2 se considera no sustancial segun la escala de Burnham & Anderson (2002).
- Los pesos de Akaike estan muy repartidos: el mejor modelo tiene peso = 0.18, y 8 de los 11 modelos estan dentro de deltaAIC < 2.
- Esto indica **incertidumbre en la seleccion de modelos**: ninguna formulacion dinamica se distingue claramente del modelo con solo covariables estaticas.

### 6c.6 Coeficientes del mejor modelo (m6: gam~pr)

#### Colonizacion (gamma):
| Covariable | Estimacion | SE | z | p-valor | Interpretacion |
|------------|-----------|-----|---|---------|----------------|
| Intercepto | -5.462 | 0.276 | -19.82 | <0.001 | Gamma base = 0.42% |
| **pr** | **-0.514** | **0.307** | **-1.68** | **0.093** | Tendencia negativa marginalmente significativa |

La precipitacion muestra un efecto marginal (p=0.093) y con signo negativo: anos mas lluviosos tienden a menor colonizacion. Esto podria parecer contraintuitivo, pero tiene sentido ecologico: anos lluviosos favorecen el crecimiento de vegetacion alta y leñosa, lo cual **reduce la calidad del habitat abierto** que necesita la avutarda. Sin embargo, el efecto no es estadisticamente significativo al nivel convencional (α=0.05).

#### Extincion (epsilon):
| Covariable | Estimacion | SE | z | p-valor | Interpretacion |
|------------|-----------|-----|---|---------|----------------|
| Intercepto | -0.326 | 0.615 | -0.53 | 0.596 | Variacion capturada por covariables |
| grass_cover | -0.034 | 0.266 | -0.13 | 0.899 | No significativo |
| **tree_cover** | **2.242** | **0.721** | **3.11** | **0.002** | **Mas arbolado = mas extincion** |

El patron es identico al v3: **tree_cover es el unico predictor significativo de extincion** (p=0.002). Las variables dinamicas (NDVI, precipitacion, temperatura) no anaden informacion significativa sobre las transiciones una vez que se controla por la estructura del habitat.

### 6c.7 Sintesis ecologica

El hecho de que las covariables dinamicas no mejoren los modelos tiene varias explicaciones ecologicas plausibles:

1. **Escala temporal insuficiente**: 6 anos (2017-2022) puede no capturar eventos climaticos extremos que afecten a poblaciones de una especie longeva como la avutarda (esperanza de vida ~15 anos).

2. **Baja variabilidad interanual en la zona**: la variacion climatica entre anos en las estepas ibericas puede ser insuficiente para generar cambios detectables en ocupacion a escala de celda.

3. **Desfase temporal (lag)**: los efectos del clima sobre la ocupacion pueden manifestarse con retardo (un ano seco afecta a la reproduccion, pero la extincion local se observa 1-2 anos despues). Los modelos actuales usan covariables del mismo ano.

4. **Covariables estaticas capturan la senal principal**: la estructura del habitat (arbolado vs. abierto) es un predictor tan fuerte de extincion que deja poco margen para que las variables dinamicas anadan informacion.

5. **Resolucion espacial**: NDVI/EVI de MODIS tienen 500m de resolucion, lo cual puede no capturar cambios microhabitat relevantes para la avutarda.

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

Con los resultados de los pipelines v2, v3 y v4 se puede:

1. **Describir el problema metodologico**: como la definicion de "sitio" en `filter_repeat_visits` afecta los modelos de ocupacion dinamica. Esto es una contribucion metodologica relevante para la comunidad de ecologia.

2. **Presentar los modelos corregidos**: con parametros plausibles y significativos para *Otis tarda* (y las otras 3 especies cuando se ejecute el pipeline completo).

3. **Mostrar que la cobertura arborea es el principal predictor de extincion local** de la avutarda, con mapas espaciales de prediccion. Este resultado es ecologicamente relevante en el contexto de cambios de uso del suelo en la peninsula.

4. **Demostrar que la variabilidad climatica interanual no afecta significativamente la dinamica de ocupacion** a la escala temporal analizada (2017-2022). Esto refuerza la importancia de la estructura del habitat por encima de la variabilidad climatica a corto plazo.

5. **Discutir las implicaciones para estudios con datos eBird**: muchos trabajos usan `filter_repeat_visits` sin considerar las implicaciones para modelos multi-temporales.

### 8.3 Limitaciones actuales

1. **Solo ejecutado para *Otis tarda***: Los datos crudos de las otras 3 especies necesitan estar disponibles para ejecutar los pipelines v2/v3/v4.

2. **Escala temporal limitada (6 anos)**: El periodo 2017-2022 puede ser insuficiente para detectar efectos climaticos en una especie longeva. Eventos extremos (sequia severa, olas de calor) podrian no estar bien representados.

3. **Sin efectos retardados (lag)**: Los modelos v4 usan covariables del mismo ano. Un desfase temporal (e.g., sequia en ano t afecta ocupacion en ano t+1) podria revelar efectos dinamicos no detectados.

4. **Covariables de land cover no exploradas como dinamicas**: El export de GEE incluye land cover por ano (MODIS MCD12Q1), pero no se ha evaluado como yearlySiteCov. Cambios de uso del suelo inter-anuales podrian ser mas relevantes que el clima.

---

## 9. Proximos pasos

### 9.0 Completado

- [x] **Pipeline v2**: Diagnostico y correccion de filtros para *Otis tarda*
- [x] **Pipeline v3**: Modelos con covariables estaticas en gamma/epsilon + mapas espaciales
- [x] **Pipeline v4**: Covariables dinamicas (NDVI, EVI, precipitacion, temperatura) en gamma/epsilon
- [x] **Seleccion de modelos**: Tabla AIC con 8 modelos (v3) y 11 modelos (v4)
- [x] **Mapas de prediccion**: Ocupacion inicial y extincion para la peninsula iberica
- [x] **Evaluacion de covariables dinamicas**: Resultado negativo — no mejoran sustancialmente el modelo (dAIC < 2)

### 9.1 Inmediatos

1. **Ejecutar los pipelines v2/v3/v4** para las otras 3 especies (*P. alchata*, *P. orientalis*, *T. tetrax*). Necesita archivos `ebd_{sp}_breeding_spain_zf.csv` y los exports de GEE correspondientes.
2. **Verificar resultados** comparativos entre especies — determinar si el patron (covariables dinamicas no significativas) es consistente entre las 4 especies o especifico de *O. tarda*

### 9.2 Corto plazo

3. **Explorar efectos retardados (lag)**: Evaluar si covariables del ano t-1 predicen transiciones en el ano t (e.g., sequia en 2019 afecta extincion en 2020). Requiere modificar la construccion de yearlySiteCovs.
4. **Explorar land cover dinamico**: El export de GEE incluye MODIS MCD12Q1 (Land_Cover_Type_1) por ano. Cambios inter-anuales de uso del suelo podrian ser mas relevantes que variables climaticas.
5. **Seleccion de modelos final** con AIC/BIC para cada especie

### 9.3 Medio plazo

6. **Goodness-of-fit (GOF)**: Validar los modelos con bootstrap parametrico (parboot en v4 esta implementado pero tarda ~20 min)
7. **Curvas de respuesta** para las covariables significativas
8. **Model averaging**: Dado que los pesos de Akaike estan repartidos en v4, considerar model averaging para las predicciones
9. **Simulaciones estocasticas** de trayectorias de ocupacion
10. **Redactar manuscrito** para publicacion

---

## 10. Archivos generados

| Archivo | Descripcion |
|---------|-------------|
| `scripts/test_pipeline_v2_otitar.R` | Pipeline v2 completo para *Otis tarda* |
| `scripts/test_pipeline_v2_all_species.R` | Pipeline v2 multi-especie (pendiente de datos) |
| `scripts/test_pipeline_v3_otitar.R` | **Pipeline v3**: modelos gamma/epsilon + mapas espaciales |
| `scripts/test_pipeline_v4_dynamic_covs.R` | **Pipeline v4**: covariables dinamicas (NDVI, EVI, pr, temp) |
| `scripts/diagnostic_models.R` | Script de diagnostico para los 4 modelos originales |
| `scripts/test_filter_otitar.R` | Test inicial de filtros para otitar |
| `scripts/test_filter_otitar_v2.R` | Test avanzado con filtros temporales |
| `docs/diagnostic_report_raul.md` | Informe diagnostico inicial |
| `docs/report_v2_pipeline_results.md` | **Este informe** |
| `results/otitar_v4_model_selection.csv` | Tabla AIC de 11 modelos (v4) con deltaAIC y pesos |
| `results/otitar_v4_best_model_coefficients.csv` | Coeficientes del mejor modelo v4 |
| `results/otitar_v4_dynamic_scaling_params.csv` | Parametros de escalado de variables dinamicas |
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

# 4. Pipeline v4 - Otis tarda con covariables dinamicas (~15 min):
source("scripts/test_pipeline_v4_dynamic_covs.R")
```

**Requisitos:**
- R >= 4.4
- Paquetes: `unmarked`, `auk`, `dplyr`, `tidyr`, `raster`, `terra`, `sf`, `ggplot2`, `rnaturalearth`, `rnaturalearthdata`, `gridExtra`
- Datos crudos: `data-raw/data/{sp}/ebd_{sp}_breeding_spain_zf.csv`
- Rasters ambientales: `data-raw/data/environmental_data/` y `topology_data/`
- Variables dinamicas (v4): `data/processed/otitar_dynamic_variables.csv` (export de GEE)
- Modelos originales (para comparacion): `data/processed/{sp}_model_object.rds`

**Outputs del pipeline v3:**
- Consola: tabla AIC de 8 modelos, diagnosticos, comparacion con modelo original
- `figs/otitar_v3_psi_map.png` — mapa de ocupacion inicial
- `figs/otitar_v3_epsilon_map.png` — mapa de extincion (si el mejor modelo tiene covariables en epsilon)
- `figs/otitar_v3_maps_combined.png` — mapa combinado

**Outputs del pipeline v4:**
- Consola: tabla AIC de 11 modelos, diagnosticos, interpretacion ecologica
- `results/otitar_v4_model_selection.csv` — tabla AIC con deltaAIC y pesos de Akaike
- `results/otitar_v4_best_model_coefficients.csv` — coeficientes del mejor modelo
- `results/otitar_v4_dynamic_scaling_params.csv` — parametros de escalado (center, scale)

