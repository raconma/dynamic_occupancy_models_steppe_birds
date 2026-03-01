# Informe Diagnostico: Modelos de Ocupacion Dinamica para Aves Esteparias

**Proyecto:** Dynamic occupancy models for steppe birds in Spain
**Para:** Raul
**Autor:** Guillermo Fandos
**Fecha:** 28 de febrero de 2026

---

## 1. Resumen ejecutivo

Se ajustaron 4 modelos de ocupacion dinamica (`unmarked::colext`) para
*Otis tarda*, *Pterocles alchata*, *Pterocles orientalis* y *Tetrax tetrax*
usando datos de eBird (2017--2022) en celdas de 2.5 km en Espana. Los
resultados son **biologicamente implausibles**: las tasas de extincion
estimadas oscilan entre el 92 y el 100 %, la colonizacion es inferior al
0.05 % y la ocupacion inicial es practicamente nula. Esto implica un
colapso demografico total que contradice los datos de atlas y censos
nacionales.

Tras un diagnostico exhaustivo, la **causa raiz no es un error de
programacion ni un problema de deteccion baja**. El problema es
**estructural**: la combinacion de (i) una matriz de deteccion con un
86 % de valores ausentes, (ii) el 69 % de los sitios con datos en un solo
ano y (iii) un crecimiento temporal del esfuerzo de muestreo de eBird hace
que el modelo sea incapaz de separar las transiciones ecologicas reales
del ruido generado por la cobertura de muestreo desigual.

El presente informe detalla la evidencia, identifica cuatro problemas
interconectados y propone un plan de accion concreto para reconfigurar el
analisis de cara a una publicacion.

---

## 2. Resumen de resultados actuales

### 2.1 Parametros estimados (interceptos en escala de probabilidad)

| Especie | psi1 | gamma | epsilon | p | Convergencia |
|---------|-----:|------:|--------:|----:|:------------:|
| *O. tarda*       | 0.06 % | 0.02 % | 100 %  | 20 % | OK  |
| *P. alchata*     | 0.05 % | 0.00 % | 100 %  | 28 % | **NO** |
| *P. orientalis*  | 0.08 % | 0.01 % | 100 %  | 17 % | OK  |
| *T. tetrax*      | 0.53 % | 0.01 % |  92 %  | 40 % | **NO** |

> **Lectura:** el modelo estima que la probabilidad de que un sitio medio
> estuviera ocupado en 2017 es del 0.05--0.5 %, que cada ano se extingue
> practicamente de todos los sitios (92--100 %) y que la colonizacion es
> despreciable. Esto es incompatible con la persistencia de estas especies
> en Espana.

### 2.2 Problemas de identificabilidad

| Especie | Params en boundary | SE > 10 | Correlaciones > 0.9 |
|---------|:---:|:---:|:---:|
| *O. tarda*       | 4 | 2 | 3 (una r = 0.9999) |
| *P. alchata*     | 9 | 8 | 7 (multiples r > 0.99) |
| *P. orientalis*  | 4 | 2 | 2 (una r = 0.9999) |
| *T. tetrax*      | 2 | 0 | 2 |

En *P. alchata* el optimizador no convergio y las estimaciones de
colonizacion y extincion tienen errores estandar del orden de
100--1 250, lo que indica que estos parametros carecen de contenido
informativo.

---

## 3. Diagnostico detallado

### 3.1 Esparsidad extrema de datos

La caracteristica mas llamativa de los datos es la proporcion de celdas
ausentes (NA) en la matriz de deteccion `y`:

| Especie | Sitios | Detecciones totales | Celdas NA | Sitios detectados |
|---------|-------:|--------------------:|----------:|------------------:|
| *O. tarda*       | 3 135 | 318 | 86.4 % | 118 (3.8 %) |
| *P. alchata*     | 3 489 | 385 | 85.3 % | 109 (3.1 %) |
| *P. orientalis*  | 3 489 | 238 | 85.3 % | 109 (3.1 %) |
| *T. tetrax*      | 3 135 | 240 | 86.4 % | 87 (2.8 %)  |

Con un 85--86 % de la matriz como NA, el modelo dispone de muy poca
informacion para estimar cuatro submodelos simultaneamente.

**Dato clave:** la ratio de detecciones frente a no-detecciones es del
orden de 0.008--0.013. Es decir, por cada deteccion hay entre 80 y
130 visitas sin observar la especie. Esto dificulta la estimacion de
cualquier covariable.

### 3.2 Cobertura temporal desigual

La mayoria de sitios carecen de datos interanuales, que es precisamente
la informacion que un modelo **dinamico** necesita:

| Numero de anos con datos | otitar / tettet | ptealc / pteori |
|:------------------------:|:---------------:|:---------------:|
| 1 ano  | 69.3 % | 67.2 % |
| 2 anos | 15.7 % | 15.9 % |
| 3 anos |  7.1 % |  6.9 % |
| 4 anos |  3.6 % |  4.4 % |
| 5 anos |  2.7 % |  3.2 % |
| 6 anos |  1.6 % |  2.5 % |

> **El 69 % de los sitios solo aportan datos de un unico ano.** Para
> estos sitios el modelo no puede estimar colonizacion ni extincion;
> simplemente observa presencia o ausencia puntual sin contexto temporal.

### 3.3 Desbalance temporal del esfuerzo de muestreo

El esfuerzo de eBird crece exponencialmente con el tiempo. En otitar /
tettet, el numero de sitios con datos pasa de 343 (11 %) en 2017 a
1 200 (38 %) en 2021:

| Ano  | Sitios con datos (otitar) | Replicas medias |
|------|:-------------------------:|:---------------:|
| 2017 | 343 (11 %)  | 0.5 |
| 2018 | 616 (20 %)  | 1.0 |
| 2019 | 776 (25 %)  | 1.3 |
| 2020 | 947 (30 %)  | 1.7 |
| 2021 | 1 200 (38 %) | 1.9 |
| 2022 | 1 120 (36 %) | 1.7 |

El resultado es que **la mitad tardio del estudio tiene el doble de
observaciones validas que la mitad temprana** (ratio = 1.9--2.0x). Esto
genera un artefacto critico:

1. Un sitio no tiene datos en 2017--2019 (ausencia de observador).
2. En 2020 aparece un observador y detecta la especie.
3. En 2021 el observador no vuelve -> dato ausente.
4. El modelo interpreta: colonizacion en 2020, posible extincion en 2021.
5. Resultado: gamma bajo pero positivo, epsilon muy alto.

### 3.4 Aparicion tardia de sitios detectados

La consecuencia directa del desbalance temporal es visible en la
distribucion de primeras detecciones:

| Primera deteccion | otitar | ptealc | pteori | tettet |
|:-----------------:|:------:|:------:|:------:|:------:|
| Anos 1-3 (2017-2019) | 40.7 % | 28.4 % | 32.1 % | 39.1 % |
| **Anos 4-6 (2020-2022)** | **59.3 %** | **71.6 %** | **67.9 %** | **60.9 %** |

> **Entre el 59 y el 72 % de los sitios donde se detecto la especie solo
> aparecen por primera vez en los anos tardios.** Esto no es colonizacion
> real: es un artefacto del aumento del esfuerzo de muestreo.

### 3.5 Alta proporcion de singletons

Los "singletons" son sitios con una sola deteccion en todo el periodo de
estudio (2017--2022). Representan entre el 35 y el 53 % de todos los
sitios detectados:

| Especie | Singletons (% de sitios detectados) | Mediana det/sitio |
|---------|:-----------------------------------:|:-----------------:|
| *O. tarda*       | 48.3 % | 2.0 |
| *P. alchata*     | 34.9 % | 2.0 |
| *P. orientalis*  | **53.2 %** | **1.0** |
| *T. tetrax*      | 39.1 % | 2.0 |

Un singleton no proporciona informacion para estimar persistencia o
transicion. Ademas, no es posible separar una deteccion unica genuina de
un falso positivo o un individuo en paso.

### 3.6 La deteccion NO es el problema principal

Contra lo que cabria esperar, la probabilidad de deteccion es razonable
en las cuatro especies:

- *O. tarda*: p = 0.20 (moderada)
- *P. alchata*: p = 0.28 (buena)
- *P. orientalis*: p = 0.17 (moderada-baja)
- *T. tetrax*: p = 0.40 (buena)

Los efectos de las covariables de deteccion son biologicamente coherentes
en las cuatro especies: mas duracion de la lista y mas esfuerzo aumentan
la deteccion; la lluvia la disminuye. **El problema no es que p sea baja,
sino que los datos carecen de la estructura temporal repetida que el
modelo dinamico necesita.**

---

## 4. Diagnostico causal

### Causa principal: datos sin estructura temporal para un modelo dinamico

El modelo `colext` asume que cada sitio se observa repetidamente a lo
largo de multiples anos (periodos primarios) y que dentro de cada ano se
realizan multiples visitas (periodos secundarios). Con la estructura
actual:

- **69 % de sitios** con datos en un solo ano: el modelo ve "ocupado en
  2020, NA en todos los demas anos". No puede determinar si el sitio
  estaba ocupado antes (psi1 subestimado) ni si sigue ocupado despues
  (epsilon sobreestimado).
- **86 % de celdas son NA**: la inmensa mayoria de las celdas sitio-ano
  no existen, asi que el modelo rellena con su prior implicito (en MLE,
  el que maximiza la verosimilitud), que en este caso es epsilon = 1.

### Causas secundarias

1. **Sobreparametrizacion de gamma y epsilon.** Con 87--118 sitios
   detectados, los submodelos de extincion tienen 4--6 parametros y los
   de colonizacion 3--5. Varios Land Cover classes estan correlacionados
   entre si (r > 0.99). Los SEs resultantes son catastróficos.

2. **Definicion de "sitio" demasiado granular.** El script
   `filter_repeat_visits()` define sitio como la combinacion
   `locality_id + observer_id`. Esto fragmenta las observaciones: la
   misma localidad visitada por dos observadores se trata como dos
   sitios independientes, cada uno con la mitad de las visitas.

3. **Escala espacial de 2.5 km.** Para especies con home ranges de
   10--50 km^2 (avutarda), una celda de 2.5 km (6.25 km^2) es mas
   pequena que el area vital del animal. El movimiento normal entre
   celdas se registra como extincion + colonizacion aparente.

---

## 5. Plan de accion para publicacion

### Fase 1: Restructuracion de datos (semanas 1-2)

#### 1.1 Cambiar la definicion de sitio

**Que hacer:** en `2_prepare_static_variables.R`, modificar
`filter_repeat_visits()`:

```r
# ACTUAL (demasiado restrictivo):
occ <- filter_repeat_visits(occ_raw,
  min_obs = 3, max_obs = 10,
  site_vars = c("locality_id", "observer_id"),
  annual_closure = TRUE)

# PROPUESTO (agrupar por celda de la grilla):
occ <- filter_repeat_visits(occ_raw,
  min_obs = 2, max_obs = 10,
  site_vars = c("grid_cell_id"),
  annual_closure = TRUE)
```

**Por que:** eliminar `observer_id` de la definicion de sitio consolida
todas las visitas a la misma localidad independientemente del observador,
aumentando las replicas por sitio.

#### 1.2 Aumentar la escala espacial

**Que hacer:** Pasar de celdas de 2.5 km a **10 km** (100 km^2).

**Por que:**
- 10 km^2 es mas acorde con el home range de las especies objetivo.
- Celdas mas grandes agregan mas listas por celda, reduciendo NAs.
- Se reduce la "extincion aparente" causada por movimiento entre celdas.
- Las celdas de 10 km son el estandar habitual en atlas de aves en
  Europa (EBBA2, Atlas de Aves Reproductoras de Espana).

**Impacto esperado:** con celdas de 10 km, cada celda contendra ~16 x
mas listas que una celda de 2.5 km (16 celdas pequenas por celda
grande), reduciendo drasticamente la proporcion de NAs y aumentando el
numero de anos con datos por sitio.

#### 1.3 Filtrar sitios por cobertura temporal minima

**Que hacer:** Anadir un filtro posterior que exija:

```r
# Retener solo sitios con datos en >= 3 anos,
# incluyendo al menos 1 ano temprano y 1 tardio
sites_keep <- site_summary %>%
  filter(
    n_years_with_data >= 3,
    has_early_year == TRUE,   # al menos 1 ano en 2017-2019
    has_late_year == TRUE     # al menos 1 ano en 2020-2022
  )
```

**Por que:** Un modelo dinamico necesita observaciones repetidas a traves
del tiempo. Un sitio con datos en un solo ano no aporta informacion sobre
transiciones. Exigir al menos 3 anos con al menos un dato temprano y uno
tardio garantiza que el modelo pueda estimar persistencia y transicion.

**Trade-off:** Se perdera la mayoria de sitios (actualmente ~69 % solo
tienen 1 ano). Es aceptable: es preferible tener 300 sitios con datos
temporales solidos que 3 000 sitios con un ano cada uno.

### Fase 2: Simplificacion del modelo (semanas 2-3)

#### 2.1 Modelo null como linea base

**Que hacer:** Antes del modelo completo, ajustar el modelo minimal:

```r
mod_null <- colext(
  psiformula     = ~ 1,
  gammaformula   = ~ 1,
  epsilonformula = ~ 1,
  pformula       = ~ 1,
  data = umf_filtrado
)
```

**Por que:** Si el modelo null ya da epsilon ~ 1.0, el problema sigue
siendo estructural en los datos. Si el null da parametros razonables,
el problema estaba en la sobreparametrizacion.

**Criterio de exito:** psi1 entre 5--30 %, epsilon < 30 %, gamma < 20 %
para las cuatro especies.

#### 2.2 Reducir el numero de covariables en gamma y epsilon

**Que hacer:** Limitar gamma y epsilon a un **maximo de 2 covariables**
cada uno (ademas del intercepto):

```r
# ACTUAL (sobreparametrizado para los datos disponibles):
epsilon_formula = ~ LC_Class_7 + LC_Class_12 + LC_Class_13 + NDVI + pr
# 6 parametros para ~109 sitios detectados

# PROPUESTO:
epsilon_formula = ~ NDVI + tmmx
# 3 parametros - mas razonable
```

**Criterio de seleccion:** elegir covariables con base biologica clara y
sin colinealidad. Evitar multiples clases de Land Cover en el mismo
submodelo (correlaciones > 0.95 entre ellas).

#### 2.3 Estrategia de seleccion de modelo por etapas

1. Modelo null para las cuatro especies -> verificar interceptos.
2. Anadir covariables de deteccion (p) con el resto null.
3. Anadir psi1 manteniendo gamma y epsilon null.
4. Probar covariables de gamma y epsilon una por una.
5. Seleccion por AIC entre modelos con 0--2 covariables en gamma/epsilon.

### Fase 3: Validacion (semanas 3-4)

#### 3.1 Comparacion con atlas

**Que hacer:** Comparar ocupacion estimada con el Atlas de Aves
Reproductoras de Espana (SEO/BirdLife) y EBBA2. El script
`5_validation.R` ya existe y puede adaptarse.

**Criterio:** correlacion positiva significativa (r > 0.3) entre
ocupacion predicha y presencia en atlas. Si no se alcanza, el modelo no
es util para inferencia espacial.

#### 3.2 Goodness-of-fit

Mantener el test parametric bootstrap y MacKenzie-Bailey del script
actual, pero ajustar `nsim` segun el tiempo disponible (nsim = 100 para
exploracion, nsim = 1 000 para publicacion).

#### 3.3 Analisis de sensibilidad

Documentar como varian los resultados segun:
- Escala espacial (5 vs 10 km).
- Filtro temporal minimo (2 vs 3 vs 4 anos con datos).
- Complejidad del modelo (null vs 1 cov vs 2 cov en gamma/epsilon).

Esto fortalecera el articulo mostrando robustez (o falta de ella).

### Fase 4: Alternativas si el modelo colext sigue sin funcionar (semanas 4-5)

Si tras los ajustes de las fases 1--3 los parametros siguen siendo
implausibles, considerar las siguientes alternativas:

#### 4.1 Modelos de ocupacion estaticos por ano

Ajustar modelos `occu()` independientes por ano y analizar la tendencia
de la ocupacion estimada a traves del tiempo. Es menos elegante que
`colext`, pero evita los problemas de estimacion de gamma/epsilon.

#### 4.2 Modelos Bayesianos con priors informativos

Usar `spOccupancy` o `ubms` con priors informativos sobre gamma y
epsilon basados en:
- Tasas de rotacion conocidas de censos nacionales (~5--15 % anual).
- Rangos biologicamente plausibles (epsilon < 30 %, gamma < 20 %).

Los priors informative regularizan los parametros y evitan boundary
estimates.

#### 4.3 Modelo de ocupacion multi-estacion sin estimar transiciones

Usar `unmarked::occuMS()` o un modelo de ocupacion multi-estacion que
estime ocupacion por ano sin asumir un proceso markoviano de
colonizacion/extincion.

---

## 6. Argumentacion para publicacion

### Puntos fuertes del estudio

1. **Primer analisis de ocupacion dinamica de aves esteparias a escala
   nacional en Espana** usando datos de citizen science.
2. **Cuatro especies con diferente ecologia** cubren un gradiente de
   tamano corporal y movilidad.
3. **Integracion de covariables ambientales de GEE** (NDVI, Land Cover,
   clima) a escala anual.
4. **Pipeline completamente reproducible** con R, renv y seed fijado.

### Narrativa sugerida para el articulo

El articulo deberia presentarse como un estudio que:

1. Evalua la viabilidad de usar datos de eBird para modelar dinamica de
   ocupacion de aves esteparias.
2. Identifica las limitaciones de los datos de citizen science para
   modelos dinamicos (cobertura temporal desigual, datos escasos).
3. Propone soluciones metodologicas (escala espacial, filtrado temporal,
   simplificacion de modelos).
4. Presenta los resultados corregidos como la mejor estimacion disponible
   de tendencias de ocupacion.

Este enfoque es honesto, util para la comunidad y publicable en revistas
como *Ecological Indicators*, *Biological Conservation*, *Diversity and
Distributions* o *Bird Conservation International*.

### Revistas objetivo

| Revista | IF | Encaje |
|---------|:---:|--------|
| *Biological Conservation* | 5.9 | Conservacion de esteparias + eBird |
| *Diversity and Distributions* | 5.0 | Macroecologia, modelos de distribucion |
| *Ecological Indicators* | 6.9 | Indicadores basados en citizen science |
| *Bird Conservation International* | 2.2 | Aves amenazadas, enfoque aplicado |
| *Ibis* | 2.3 | Ornitologia europea, metodologia |

---

## 7. Cronograma propuesto

| Semana | Actividad | Entregable |
|:------:|-----------|------------|
| 1 | Cambiar escala a 10 km, redefinir sitios | Nuevo dataset |
| 2 | Filtrado temporal, modelo null | Resultados preliminares |
| 3 | Seleccion de covariables, modelos finales | Tabla de modelos |
| 4 | Validacion con atlas, GOF, sensibilidad | Figuras de validacion |
| 5 | Redaccion del manuscrito | Borrador completo |
| 6 | Revision interna (Raul) | Manuscrito revisado |
| 7 | Envio a revista | Submission |

---

## 8. Resumen de cambios necesarios en el codigo

| Script | Cambio | Prioridad |
|--------|--------|:---------:|
| `2_prepare_static_variables.R` | Cambiar `site_vars`, grilla 10 km | CRITICA |
| `2_prepare_static_variables.R` | Filtro de cobertura temporal >= 3 anos | CRITICA |
| `R/model_configs.R` | Reducir covariables de gamma/epsilon a <= 2 | ALTA |
| `4_occupancy_models.R` | Anadir modelo null como referencia | ALTA |
| `4_occupancy_models.R` | Anadir seleccion de modelos por etapas | MEDIA |
| `5_validation.R` | Adaptacion a nueva escala | MEDIA |
| Nuevo script | Analisis de sensibilidad a escala y filtrado | MEDIA |

---

## Anexo A: Evidencia detallada por especie

### A.1 Otis tarda

```
Convergencia:     OK (code 0)
psi1 intercepto:  logit = -7.36  -> prob = 0.06 %
gamma intercepto: logit = -8.68  -> prob = 0.02 %
epsilon intercepto: logit = 60.4 -> prob = 100 % (BOUNDARY)
p intercepto:     logit = -1.38  -> prob = 20 %
Sitios detectados: 118 / 3 135 (3.8 %)
Singletons: 57 (48.3 % de detectados)
Sitios solo en anos tardios: 59.3 %
Correlacion critica: ext(Int) <-> ext(LC_13) r = 0.9999
SE extremos: ext(Int) = 54, ext(LC_13) = 142
```

### A.2 Pterocles alchata

```
Convergencia:     NO (code 1) *** NO CONVERGIO ***
psi1 intercepto:  logit = -7.62  -> prob = 0.05 %
gamma intercepto: logit = -169.4 -> prob = 0.00 % (BOUNDARY)
epsilon intercepto: logit = 90.8 -> prob = 100 % (BOUNDARY)
p intercepto:     logit = -0.95  -> prob = 28 %
Sitios detectados: 109 / 3 489 (3.1 %)
Singletons: 38 (34.9 % de detectados)
Sitios solo en anos tardios: 71.6 % (el peor caso)
SE extremos: col(Int) = 571, col(LC_0) = 736, col(LC_13) = 1 255
  ext(pr) = 77, ext(tmmx) = 211
Num. parametros con SE > 10: 8
```

### A.3 Pterocles orientalis

```
Convergencia:     OK (code 0)
psi1 intercepto:  logit = -7.08  -> prob = 0.08 %
gamma intercepto: logit = -8.94  -> prob = 0.01 %
epsilon intercepto: logit = 73.2 -> prob = 100 % (BOUNDARY)
p intercepto:     logit = -1.58  -> prob = 17 %
Sitios detectados: 109 / 3 489 (3.1 %)
Singletons: 58 (53.2 % de detectados, el peor caso)
Sitios solo en anos tardios: 67.9 %
Mediana de detecciones por sitio detectado: 1 (el mas bajo)
epsilon tiene 6 parametros para ~50 eventos de extincion aparente
```

### A.4 Tetrax tetrax

```
Convergencia:     NO (code 1) *** NO CONVERGIO ***
psi1 intercepto:  logit = -5.24  -> prob = 0.53 %
gamma intercepto: logit = -9.87  -> prob = 0.01 %
epsilon intercepto: logit = 2.47 -> prob = 92 % (alto, no extremo)
p intercepto:     logit = -0.41  -> prob = 40 % (la mejor deteccion)
Sitios detectados: 87 / 3 135 (2.8 %, el menor)
Singletons: 34 (39.1 % de detectados)
Sitios solo en anos tardios: 60.9 %
Nota: el "menos malo" de los 4, pero aun implausible y no converge
```

---

## Anexo B: Tabla comparativa de las 4 especies

| Metrica | otitar | ptealc | pteori | tettet |
|---------|-------:|-------:|-------:|-------:|
| Sitios totales | 3 135 | 3 489 | 3 489 | 3 135 |
| Sitios detectados | 118 | 109 | 109 | 87 |
| % sitios detectados | 3.8 | 3.1 | 3.1 | 2.8 |
| Detecciones totales | 318 | 385 | 238 | 240 |
| % NAs en matriz y | 86.4 | 85.3 | 85.3 | 86.4 |
| Ratio det./no-det. | 0.013 | 0.013 | 0.008 | 0.010 |
| Sitios con 1 ano datos | 69.3 % | 67.2 % | 67.2 % | 69.3 % |
| Sitios con 6 anos datos | 1.6 % | 2.5 % | 2.5 % | 1.6 % |
| Detectados solo anos 4-6 | 59.3 % | 71.6 % | 67.9 % | 60.9 % |
| Singletons (% det.) | 48.3 % | 34.9 % | 53.2 % | 39.1 % |
| p intercepto (prob) | 0.20 | 0.28 | 0.17 | 0.40 |
| psi1 (prob) | 0.0006 | 0.0005 | 0.0008 | 0.0053 |
| epsilon (prob) | 1.00 | 1.00 | 1.00 | 0.92 |
| gamma (prob) | 0.0002 | 0.0000 | 0.0001 | 0.0001 |
| Convergencia | OK | NO | OK | NO |
| Params boundary | 4 | 9 | 4 | 2 |
