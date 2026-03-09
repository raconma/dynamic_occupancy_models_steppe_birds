# Reviewer Comments - Dynamic Occupancy Models Steppe Birds

Valoración crítica como revisor de revista top-tier (e.g., Diversity & Distributions, Biological Conservation).

**Veredicto general:** Major Revision

---

## 🔴 Críticos (5)

### 1. Ocupación muy baja en algunas especies
- P. orientalis y P. alchata muestran ocupaciones iniciales muy bajas (<5%)
- Puede afectar la estimabilidad de parámetros de colonización/extinción
- **Acción:** Discutir limitaciones, justificar que el modelo colext tolera baja prevalencia con suficientes sitios

### 2. No hay test de bondad de ajuste (GOF)
- No se presenta ningún test de GOF (e.g., parametric bootstrap con parboot)
- Imprescindible para validar que el modelo se ajusta a los datos
- **Acción:** Realizar parboot con 200-500 simulaciones, reportar c-hat y p-value

### 3. P. orientalis: modelo baseline como mejor modelo
- El mejor modelo para P. orientalis es m0 (baseline estático, sin covariables dinámicas)
- Sugiere que no hay suficiente señal en los datos o que la especie es demasiado rara
- **Acción:** Discutir honestamente, considerar si merece la pena incluir esta especie o tratarla como caso de estudio negativo

### 4. Serie temporal corta (6 años, 2017-2022)
- 6 años es el mínimo para modelos de ocupación dinámica
- Limita la capacidad de detectar tendencias reales vs. fluctuaciones estocásticas
- **Acción:** Reconocer como limitación, argumentar que eBird permite esta resolución temporal y que los resultados son coherentes ecológicamente. Explorar extensión a 2023-2025

### 5. No hay model averaging
- Se selecciona un único mejor modelo por especie sin considerar incertidumbre en la selección
- Modelos competitivos (deltaAIC < 2) deberían promediarse
- **Acción:** Calcular AIC weights, identificar modelos competitivos, promediar coeficientes con pesos AIC

---

## 🟡 Importantes (5)

### 6. Autocorrelación espacial
- Los modelos colext asumen independencia entre sitios
- Con datos eBird en grid, es probable que exista autocorrelación espacial
- **Acción:** Test de Moran's I sobre residuos/OccuMap, discutir implicaciones

### 7. Sesgo de muestreo eBird
- eBird tiene sesgos conocidos: geográfico (más cerca de ciudades/carreteras), temporal (fines de semana), observador
- No se discute ni se corrige este sesgo
- **Acción:** Mapas de esfuerzo de muestreo, discutir cómo el diseño del estudio (grid, filtros de esfuerzo) mitiga estos sesgos

### 8. Validación solo para T. tetrax
- El script de validación con datos del censo nacional solo cubre una especie
- **Acción:** Buscar datos de validación independientes para las otras especies (censos de avutarda, atlas, etc.)

### 9. Escala estandarizada en curvas de respuesta
- Las curvas de respuesta se presentan en escala estandarizada (z-scores), difícil de interpretar ecológicamente
- **Acción:** Presentar curvas también en escala original (%, °C, metros)

### 10. Procedencia de OccuMap
- No queda claro cómo se generan los mapas de ocupación (OccuMap)
- **Acción:** Documentar el método de generación de mapas espaciales

---

## 🟢 Menores (4)

### 11. Narrativa ecológica
- Los resultados se presentan de forma descriptiva, falta una narrativa ecológica que conecte los hallazgos con la biología de las especies
- **Acción:** Discutir por qué tree_cover es predictor universal de extinción (mecanismo: pérdida de hábitat estepario)

### 12. Comparación con literatura
- No se comparan las tasas de ocupación/extinción/colonización con estudios previos
- **Acción:** Buscar estudios de referencia para estas especies en la Península Ibérica

### 13. Implicaciones para la gestión
- Faltan recomendaciones concretas de gestión/conservación
- **Acción:** Traducir los resultados en recomendaciones prácticas (e.g., umbrales de cobertura arbórea, protección de cultivos cerealistas)

### 14. Tabla de detección
- No se presenta una tabla detallada de tasas de detección por año y especie
- **Acción:** Tabla con detección naive, detección modelada, esfuerzo por año

---

## Análisis implementados (script `scripts/reviewer_analyses.R`)

El script aborda los puntos: 2, 5, 6, 7, 9 y 14.

| Sección | Punto | Descripción | Estado |
|---------|-------|-------------|--------|
| 1 | #2 | GOF (parboot) | Requiere refit de modelos |
| 2 | #5 | Model averaging (AIC weights) | Requiere implementación manual |
| 3 | #6 | Autocorrelación espacial (Moran's I) | Listo |
| 4 | #7 | Mapas de esfuerzo muestreo | Listo |
| 5 | #9 | Curvas respuesta escala original | Listo |
| 6 | #14 | Tabla detección extendida | Listo |

**Nota:** Los puntos 3 (P. orientalis) y 4 (serie temporal corta) no se abordan computacionalmente — requieren discusión en el manuscrito. La extensión temporal a 2023+ está pendiente de nuevos datos.

---

## Notas técnicas

- Los modelos guardados (.rds) tienen incompatibilidad de versión con unmarked v1.5.1
- `fitted()`, `residuals()`, `parboot()`, `fitList()` fallan con error "no hay un slot de nombre 'formlist'"
- **Solución:** Refitear modelos desde UMF + fórmulas almacenadas → funciona correctamente
- `coef()`, `vcov()`, `AIC`, `projected()`, `smoothed()` funcionan sin refitear
