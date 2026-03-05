# Dynamic Occupancy Models for Iberian Steppe Birds (2017--2023)

**Analysis Report -- March 2026**
**Branch:** `analysis-2023` | **Pipeline run:** 2026-03-05

---

## 1. Overview

This document summarises the results of fitting multi-season dynamic occupancy models to seven years of eBird citizen-science data (2017--2023) for four threatened steppe-dependent bird species in the Iberian Peninsula. The analysis pipeline consists of six steps: (1) eBird data filtering, (2) static covariate extraction, (3) dynamic covariate extraction, (4) dynamic occupancy model fitting via `unmarked::colext`, (5) independent validation against the Spanish Biodiversity Atlas, and (6) a spatial autocorrelation assessment using Bayesian occupancy models in `spOccupancy`.

### Study species

| Code | Species | Common name | IUCN | Notes |
|------|---------|-------------|------|-------|
| `otitar` | *Otis tarda* | Great Bustard | VU | Flagship steppe species |
| `ptealc` | *Pterocles alchata* | Pin-tailed Sandgrouse | LC | Arid specialist, rare in eBird |
| `pteori` | *Pterocles orientalis* | Black-bellied Sandgrouse | LC | Similar ecology, wider range |
| `tettet` | *Tetrax tetrax* | Little Bustard | NT | Severe decline across Europe |

### Data sources

- **Occurrence data:** eBird Basic Dataset, filtered with `auk` (complete checklists, stationary/travelling, duration 5--300 min, distance < 10 km, Jan--Dec each year).
- **Static covariates:** WorldClim bioclimatic variables (bio1, bio2), MODIS tree/grass cover, topographic variables (elevation, aspect, slope).
- **Dynamic covariates (yearly):** NDVI (MODIS), precipitation and temperature (TerraClimate), MODIS MCD12Q1 land-cover proportions (Classes 0, 6, 7, 12, 13) -- all extracted via Google Earth Engine.
- **Validation:** Spanish Biodiversity Atlas (atlas de biodiversidad, breeding-season distribution polygons for all species).

---

## 2. Dynamic Occupancy Models (`colext`)

We fitted MacKenzie et al. (2003) dynamic occupancy models using `unmarked::colext()`, estimating four parameter sets: initial occupancy ($\psi_1$), colonisation ($\gamma$), extinction ($\varepsilon$), and detection ($p$). All continuous covariates were standardised (z-scored) prior to fitting. Model selection followed AIC-based comparison; the final formulas are centralised in `R/model_configs.R`.

### 2.1 Model specifications

| Species | $\psi_1$ | $\gamma$ | $\varepsilon$ | $p$ | AIC |
|---------|----------|----------|----------------|-----|-----|
| `otitar` | bio1 + bio2 + tree_cover + grass_cover + topo_elev | NDVI + pr + tmmn + tmmx | LC6 + LC13 + tmmx | effort + observers + NDVI_obs + pr_obs + topo_aspect_obs | **2267.9** |
| `ptealc` | bio1 + bio2 + tree_cover + grass_cover + topo_aspect | NDVI + pr | pr + tmmx | time + duration + effort + observers + NDVI_obs + pr_obs | **1973.9** |
| `pteori` | bio2 + tree_cover + grass_cover | LC7 + NDVI + tmmn + tmmx | LC12 + NDVI + pr | time + duration + observers + pr_obs + topo_aspect_obs | **2060.7** |
| `tettet` | bio2 + tree_cover + grass_cover + topo_elev | LC12 | LC12 | effort + observers + time | **1780.1** |

*LC = Land_Cover_Type_1_Percent_Class (MODIS MCD12Q1)*

> [!NOTE]
> **For Raul:** The formulas for `ptealc`, `pteori`, and especially `tettet` had to be substantially simplified from the original specifications to resolve convergence failures caused by complete/quasi-complete separation. Details in Section 7.

### 2.2 Parameter estimates

#### *Otis tarda* (Great Bustard) -- AIC = 2267.9

| Sub-model | Parameter | Estimate | SE | *z* | *P* |
|-----------|-----------|----------|----|-----|-----|
| Initial ($\psi_1$) | Intercept | -8.085 | 0.867 | -9.33 | <0.001 |
| | bio1 | 0.727 | 0.372 | 1.95 | 0.051 |
| | bio2 | 0.618 | 0.271 | 2.28 | **0.023** |
| | tree_cover | -6.449 | 0.884 | -7.30 | **<0.001** |
| | grass_cover | 0.355 | 0.151 | 2.35 | **0.019** |
| | topo_elev | 1.213 | 0.442 | 2.74 | **0.006** |
| Colonisation ($\gamma$) | Intercept | -10.00 | 2.015 | -4.96 | <0.001 |
| | NDVI | -1.45 | 0.583 | -2.49 | **0.013** |
| | pr | -2.30 | 1.475 | -1.56 | 0.119 |
| | tmmn | -4.52 | 1.949 | -2.32 | **0.020** |
| | tmmx | 2.59 | 1.684 | 1.54 | 0.124 |
| Extinction ($\varepsilon$) | Intercept | -1.25 | 0.426 | -2.94 | **0.003** |
| | LC6 | 1.86 | 1.537 | 1.21 | 0.225 |
| | LC13 | 2.80 | 0.833 | 3.36 | **<0.001** |
| | tmmx | -1.51 | 0.761 | -1.99 | **0.047** |
| Detection ($p$) | Intercept | -1.335 | 0.152 | -8.81 | <0.001 |
| | effort | 0.225 | 0.060 | 3.78 | **<0.001** |
| | observers | 0.223 | 0.053 | 4.23 | **<0.001** |
| | NDVI_obs | 0.428 | 0.119 | 3.61 | **<0.001** |
| | pr_obs | -0.758 | 0.200 | -3.80 | **<0.001** |
| | topo_aspect_obs | 0.204 | 0.079 | 2.57 | **0.010** |

**Interpretation:** Initial occupancy strongly driven by low tree cover and open grassland habitats at moderate-to-high elevations. Urban land cover (LC13) substantially increases local extinction risk ($\beta$ = 2.80, *P* < 0.001). Colonisation appears mediated by climatic conditions, with lower minimum temperatures favouring colonisation. Detection increases with survey effort, number of observers, and NDVI (greener vegetation = easier to detect in open fields), and decreases with precipitation.

#### *Pterocles alchata* (Pin-tailed Sandgrouse) -- AIC = 1973.9

| Sub-model | Parameter | Estimate | SE | *z* | *P* |
|-----------|-----------|----------|----|-----|-----|
| Initial ($\psi_1$) | Intercept | -6.705 | 0.748 | -8.96 | <0.001 |
| | bio1 | -0.615 | 0.201 | -3.06 | **0.002** |
| | bio2 | 1.517 | 0.293 | 5.18 | **<0.001** |
| | tree_cover | -4.538 | 0.749 | -6.06 | **<0.001** |
| | grass_cover | 0.519 | 0.150 | 3.45 | **<0.001** |
| | topo_aspect | -0.066 | 0.132 | -0.50 | 0.617 |
| Colonisation ($\gamma$) | Intercept | -6.522 | 0.413 | -15.79 | <0.001 |
| | NDVI | -1.081 | 0.296 | -3.65 | **<0.001** |
| | pr | 0.162 | 0.359 | 0.45 | 0.652 |
| Extinction ($\varepsilon$) | Intercept | -3.46 | 4.45 | -0.78 | 0.437 |
| | pr | -52.07 | 30.67 | -1.70 | 0.090 |
| | tmmx | -70.87 | 40.79 | -1.74 | 0.082 |
| Detection ($p$) | Intercept | -1.131 | 0.210 | -5.40 | <0.001 |
| | time | -0.364 | 0.071 | -5.13 | **<0.001** |
| | duration | 0.900 | 0.094 | 9.61 | **<0.001** |
| | effort | 0.197 | 0.075 | 2.64 | **0.008** |
| | observers | 0.114 | 0.077 | 1.48 | 0.138 |

> [!WARNING]
> **For Raul -- Extinction coefficients `ptealc`:** The extinction estimates for `pr` (-52.1) and `tmmx` (-70.9) are extremely large in absolute value despite borderline-significant *P*-values. This is a classic symptom of sparse data in the extinction sub-model (only ~18 extinction events across 7 years at 4131 sites). The SEs are also very large (30--40), so these estimates are highly uncertain. We should discuss whether to simplify this sub-model to intercept-only, as was done for `tettet`. Despite this, the overall model AIC (1973.9) and validation metrics (AUC = 0.844) are acceptable.

**Interpretation:** Temperature seasonality (bio2) is the strongest positive predictor of initial occupancy, together with grass cover and the absence of trees. Colonisation is negatively associated with NDVI, consistent with the species' preference for bare/sparse ground. The extinction parameters are unreliable due to data sparsity (see warning above).

#### *Pterocles orientalis* (Black-bellied Sandgrouse) -- AIC = 2060.7

| Sub-model | Parameter | Estimate | SE | *z* | *P* |
|-----------|-----------|----------|----|-----|-----|
| Initial ($\psi_1$) | Intercept | -5.632 | 0.573 | -9.83 | <0.001 |
| | bio2 | 0.359 | 0.205 | 1.76 | 0.079 |
| | tree_cover | -3.954 | 0.603 | -6.55 | **<0.001** |
| | grass_cover | 0.885 | 0.144 | 6.16 | **<0.001** |
| Colonisation ($\gamma$) | Intercept | -14.78 | 3.805 | -3.88 | **<0.001** |
| | LC7 | -9.80 | 6.463 | -1.52 | 0.129 |
| | NDVI | -0.563 | 0.492 | -1.14 | 0.253 |
| | tmmn | -8.105 | 2.670 | -3.04 | **0.002** |
| | tmmx | 9.135 | 2.682 | 3.41 | **<0.001** |
| Extinction ($\varepsilon$) | Intercept | -0.239 | 0.349 | -0.68 | 0.494 |
| | LC12 | -0.831 | 0.263 | -3.16 | **0.002** |
| | NDVI | 2.047 | 0.557 | 3.68 | **<0.001** |
| | pr | 1.687 | 0.734 | 2.30 | **0.022** |
| Detection ($p$) | Intercept | -1.585 | 0.157 | -10.07 | <0.001 |
| | time | -0.244 | 0.071 | -3.44 | **<0.001** |
| | duration | 0.423 | 0.083 | 5.11 | **<0.001** |

**Interpretation:** The cleanest extinction sub-model among the four species. Cropland cover (LC12) significantly reduces local extinction risk ($\beta$ = -0.831, *P* = 0.002), consistent with the species' reliance on cereal pseudo-steppe. Higher NDVI and precipitation increase extinction probability, possibly reflecting encroachment of denser vegetation into previously suitable open habitats. Colonisation is strongly modulated by temperature: the combination of low tmmn and high tmmx (i.e., continental climate with large diurnal/seasonal amplitude) favours colonisation.

#### *Tetrax tetrax* (Little Bustard) -- AIC = 1780.1

| Sub-model | Parameter | Estimate | SE | *z* | *P* |
|-----------|-----------|----------|----|-----|-----|
| Initial ($\psi_1$) | Intercept | -6.224 | 0.588 | -10.59 | <0.001 |
| | bio2 | 1.149 | 0.256 | 4.49 | **<0.001** |
| | tree_cover | -3.724 | 0.597 | -6.24 | **<0.001** |
| | grass_cover | 0.772 | 0.159 | 4.87 | **<0.001** |
| | topo_elev | -0.106 | 0.165 | -0.64 | 0.520 |
| Colonisation ($\gamma$) | Intercept | -7.294 | 0.594 | -12.28 | <0.001 |
| | LC12 | 0.902 | 0.437 | 2.06 | **0.039** |
| Extinction ($\varepsilon$) | Intercept | -1.836 | 0.306 | -6.01 | <0.001 |
| | LC12 | -0.565 | 0.278 | -2.03 | **0.042** |
| Detection ($p$) | Intercept | -0.621 | 0.105 | -5.93 | <0.001 |
| | effort | 0.218 | 0.075 | 2.91 | **0.004** |
| | observers | 0.069 | 0.052 | 1.32 | 0.186 |
| | time | -0.219 | 0.072 | -3.02 | **0.003** |

> [!IMPORTANT]
> **For Raul -- tettet model simplification:** This species required the most aggressive formula simplification. The original model with 5 covariates in $\varepsilon$ and `pr` in $\gamma$ failed entirely (only 25 transitions from occupied sites: 9 extinctions, 16 persistences). Progressive model building revealed that **any model containing `pr` (precipitation) in colonisation or extinction diverges**, likely due to complete separation. Only `Land_Cover_Type_1_Percent_Class_12` (cropland) could be supported in both dynamic sub-models. Despite this, the model has the **lowest AIC (1780.1)** and the **highest validation Spearman correlation (0.521)** among all four species. The ecological signal is clear and consistent: cropland promotes colonisation and reduces extinction.

**Interpretation:** Cropland (LC12) is the single most important dynamic driver for *T. tetrax*: it significantly promotes colonisation ($\beta$ = +0.90, *P* = 0.039) and reduces local extinction ($\beta$ = -0.57, *P* = 0.042). Initial occupancy follows the expected steppe pattern (open grassland, low tree cover, high temperature seasonality). Detection is highest with greater survey effort and during earlier times of day.

---

## 3. Occupancy Maps (Initial Occupancy Predictions)

The maps below show predicted initial occupancy probability ($\psi_1$) across the study area, based on static bioclimatic and land-cover covariates.

### *Otis tarda*
![otitar occupancy map](figures/otitar_occupancy_map.png)

### *Pterocles alchata*
![ptealc occupancy map](figures/ptealc_occupancy_map.png)

### *Pterocles orientalis*
![pteori occupancy map](figures/pteori_occupancy_map.png)

### *Tetrax tetrax*
![tettet occupancy map](figures/tettet_occupancy_map.png)

> [!NOTE]
> **For Raul:** The predicted $\psi_1$ probabilities are very low in absolute terms (< 1% for most sites). This is partly a consequence of the eBird sampling framework: occupancy is estimated at the checklist-grid-cell level (~10 km), not at the territory level, so low absolute values are expected. What matters is the relative spatial pattern and whether it aligns with known distribution -- see validation (Section 5). We should also check whether the double-scaling issue (data already z-scored in GEE, then re-scaled in R) is inflating the intercept.

---

## 4. Response Curves

Marginal response curves for initial occupancy ($\psi_1$), holding other covariates at their mean value.

| *O. tarda* | *P. alchata* |
|:-:|:-:|
| ![otitar responses](figures/otitar_response_occupancy.png) | ![ptealc responses](figures/ptealc_response_occupancy.png) |

| *P. orientalis* | *T. tetrax* |
|:-:|:-:|
| ![pteori responses](figures/pteori_response_occupancy.png) | ![tettet responses](figures/tettet_response_occupancy.png) |

For `tettet`, we also obtained response curves for the dynamic sub-models (colonisation, extinction, detection), which failed for the other three species due to numerical issues in univariate `colext` refitting:

| Colonisation | Extinction | Detection |
|:-:|:-:|:-:|
| ![tettet colonisation](figures/tettet_response_colonization.png) | ![tettet extinction](figures/tettet_response_extinction.png) | ![tettet detection](figures/tettet_response_detection.png) |

---

## 5. Independent Validation (Spanish Biodiversity Atlas)

Models were validated against the Spanish Biodiversity Atlas (5,202 UTM 10x10 km polygons with presence/absence records) using spatially-blocked 5-fold cross-validation (`blockCV`). Predicted occupancy probabilities from the `colext` model were aggregated to atlas polygons by taking the median predicted value across all grid cells within each polygon.

### 5.1 Summary metrics

| Species | Atlas prev. | AUC | TSS | RMSE | Spearman $\rho$ | Threshold |
|---------|-------------|-----|-----|------|-----------------|-----------|
| *O. tarda* | 10.1% | **0.884** | **0.612** | 0.275 | 0.405 | 0.0028 |
| *P. alchata* | 8.4% | **0.844** | **0.569** | 0.261 | 0.334 | 0.0014 |
| *P. orientalis* | 16.2% | **0.823** | **0.536** | 0.366 | 0.410 | 0.0031 |
| *T. tetrax* | 24.8% | **0.850** | **0.546** | 0.460 | **0.521** | 0.014 |

All four species show good-to-excellent discriminative ability (AUC > 0.82) and positive, highly significant continuous calibration (Spearman $\rho$ all *P* < 0.001). The highest AUC (0.884) is for *O. tarda*, the most data-rich species; the highest calibration ($\rho$ = 0.521) is for *T. tetrax*, despite having the simplest dynamic sub-models.

> [!IMPORTANT]
> **For Raul:** These validation results are strong. AUC > 0.80 for all species, with TSS > 0.50, indicates that the eBird-derived dynamic occupancy models capture the broad-scale distribution patterns well. The Spearman correlations confirm that the predicted probabilities are not just discriminative but also reasonably well-calibrated in rank order. I'd highlight the *O. tarda* AUC = 0.884 and the *T. tetrax* Spearman = 0.521 as key results for the paper.

### 5.2 Validation maps

| *O. tarda* | *P. alchata* |
|:-:|:-:|
| ![otitar validation](figures/otitar_validation_maps.png) | ![ptealc validation](figures/ptealc_validation_maps.png) |

| *P. orientalis* | *T. tetrax* |
|:-:|:-:|
| ![pteori validation](figures/pteori_validation_maps.png) | ![tettet validation](figures/tettet_validation_maps.png) |

### 5.3 Calibration plots

| *O. tarda* | *P. alchata* |
|:-:|:-:|
| ![otitar calibration](figures/otitar_validation_calibration.png) | ![ptealc calibration](figures/ptealc_validation_calibration.png) |

| *P. orientalis* | *T. tetrax* |
|:-:|:-:|
| ![pteori calibration](figures/pteori_validation_calibration.png) | ![tettet calibration](figures/tettet_validation_calibration.png) |

### 5.4 Spatial residuals

Moran's I on model residuals indicates significant positive spatial autocorrelation for all species:

| Species | Moran's I | *P* |
|---------|-----------|-----|
| *O. tarda* | 0.358 | <0.001 |
| *P. alchata* | 0.491 | <0.001 |
| *P. orientalis* | 0.374 | <0.001 |
| *T. tetrax* | 0.423 | <0.001 |

This is expected for environment-only models that do not account for dispersal limitation, historical biogeographic processes, or spatial dependence in detectability. The spatial occupancy models (Section 6) address this explicitly.

| *O. tarda* | *P. alchata* |
|:-:|:-:|
| ![otitar residuals](figures/otitar_validation_residuals_map.png) | ![ptealc residuals](figures/ptealc_validation_residuals_map.png) |

| *P. orientalis* | *T. tetrax* |
|:-:|:-:|
| ![pteori residuals](figures/pteori_validation_residuals_map.png) | ![tettet residuals](figures/tettet_validation_residuals_map.png) |

---

## 6. Spatial Autocorrelation Assessment (`spOccupancy`)

To evaluate the importance of spatial autocorrelation and quantify the spatial range of dependence, we fitted two Bayesian occupancy models per species using `spOccupancy`:

1. **tPGOcc** -- Multi-season Polya-Gamma occupancy (temporal, non-spatial) across all 7 years
2. **spPGOcc** -- Single-season spatial occupancy with NNGP (nearest-neighbour Gaussian process) for the most recent year (2023)

Models used the same occupancy covariates as the `colext` initial occupancy formula plus NDVI and tmmx as time-varying covariates.

### 6.1 Model comparison

| Species | tPGOcc WAIC | spPGOcc WAIC | $\Delta$WAIC | Moran's I (tPGOcc) | Spatial range (km) |
|---------|-------------|--------------|------------|--------------------|--------------------|
| *O. tarda* | 2496.9 | 421.2 | **2075.7** | 0.145 | **182** |
| *P. alchata* | 2103.8 | 304.8 | **1799.0** | 0.236 | **334** |
| *P. orientalis* | 2273.7 | 445.1 | **1828.6** | 0.173 | **263** |
| *T. tetrax* | 1890.4 | 374.2 | **1516.2** | 0.144 | **23** |

The spatial models (spPGOcc) show dramatically lower WAIC across all species, confirming that spatial autocorrelation is a substantial component of the residual variation. All tPGOcc residuals exhibit significant Moran's I (*P* < 10$^{-79}$).

> [!IMPORTANT]
> **For Raul -- Key finding for the paper:** The spatial range estimates are ecologically informative:
> - *P. alchata* has the largest spatial range (**334 km**), consistent with its nomadic/semi-nomadic movements tracking ephemeral food resources across arid landscapes.
> - *O. tarda* and *P. orientalis* show intermediate ranges (182--263 km), reflecting regional-scale habitat connectivity.
> - *T. tetrax* has a strikingly short range (**23 km**), suggesting highly localised spatial dependence -- likely reflecting its strong site fidelity and limited dispersal in the breeding season.
>
> These spatial range parameters could form a central result in the paper: they link directly to species ecology (dispersal, site fidelity) and have implications for conservation planning (connectivity corridors, protected area design). However, we should note that these are from single-year (2023) snapshots and the MCMC convergence for spPGOcc was uneven (see Section 6.2).

### 6.2 MCMC diagnostics

| Species | tPGOcc max $\hat{R}$ | spPGOcc max $\hat{R}$ | spPGOcc min ESS | Notes |
|---------|----------------------|-----------------------|-----------------|-------|
| *O. tarda* | 1.19 ($\sigma^2_t$) | **3.95** ($\sigma^2$) | **10.8** ($\sigma^2$) | Spatial variance poorly mixed |
| *P. alchata* | 1.09 | **2.05** ($\sigma^2$) | **20** ($\sigma^2$) | Same issue |
| *P. orientalis* | 1.04 | **1.69** ($\sigma^2$) | **26** ($\sigma^2$) | Marginal |
| *T. tetrax* | 1.08 | **2.00** ($\phi$) | **24** ($\phi$) | Range parameter poorly mixed |

> [!CAUTION]
> **For Raul -- MCMC convergence:** The spPGOcc spatial variance ($\sigma^2$) and range ($\phi$) parameters show $\hat{R}$ > 1.5 and ESS < 50 for all species. This means the spatial range estimates in Section 6.1, while informative, should be treated with caution. For the paper, we should either: (a) increase the number of MCMC iterations (currently 3000 post-burnin) substantially (e.g., 25,000+); (b) use informative priors on the spatial range based on ecological knowledge; or (c) present these results as exploratory/supplementary while flagging the convergence issue. The occupancy coefficients ($\beta$) and detection coefficients ($\alpha$) converge well for all models.

### Spatial diagnostics

| *O. tarda* | *P. alchata* |
|:-:|:-:|
| ![otitar spatial](figures/otitar_spatial_diagnostics.png) | ![ptealc spatial](figures/ptealc_spatial_diagnostics.png) |

| *P. orientalis* | *T. tetrax* |
|:-:|:-:|
| ![pteori spatial](figures/pteori_spatial_diagnostics.png) | ![tettet spatial](figures/tettet_spatial_diagnostics.png) |

---

## 7. Population Trends (Stochastic Simulations)

We projected occupancy, colonisation, and extinction rates over the 7-year study period using the fitted `colext` models (500 simulations). Mean occupancy was computed by applying predicted $\gamma$ and $\varepsilon$ matrices to the estimated initial occupancy.

### 7.1 Summary

| Species | $\psi_{2017}$ | $\psi_{2023}$ | Trend | Mean $\gamma$ | Mean $\varepsilon$ |
|---------|----------------|----------------|-------|----------------|---------------------|
| *O. tarda* | ~0% | 0.67% | Increasing | 0.20% | 26.4% |
| *P. alchata* | 0.02% | 0.87% | Increasing | 0.24% | 46.6% |
| *P. orientalis* | 0.44% | 1.15% | Increasing | 0.23% | 43.0% |
| *T. tetrax* | ~0% | 0.52% | Increasing | 0.11% | 15.0% |

> [!NOTE]
> **For Raul:** All species show an apparent increasing occupancy trend over 2017--2023. However, this should be interpreted with caution:
> 1. **The initial occupancy ($\psi_{2017}$) estimates are near zero** for `otitar` and `tettet`. This likely reflects the low base-rate problem in eBird data rather than genuine near-absence -- these species were present in Spain throughout.
> 2. **The increasing trend may partly reflect increasing eBird coverage** over this period rather than genuine population expansion. The number of eBird checklists in Spain has grown substantially each year.
> 3. **Extinction rates are high** for the sandgrouse species (43--47%), which seems ecologically unrealistic for site-level occupancy. This may reflect high turnover at the checklist-grid level driven by detectability variation rather than true local extinction.
> 4. Despite these caveats, the *relative* differences between species are informative: *T. tetrax* has the lowest extinction rate (15%), suggesting the most stable site occupancy.

### 7.2 Prevalence plots

| *O. tarda* | *P. alchata* |
|:-:|:-:|
| ![otitar trend](figures/otitar_prevalence_over_time.png) | ![ptealc trend](figures/ptealc_prevalence_over_time.png) |

| *P. orientalis* | *T. tetrax* |
|:-:|:-:|
| ![pteori trend](figures/pteori_prevalence_over_time.png) | ![tettet trend](figures/tettet_prevalence_over_time.png) |

---

## 8. Technical Notes

### 8.1 Complete separation issues (model formula modifications)

Several species exhibited complete or quasi-complete separation in their colonisation and extinction sub-models, caused by zero-inflated land-cover variables that had no variance at occupied sites. This is a well-known problem in logistic regression with sparse binary data (Albert & Anderson, 1984; Heinze & Schemper, 2002).

**Species-specific fixes:**

| Species | Sub-model | Removed covariates | Reason |
|---------|-----------|--------------------|--------|
| `ptealc` | $\gamma$ | LC0, LC13 | LC0 = 0 at 98% of occupied sites; LC13: only 2/16 colonisation events had >0 |
| `ptealc` | $\varepsilon$ | LC0, tmmn | LC0: same as above; tmmn: $r$ = 0.81 with tmmx, only 18 extinction events |
| `pteori` | $\varepsilon$ | LC7, LC13 | LC13 = 0 in 43/44 transitions (estimate diverged to +109); LC7 near-zero at most sites (estimate = -5.4, SE = NaN) |
| `tettet` | $\gamma$ | pr | Complete separation -- any model with `pr` in $\gamma$ or $\varepsilon$ fails to converge |
| `tettet` | $\varepsilon$ | LC0, LC13, NDVI, tmmx | LC0 and LC13: zero variance; NDVI and tmmx: only 25 transitions total (9 ext, 16 persist) cannot support >2 parameters |

### 8.2 Goodness-of-fit

- **Parametric bootstrap** (`parboot`, nsim = 10): Completed for all 4 species but sample size too small for reliable *P*-values. Should be increased to nsim = 500+ for publication.
- **MacKenzie-Bailey GOF** (`mb.gof.test`): Failed for `otitar`, `ptealc`, and `pteori` due to a row-count mismatch (off-by-one, likely a bug in the `AICcmodavg` package with 7-year data). Succeeded for `tettet`.
- **Response curves:** Generated successfully only for `tettet` (simplest model). The other three species fail with "initial value in vmmin not finite" during univariate `colext` refitting -- this occurs because fixing some covariates at their mean while varying one at a time pushes the optimizer into a flat likelihood surface.

### 8.3 Known issues and recommendations

1. **Increase GOF simulations** for publication (parboot nsim = 500, MB GOF if row-count bug is resolved).
2. **Increase spPGOcc MCMC iterations** to >=25,000 post-burnin for reliable spatial range estimates.
3. **Consider Firth-type penalised likelihood** or Bayesian priors for the extinction sub-models of rare species to handle separation more principally.
4. **Investigate the double-scaling issue**: some covariates may have been pre-standardised in GEE and then re-standardised in R, which inflates the intercept while preserving slope ratios. This does not affect model selection or prediction but complicates direct interpretation of coefficient magnitudes.
5. **Response curves for publication** should be generated using predict() with a newdata grid rather than refitting univariate models.

---

## 9. File inventory

### Scripts (tracked in git)
| File | Description |
|------|-------------|
| `scripts/1_prepare_ebird_data.R` | eBird filtering with auk |
| `scripts/2_prepare_static_variables.R` | Static covariate extraction (WorldClim, MODIS, topography) |
| `scripts/3_prepare_dynamic_variables.R` | Dynamic covariate extraction (GEE exports) |
| `scripts/4_occupancy_models.R` | colext model fitting, GOF, maps, simulations |
| `scripts/5_validation.R` | Atlas validation with spatial block CV |
| `scripts/6_spatial_occupancy_test.R` | spOccupancy (tPGOcc + spPGOcc) |
| `R/model_configs.R` | Centralised model formulas |

### Key outputs (NOT tracked -- regenerated by pipeline)
| Directory | Contents | Approx. size |
|-----------|----------|-------------|
| `results/` | Model objects (.rds), summaries, simulations, MCMC diagnostics | ~2.8 GB |
| `figs/` | All generated figures (PNG) | ~9 MB |
| `data/processed_2023/` | Processed data per species | ~500 MB |

---

## 10. Next Steps

### 10.1 Spatial dynamic occupancy models (`stPGOcc`)

The single-year spatial models (spPGOcc, Section 6) demonstrate that spatial autocorrelation substantially improves model fit ($\Delta$WAIC > 1,500 for all species). The logical next step is to fit the **full spatio-temporal model** using `spOccupancy::stPGOcc`, which combines:

- **NNGP spatial random effects** (as in spPGOcc) to account for spatial dependence
- **AR(1) temporal autocorrelation** across the 7 study years (2017--2023)
- **Detection sub-model** with observation-level covariates

#### What `stPGOcc` does and does not do

| | `colext` (current) | `stPGOcc` (proposed) |
|---|---|---|
| Colonisation/extinction rates | Explicit $\gamma$ and $\varepsilon$ with covariates | Not estimated (no Markov transitions) |
| Spatial dependence | None (sites independent) | NNGP with estimated spatial range |
| Temporal dynamics | Mechanistic (Markov chain) | Phenomenological (AR(1) random effects) |
| Ecological interpretation | Demographic rates (colonisation, extinction) | Spatio-temporal occurrence patterns |

**These two approaches are complementary, not competing.** `colext` answers *"what environmental drivers affect colonisation and extinction?"* while `stPGOcc` answers *"how much residual spatial structure remains unexplained, and at what spatial scale does it operate?"*

#### Why it is worth running

1. The Moran's I values from tPGOcc residuals (0.14--0.24, all *P* < 10$^{-79}$) confirm genuine spatial autocorrelation
2. The spPGOcc spatial range estimates (23--334 km) are ecologically meaningful and could be a central result linking to dispersal ecology and conservation planning
3. Reviewers will likely question the independence assumption in `colext` -- having `stPGOcc` results pre-empts this

#### Recommended MCMC settings

The current spPGOcc settings (800 batches, 10,000 burn-in) yielded poor convergence for spatial parameters ($\hat{R}$ > 1.5, ESS < 50). For `stPGOcc` with 7 years of data:

```r
stPGOcc(
  n.batch    = 1600,      # 1600 × 25 = 40,000 total samples
  batch.length = 25,
  n.burn     = 20000,     # discard first 20,000
  n.thin     = 20,        # keep every 20th → 1,000 post-burnin per chain
  n.chains   = 3,
  NNGP       = TRUE,
  n.neighbors = 5,
  ar1        = TRUE       # AR(1) temporal autocorrelation
)
```

Expected runtime: **6--24 hours per species** (total ~1--4 days for all four).

Target convergence: $\hat{R}$ < 1.1 for $\beta$ and $\alpha$; $\hat{R}$ < 1.3 for $\sigma^2$ and $\phi$.

### 10.2 Goodness-of-fit

- Increase parametric bootstrap simulations from nsim = 10 (testing) to **nsim = 500+** for publication-quality *P*-values.
- The MacKenzie-Bailey GOF test (`mb.gof.test`) fails due to a row-count off-by-one error with 7-year data. Consider filing a bug report or using an alternative implementation.

### 10.3 Response curves

The current response curves are generated by refitting univariate `colext` models, which fails for 3/4 species due to flat likelihood surfaces. For publication, use `predict()` with a `newdata` grid instead:

```r
# Create prediction grid varying one covariate at a time
newdata <- data.frame(covariate = seq(min, max, length.out = 100))
# Fix all other covariates at their mean
pred <- predict(model, type = "state", newdata = newdata)
```

### 10.4 Double-scaling issue

Some covariates may have been pre-standardised in Google Earth Engine and then re-standardised in R (`scale()`). This does not affect model selection or spatial predictions but it compresses the logit-scale intercept and complicates direct interpretation of coefficient magnitudes. Verify whether GEE exports are already z-scored and, if so, skip the R-side `scale()` for those variables.

---

*Report generated 2026-03-05. Pipeline executed on branch `analysis-2023`.*
*For questions or discussion, contact Guillermo Fandos.*
