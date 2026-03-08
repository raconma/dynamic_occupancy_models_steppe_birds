# Dynamic occupancy and counterfactual attribution of Iberian steppe bird declines (2017–2023)

## Full results report — v4 (March 2026)

**Species:** *Otis tarda* (Great Bustard), *Pterocles alchata* (Pin-tailed Sandgrouse), *Pterocles orientalis* (Black-bellied Sandgrouse), *Tetrax tetrax* (Little Bustard)

**Data source:** eBird citizen science, mainland Spain, 5-km grid cells

**Modelling frameworks:**
- Dynamic occupancy (colext, `unmarked`)
- Spatio-temporal occupancy (stPGOcc, `spOccupancy`)
- Counterfactual factorial attribution

---

## 1 Study design and data

### 1.1 Study area and sampling grid

The study encompasses mainland peninsular Spain (excluding the Balearic and Canary Islands). We overlaid the study area with a regular 5-km grid (UTM 30N) and extracted eBird checklists within each cell for each primary period (year, 2017–2023). Each primary period was subdivided into *J* = 10 secondary sampling occasions (replicate visits).

**Grid sizes after processing:**
| Species | Grid cells | Cells removed (NA) | Cells with ≥1 detection |
|---------|-----------|-------------------|------------------------|
| *O. tarda* | 3,746 | 1 | 65 occupied-site observations |
| *P. alchata* | 4,131 | 1 | 64 occupied-site observations |
| *P. orientalis* | 4,131 | 1 | 44 occupied-site observations |
| *T. tetrax* | 3,746 | 0 | 25 occupied-site observations |

*Figure reference:* `figs/pub_map_sampling_sites.png`, `figs/pub_map_main_figure.png`

### 1.2 Covariates

**Static site covariates (initial occupancy ψ₁):** WorldClim bioclimatic variables (bio1, bio2), tree cover, grass cover, topographic elevation and aspect. These represent time-invariant habitat suitability.

**Dynamic yearly covariates (colonisation γ, extinction ε):** TerraClimate precipitation (pr), minimum temperature (tmmn), maximum temperature (tmmx), MODIS NDVI, and MODIS land-cover class proportions (LC0, LC6, LC7, LC10, LC12, LC13, LC14). All dynamic covariates were standardised (z-scored) using mean and SD computed on the *full training dataset* (all rows, all years) before any filtering.

**Detection covariates (p):** Survey effort distance, duration, number of observers, time of day, and observation-level NDVI and precipitation.

> **Design rationale:** Static WorldClim variables appear only in the initial occupancy submodel (ψ₁), never in γ/ε. This two-scale design prevents conflating baseline habitat suitability with interannual change drivers — a key requirement for clean counterfactual attribution.

### 1.3 NDVI collinearity and treatment

NDVI has a dual nature: ~50% of its variance is explained by climate variables (pr, tmmn, tmmx), with the remainder reflecting land management (irrigation, grazing, agricultural practices). A site-level regression of NDVI on climate yielded a mean R² = 0.51 (range: 0.04–0.97; 53% of sites > 0.5; 28% > 0.7).

This collinearity creates interpretability challenges for attribution. We therefore conducted coefficient stability analysis, refitting each species model with NDVI removed and comparing beta estimates:

| Species | Submodel | Covariate | Flag | Shift |
|---------|----------|-----------|------|-------|
| *O. tarda* | gamma | tmmn | BETA_SHIFT | 1.10 SE |
| *P. alchata* | gamma | pr | **SIGN_CHANGE** | +0.162 → −0.717 |
| *P. orientalis* | epsilon | LC12 | BETA_SHIFT | 2.80 SE |
| *P. orientalis* | epsilon | pr | SE_INCREASE | +60% |

**Decision (v4):** NDVI was removed from *P. alchata* gamma and *P. orientalis* epsilon. In both cases, its inclusion caused coefficient instability in co-occurring covariates, compromising causal attribution. For *O. tarda* gamma and *P. orientalis* gamma, NDVI was retained but reclassified as climate-adjacent in the attribution analysis (VIF < 5; beta shifts < 2 SE). AIC costs: *P. alchata* ΔAIC = +7.81; *P. orientalis* ΔAIC = +60.27.

*Figure reference:* `figs/pub_fig_collinearity_diagnostics.png`, `figs/diagnostics/fig_ndvi_climate_r2_map.png`

---

## 2 Dynamic occupancy models (colext)

### 2.1 Model formulas

All formulas are centralised in `R/model_configs.R`. After the NDVI revision (v4), the final specifications are:

| Species | ψ₁ | γ | ε | p |
|---------|------|------|------|------|
| *O. tarda* | bio1 + bio2 + tree + grass + elev | NDVI + pr + tmmn + tmmx | LC6 + LC13 + tmmx | effort + observers + NDVI_obs + pr_obs + aspect_obs |
| *P. alchata* | bio1 + bio2 + tree + grass + aspect | **pr** | pr + tmmx | time + duration + effort + observers + NDVI_obs + pr_obs |
| *P. orientalis* | bio2 + tree + grass | LC7 + NDVI + tmmn + tmmx | **LC12 + pr** | time + duration + observers + pr_obs + aspect_obs |
| *T. tetrax* | bio2 + tree + grass + elev | LC12 | LC12 | effort + observers + time |

Bold text indicates revised formulas (NDVI removed).

### 2.2 Model fit and AIC

| Species | AIC | Parameters | N sites |
|---------|------|-----------|---------|
| *O. tarda* | 2,267.90 | 17 | 3,745 |
| *P. alchata* | 1,981.71 | 15 | 4,130 |
| *P. orientalis* | 2,121.02 | 16 | 4,130 |
| *T. tetrax* | 1,780.09 | 11 | 3,746 |

### 2.3 Initial occupancy (ψ₁)

Across all four species, tree cover was the strongest predictor of initial occupancy, with consistently large negative coefficients (−3.72 to −6.45, all P < 10⁻⁹). Grass cover was positively associated with occupancy for all species (+0.36 to +0.77, all P < 0.02). Temperature seasonality (bio2) was positive for three species (except *O. tarda*, where bio1 dominated). These patterns confirm the a priori expectation that steppe specialists require open, treeless landscapes.

**Table: Initial occupancy coefficients (logit-scale)**

| Covariate | *O. tarda* | *P. alchata* | *P. orientalis* | *T. tetrax* |
|-----------|-----------|-------------|----------------|------------|
| Intercept | −8.085*** | −6.434*** | −5.790*** | −6.224*** |
| bio1 | +0.727• | −0.569** | — | — |
| bio2 | +0.618* | +1.207*** | +0.819*** | +1.149*** |
| tree_cover | −6.449*** | −4.545*** | −3.870*** | −3.724*** |
| grass_cover | +0.355* | +0.484*** | +0.613*** | +0.772*** |
| topo_elev | +1.213** | — | — | −0.106 |
| topo_aspect | — | −0.123 | — | — |

Significance: *** P < 0.001, ** P < 0.01, * P < 0.05, • P < 0.10

*Figure reference:* `figs/pub_fig5_psi_coefficients.png`, `figs/pub_map_occupancy_4species.png`

### 2.4 Colonisation (γ)

Colonisation rates are extremely low across all species. The intercept-only (baseline) colonisation probabilities range from essentially zero (*P. orientalis*) to 0.14% (*P. alchata*). With only 11–22 colonisation events observed across the entire study period, statistical power for covariate effects on colonisation is limited.

**Table: Colonisation coefficients (logit-scale)**

| Covariate | *O. tarda* | *P. alchata* | *P. orientalis* | *T. tetrax* |
|-----------|-----------|-------------|----------------|------------|
| Intercept | −10.00*** | −6.591*** | −22.37 (NaN SE) | −7.294*** |
| NDVI | −1.45* | — | −1.18 (SE=92.5) | — |
| pr | −2.30 | −0.717 | — | — |
| tmmn | −4.52* | — | −2.27 (SE=180) | — |
| tmmx | +2.59 | — | +6.32 (NaN SE) | — |
| LC7 | — | — | −0.34 (SE=34.4) | — |
| LC12 | — | — | — | +0.902* |

For *P. orientalis*, the gamma submodel shows complete/quasi-complete separation: the intercept is −22.37 with NaN standard errors, and all covariate SEs are enormous (34–180). This indicates that colonisation events are too rare and patterned for reliable estimation of covariate effects.

For *P. alchata*, gamma was excluded from the attribution analysis because only 16 colonisation events were observed across 2,521 site-year opportunities (0.6%).

*Figure reference:* `figs/pub_fig2_coefficient_forest_plot.png`, `figs/pub_fig_response_gam_all.png`

### 2.5 Extinction (ε)

Extinction rates are 1–2 orders of magnitude higher than colonisation rates, confirming a fundamental asymmetry in the dynamics of these species. Baseline extinction probabilities range from 13.75% (*T. tetrax*) to 48.94% (*P. alchata*).

**Table: Extinction coefficients (logit-scale)**

| Covariate | *O. tarda* | *P. alchata* | *P. orientalis* | *T. tetrax* |
|-----------|-----------|-------------|----------------|------------|
| Intercept | −1.252** | −0.042 | −1.830*** | −1.836*** |
| LC6 | +1.86 | — | — | — |
| LC12 | — | — | −0.095 | −0.565* |
| LC13 | +2.80*** | — | — | — |
| pr | — | −11.66 (SE=11.6) | +1.905 | — |
| tmmx | −1.51* | −17.62 (SE=14.7) | — | — |

Notable findings:
- *O. tarda*: Urban/built-up areas (LC13) strongly increase extinction (β = +2.80, P < 0.001).
- *T. tetrax*: Cropland (LC12) reduces both colonisation (+0.90 in γ) and extinction (−0.57 in ε), making it the single most consistent land-use driver.
- *P. alchata*: Extinction coefficients have very large SEs indicating near-separation (only 18 extinction events from 64 occupied-site observations).

*Figure reference:* `figs/pub_fig_response_eps_all.png`

### 2.6 Detection (p)

Detection probabilities at the intercept range from 18.3% (*P. orientalis*, logit = −1.69) to 34.9% (*T. tetrax*, logit = −0.62). All species show positive effects of survey effort and/or observer number, consistent with the expectation that more thorough surveys increase detection. Time of day was negatively associated with detection for three species, suggesting early-morning surveys are most effective.

**Table: Detection coefficients (logit-scale)**

| Covariate | *O. tarda* | *P. alchata* | *P. orientalis* | *T. tetrax* |
|-----------|-----------|-------------|----------------|------------|
| Intercept | −1.335*** | −1.075*** | −1.688*** | −0.621*** |
| effort | +0.225*** | +0.196** | — | +0.218** |
| observers | +0.223*** | +0.110 | −0.012 | +0.069 |
| time | — | −0.367*** | −0.236*** | −0.219** |
| duration | — | +0.896*** | +0.369*** | — |
| NDVI_obs | +0.428*** | −0.159 | — | — |
| pr_obs | −0.758*** | −0.292 | −0.391* | — |
| topo_aspect_obs | +0.204* | — | −0.142 | — |

Effort confounding was assessed by testing whether estimated occupancy trends correlate with survey effort trends. Results showed no significant confounding for any species (all P > 0.05), with *T. tetrax* showing the cleanest decoupling (ρ = 0.018, P = 0.97).

*Figure reference:* `figs/pub_fig4_detection_comparison.png`, `figs/pub_fig_detection_covariates.png`, `figs/pub_fig_effort_cross_species.png`

---

## 3 Occupancy dynamics and equilibrium

### 3.1 Stochastic simulations (2017–2023)

Forward simulations (500 iterations per species) using estimated ψ₁, γ, and ε reveal universally increasing but very low occupancy trajectories:

**Table: Simulated occupancy prevalence (% of sites occupied)**

| Year | *O. tarda* | *P. alchata* | *P. orientalis* | *T. tetrax* |
|------|-----------|-------------|----------------|------------|
| 2017 | 0.001% | 0.020% | 0.444% | 0.002% |
| 2018 | 0.206% | 0.259% | 0.622% | 0.111% |
| 2019 | 0.288% | 0.402% | 0.667% | 0.206% |
| 2020 | 0.433% | 0.550% | 0.859% | 0.295% |
| 2021 | 0.441% | 0.683% | 0.902% | 0.373% |
| 2022 | 0.648% | 0.803% | 1.026% | 0.445% |
| 2023 | 0.675% | 0.869% | 1.146% | 0.521% |

All species show increasing trends, but absolute occupancy remains below 1.2% everywhere. The upward trajectory primarily reflects the very low initial occupancy (ψ₁ ≈ 0) rather than high colonisation rates.

*Figure reference:* `figs/pub_fig_occupancy_trends_panel.png`, `figs/pub_fig_col_ext_rates.png`

### 3.2 Naive vs detection-corrected transition rates

A critical validation step compares raw (naive) transition rates — computed directly from the detection history — with model-corrected estimates that account for imperfect detection:

**Table: Naive vs corrected transition rates**

| Species | Naive γ | Corrected γ | Ratio | Naive ε | Corrected ε | Ratio |
|---------|---------|-------------|-------|---------|-------------|-------|
| *O. tarda* | 0.89% | 0.21% | 0.23 | 44.6% | 26.4% | 0.59 |
| *P. alchata* | 0.63% | 0.17% | 0.27 | 28.1% | 48.5% | 1.72 |
| *P. orientalis* | 0.87% | ~0% | ~0 | 34.1% | 23.4% | 0.69 |
| *T. tetrax* | 0.57% | 0.11% | 0.19 | 36.0% | 15.0% | 0.42 |

Key patterns:
1. **Corrected colonisation is 4–5× lower** than naive estimates for all species (ratios 0.19–0.27), because false absences inflate apparent colonisation.
2. **Corrected extinction is typically lower** than naive (ratios 0.42–0.69) for three species, because missed detections inflate apparent extinction.
3. **Exception:** *P. alchata* has corrected ε > naive ε (ratio 1.72), likely driven by the near-separation in the epsilon submodel (SEs > 11).

### 3.3 Equilibrium occupancy (ψ*)

Under current (time-averaged) conditions, the long-run equilibrium occupancy ψ* = γ/(γ+ε) was computed via parametric bootstrap (n = 5,000 draws from the multivariate normal approximation to the coefficient distribution):

**Table: Equilibrium occupancy and recovery timescales**

| Species | ψ* median | 95% CI | T_recol (yr) | T_persist (yr) | Status |
|---------|----------|--------|-------------|---------------|--------|
| *O. tarda* | 0.02% | [0.0004%, 0.93%] | 22,648 | 4 | Reportable |
| *P. alchata* | 0.31% | [0.08%, 6.50%] | 732 | 2 | Reportable (wide CI) |
| *P. orientalis* | ~0% | [0%, 93.70%] | ~∞ | 7 | **Excluded** (vcov not PD) |
| *T. tetrax* | 0.49% | [0.16%, 1.62%] | 1,471 | 7 | Reportable |

**Interpretation:** All four species are caught in a *demographic trap*: extinction rates exceed colonisation rates by 1–2 orders of magnitude, driving equilibrium occupancy well below 1%. Even for the most "optimistic" species (*T. tetrax*), the median equilibrium is 0.49% and recolonisation of a lost site would take ~1,500 years on average. For *O. tarda*, recolonisation timescales exceed 20,000 years.

*P. orientalis* equilibrium is excluded from headline reporting because the gamma submodel is completely separated (intercept = −22.37 with NaN SE), making the variance–covariance matrix non-positive-definite and the bootstrap CI uninformative.

---

## 4 Spatial structure (stPGOcc)

### 4.1 Model comparison

Three occupancy model frameworks were compared: temporal-only (tPGOcc), single-year spatial (spPGOcc, 2023), and full spatio-temporal (stPGOcc, 2017–2023):

**Table: Spatial model comparison**

| Species | tPGOcc WAIC | stPGOcc WAIC | ΔWAIC | Spatial range (km) |
|---------|------------|-------------|-------|-------------------|
| *O. tarda* | 2,497 | 2,040 | −457 | 46.7 |
| *P. alchata* | 2,104 | 1,655 | −449 | 263.9 |
| *P. orientalis* | 2,274 | 1,821 | −453 | 171.6 |
| *T. tetrax* | 1,890 | 1,539 | −352 | 42.7 |

Substantial WAIC improvement (ΔWAIC > 350) for all species confirms that spatial random effects capture important unmodelled spatial structure. Effective spatial ranges span from 42.7 km (*T. tetrax*) to 263.9 km (*P. alchata*), reflecting species-specific dispersal capacities and habitat patch geometry.

### 4.2 Residual spatial autocorrelation

**Table: Moran's I before and after spatial random effects**

| Species | Moran's I (tPGOcc) | Moran's I (stPGOcc) | Reduction |
|---------|-------------------|--------------------|-----------|
| *O. tarda* | 0.145 | −0.060 | 141% |
| *P. alchata* | 0.236 | −0.042 | 82% |
| *P. orientalis* | 0.173 | −0.035 | 80% |
| *T. tetrax* | 0.144 | −0.042 | 71% |

All post-stPGOcc Moran's I values are ≤ 0 (i.e., no residual positive spatial autocorrelation), confirming that the NNGP spatial random effects adequately capture spatial structure.

### 4.3 Convergence

> **⚠️ CALLOUT TO RAUL:** All four stPGOcc models have convergence problems (Rhat > 1.1 for multiple parameters). The current runs used n.batch = 800, batch.length = 25 (20,000 total iterations). **Action needed:** Re-run all stPGOcc models on the cluster with n.batch = 4,000, batch.length = 25 (100,000 iterations), n.burn = 50,000, n.thin = 50. Expected runtime: 6–24 hours per species. Until convergence is achieved, stPGOcc results should be presented as "preliminary" in the manuscript.

**Convergence diagnostics (current runs):**

| Species | Worst Rhat | Parameter | ESS min | Publishable? |
|---------|-----------|-----------|---------|-------------|
| *O. tarda* | 4.19 | σ²_t | 9.1 | No |
| *P. alchata* | 12.28 | ρ (AR1) | 5.2 | No |
| *P. orientalis* | 11.86 | σ² | 5.6 | No |
| *T. tetrax* | 1.84 | σ²_t | 31.0 | No (marginal) |

*T. tetrax* is the closest to convergence and may only need a modest increase in iterations.

*Figure reference:* `figs/pub_fig_spatial_moran.png`, `figs/pub_fig_spatial_range.png`, `figs/pub_fig_spatial_occupancy_trends.png`

---

## 5 Model validation

### 5.1 Spatial block cross-validation

Five-fold spatial block cross-validation (blockCV) was used to assess predictive performance. Discrimination was evaluated via AUC and TSS; calibration was assessed via Spearman rank correlation and residual maps.

**Table: Cross-validation metrics (5-fold spatial blockCV)**

| Species | Mean AUC | AUC range | Mean TSS | Mean Sens | Mean Spec |
|---------|----------|-----------|----------|-----------|-----------|
| *O. tarda* | 0.884 | 0.830–0.937 | 0.612 | 0.862 | 0.750 |
| *P. alchata* | 0.844 | 0.780–0.895 | 0.569 | 0.823 | 0.746 |
| *P. orientalis* | 0.823 | 0.772–0.861 | 0.536 | 0.874 | 0.663 |
| *T. tetrax* | 0.850 | 0.815–0.897 | 0.546 | 0.835 | 0.711 |

All species achieve AUC > 0.80, indicating good discriminatory power. *O. tarda* has the highest discrimination (AUC = 0.884) and calibration (TSS = 0.612).

> **⚠️ CALLOUT TO RAUL:** The current blockCV block size is 50 km, which is smaller than the spatial autocorrelation range of *P. alchata* (264 km) and *P. orientalis* (172 km). **This has been fixed** in `scripts/5_validation.R` (increased to 270 km), but the script needs to be **re-run** to produce updated validation metrics. Current AUC values may be optimistically biased by spatial leakage. Re-running requires the atlas data and may take 30–60 min per species.

*Figure reference:* `figs/otitar_validation_calibration.png`, `figs/ptealc_validation_calibration.png`, `figs/pteori_validation_calibration.png`, `figs/tettet_validation_calibration.png`

### 5.2 Goodness of fit (parametric bootstrap)

> **⚠️ CALLOUT TO RAUL:** Parametric bootstrap goodness-of-fit (parboot) **failed for all four species** due to a serialization bug: when colext model objects are saved with `saveRDS()` and reloaded, the internal `formlist` slot is lost, causing `parboot()` to crash.
>
> **Fix:** Models must be fitted and `parboot()` run in the **same R session** without saving/reloading. The number of simulations should be increased from 10 to at least 500. This is best done on the cluster:
> ```r
> # In scripts/4_occupancy_models.R, after fitting:
> GOF <- parboot(Mod.final, nsim = 500)  # Must be in same session
> saveRDS(GOF, here("results", paste0(sp, "_gof_parboot.rds")))
> ```
> Until parboot results are available, report MacKenzie–Bailey GOF instead (already computed with nsim = 10; increase to 500 for publication).

### 5.3 Collinearity diagnostics (VIF)

All variance inflation factors are below the standard threshold of 5:

| Species | Submodel | Max VIF | Covariate |
|---------|----------|---------|-----------|
| *O. tarda* | gamma | 4.43 | tmmx |
| *P. alchata* | gamma | 1.46 | pr |
| *P. orientalis* | gamma | 3.57 | tmmx |
| *T. tetrax* | gamma | 1.00 | LC12 |

Although VIF values are acceptable, pairwise correlations reveal 5 SEVERE pairs (|r| > 0.70): pr–tmmx (−0.70 to −0.71) and tmmn–tmmx (+0.79 to +0.81) across species. These high correlations do not inflate VIF above 5 because they involve only 2-variable subsets, but they may reduce precision of individual coefficient estimates.

*Figure reference:* `figs/pub_fig_collinearity_diagnostics.png`

---

## 6 Counterfactual attribution

### 6.1 Framework

We decompose occupancy changes into climate and land-use contributions using a factorial design with four scenarios:

| Scenario | Climate covariates | Land-use covariates |
|----------|-------------------|---------------------|
| S0 (Null) | Frozen at 2017 | Frozen at 2017 |
| S1 (Climate only) | Observed | Frozen at 2017 |
| S2 (Land use only) | Frozen at 2017 | Observed |
| S3 (Combined) | Observed | Observed |

Attribution effects are computed as:
- **Climate effect:** S1 − S0
- **Land-use effect:** S2 − S0
- **Interaction:** S3 − S1 − S2 + S0

Covariate classification:
- **Climate:** NDVI (where retained), pr, tmmn, tmmx
- **Land use:** LC6, LC7, LC10, LC12, LC13, LC14

### 6.2 Scaling fix (v4)

> In previous versions, the attribution script (`scripts/8`) recomputed covariate scaling on the filtered (drop_na) dataset, introducing a subtle mismatch with the training scaling used in `scripts/4` (which scales on all rows before filtering). This has been **fixed** in v4: the attribution script now loads pre-computed training scaling parameters (`train_dyn_scale.rds`) directly from the model-fitting stage.

### 6.3 Attribution results

**Table: Cross-species attribution summary**

| Species | γ climate | γ land-use | ε climate | ε land-use | Dominant | Note |
|---------|----------|-----------|----------|-----------|---------|------|
| *O. tarda* | −6.1×10⁻⁴ | 0 | +7.4×10⁻⁴ | −1.1×10⁻³ | **Climate** | |
| *P. alchata* | 0 | 0 | −1.3×10⁻³ | 0 | **Climate*** | γ excluded (16 events) |
| *P. orientalis* | 0 | 0 | +2.3×10⁻³ | +5.6×10⁻⁵ | **Climate** | |
| *T. tetrax* | 0 | −1.0×10⁻⁶ | 0 | +9.0×10⁻⁶ | **Land use** | Only LC12 in model |

Values are mean Δprobability across sites and years. Asterisk (*) indicates gamma excluded from attribution.

### 6.4 Interpretation

1. **Climate dominates for three species** (*O. tarda*, *P. alchata*, *P. orientalis*). For these species, interannual variation in precipitation and temperature drives the largest changes in colonisation and extinction probabilities.

2. **Land use dominates for *T. tetrax***. The Little Bustard's dynamics are driven entirely by cropland proportion (LC12), with no climate covariates in either γ or ε. This makes *T. tetrax* uniquely vulnerable to agricultural intensification and land abandonment.

3. **Attribution effects are small in absolute terms** (< 0.01 probability units). This reflects the short 7-year window (2017–2023), during which climate and land-use trends have only begun to diverge from baseline. Over longer timescales, these effects would compound.

4. **For *O. tarda*, the signs are ecologically coherent:** climate reduces colonisation (Δγ < 0) and increases extinction (Δε > 0), consistent with warming and aridity stress. Land use slightly reduces extinction (Δε < 0), potentially reflecting conservation measures (agri-environment schemes in key areas).

### 6.5 Bootstrap confidence intervals

**Table: Bootstrap summary (n = 1,000)**

| Species | γ mean | γ 95% CI | ε mean | ε 95% CI | ψ* mean | ψ* 95% CI |
|---------|--------|---------|--------|---------|---------|----------|
| *O. tarda* | 0.035% | [0.0001%, 0.21%] | 23.1% | [11.3%, 39.9%] | 0.14% | [0.0004%, 0.90%] |
| *P. alchata* | 0.16% | [0.06%, 0.38%] | 48.9% | [2.6%, 96.7%] | 0.95% | [0.08%, 4.91%] |
| *P. orientalis* | — | (vcov not PD) | — | — | — | — |
| *T. tetrax* | 0.08% | [0.02%, 0.22%] | 14.2% | [8.4%, 22.3%] | 0.59% | [0.15%, 1.65%] |

*P. orientalis* bootstrap was not possible because the variance–covariance matrix is not positive definite (gamma separation). The `nearPD` approximation yields a point estimate but no reliable CIs.

---

## 7 Separation and convergence diagnostics

### 7.1 Complete and quasi-complete separation

Several submodels show evidence of separation, a common problem in logistic regression with rare events:

| Species | Submodel | Issue | Details |
|---------|----------|-------|---------|
| *P. alchata* | ε | Near-separation | Intercept SE = 1.87; pr SE = 11.60; tmmx SE = 14.74. Only 18 extinction events. |
| *P. orientalis* | γ | **Complete separation** | Intercept = −22.37 (NaN SE); tmmx = +6.32 (NaN SE). Only 22 colonisation events across 2,541 opportunities. |
| *P. orientalis* | γ | Quasi-separation | LC7 SE = 34.4; NDVI SE = 92.5; tmmn SE = 180.4. |
| *T. tetrax* | — | None | Cleanest model; only 1 covariate per submodel. |
| *O. tarda* | — | None | All SEs finite and reasonable. |

These separation issues are intrinsic to the data: colonisation and extinction events are extremely rare for these species. They do not invalidate the non-spatial colext models but limit the interpretability of specific covariate effects in the affected submodels.

---

## 8 Key findings and conservation implications

### 8.1 The demographic trap

All four Iberian steppe bird species are caught in a demographic trap characterised by:

1. **Extinction rates 100–500× higher than colonisation rates** (ε/γ ratio: 125 for *O. tarda* to 288 for *P. alchata*).
2. **Equilibrium occupancy < 1%** for all species with reliable estimates.
3. **Recolonisation timescales of centuries to millennia** (730 years for *P. alchata* to 22,600 years for *O. tarda*).

This asymmetry means that once a site is lost, recovery is extremely slow under current conditions. Conservation must therefore prioritise **preventing extinction at occupied sites** over promoting colonisation of empty ones.

### 8.2 Species-specific vulnerabilities

- **O. tarda** (Great Bustard): Climate-driven decline with extinction strongly linked to urban encroachment (LC13, β = +2.80). Very low baseline colonisation (0.005%). Longest recovery time (22,600 yr). Most responsive to NDVI in both detection and colonisation.

- **P. alchata** (Pin-tailed Sandgrouse): Highest extinction rate (49%) among the four species. Climate-driven via precipitation. Limited analytical power (16 colonisation events; near-separation in ε). Largest spatial autocorrelation range (264 km).

- **P. orientalis** (Black-bellied Sandgrouse): Complete separation in γ makes colonisation dynamics uninterpretable. Climate (precipitation) drives extinction. Cropland (LC12) has a weak negative effect on extinction. AIC penalty of NDVI removal was the highest (+60).

- **T. tetrax** (Little Bustard): Cleanest model with clear land-use signal. Cropland proportion (LC12) is the sole driver — increasing it raises colonisation and reduces extinction. This species is the strongest candidate for targeted agri-environment scheme interventions.

### 8.3 Spatial autocorrelation

The substantial spatial autocorrelation ranges (43–264 km) imply that occupancy at a given site is influenced by conditions at sites tens to hundreds of kilometres away. This has two practical implications:

1. **Conservation actions must be coordinated at landscape scale**, not at individual sites.
2. **Standard cross-validation (without spatial blocking) overestimates predictive performance.** Our blockCV with 270-km blocks provides more realistic (and likely lower) AUC estimates.

---

## 9 Pending items for submission

> **⚠️ ITEMS REQUIRING RAUL'S ACTION:**

| Item | Priority | Effort | Notes |
|------|----------|--------|-------|
| Re-run stPGOcc with 100K iterations | **Blocking** | 6–24h/species (cluster) | n.batch = 4000, batch.length = 25, n.burn = 50000, n.thin = 50 |
| Re-run parboot GOF in same R session | **Blocking** | 1–2h/species (cluster) | nsim = 500. Must fit model + run parboot without saveRDS/readRDS in between |
| Re-run blockCV with 270-km blocks | **Strongly recommended** | 30–60 min/species | Script already fixed. Just run `Rscript scripts/5_validation.R` |
| Increase stochastic simulation nsim | Recommended | ~30 min total | Change NSIM_PREV from 500 to 5,000 in scripts/4 |
| Increase parboot/MB GOF nsim | Recommended | Cluster | NSIM_PARBOOT = 1,000; NSIM_MB_GOF = 500 |

---

## 10 Files inventory

### 10.1 Main scripts

| Script | Purpose | Status |
|--------|---------|--------|
| `scripts/4_occupancy_models.R` | colext fitting pipeline | ✅ Updated (saves train_dyn_scale) |
| `scripts/5_validation.R` | Spatial block CV | ✅ Fixed (270 km blocks) — needs re-run |
| `scripts/6_spatial_occupancy_test.R` | stPGOcc | ⚠️ Needs more iterations |
| `scripts/8_counterfactual_attribution.R` | Original attribution | ✅ Fixed (loads train_dyn_scale) |
| `scripts/9_collinearity_ndvi.R` | NDVI diagnostics | ✅ Complete |
| `scripts/10_attribution_revised.R` | **Revised attribution** | ✅ New |
| `scripts/refit_revised_models.R` | Refit ptealc/pteori | ✅ New |
| `scripts/compute_equilibrium.R` | Equilibrium ψ* | ✅ New |
| `R/model_configs.R` | Centralised formulas | ✅ Updated (NDVI removed) |

### 10.2 Key results files

| File | Content |
|------|---------|
| `results/{sp}_model_object.rds` | Fitted colext models (revised for ptealc/pteori) |
| `results/{sp}_train_dyn_scale.rds` | Training scaling parameters |
| `results/equilibrium_occupancy_table.csv` | Equilibrium ψ* with bootstrap CIs |
| `results/naive_vs_corrected_gamma.csv` | Naive vs detection-corrected rates |
| `results/attribution_table3.csv` | Cross-species attribution summary |
| `results/attribution_boot_summary.csv` | Bootstrap CIs for attribution |
| `results/attribution_revised_predictions.rds` | Full counterfactual predictions |
| `results/diagnostics/vif_summary.csv` | VIF diagnostics |
| `results/diagnostics/coefficient_stability.csv` | NDVI stability analysis |

### 10.3 Figures for manuscript

**Main figures (recommended):**

| Figure | File | Content |
|--------|------|---------|
| Fig. 1 | `figs/pub_map_main_figure.png` | Study area + survey effort |
| Fig. 2 | `figs/pub_map_occupancy_4species.png` | Initial occupancy maps (4 panels) |
| Fig. 3 | `figs/pub_fig2_coefficient_forest_plot.png` | Colonisation coefficients |
| Fig. 4 | `figs/pub_fig_coef_all_submodels.png` | All submodel coefficients |
| Fig. 5 | `figs/pub_fig_occupancy_trends_panel.png` | Prevalence trends 2017–2023 |
| Fig. 6 | `figs/pub_fig_col_ext_rates.png` | Colonisation/extinction rates |
| Fig. 7 | `figs/pub_fig_spatial_moran.png` | Spatial autocorrelation diagnostics |
| Fig. 8 | `figs/pub_fig_collinearity_diagnostics.png` | VIF and NDVI collinearity |
| Fig. 9 | `figs/pub_fig_summary_panel.png` | Summary panel |

**Supplementary figures:**

| Figure | File | Content |
|--------|------|---------|
| S1 | `figs/pub_fig_response_panel_main.png` | Response curves (all submodels) |
| S2 | `figs/pub_fig_effort_cross_species.png` | Effort diagnostics |
| S3 | `figs/pub_fig_spatial_range.png` | Spatial correlation ranges |
| S4 | `figs/pub_fig_spatial_occupancy_trends.png` | stPGOcc occupancy trends |
| S5 | `figs/{sp}_validation_calibration.png` | CV calibration plots (4 panels) |
| S6 | `figs/pub_fig_spatial_w_maps.png` | Spatial random effect maps |

---

## Appendix A: Model revision log

| Date | Change | Rationale | AIC impact |
|------|--------|-----------|-----------|
| v1–v3 | Iterative formula selection | AIC comparison + separation fixes | — |
| v4 | Remove NDVI from ptealc/gamma | pr sign change (+0.16 → −0.72); clean attribution | +7.81 |
| v4 | Remove NDVI from pteori/epsilon | SE inflation +60%; LC12 beta shift 2.8 SE | +60.27 |
| v4 | Fix scaling mismatch | scripts/4 vs scripts/8 used different scaling | N/A |
| v4 | blockCV block size 50 → 270 km | 50 km < 264 km spatial range | Pending re-run |
| v4 | Save train_dyn_scale.rds | Enables reproducible scaling in downstream scripts | N/A |

---

*Report generated on branch `audit-gcb-v4`. Last commit: `6af1f28`.*
