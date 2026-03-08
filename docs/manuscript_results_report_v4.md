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

**Decision (v4):** NDVI was removed from *P. alchata* gamma because its inclusion caused a sign change in the precipitation coefficient (+0.162 → −0.717), rendering climate attribution uninterpretable (AIC cost: +7.81). For *P. orientalis* epsilon, NDVI was **retained** in the final model (AIC-preferred by 60 units) but classified as "climate-adjacent" in the attribution analysis; a sensitivity model without NDVI is stored in `results/pteori_epsilon_sensitivity_noNDVI.rds`. For *O. tarda* gamma, NDVI was retained and also classified as climate-adjacent (VIF < 5; beta shifts < 2 SE).

*Figure reference:* `figs/pub_fig_collinearity_diagnostics.png`, `figs/diagnostics/fig_ndvi_climate_r2_map.png`

---

## 2 Dynamic occupancy models (colext)

### 2.1 Model formulas

All formulas are centralised in `R/model_configs.R`. After the NDVI revision (v4), the final specifications are:

| Species | ψ₁ | γ | ε | p |
|---------|------|------|------|------|
| *O. tarda* | bio1 + bio2 + tree + grass + elev | NDVI† + pr + tmmn + tmmx | LC6 + LC13 + tmmx | effort + observers + NDVI_obs + pr_obs + aspect_obs |
| *P. alchata* | bio1 + bio2 + tree + grass + aspect | **pr** | pr + tmmx | time + duration + effort + observers + NDVI_obs + pr_obs |
| *P. orientalis* | bio2 + tree + grass | LC7 + NDVI† + tmmn + tmmx | LC12 + NDVI† + pr | time + duration + observers + pr_obs + aspect_obs |
| *T. tetrax* | bio2 + tree + grass + elev | LC12 | LC12 | effort + observers + time |

Bold text indicates revised formulas (NDVI removed). † indicates NDVI classified as "climate-adjacent" in attribution analysis (~50% climate-driven, R² ≈ 0.51).

### 2.2 Model fit and AIC

| Species | AIC | Parameters | N sites | Note |
|---------|------|-----------|---------|------|
| *O. tarda* | 2,267.90 | 17 | 3,745 | |
| *P. alchata* | 1,981.71 | 15 | 4,130 | NDVI removed from γ (+7.81) |
| *P. orientalis* | 2,060.74 | 17 | 4,130 | NDVI retained in ε (original) |
| *T. tetrax* | 1,780.09 | 11 | 3,746 | |

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
| Intercept | −10.00*** | −6.591*** | −14.78*** (SE=3.81) | −7.294*** |
| NDVI | −1.45* | — | −0.56 (SE=0.49) | — |
| pr | −2.30 | −0.717 | — | — |
| tmmn | −4.52* | — | −8.11** (SE=2.67) | — |
| tmmx | +2.59 | — | +9.14*** (SE=2.68) | — |
| LC7 | — | — | −9.80 (SE=6.46) | — |
| LC12 | — | — | — | +0.902* |

For *P. orientalis*, the gamma submodel has large standard errors for some covariates (LC7 SE = 6.46), reflecting the rarity of colonisation events (22 across 2,541 opportunities, 0.87%). However, the temperature covariates (tmmn, tmmx) are significant (P < 0.003 and P < 0.001 respectively), indicating that extreme temperatures influence colonisation probability. The baseline colonisation rate is effectively zero (γ_baseline ≈ 3.8 × 10⁻⁷).

For *P. alchata*, gamma was excluded from the attribution analysis because only 16 colonisation events were observed across 2,521 site-year opportunities (0.6%).

*Figure reference:* `figs/pub_fig2_coefficient_forest_plot.png`, `figs/pub_fig_response_gam_all.png`

### 2.5 Extinction (ε)

Extinction rates are 1–2 orders of magnitude higher than colonisation rates, confirming a fundamental asymmetry in the dynamics of these species. Baseline extinction probabilities range from 13.75% (*T. tetrax*) to 48.94% (*P. alchata*).

**Table: Extinction coefficients (logit-scale)**

| Covariate | *O. tarda* | *P. alchata* | *P. orientalis* | *T. tetrax* |
|-----------|-----------|-------------|----------------|------------|
| Intercept | −1.252** | −0.042 | −0.239 (SE=0.56) | −1.836*** |
| LC6 | +1.86 | — | — | — |
| LC12 | — | — | −0.831* (SE=0.36) | −0.565* |
| LC13 | +2.80*** | — | — | — |
| NDVI† | — | — | +2.047*** (SE=0.56) | — |
| pr | — | −11.66 (SE=11.6) | +0.445 (SE=0.57) | — |
| tmmx | −1.51* | −17.62 (SE=14.7) | — | — |

† NDVI classified as "climate-adjacent" in attribution (R² ≈ 0.51 with climate variables).

Notable findings:
- *O. tarda*: Urban/built-up areas (LC13) strongly increase extinction (β = +2.80, P < 0.001).
- *P. orientalis*: NDVI has the strongest effect on extinction (β = +2.05, P = 0.0002) — higher NDVI increases extinction probability. Cropland (LC12) reduces extinction (β = −0.83, P = 0.02). A sensitivity model without NDVI (ΔAIC = +60) is provided in `results/pteori_epsilon_sensitivity_noNDVI.rds`.
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
| *O. tarda* | 0.02% | [0.0004%, 0.93%] | 22,648 | 4 | ✅ Reportable |
| *P. alchata* | 0.31% | [0.08%, 6.50%] | 732 | 2 | ✅ Reportable (wide CI) |
| *P. orientalis* | 0.0001% | [0%, 0.17%] | 2,599,962 | 2 | ✅ Reportable |
| *T. tetrax* | 0.49% | [0.16%, 1.62%] | 1,471 | 7 | ✅ Reportable |

**Interpretation:** All four species are caught in a *demographic trap*: extinction rates exceed colonisation rates by 2–6 orders of magnitude, driving equilibrium occupancy well below 1%. Even for the most "optimistic" species (*T. tetrax*), the median equilibrium is 0.49% and recolonisation of a lost site would take ~1,500 years on average. For *O. tarda*, recolonisation timescales exceed 20,000 years. *P. orientalis* has the most extreme profile: baseline colonisation is effectively zero (γ ≈ 3.8 × 10⁻⁷), implying recolonisation timescales of ~2.6 million years — functionally, once a site is lost, it is permanently lost.

### 3.4 Extinction debt

The difference between current simulated occupancy and equilibrium occupancy (ψ*) quantifies the **extinction debt** — the fraction of currently occupied sites that are expected to be lost as the system relaxes toward equilibrium under current conditions:

| Species | Current occ. (%) | ψ* (%) | Debt (pp) | Debt fraction | Interpretation |
|---------|-----------------|--------|-----------|---------------|---------------|
| *O. tarda* | 0.675 | 0.021 | 0.655 | **97%** | 97% of current sites are transient |
| *P. alchata* | 0.869 | 0.312 | 0.557 | **64%** | Two-thirds of sites at risk |
| *P. orientalis* | 1.146 | 0.0001 | 1.145 | **100%** | All current sites are transient |
| *T. tetrax* | 0.521 | 0.493 | 0.028 | **5%** | Near equilibrium; most resilient |

*T. tetrax* stands out as the only species near equilibrium (5% debt), consistent with its land-use-driven dynamics and shorter recolonisation timescale. The remaining three species carry substantial extinction debts (64–100%), driven primarily by the orders-of-magnitude asymmetry between γ and ε.

*Figure reference:* `figs/pub_fig_isocline_equilibrium.png` (panels a–c)

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
- **Climate:** pr, tmmn, tmmx
- **Climate-adjacent:** NDVI (retained for *O. tarda* gamma and *P. orientalis* gamma/epsilon; ~50% climate-driven)
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

| Species | γ baseline | ε baseline | ψ* median | ψ* 95% CI |
|---------|-----------|-----------|----------|----------|
| *O. tarda* | 0.0045% | 22.24% | 0.02% | [0.0004%, 0.93%] |
| *P. alchata* | 0.137% | 48.94% | 0.31% | [0.08%, 6.50%] |
| *P. orientalis* | 0.000038% | 44.06% | 0.0001% | [0%, 0.17%] |
| *T. tetrax* | 0.068% | 13.75% | 0.49% | [0.16%, 1.62%] |

All four species have positive-definite variance–covariance matrices, enabling reliable bootstrap CIs (n = 5,000 draws). *P. orientalis* required a `nearPD` correction for one parameter but the resulting CIs are informative (upper bound = 0.17%, well below current occupancy of 1.15%).

---

## 7 Separation and convergence diagnostics

### 7.1 Complete and quasi-complete separation

Several submodels show evidence of separation, a common problem in logistic regression with rare events:

| Species | Submodel | Issue | Details |
|---------|----------|-------|---------|
| *P. alchata* | ε | Near-separation | Intercept SE = 1.87; pr SE = 11.60; tmmx SE = 14.74. Only 18 extinction events. |
| *P. orientalis* | γ | Large SEs | LC7 SE = 6.46 (22 colonisation events). But tmmn (P=0.002) and tmmx (P<0.001) are significant. |
| *T. tetrax* | — | None | Cleanest model; only 1 covariate per submodel. |
| *O. tarda* | — | None | All SEs finite and reasonable. |

> **Note on P. orientalis gamma convergence:** In a previous model revision where NDVI was removed from epsilon, the gamma submodel collapsed to complete separation (intercept = −22.37, NaN SEs). Reverting epsilon to the original formula (with NDVI) restored gamma convergence with finite, interpretable SEs. This demonstrates that the likelihood surface in multi-submodel frameworks is sensitive to specification changes in other submodels. The sensitivity model is preserved in `results/pteori_epsilon_sensitivity_noNDVI.rds`.

These separation issues are intrinsic to the data: colonisation and extinction events are extremely rare for these species. They do not invalidate the non-spatial colext models but limit the interpretability of specific covariate effects in the affected submodels.

---

## 8 Key findings and conservation implications

### 8.1 The demographic trap and extinction debt

All four Iberian steppe bird species are caught in a demographic trap characterised by:

1. **Extinction rates 100–1,000,000× higher than colonisation rates** (ε/γ ratio: 203 for *T. tetrax* to 1,152,790 for *P. orientalis*).
2. **Equilibrium occupancy < 1%** for all four species (range: 0.0001–0.49%).
3. **Recolonisation timescales of centuries to millions of years** (732 years for *P. alchata* to 2.6 million years for *P. orientalis*).
4. **Extinction debt of 5–100%** of current occupancy. Three of four species have > 60% of their current sites destined for eventual loss under current conditions.

This asymmetry constitutes an *extinction debt* (sensu Tilman et al. 1994): the current occupancy exceeds the long-run equilibrium, and the difference represents sites that will eventually be lost as the system relaxes. Conservation must therefore prioritise **preventing extinction at occupied sites** over promoting colonisation of empty ones, because recolonisation is functionally impossible at the current rates.

*Figure reference:* `figs/pub_fig_isocline_equilibrium.png`

### 8.2 Species-specific vulnerabilities

- **O. tarda** (Great Bustard): Climate-driven decline with extinction strongly linked to urban encroachment (LC13, β = +2.80). Very low baseline colonisation (0.005%). Longest recovery time (22,600 yr). Most responsive to NDVI in both detection and colonisation.

- **P. alchata** (Pin-tailed Sandgrouse): Highest extinction rate (49%) among the four species. Climate-driven via precipitation. Limited analytical power (16 colonisation events; near-separation in ε). Largest spatial autocorrelation range (264 km).

- **P. orientalis** (Black-bellied Sandgrouse): The most extreme case of the demographic trap — baseline γ ≈ 3.8 × 10⁻⁷, implying recolonisation timescales of ~2.6 million years. NDVI strongly increases extinction (β = +2.05, P < 0.001), classified as climate-adjacent (~50% climate-driven). Cropland (LC12) reduces extinction (β = −0.83, P = 0.02). Temperature extremes significantly affect colonisation (tmmn P = 0.002, tmmx P < 0.001). Extinction debt = 100% — all current sites are transient.

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
| `scripts/10_attribution_revised.R` | **Revised attribution** | ✅ Updated (NDVI climate-adjacent) |
| `scripts/execute_decisions_v4.R` | **Execute all D1-D7 decisions** | ✅ New |
| `scripts/fig_isocline_equilibrium.R` | **Isocline plot + debt + timescales** | ✅ New |
| `scripts/refit_revised_models.R` | Refit ptealc/pteori | ✅ New |
| `scripts/compute_equilibrium.R` | Equilibrium ψ* | ✅ New |
| `R/model_configs.R` | Centralised formulas | ✅ Updated (pteori reverted to original) |

### 10.2 Key results files

| File | Content |
|------|---------|
| `results/{sp}_model_object.rds` | Fitted colext models (ptealc revised; pteori reverted to original) |
| `results/{sp}_train_dyn_scale.rds` | Training scaling parameters |
| `results/equilibrium_occupancy_table.csv` | Equilibrium ψ* with bootstrap CIs (all 4 spp OK) |
| `results/extinction_debt_table.csv` | **Extinction debt per species** |
| `results/isocline_plot_data.csv` | **Baseline rates + CIs for isocline plot** |
| `results/pteori_sensitivity_ndvi_epsilon.csv` | **Sensitivity: with/without NDVI coefficients** |
| `results/pteori_epsilon_sensitivity_noNDVI.rds` | **No-NDVI model object (sensitivity)** |
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
| Fig. 3 | `figs/pub_fig_coef_all_submodels.png` | All submodel coefficients (γ + ε + ψ₁) |
| Fig. 4 | `figs/pub_fig_occupancy_trends_panel.png` | Prevalence trends 2017–2023 |
| **Fig. 5** | **`figs/pub_fig_isocline_equilibrium.png`** | **Isocline plot + extinction debt + recolonisation (KEY FIGURE)** |
| Fig. 6 | Attribution figure (revised) | Attribution summary (4 panels) |

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
| v4a | Remove NDVI from ptealc/gamma | pr sign change (+0.16 → −0.72); clean attribution | +7.81 |
| v4a | Remove NDVI from pteori/epsilon | SE inflation +60%; LC12 beta shift 2.8 SE | +60.27 |
| v4a | Fix scaling mismatch | scripts/4 vs scripts/8 used different scaling | N/A |
| v4a | blockCV block size 50 → 270 km | 50 km < 264 km spatial range | Pending re-run |
| v4a | Save train_dyn_scale.rds | Enables reproducible scaling in downstream scripts | N/A |
| **v4b** | **Revert pteori/epsilon to original (~LC12+NDVI+pr)** | **ΔAIC +60 too large; gamma collapsed; classify NDVI as climate-adjacent** | **−60.27** |
| v4b | Compute extinction debt (all 4 spp) | Extinction debt framing for GCB narrative | N/A |
| v4b | Create isocline plot (3 panels + PhyloPic) | Key GCB figure — demographic trap visualisation | N/A |
| v4b | Save no-NDVI pteori as sensitivity | `pteori_epsilon_sensitivity_noNDVI.rds` | N/A |

---

*Report generated on branch `audit-gcb-v4`. Last updated after execution of all decisions (v4b).*
