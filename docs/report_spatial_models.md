# Bayesian Spatial Occupancy Models for Iberian Steppe Birds (2017--2023)

**Supplementary Analysis Report -- March 2026**
**Branch:** `analysis-2023` | **Script:** `scripts/7_spatial_results.R`

---

## 1. Introduction

This report documents the results of fitting Bayesian spatial occupancy models using the `spOccupancy` package as a complement to the frequentist dynamic occupancy models (`unmarked::colext`) presented in the main analysis report (`docs/analysis_report_2023.md`).

### 1.1 Motivation

The colext models assume independence among sites conditional on covariates. Validation of the colext predictions against the Spanish Biodiversity Atlas revealed significant positive spatial autocorrelation in residuals (Moran's I = 0.14--0.24, *P* < 10^-79^) for all four species. This indicates that occupancy at a given site is partially predicted by occupancy at nearby sites, even after accounting for environmental covariates.

Spatial autocorrelation can arise from:
- **Unmeasured spatially-structured covariates** (e.g., soil type, agricultural practices)
- **Dispersal-mediated spatial dependence** (conspecific attraction, metapopulation connectivity)
- **Spatial mismatch** between sampling resolution and the scale of ecological processes

Including an explicit spatial random effect addresses this by absorbing residual spatial structure, leading to more accurate occupancy estimates and providing ecologically informative spatial range parameters.

### 1.2 What each framework contributes

The two modelling frameworks answer **different ecological questions** and neither can replace the other:

| | colext (unmarked) | stPGOcc (spOccupancy) |
|---|---|---|
| **What it estimates** | Explicit colonisation ($\gamma$) and extinction ($\varepsilon$) as Markov transition probabilities, with environmental covariates on each | Year-specific occupancy ($\psi_{j,t}$) with a spatial random effect ($w_j$); colonisation and extinction are *derived* from the latent states, not explicitly modelled |
| **Unique contribution** | Identifies **which environmental variables drive colonisation and extinction** (e.g., precipitation promotes *P. alchata* colonisation; urban land cover increases *O. tarda* extinction). This mechanistic insight is **not available from spOccupancy** | Provides **spatially-explicit maps** of occupancy, colonisation hotspots, extinction hotspots, and temporal trends, all accounting for spatial autocorrelation. Also estimates the **effective spatial range** of each species |
| **Key limitation** | Assumes site independence (violated: Moran's I significant for all species) | Does not parameterise $\gamma$ and $\varepsilon$ as functions of environmental covariates; cannot tell us *why* colonisation or extinction occurs |

> [!IMPORTANT]
> **Summary of the division of labour:**
> - **colext** answers: *"What environmental factors cause a site to be colonised or go locally extinct?"* (e.g., cropland proportion drives *T. tetrax* colonisation and extinction)
> - **stPGOcc** answers: *"Where and when does colonisation/extinction happen, at what spatial scale, and what do the spatial patterns look like?"* (maps, trends, hotspots, spatial range)
>
> The two frameworks are **complementary, not competing**. All occupancy maps, colonisation/extinction maps, and temporal trends in this report come from the Bayesian spatial model (stPGOcc). The environmental drivers of turnover come exclusively from colext.

---

## 2. Methods

### 2.1 spOccupancy framework

We used the `spOccupancy` R package (Doser et al. 2022) to fit three model variants per species in a progressive strategy:

1. **tPGOcc** -- Multi-season Polya-Gamma occupancy model *without* spatial structure (non-spatial baseline). Includes AR(1) temporal correlation.
2. **spPGOcc** -- Single-season spatial model for 2023 only, with a Nearest Neighbor Gaussian Process (NNGP, *m* = 5 neighbours) spatial random effect. Proof of concept for spatial structure.
3. **stPGOcc** -- Full spatio-temporal model (2017--2023) combining NNGP spatial structure with AR(1) temporal correlation. This is the primary model for inference.

### 2.2 Model specification

**Occupancy sub-model:**

$$\text{logit}(\psi_{j,t}) = \mathbf{x}_j'\boldsymbol{\beta} + w_j + \eta_t$$

where $w_j$ is the spatial random effect (Gaussian Process) and $\eta_t$ captures temporal autocorrelation via AR(1).

**Detection sub-model:**

$$\text{logit}(p_{j,t,k}) = \mathbf{v}_{j,t,k}'\boldsymbol{\alpha}$$

**Covariates:**

The occupancy covariates include the species-specific `psi_vars` from the colext model configuration (see `R/model_configs.R`) plus dynamic variables NDVI and maximum temperature (tmmx). Detection covariates are duration and number of observers.

| Species | Occupancy covariates | Detection covariates |
|---------|---------------------|---------------------|
| *O. tarda* | bio1, bio2, tree_cover, grass_cover, topo_elev, NDVI, tmmx | duration, observers |
| *P. alchata* | bio1, bio2, tree_cover, grass_cover, topo_aspect, NDVI, tmmx | duration, observers |
| *P. orientalis* | bio2, tree_cover, grass_cover, NDVI, tmmx | duration, observers |
| *T. tetrax* | bio2, tree_cover, grass_cover, topo_elev, NDVI, tmmx | duration, observers |

### 2.3 MCMC settings

| Setting | tPGOcc | spPGOcc | stPGOcc |
|---------|--------|---------|---------|
| Samples/chain | 10,000 | 10,000 | 20,000 |
| Burn-in | 5,000 | 5,000 | 10,000 |
| Thin | 5 | 5 | 20 |
| Chains | 3 | 3 | 3 |
| Post-burn samples/chain | 1,000 | 1,000 | 500 |
| Total post-burn samples | 3,000 | 3,000 | 1,500 |
| NNGP neighbours | -- | 5 | 5 |
| AR(1) | Yes | -- | Yes |
| Covariance model | -- | Exponential | Exponential |

### 2.4 Spatial parameters

- **Spatial variance ($\sigma^2$):** Inverse-Gamma(2, 1) prior
- **Decay parameter ($\phi$):** Uniform prior bounded by 3/max_distance (long-range) to 3/(0.01 * max_distance) (short-range)
- **Effective spatial range:** For the exponential covariance function, the distance at which spatial correlation drops to ~5% is $3/\phi$ (in metres, converted to km)

### 2.5 Derived turnover rates from stPGOcc

> [!NOTE]
> **Key distinction:** In colext, colonisation ($\gamma$) and extinction ($\varepsilon$) are *explicit model parameters* with environmental covariates attached to each (e.g., $\gamma \sim \text{NDVI} + \text{precipitation}$). In stPGOcc, there are no $\gamma$ or $\varepsilon$ parameters in the model. Instead, we **derive** turnover rates *post hoc* from the latent occupancy states $z_{j,t}$.

For each MCMC sample, we compute:

$$\hat{\gamma}_t = \frac{\sum_j z_{j,t}(1 - z_{j,t-1})}{\sum_j (1 - z_{j,t-1})} \qquad \hat{\varepsilon}_t = \frac{\sum_j (1 - z_{j,t})z_{j,t-1}}{\sum_j z_{j,t-1}}$$

This gives us posterior distributions of annual colonisation and extinction rates, as well as site-specific colonisation/extinction hotspot maps. These derived rates tell us **where and when** turnover occurs, but not **why** (for the *why*, see the colext gamma/epsilon covariate results in the main report).

### 2.6 Diagnostics

Model convergence assessed via:
- **Rhat** (potential scale reduction factor; target < 1.1)
- **ESS** (effective sample size; target > 100)
- **WAIC** (Widely Applicable Information Criterion; lower = better)
- **Moran's I** on model residuals (target: non-significant after spatial correction)

---

## 3. Results

### 3.1 Model comparison (WAIC)

All four species show **strong support** for the spatial model (stPGOcc) over the non-spatial baseline (tPGOcc), with $\Delta$WAIC ranging from 352 to 457.

| Species | tPGOcc WAIC | stPGOcc WAIC | $\Delta$WAIC | Interpretation |
|---------|-------------|--------------|--------------|----------------|
| *O. tarda* | 2496.9 | 2040.4 | **456.5** | Very strong |
| *P. alchata* | 2103.8 | 1655.2 | **448.6** | Very strong |
| *P. orientalis* | 2273.7 | 1821.3 | **452.4** | Very strong |
| *T. tetrax* | 1890.4 | 1538.8 | **351.6** | Very strong |

![WAIC comparison](../figs/pub_fig_spatial_waic.png)

### 3.2 Spatial autocorrelation reduction

The spatial random effect in stPGOcc effectively eliminates the residual spatial autocorrelation that was present in the non-spatial model. In all cases, Moran's I drops from significantly positive to slightly negative (non-significant), indicating slight over-correction -- expected behaviour for spatial models.

| Species | Moran's I (tPGOcc) | *P*-value | Moran's I (stPGOcc) | *P*-value | Reduction |
|---------|--------------------|---------|-----------------------|---------|-----------|
| *O. tarda* | 0.145 | 1.3 x 10^-79^ | -0.060 | 1.00 | 141% |
| *P. alchata* | 0.236 | 1.4 x 10^-227^ | -0.042 | 1.00 | 118% |
| *P. orientalis* | 0.173 | 2.2 x 10^-123^ | -0.035 | 1.00 | 120% |
| *T. tetrax* | 0.144 | 6.9 x 10^-79^ | -0.042 | 1.00 | 129% |

![Moran's I comparison](../figs/pub_fig_spatial_moran.png)

### 3.3 Effective spatial ranges

The estimated spatial ranges vary substantially across species and align with known movement ecology:

| Species | Range (km) | 95% CI | Ecological interpretation |
|---------|-----------|--------|-------------------------|
| *T. tetrax* | **42.7** | See figure | Local: strong breeding-season site fidelity, limited natal dispersal (~3--5 km). Range reflects the spatial scale of suitable habitat patches (cereal pseudo-steppe) rather than individual movement. |
| *O. tarda* | **46.7** | See figure | Mesoscale: matches inter-lek connectivity in the Iberian lek network. Males disperse 20--50 km between leks; females show high breeding site fidelity. |
| *P. orientalis* | **171.6** | See figure | Intermediate: broader habitat connectivity across the cereal pseudo-steppe belt. More sedentary than *P. alchata* but occupies large landscape units. |
| *P. alchata* | **263.9** | See figure | Regional/nomadic: consistent with semi-nomadic ecology. Tracks ephemeral food resources (seeds) and water sources across large semi-arid landscapes of central-eastern Spain. |

![Spatial range comparison](../figs/pub_fig_spatial_range.png)

**Conservation implication:** Protected areas or habitat management zones should encompass areas with diameters at least matching the species' spatial range to capture the relevant spatial dependence structure. For *P. alchata*, this implies landscape-scale conservation strategies spanning >250 km, while *T. tetrax* conservation can be effective at local scales (~40--50 km).

### 3.4 Coefficient estimates

The stPGOcc occupancy coefficients are broadly consistent with the colext initial occupancy ($\psi_1$) estimates, validating the colext results despite the violated independence assumption.

Key patterns conserved across both frameworks:
- **Tree cover** is strongly negative for all species (steppe specialists avoid forest)
- **Grass cover** is generally positive (open grassland habitat)
- **bio2** (temperature seasonality) is positive for most species (continental climate preference)

However, some coefficients change when spatial structure is accounted for:
- The **intercept** shifts substantially (more negative in stPGOcc), reflecting the baseline occupancy adjustment
- Some covariate effects are attenuated by the spatial random effect, suggesting partial confounding with spatially-structured unmeasured variables

![Coefficient forest plot](../figs/pub_fig_spatial_forest_betas.png)

![colext vs spOccupancy comparison](../figs/pub_fig_spatial_colext_comparison.png)

### 3.5 Spatial random effects

The spatial random effect ($w_j$) maps reveal the residual spatial structure not captured by environmental covariates.

![Spatial random effect maps](../figs/pub_fig_spatial_w_maps.png)

Interpretation by species:

- **O. tarda:** Positive spatial residuals concentrated in Castilla y Leon and Extremadura, corresponding to the core lek network. Negative residuals in Mediterranean coast and Galicia.
- **P. alchata:** Strong positive residuals in the Ebro Valley and south-eastern arid belt, reflecting the nomadic network. Broad spatial scale matches the 264 km range estimate.
- **P. orientalis:** Similar to *P. alchata* but more diffuse. Positive residuals across the central-eastern meseta.
- **T. tetrax:** Highly localised positive residuals, consistent with the 43 km spatial range. Clustered in traditional breeding areas (Lleida, Castilla-La Mancha, Alentejo frontier).

### 3.6 Occupancy trends

Derived from the posterior distribution of latent occupancy states ($z_{j,t}$), the stPGOcc occupancy trends provide a spatially-explicit estimate of the proportion of sites occupied over 2017--2023.

![Occupancy trends](../figs/pub_fig_spatial_occupancy_trends.png)

> [!NOTE]
> The occupancy estimates from stPGOcc are higher than from colext because the spatial random effect captures additional spatial structure that inflates site-level occupancy probability relative to the environment-only model. Trends should be interpreted in relative terms (direction and magnitude of change), not as absolute occupancy levels. Both frameworks show increasing trends, likely influenced by increasing eBird survey effort over the study period.

### 3.7 Colonisation and extinction dynamics

Turnover rates derived from the stPGOcc latent occupancy states ($z_{j,t}$) show the balance between colonisation and extinction events across the study period. These rates tell us **where and when** turnover is happening spatially, while the colext results (main report) tell us **which environmental variables** drive that turnover.

![Colonization and extinction rates](../figs/pub_fig_spatial_col_ext_rates.png)

The spatial distribution of colonisation and extinction events reveals hotspots of demographic activity:

![Colonization and extinction hotspot maps](../figs/pub_fig_spatial_col_ext_maps.png)

Key patterns:
- **Colonisation hotspots** tend to be located at range margins, where new sites are being occupied
- **Extinction hotspots** overlap with areas of land-use change and habitat fragmentation
- *T. tetrax* shows the most spatially concentrated turnover, consistent with its localised ecology and the colext finding that cropland (LC12) is the sole driver of both colonisation and extinction
- *P. alchata* shows the most diffuse turnover, reflecting nomadic range-wide dynamics; colext identifies NDVI and precipitation as the drivers of these colonisation events

### 3.8 Occupancy change maps

The net change in occupancy probability ($\Delta\psi$ = $\psi_{2023}$ - $\psi_{2017}$) maps provide a spatial perspective on population trends.

![Occupancy change maps](../figs/pub_fig_spatial_delta_psi_maps.png)

### 3.9 Cross-species synthesis

**Model support:** All four species show overwhelming support for including spatial structure ($\Delta$WAIC > 350). This is consistent with the significant Moran's I values found in the colext validation.

**Spatial ecology gradient:** The four species span a gradient from localised (43 km) to regional/nomadic (264 km) spatial structure, which directly maps onto their known movement ecology:

```
Local (< 50 km)          Mesoscale (50-100 km)       Regional (> 100 km)
T. tetrax (43 km)        O. tarda (47 km)            P. orientalis (172 km)
                                                      P. alchata (264 km)
```

**Coefficient consistency:** The direction and relative magnitude of environmental effects are largely preserved between colext and stPGOcc, which validates the colext inference despite the violated site-independence assumption.

---

## 4. Integration with colext results

### 4.1 What colext tells us that stPGOcc cannot

The colext model parameterises colonisation ($\gamma$) and extinction ($\varepsilon$) as explicit functions of environmental covariates. This is its unique and irreplaceable contribution. The stPGOcc model has no equivalent -- it derives turnover from latent states but cannot attribute it to specific variables.

**Environmental drivers of colonisation ($\gamma$) from colext:**

| Species | Key drivers | Ecological interpretation |
|---------|------------|--------------------------|
| *O. tarda* | NDVI (+), pr (+), tmmn (-), tmmx (+) | Colonisation favoured by productive sites with cold winters |
| *P. alchata* | NDVI (+), pr (+) | Tracks vegetation greenness and rainfall -- consistent with nomadic seed-tracking |
| *P. orientalis* | LC7 (+), NDVI (+), tmmn (-), tmmx (+) | Open shrubland and NDVI promote colonisation |
| *T. tetrax* | LC12 (+) | Cropland (cereals) is the **sole** driver -- highly simplified due to sparse transitions |

**Environmental drivers of extinction ($\varepsilon$) from colext:**

| Species | Key drivers | Ecological interpretation |
|---------|------------|--------------------------|
| *O. tarda* | LC6 (+), LC13 (+), tmmx (+) | Urban/artificial surfaces and heat increase local extinction |
| *P. alchata* | pr (-), tmmx (+) | Drought and heat stress drive local extinction |
| *P. orientalis* | LC12 (-), NDVI (-), pr (-) | Cropland and vegetation **reduce** extinction (habitat stability) |
| *T. tetrax* | LC12 (-) | Cropland **reduces** extinction -- cereal pseudo-steppe is critical habitat |

> [!NOTE]
> These environmental drivers of turnover are **exclusively from colext**. The stPGOcc model cannot identify which variables cause colonisation or extinction -- it only tells us where and when these events occur.

### 4.2 What stPGOcc tells us that colext cannot

The stPGOcc model provides three types of information that colext cannot:

1. **Spatially-explicit occupancy, colonisation, and extinction maps** -- colext assumes site independence and produces predictions that ignore spatial autocorrelation. The stPGOcc maps (Figures A, I, J) account for this structure and are therefore more spatially coherent.

2. **Effective spatial range** -- a model-based estimate of the spatial scale of occupancy dependence (43--264 km), which links directly to species' dispersal ecology and metapopulation connectivity.

3. **Spatial random effect** ($w_j$) -- reveals where occupancy is higher or lower than predicted by environment alone, exposing unmeasured spatial drivers (e.g., traditional land management, local conservation efforts, dispersal corridors).

### 4.3 Shared covariate comparison (occupancy $\psi$)

Both frameworks estimate the effect of environmental covariates on occupancy probability. For the three covariates shared across all species (bio2, tree_cover, grass_cover):

- **tree_cover**: Consistently strongly negative in both frameworks (all species). The Bayesian estimates are generally attenuated relative to colext, likely because the spatial random effect absorbs some spatially-structured variation that colext attributed to tree cover.
- **grass_cover**: Positive in both frameworks, though some attenuation in stPGOcc.
- **bio2**: Positive in both frameworks for most species; strongest for *O. tarda*.

The broad consistency of covariate effects across frameworks **validates the colext occupancy estimates** despite the violated site-independence assumption. The environmental signal is real; colext simply lacked a spatial component to absorb the residual structure.

### 4.4 Integrated ecological narrative

Combining both frameworks yields a complete picture for each species:

**Example -- *Tetrax tetrax* (Little Bustard):**
- *Why* does it colonise new sites? Because cropland (cereals) increases locally (colext: $\gamma \sim$ LC12)
- *Why* does it go locally extinct? When cropland decreases (colext: $\varepsilon \sim$ LC12)
- *Where* are the colonisation and extinction hotspots? See stPGOcc maps (Figure I) -- concentrated in Castilla-La Mancha and Lleida
- *At what spatial scale* does occupancy operate? 43 km (stPGOcc) -- highly localised, matching known site fidelity
- *What is the overall trend?* Increasing occupancy 2017--2023 (stPGOcc Figure G), though likely inflated by increasing eBird effort

---

## 5. Discussion

### 5.1 Why spatial models matter for steppe birds

Steppe birds in the Iberian Peninsula occupy a fragmented landscape of cereal pseudo-steppe embedded in a matrix of forest, scrubland, and urbanised areas. This spatial configuration creates occupancy patterns that violate the independence assumption of standard occupancy models. Our results demonstrate that including spatial structure:
- Dramatically improves model fit ($\Delta$WAIC > 350 for all species)
- Eliminates residual spatial autocorrelation
- Provides ecologically meaningful spatial range estimates
- Does not fundamentally alter the environmental covariate effects identified by colext

### 5.2 The value of combining both frameworks

Neither framework alone tells the full story:

- **colext alone** identifies the environmental drivers of turnover (e.g., "cropland promotes colonisation") but cannot map where those dynamics play out spatially, and its occupancy maps suffer from unmodelled spatial autocorrelation.
- **stPGOcc alone** produces high-quality spatial maps and trend estimates but cannot explain *why* turnover occurs at specific sites -- it only shows *that* it occurs.

By combining both:
- We know *what drives* colonisation and extinction (colext), and we can map *where and when* those events concentrate spatially (stPGOcc)
- We can validate the colext environmental effects against a framework that accounts for spatial autocorrelation -- the consistency we observe (Section 4.3) strengthens confidence in both sets of results
- Conservation recommendations can be both **mechanistic** (manage cropland proportion to benefit *T. tetrax*) and **spatially targeted** (focus on the colonisation hotspots identified in Figure I)

### 5.3 Conservation implications

1. **Protected area design:** The spatial range parameter directly informs the minimum spatial extent that conservation plans should encompass. For *T. tetrax* and *O. tarda* (~43--47 km), SPAs and Natura 2000 sites at the landscape scale are sufficient. For *P. alchata* (~264 km), regional-scale habitat networks are required.

2. **Connectivity corridors:** The spatial random effect maps identify areas where occupancy is higher (or lower) than predicted by environment alone. Positive spatial residuals in the absence of suitable habitat suggest stepping-stone connectivity; negative residuals in suitable habitat suggest barrier effects.

3. **Priority areas:** The colonisation/extinction hotspot maps (stPGOcc, Figure I) identify *where* to intervene, while the colext covariate effects identify *what* to manage (e.g., maintain cereal pseudo-steppe for *T. tetrax*; reduce urban expansion near *O. tarda* leks).

---

## 6. Caveats and Limitations

### 6.1 MCMC convergence

The spatial hyperparameters ($\sigma^2$, $\phi$) and temporal parameters ($\sigma^2_t$, $\rho$) show poor convergence in several stPGOcc models:

| Parameter | Typical Rhat | Typical ESS | Target |
|-----------|-------------|-------------|--------|
| $\beta$ coefficients | 1.00--1.07 | 185--870 | < 1.1 / > 100 |
| $\alpha$ coefficients | 1.00--1.01 | 1,300--3,000 | < 1.1 / > 100 |
| $\sigma^2$ (spatial) | 1.01--1.07 | 31--48 | **< 1.1 / > 100** |
| $\phi$ (decay) | 1.25--1.31 | 44--56 | **< 1.1 / > 100** |
| $\sigma^2_t$ (temporal) | 1.84--4.19 | 9--108 | **< 1.1 / > 100** |

> [!IMPORTANT]
> The spatial variance and range parameters have not fully converged. The coefficient estimates ($\beta$, $\alpha$) are reliable, but the spatial range estimates should be treated as preliminary. **Recommendation:** Increase MCMC iterations to at least 25,000 post-burnin samples (currently 1,500) for publication-grade inference on spatial hyperparameters.

### 6.2 Derived vs explicit turnover rates

The colonisation and extinction rates derived from stPGOcc are **not directly comparable** to the colext $\gamma$ and $\varepsilon$ estimates because:
- colext conditions turnover on the previous year's occupancy state via Markov transitions
- stPGOcc derives turnover from the posterior of latent states, which are themselves influenced by the spatial random effect
- The derived rates may be more spatially coherent but less mechanistically interpretable

### 6.3 eBird sampling bias

As noted in the main report, increasing eBird survey effort over 2017--2023 likely inflates apparent colonisation rates. This affects both frameworks equally but is particularly important when interpreting the occupancy trend and colonisation hotspot maps.

---

## 7. Figure and table inventory

### Tables

| File | Description |
|------|-------------|
| `results/pub_spatial_table1_model_comparison.csv` | WAIC, Moran's I, spatial range for all species |
| `results/pub_spatial_table2_stpgocc_coefficients.csv` | stPGOcc posterior summaries (mean, SD, 95% CI, Rhat, ESS) |

### Figures

| File | Description |
|------|-------------|
| `figs/pub_fig_spatial_w_maps.png` | Spatial random effect (w) maps, 4 species |
| `figs/pub_fig_spatial_forest_betas.png` | Coefficient forest plots (tPGOcc vs stPGOcc) |
| `figs/pub_fig_spatial_range.png` | Effective spatial range comparison with ecological bands |
| `figs/pub_fig_spatial_waic.png` | WAIC model comparison |
| `figs/pub_fig_spatial_moran.png` | Moran's I before/after spatial model |
| `figs/pub_fig_spatial_colext_comparison.png` | colext vs stPGOcc shared covariate comparison |
| `figs/pub_fig_spatial_occupancy_trends.png` | Occupancy trends 2017--2023 |
| `figs/pub_fig_spatial_col_ext_rates.png` | Derived colonisation/extinction rates |
| `figs/pub_fig_spatial_col_ext_maps.png` | Colonisation/extinction hotspot maps |
| `figs/pub_fig_spatial_delta_psi_maps.png` | Occupancy change maps ($\Delta\psi$) |

---

*Generated by `scripts/7_spatial_results.R` on the `analysis-2023` branch.*
