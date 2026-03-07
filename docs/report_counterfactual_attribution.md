# Counterfactual Attribution Analysis: Climate vs Land-Use Drivers of Steppe Bird Occupancy Dynamics (2017--2023)

**Analysis Report -- March 2026**
**Branch:** `counterfactual-attribution` | **Script:** `scripts/8_counterfactual_attribution.R`

---

## 1. Rationale and Conceptual Framework

### 1.1 Why counterfactual attribution?

The dynamic occupancy models fitted in the main analysis (Section 2 of the analysis report) establish that colonisation ($\gamma$) and extinction ($\varepsilon$) rates vary as a function of both climate variables (temperature, precipitation, NDVI) and land-use variables (MODIS land-cover proportions). A natural follow-up question is: **over the 2017--2023 study period, how much of the observed variation in turnover rates is attributable to each driver category?**

This question cannot be answered simply by examining coefficient magnitudes, because the *effect* of each covariate on predicted turnover depends not only on its coefficient but also on *how much it has actually changed* over the study period. A covariate with a large coefficient but no temporal trend will contribute nothing to observed changes; conversely, a covariate with a modest coefficient but a strong, widespread trend can dominate the dynamics.

### 1.2 Factorial scenario design

We decompose the attribution using a $2 \times 2$ factorial design, analogous to a two-factor ANOVA:

| Scenario | Climate covariates | Land-use covariates | Interpretation |
|----------|--------------------|---------------------|----------------|
| **S0** (Null baseline) | Frozen at 2017 | Frozen at 2017 | No environmental change |
| **S1** (Climate only) | Observed 2017--2023 | Frozen at 2017 | Climate trends in isolation |
| **S2** (Land use only) | Frozen at 2017 | Observed 2017--2023 | Land-use change in isolation |
| **S3** (Combined) | Observed 2017--2023 | Observed 2017--2023 | Actual observed conditions |

**"Frozen at 2017"** means each site retains its own 2017 value for all subsequent years. This preserves spatial heterogeneity but removes temporal change.

**Covariate classification:**
- *Climate:* NDVI, precipitation (pr), minimum temperature (tmmn), maximum temperature (tmmx). NDVI is classified as climate-adjacent given its strong coupling to precipitation and temperature at annual timescales.
- *Land use:* All MODIS MCD12Q1 land-cover proportions (Classes 6, 7, 12, 13).

### 1.3 Attribution decomposition

For each site $i$, year $t$, and turnover parameter $\theta \in \{\gamma, \varepsilon\}$:

$$\Delta\theta_{\text{climate}} = \theta_{S1} - \theta_{S0}$$

$$\Delta\theta_{\text{landuse}} = \theta_{S2} - \theta_{S0}$$

$$\Delta\theta_{\text{interaction}} = \theta_{S3} - \theta_{S1} - \theta_{S2} + \theta_{S0}$$

$$\Delta\theta_{\text{total}} = \theta_{S3} - \theta_{S0}$$

The interaction term captures non-additive effects -- i.e., cases where the combined impact of climate and land-use change differs from the sum of their individual effects (due to nonlinearities in the logistic link function).

---

## 2. Covariate Trends (2017--2023)

Before examining attribution, we quantify whether the environmental covariates have actually changed systematically over the study period. For each covariate and site, we fit a linear regression (value ~ year) and extract the slope.

### 2.1 Climate covariates

| Covariate | Mean trend (per year) | SD | Sig. positive | Sig. negative |
|-----------|----------------------|-----|---------------|---------------|
| NDVI | -22.7 NDVI units | 82.1 | 1.3% | 2.9% |
| Precipitation (pr) | -0.27 mm | 2.11 | 6.2% | 0% |
| Min. temperature (tmmn) | +1.23 $^\circ$C | 0.69 | 0.1% | 0% |
| Max. temperature (tmmx) | +0.51 $^\circ$C | 0.64 | 0% | 0% |

*Values shown for otitar sites; ptealc/pteori sites show similar patterns.*

**Interpretation:** Minimum temperatures show a consistent warming trend across virtually all sites (+1.23 $^\circ$C/yr mean), though few individual sites reach significance over the 7-year window. Precipitation shows a weak declining trend. NDVI declines are spatially heterogeneous (large SD relative to mean). Maximum temperature trends are modest.

### 2.2 Land-use covariates

| Covariate | Mean trend (per year) | SD | Sig. positive | Sig. negative |
|-----------|----------------------|-----|---------------|---------------|
| LC6 (closed shrubland) | -0.011 pp | 0.12 | 0.1% | 1.2% |
| LC7 (open shrubland) | +0.016 pp | 1.17 | 5.8% | 8.1% |
| LC12 (cropland) | +0.046 pp | 1.09 | 9.3--9.6% | 6.0--6.4% |
| LC13 (urban/built-up) | +0.004 pp | 0.04 | 1.4% | 0.1% |

*pp = percentage points of land cover proportion per year.*

**Interpretation:** Cropland (LC12) shows the most dynamic land-cover change, with both gains and losses occurring at comparable frequencies (9.3--9.6% of sites gaining, 6.0--6.4% losing, $P < 0.05$). Open shrubland (LC7) also shows bidirectional change. Urban expansion (LC13) is slow but consistently positive. Closed shrubland (LC6) is stable or declining slightly.

### 2.3 Key finding

The relatively weak and spatially heterogeneous covariate trends mean that **the counterfactual attribution effects will be small in absolute terms** -- typically on the order of $10^{-3}$ to $10^{-4}$ change in turnover probability. This is expected: a 7-year window is short relative to the timescales of climate change and land-use transformation. The value of the analysis lies in identifying the *relative* contributions and *spatial patterns* of each driver, not in detecting large absolute effects.

---

## 3. Attribution Results

### 3.1 Cross-species summary

| Species | $\Delta\gamma_{\text{climate}}$ | $\Delta\gamma_{\text{landuse}}$ | $\Delta\varepsilon_{\text{climate}}$ | $\Delta\varepsilon_{\text{landuse}}$ | Dominant driver |
|---------|------|------|------|------|------|
| *Otis tarda* | -0.00061 | 0 | +0.00074 | -0.00110 | **Climate** |
| *Pterocles alchata* | +0.00007 | 0 | -0.00734 | 0 | **Climate** |
| *Pterocles orientalis* | -0.00034 | -0.00018 | +0.00573 | -0.00021 | **Climate** |
| *Tetrax tetrax* | 0 | -0.00000 | 0 | +0.00001 | **Land use** |

*Values are mean $\pm$ SD across all sites and years. Interaction terms omitted for clarity (all $< 5 \times 10^{-4}$).*

### 3.2 Species-specific results

#### 3.2.1 *Otis tarda* (Great Bustard)

![Attribution summary -- Otis tarda](../figs/fig_attribution_summary_otitar.png)

For *Otis tarda*, climate trends were the dominant driver of colonisation dynamics, while land-use change had a larger effect on extinction than climate. The climate-driven reduction in colonisation probability ($\Delta\gamma_{\text{climate}} = -0.00061$) reflects the negative coefficients for NDVI ($\beta = -1.45$) and tmmn ($\beta = -4.52$) in the colonisation submodel, combined with declining NDVI and rising minimum temperatures. These trends reduce the probability of new site colonisation.

For extinction, land-use change contributed more strongly ($\Delta\varepsilon_{\text{landuse}} = -0.00110$) than climate ($\Delta\varepsilon_{\text{climate}} = +0.00074$). The land-use effect is driven by changes in closed shrubland (LC6, $\beta_\varepsilon = 1.86$) and urban cover (LC13, $\beta_\varepsilon = 2.80$). The slight urban expansion (LC13 trend = +0.004 pp/yr) amplifies local extinction risk at affected sites, consistent with the species' known sensitivity to habitat fragmentation and human disturbance.

The interaction term is negligible ($< 5 \times 10^{-5}$), indicating that climate and land-use effects are approximately additive for this species.

> **Caveat:** NDVI is classified as climate-adjacent in this analysis. Its partial coupling to land management practices (irrigation, crop rotation) means that the climate attribution may partially capture land-use effects mediated through vegetation productivity.

#### 3.2.2 *Pterocles alchata* (Pin-tailed Sandgrouse)

![Attribution summary -- Pterocles alchata](../figs/fig_attribution_summary_ptealc.png)

For *P. alchata*, the attribution is entirely climate-driven by construction: the fitted model contains no land-use covariates in either the colonisation ($\gamma \sim$ NDVI + pr) or extinction ($\varepsilon \sim$ pr + tmmx) submodels. This reflects the convergence constraints that forced removal of land-cover variables during model selection (complete separation with $<$ 2% variation in LC classes at occupied sites).

The climate effect on extinction ($\Delta\varepsilon_{\text{climate}} = -0.00734$) is the largest among all species, but with very high spatial variance (SD = 0.296). This reflects the extreme coefficient magnitudes in the extinction submodel ($\beta_{\text{pr}} = -52.1$; $\beta_{\text{tmmx}} = -70.9$), which amplify even small covariate changes into large predicted shifts at some sites.

> **Caveat:** The large extinction coefficients suggest potential estimation instability due to small sample sizes (only 18 extinction events across all years). The climate attribution for this species should be interpreted with caution; the high SD reflects a few sites with extreme predicted changes rather than a consistent landscape-level signal.

#### 3.2.3 *Pterocles orientalis* (Black-bellied Sandgrouse)

![Attribution summary -- Pterocles orientalis](../figs/fig_attribution_summary_pteori.png)

*P. orientalis* is the only species with substantive attribution to both climate and land-use drivers. For colonisation, climate accounts for approximately 75% of the total change ($\Delta\gamma_{\text{climate}} = -0.00034$; $\Delta\gamma_{\text{landuse}} = -0.00018$). This reflects the strong temperature effects in the colonisation submodel ($\beta_{\text{tmmn}} = -8.10$; $\beta_{\text{tmmx}} = 9.14$), combined with the warming trends over 2017--2023. The opposing signs of tmmn and tmmx create a nonlinear temperature response: warming minimum temperatures reduce colonisation, while warming maximum temperatures increase it, with the net effect being a slight decrease.

For extinction, climate dominates ($\Delta\varepsilon_{\text{climate}} = +0.00573$), driven primarily by NDVI ($\beta_\varepsilon = 2.05$) and precipitation ($\beta_\varepsilon = 1.69$). Declining NDVI and precipitation trends increase extinction rates. Land-use change via cropland proportion (LC12, $\beta_\varepsilon = -0.83$) contributes a small protective effect ($\Delta\varepsilon_{\text{landuse}} = -0.00021$), consistent with this species' use of dry cereal croplands as foraging habitat.

The interaction term ($\Delta\varepsilon_{\text{interaction}} = -0.00049$) is small but non-zero, indicating mild synergistic effects between climate and land-use change on extinction dynamics.

#### 3.2.4 *Tetrax tetrax* (Little Bustard)

![Attribution summary -- Tetrax tetrax](../figs/fig_attribution_summary_tettet.png)

For *T. tetrax*, the attribution is entirely to land-use change by construction: the model contains only cropland proportion (LC12) in both the colonisation ($\gamma \sim$ LC12) and extinction ($\varepsilon \sim$ LC12) submodels, with no climate covariates. This reflects the severe convergence constraints for this species (only 25 occupied-site transitions, of which 9 were extinctions), which prevented the inclusion of any additional covariates.

The absolute effects are negligible ($\Delta\gamma_{\text{landuse}} = -1 \times 10^{-6}$; $\Delta\varepsilon_{\text{landuse}} = 9 \times 10^{-6}$), reflecting the near-zero net trend in LC12 at the landscape level (mean = +0.027 pp/yr, but with bidirectional changes cancelling spatially). The spatial maps reveal that the few sites with detectable attribution are concentrated in agricultural transition zones where cropland proportion has changed most.

The key ecological signal is that cropland proportion has a protective effect on this species: higher LC12 reduces extinction ($\beta_\varepsilon = -0.57$, $P = 0.042$) and increases colonisation ($\beta_\gamma = 0.90$, $P = 0.039$). Any future conversion of cropland to other land covers at occupied sites would increase extinction risk.

---

## 4. Spatial Patterns

### 4.1 Attribution maps

Six attribution maps per species are available in `figs/attribution/`:
- `fig_attribution_delta_gamma_climate_{sp}.png` -- Climate effect on colonisation
- `fig_attribution_delta_gamma_landuse_{sp}.png` -- Land-use effect on colonisation
- `fig_attribution_delta_epsilon_climate_{sp}.png` -- Climate effect on extinction
- `fig_attribution_delta_epsilon_landuse_{sp}.png` -- Land-use effect on extinction
- `fig_attribution_delta_gamma_total_{sp}.png` -- Total colonisation change
- `fig_attribution_delta_epsilon_total_{sp}.png` -- Total extinction change

The four-panel summary figures (`figs/fig_attribution_summary_{sp}.png`) provide a compact comparison of the main attribution maps per species.

### 4.2 Spatial structure

For *O. tarda*, the strongest land-use effects on extinction are concentrated in central and southern Spain, where urban expansion and shrubland loss are most pronounced. Climate effects on colonisation are spatially diffuse, consistent with the broad-scale nature of temperature and precipitation trends.

For *P. orientalis*, the climate-driven extinction increase is spatially concentrated in the Meseta Central and Ebro Valley, where NDVI declines and warming have been most intense. The land-use protective effect (via LC12) is scattered across the major cereal-growing regions.

For *T. tetrax*, the spatial pattern is dominated by isolated sites where LC12 has changed substantially, embedded within a matrix of negligible effects.

---

## 5. Cross-Framework Validation (colext vs stPGOcc)

### 5.1 stPGOcc turnover hotspots

We derived colonisation and extinction hotspot rasters from the fitted stPGOcc (spatio-temporal occupancy) models by computing posterior mean transition rates from the latent occupancy states ($z$). For each MCMC sample, colonisation at site $i$, year $t$ was defined as $z_{i,t} \times (1 - z_{i,t-1})$, and extinction as $(1 - z_{i,t}) \times z_{i,t-1}$, averaged across posterior samples and years.

### 5.2 Correlations between attribution and hotspots

| Species | $r(\Delta\gamma_{\text{climate}}, \text{col hotspot})$ | $r(\Delta\varepsilon_{\text{climate}}, \text{ext hotspot})$ | $r(\Delta\gamma_{\text{landuse}}, \text{col hotspot})$ | $r(\Delta\varepsilon_{\text{landuse}}, \text{ext hotspot})$ |
|---------|---|---|---|---|
| *O. tarda* | -0.103 | -0.146 | NA | 0.002 |
| *P. alchata* | -0.023 | 0.020 | NA | NA |
| *P. orientalis* | -0.038 | 0.016 | -0.014 | 0.016 |
| *T. tetrax* | NA | NA | 0.081 | -0.079 |

*NA = zero variance in the attribution variable (no covariates of that type in model).*

### 5.3 Interpretation

The correlations between colext attribution maps and stPGOcc turnover hotspots are weak ($|r| < 0.15$ for all comparisons). This is expected for two methodological reasons:

1. **Different model frameworks:** The colext model estimates explicit Markov transition parameters ($\gamma$, $\varepsilon$) as functions of covariates, while stPGOcc estimates year-specific latent occupancy states ($z$) with spatial random effects. The two frameworks answer different questions: colext asks "what drives transitions?" while stPGOcc asks "where do transitions occur?"

2. **Spatial resolution mismatch:** The colext attribution is computed at the site level using covariate values, while the stPGOcc hotspots are derived from posterior occupancy states that additionally absorb spatial random effects ($w$). The spatial random effects capture unmeasured drivers that are by definition orthogonal to the measured covariates.

The strongest correlation ($r = -0.146$) is for *O. tarda* climate attribution vs extinction hotspots, suggesting that sites where climate drives extinction increases (positive $\Delta\varepsilon_{\text{climate}}$) tend to overlap weakly with stPGOcc extinction coldspots. This is consistent with the spatial random effects absorbing some of the climate-driven spatial signal.

---

## 6. Methodological Considerations

### 6.1 Sensitivity to NDVI classification

NDVI was classified as a climate-adjacent covariate throughout this analysis, based on its strong coupling to precipitation and temperature at annual timescales in semi-arid Mediterranean environments. However, NDVI also responds to land management (irrigation, crop rotation, grazing intensity), which creates ambiguity in the attribution.

For *O. tarda* and *P. orientalis*, where NDVI enters the colonisation submodel, reclassifying NDVI as a land-use covariate would shift some of the climate attribution to land use. The magnitude of this shift depends on the covariance between NDVI and the other climate covariates (pr, tmmn, tmmx) at each site.

### 6.2 Scaling approach

All counterfactual predictions used the same scaling parameters as the training data (per-column z-scores of the site $\times$ year matrix). For the "frozen" scenarios, 2017 site-level values were standardised using the 2017 scaling parameters, ensuring that the predicted values are on the same scale as the model's training data. This avoids the artifact of applying different standardisation to counterfactual vs observed data.

### 6.3 Absolute vs relative effects

The absolute attribution effects are small ($10^{-3}$ to $10^{-6}$ probability units) because:
1. The 7-year study window captures only the early stages of multi-decadal environmental trends.
2. Colonisation rates are inherently very low for rare species ($\gamma \approx 0.001$--$0.003$), limiting the dynamic range of climate/land-use effects.
3. Within-site covariate changes over 7 years are small relative to the cross-site variation used to estimate model coefficients.

These small absolute effects do not imply that climate and land-use change are unimportant. Rather, they indicate that the *cumulative* impact will be appreciable only over longer time horizons (decades), consistent with the long generation times and slow population dynamics of these large-bodied steppe birds.

### 6.4 Limitations

- **Model structure:** The colext model assumes a first-order Markov process, which may not capture lagged or threshold effects of environmental change.
- **Sample size:** For *P. alchata* and *T. tetrax*, the small number of observed transitions limits coefficient precision and, consequently, the reliability of attribution estimates.
- **Additive decomposition:** The factorial decomposition assumes that the four scenarios span the relevant space of environmental variation. Unmodelled drivers (e.g., agricultural policy changes, predation pressure) are not captured.

---

## 7. File Outputs

| File | Description |
|------|-------------|
| `results/counterfactual_predictions.rds` | Full prediction arrays (4 species $\times$ 4 scenarios $\times$ 7 years $\times$ all sites) |
| `results/pub_table_attribution_summary.csv` | Cross-species summary table |
| `results/attribution_interpretations.txt` | Auto-generated species paragraphs |
| `results/attribution_hotspot_correlation.csv` | colext--stPGOcc correlation table |
| `results/covariate_trends/trend_*.tif` | Per-covariate, per-species trend rasters (GeoTIFF) |
| `results/stpgocc_col_hotspots_*.tif` | stPGOcc colonisation hotspot rasters |
| `results/stpgocc_ext_hotspots_*.tif` | stPGOcc extinction hotspot rasters |
| `figs/attribution/fig_attribution_*.png` | Individual attribution maps (24 total, 300 dpi) |
| `figs/fig_attribution_summary_*.png` | Four-panel summary figures per species |

---

## 8. Reproducibility

The analysis is fully reproducible from `scripts/8_counterfactual_attribution.R`, which:
1. Loads the fitted colext model objects from `results/{sp}_model_object.rds`
2. Loads the raw covariate data from `data/processed_2023/{sp}/{sp}_occ_wide_dynamic.csv`
3. Recomputes the training scaling internally (matching `scripts/4_occupancy_models.R`)
4. Generates all predictions, attribution decompositions, maps, and tables

The stPGOcc hotspot rasters are generated separately by loading the posterior samples from `results/results_spatial/{sp}_spatial_stPGOcc.rds` and computing transition rates from the latent occupancy states.

**R session requirements:** `unmarked`, `terra`, `sf`, `ggplot2`, `dplyr`, `tidyr`, `here`, `rnaturalearth`, `rnaturalearthdata`. Optional: `rmapshaper` (for island filtering; falls back to bounding-box crop if unavailable).
