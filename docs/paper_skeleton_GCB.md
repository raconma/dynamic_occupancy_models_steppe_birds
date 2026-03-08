# Paper Skeleton: Steppe Bird Occupancy Dynamics (GCB/Ecography)

## Working Title Options

**Option A (driver-attribution focus):**
"Climate and land-use change contribute asymmetrically to occupancy dynamics of Iberian steppe birds: a multi-framework attribution analysis using citizen science data"

**Option B (spatial-temporal focus):**
"Spatial range, temporal dynamics, and driver attribution in four declining steppe bird species: integrating dynamic occupancy models with citizen science at scale"

**Option C (conservation focus):**
"Low colonisation rates constrain range recovery in Iberian steppe birds despite widespread habitat availability: insights from seven years of citizen science monitoring"

---

## Narrative Arc

The core story in one paragraph: **Steppe birds in the Iberian Peninsula are declining, and understanding whether climate or land-use change drives their occupancy dynamics is critical for conservation. Using seven years of eBird citizen science data (2017-2023) across >4000 5-km cells for four species, we fit complementary dynamic occupancy frameworks -- one estimating explicit colonisation/extinction as functions of environmental covariates (colext), the other resolving spatial autocorrelation via nearest-neighbor Gaussian processes (stPGOcc). A factorial counterfactual analysis attributes the observed dynamics to climate vs land-use drivers. We find that colonisation rates are uniformly low (<0.3% per site-year) across all species, extinction rates are spatially structured but driven by different factors per species, and the 7-year climate/land-use trends produce only marginal changes in turnover. The dominant message is that these species' ranges are inertial: habitat loss increases extinction locally, but recolonisation is too slow to compensate, regardless of environmental conditions.**

---

## Structure

### TITLE PAGE
- Title, authors, affiliations
- Running head (~45 characters)
- Keywords (8-10): citizen science, dynamic occupancy, colonisation-extinction, steppe birds, land-use change, climate attribution, spatial autocorrelation, Iberian Peninsula, eBird, conservation biogeography

---

### ABSTRACT (250 words max)

Structure: Context > Gap > Approach > Key Results > Implications

1. **Context** (2 sentences): Steppe habitats in the Iberian Peninsula support globally important populations of threatened bird species, but face dual pressures from agricultural intensification and climate change. Quantifying the relative contribution of each driver to population-level occupancy dynamics is essential for targeted conservation.

2. **Gap** (1 sentence): Existing assessments of steppe bird declines rely on trend estimates that identify *where* and *how fast* populations change, but not *why* -- i.e., whether climate trends or land-use change are the proximate drivers of colonisation and extinction at the landscape scale.

3. **Approach** (3 sentences): We analysed seven years (2017-2023) of eBird citizen science data for four steppe bird species (*Otis tarda*, *Pterocles alchata*, *P. orientalis*, *Tetrax tetrax*) across >4000 5-km grid cells spanning the Iberian Peninsula. We integrated two complementary occupancy modelling frameworks: dynamic occupancy models (colext) estimating colonisation and extinction as explicit functions of climate and land-use covariates, and spatio-temporal Bayesian occupancy models (stPGOcc) incorporating nearest-neighbour Gaussian process spatial random effects. A factorial counterfactual analysis decomposed the observed occupancy dynamics into contributions from climate trends, land-use change, and their interaction.

4. **Key results** (3 sentences): Extinction rates exceeded colonisation rates by 200–1,000,000× (bootstrap ε/γ medians; P(ε/γ > 100) = 81–99%), creating an extinction debt of 5–100% of current occupancy and recolonisation timescales of centuries to millions of years. Climate trends dominated attribution for three species while land-use change drove *T. tetrax*; NDVI decomposition sensitivity confirmed these classifications, with the climate-driven NDVI component amplifying extinction effects for *P. orientalis*. Detection correction revealed that naive colonisation estimates overestimate true rates by 5–22,000×, while calibration slopes (0.58–0.62) indicate moderate prediction overdispersion validated against an independent atlas.

5. **Implications** (2 sentences): The extreme asymmetry between colonisation and extinction implies that range contractions are effectively irreversible: even achieving 10% equilibrium occupancy would require 22–117,000× increases in colonisation, far beyond any management intervention. Conservation must prioritise preventing extinctions at currently occupied sites -- particularly by maintaining cropland mosaics and limiting urban expansion -- rather than relying on natural recolonisation of lost habitat.

---

### 1. INTRODUCTION

**P1 -- The global steppe bird crisis** (~150 words)
- Steppe/grassland birds among most declining bird groups globally
- European farmland bird index: dramatic declines since 1980s
- Iberian Peninsula as last stronghold for several species (O. tarda, T. tetrax)
- Dual threats: agricultural intensification + climate change
- IUCN status of study species

**P2 -- The attribution gap** (~150 words)
- Trend monitoring (atlas data, BBS indices) documents *that* populations decline
- But conservation action requires knowing *why*: climate or land use?
- Standard approaches: correlative SDMs, trend analysis
- Limitation: correlative models cannot decompose observed change into driver contributions
- Need for mechanistic, counterfactual framework

**P3 -- Dynamic occupancy as a solution** (~150 words)
- MacKenzie et al. (2003) colext framework: explicit Markov transitions (gamma, epsilon)
- Covariates on transitions -> mechanistic link between drivers and demographic processes
- Advantage over static SDMs: captures temporal dynamics, accounts for imperfect detection
- But: no spatial structure -> residual spatial autocorrelation biases inference
- Recent advances: spatio-temporal Bayesian models (spOccupancy; Doser et al. 2022)

**P4 -- Citizen science enables landscape-scale inference** (~150 words)
- Traditional surveys cannot cover >4000 sites x 7 years x 4 species
- eBird provides unprecedented spatial and temporal coverage
- Key challenge: variable effort, detection heterogeneity
- Dynamic occupancy models handle this explicitly via detection submodel
- Previous applications of eBird to occupancy modelling (cite relevant papers)
- Our contribution: first multi-species counterfactual attribution analysis using citizen science occupancy data at this scale

**P5 -- Objectives** (~100 words)
- Objective 1: Fit dynamic occupancy models (colext) for four steppe bird species, identifying the environmental covariates that drive colonisation and extinction across the Iberian Peninsula.
- Objective 2: Quantify the spatial structure in occupancy dynamics using spatio-temporal Bayesian models (stPGOcc) with nearest-neighbour Gaussian processes.
- Objective 3: Decompose the observed 2017-2023 occupancy dynamics into contributions from climate trends, land-use change, and their interaction using a factorial counterfactual framework.
- Objective 4: Assess cross-framework consistency between the covariate-driven attribution (colext) and the spatially-explicit turnover patterns (stPGOcc).

---

### 2. METHODS

**2.1 Study area and species** (~200 words)
- Iberian Peninsula, mainland Spain
- 5-km grid cells covering the species' potential range
- Four species, ecological notes (habitat, diet, threats)
- Table 1: Species overview (IUCN status, population estimates, habitat requirements)

**2.2 eBird data and filtering** (~200 words)
- eBird Basic Dataset, 2017-2023
- Filtering protocol (auk): complete checklists, stationary/travelling, duration 5-300 min, distance <10 km
- Zero-filling protocol
- Detection history: 10 replicate visits per site per year
- Final sample sizes per species (Table S1)
- Discussion of citizen science data quality and spatial/temporal coverage
- Figure 1: Map of study area with site locations and survey effort distribution

**2.3 Environmental covariates** (~250 words)

*2.3.1 Static covariates (initial occupancy)*
- WorldClim bioclimatic variables (bio1: annual mean temperature, bio2: temperature range)
- MODIS tree cover, grass cover
- Topographic variables (elevation, aspect) from DEM
- All standardised (z-scored) prior to analysis

*2.3.2 Dynamic covariates (colonisation/extinction)*
- Climate: NDVI (MODIS MOD13A3), precipitation, min/max temperature (TerraClimate)
- Land use: MODIS MCD12Q1 land-cover proportions (Classes 6, 7, 12, 13)
- All extracted via Google Earth Engine at annual resolution, breeding season means
- Table 2: Covariate descriptions with sources and temporal resolution

**2.4 Dynamic occupancy models (colext)** (~500 words)

*2.4.1 Model formulation*
- MacKenzie et al. (2003) multi-season framework implemented via `unmarked::colext()`
- State process: z_{i,1} ~ Bernoulli(psi_i); z_{i,t} | z_{i,t-1} ~ Bernoulli(z_{i,t-1}*(1-epsilon_{i,t}) + (1-z_{i,t-1})*gamma_{i,t})
- Observation process: y_{i,j,t} | z_{i,t} ~ Bernoulli(z_{i,t} * p_{i,j,t})
- Four linked submodels estimated via maximum likelihood:
  - psi (initial occupancy): site-level static covariates (bioclimatic, land cover, topography)
  - gamma (colonisation probability): yearly site-level dynamic covariates
  - epsilon (extinction probability): yearly site-level dynamic covariates
  - p (detection probability): observation-level covariates (effort, observers, time, weather)
- All continuous covariates z-standardised per year using training-set means and SDs
- Logit link for all submodels

*2.4.2 Model selection*
- Species-specific covariate sets selected via AIC-based forward/backward selection
- Candidate covariates: climate (NDVI, pr, tmmn, tmmx) and land use (LC classes 0, 6, 7, 12, 13) for gamma/epsilon
- Convergence diagnostics: complete separation detected in P. alchata (Class_0 removed from gamma/epsilon), T. tetrax (pr, tmmx removed; only LC12 retained)
- Final formulas per species in Table 3
- Full AIC model comparison tables in Supplementary S2

*2.4.3 Derived quantities*
- Annual occupancy prevalence: psi_t = sum(z_{i,t}) / N, derived via stochastic simulation from fitted gamma/epsilon
- Turnover rate: tau_t = (col_t + ext_t) / N
- Site-level predicted colonisation and extinction maps using `predict()` with fitted model and yearly covariate surfaces
- Confidence intervals via delta method (SE from Hessian)

*2.4.4 Model validation*
- Goodness of fit: parametric bootstrap (parboot, n=500 replicates)
- MacKenzie-Bailey chi-squared test
- Spatial block cross-validation (5-fold, blockCV) against Spanish Biodiversity Atlas
- Metrics: AUC, TSS, RMSE
- Detection model diagnostics: effect of effort, observers, time-of-day on p

**2.5 Spatio-temporal occupancy models (stPGOcc)** (~250 words)

*2.5.1 Model structure*
- spOccupancy::stPGOcc with AR(1) temporal random effects
- Nearest-neighbour Gaussian Process (NNGP, m=5 neighbours)
- Exponential covariance with estimated spatial range (phi)
- Same covariates as colext for comparability

*2.5.2 MCMC settings*
- 3 chains, 20,000 iterations, burn-in 10,000, thin 20
- Convergence: Rhat, ESS, traceplots
- Spatial diagnostics: Moran's I on residuals, effective spatial range

**2.6 Counterfactual attribution analysis** (~300 words)

*2.6.1 Covariate trend quantification*
- Linear regression per site: covariate ~ year (2017-2023)
- Fraction of sites with significant positive/negative trends (P < 0.05)

*2.6.2 Factorial scenario design*
- Four scenarios (S0-S3): null baseline, climate only, land use only, combined
- "Frozen" = 2017 site-level values retained for all years
- Predictions using training scaling parameters

*2.6.3 Attribution decomposition*
- Delta_gamma_climate, delta_gamma_landuse, delta_gamma_interaction, delta_gamma_total
- Same for epsilon
- Mean and SD across sites and years
- Figure concept: "Attribution decomposition diagram" showing the factorial logic

*2.6.4 NDVI decomposition sensitivity*
- NDVI is ~50% climate-driven (mean R² = 0.51 from lm(NDVI ~ pr + tmmn + tmmx))
- Site-level regression yields NDVI_climate (fitted) and NDVI_residual (residuals)
- Each component z-scored separately using training-set statistics
- Affected submodels: O. tarda gamma, P. orientalis gamma/epsilon
- Attribution recomputed with Climate pathway = (pr, tmmn, tmmx, NDVI_climate) and Land-use pathway = (LC vars, NDVI_residual)
- Provides sensitivity check on dominant-driver classification

**2.7 Demographic asymmetry quantification** (~200 words)

*2.7.1 ε/γ ratio with bootstrap CIs*
- Parametric bootstrap (n = 5,000 draws from mvrnorm(coef, vcov))
- γ = plogis(col_intercept), ε = plogis(ext_intercept) at mean covariates
- Report median ratio with 95% percentile CIs
- Report P(ε/γ > 100) and P(ε/γ > 1000) as measures of demographic asymmetry

*2.7.2 Delta-γ analysis (colonisation multiplier)*
- γ_required = ψ*_target × ε / (1 − ψ*_target)
- Multiplier = γ_required / γ_current
- Computed for ψ* targets of 5% and 10%
- Bootstrap CIs from the same 5,000 draws
- Quantifies how much colonisation must increase for population recovery

*2.7.3 Naive vs detection-corrected comparison*
- Naive γ/ε from raw detection histories (occupied = detected ≥ 1 time)
- Corrected γ/ε from fitted colext models with bootstrap CIs
- Ratio = naive/corrected quantifies detection bias

**2.8 Model validation** (~200 words)

*2.8.1 Calibration curves*
- Spatially-blocked 5-fold cross-validation (blockCV, 270 km blocks)
- Calibration slope estimated on logit scale via logistic regression
- H₀: slope = 1 (perfect calibration)
- Validated against independent Spanish Biodiversity Atlas

**2.9 Cross-framework validation** (~150 words)
- stPGOcc-derived colonisation/extinction hotspots from posterior z-samples
- Pearson correlation between attribution maps and stPGOcc hotspots
- Interpretation: high correlation = drivers explain spatial turnover; low = unmeasured factors dominate

---

### 3. RESULTS

**3.1 Survey coverage and detection** (~150 words)
- Summary statistics: sites, visits, detections per species
- Detection probabilities (mean, range)
- Effects of effort covariates on detection
- Figure S1: Detection probability by effort/observers/time

**3.2 Dynamic occupancy model results (colext)** (~600 words)

*3.2.1 Detection*
- Detection probabilities ranged from 0.08 (T. tetrax) to 0.35 (O. tarda)
- Effort (distance/duration) and number of observers consistently positive effects across all species
- Time-of-day significant for P. alchata and P. orientalis (dawn/dusk activity patterns)
- Weather covariates (NDVI_obs, pr_obs) improved detection for O. tarda and P. alchata
- Table 4: Detection model coefficients across species

*3.2.2 Initial occupancy*
- All species: strong negative effect of tree cover, positive of grassland
- Temperature and topography define range limits (bio1, bio2, elevation, aspect)
- Spatial predictions match known range boundaries
- Figure 2: Initial occupancy probability maps (4-panel, one per species)

*3.2.3 Colonisation*
- Extremely low baseline rates for all species (gamma < 0.003)
- O. tarda: negative NDVI and tmmn effects (P < 0.05); warmer/drier conditions reduce already-low colonisation
- P. alchata: weak effects of NDVI and precipitation, wide CIs due to only 16 colonisation events
- P. orientalis: strong opposing temperature effects (tmmn negative, tmmx positive); NDVI positive
- T. tetrax: only cropland proportion (LC12) significant -- purely land-use-driven
- Figure 3: Coefficient forest plot for colonisation across all 4 species (standardised betas +/- 95% CI)
- Table 5: Colonisation model coefficients with SE, z-values, and P-values

*3.2.4 Extinction*
- Moderate to high rates (0.15-0.48), contrasting sharply with near-zero colonisation
- O. tarda: urban cover (LC13) strongly increases extinction (beta = 2.80, P < 0.001); also tmmx positive
- P. alchata: precipitation reduces extinction (beta = -0.91, P < 0.05); temperature increases it
- P. orientalis: NDVI (beta = -1.21) and precipitation (beta = -0.87) reduce extinction risk; LC12 marginal
- T. tetrax: cropland reduces extinction (beta = -0.57, P = 0.042) -- protective effect of agricultural mosaics
- Figure 4: Coefficient forest plot for extinction across all 4 species
- Table 6: Extinction model coefficients with SE, z-values, and P-values

*3.2.5 Predicted colonisation and extinction maps*
- Site-level predicted gamma and epsilon for each year (2017-2023)
- Colonisation spatially homogeneous and uniformly near zero for all species
- Extinction spatially structured: concentrated at range margins and in areas with high land-use change
- Figure 5: Predicted colonisation and extinction probability maps (4x2 panel: species x process, averaged over 2017-2023)

*3.2.6 Occupancy dynamics and prevalence trends*
- Stochastic simulation from fitted gamma/epsilon -> derived annual prevalence
- All species stable or slightly declining over 2017-2023
- The colonisation-extinction asymmetry (gamma << epsilon) as the defining dynamical feature
- Figure 6: Prevalence trends 2017-2023 with confidence intervals (4-panel)

*3.2.7 Model fit and validation*
- Parametric bootstrap: adequate fit for all species (P > 0.10)
- Cross-validation AUC: 0.72-0.89 across species
- Validation against Spanish Biodiversity Atlas: TSS = 0.45-0.67
- Table 7: Model fit statistics (chi-squared, P-value, AUC, TSS, RMSE per species)

**3.3 Spatial structure (stPGOcc)** (~300 words)
- Substantial spatial autocorrelation in colext residuals (Moran's I = 0.14-0.24)
- stPGOcc reduces autocorrelation by 59-82%
- Effective spatial ranges: 23-334 km (species-dependent)
- Spatial random effects reveal unmeasured habitat quality gradients not captured by colext covariates
- Coefficient estimates broadly consistent between colext and stPGOcc (direction and significance preserved)
- Figure 7: Spatial diagnostics (residual maps, Moran's I comparison, effective spatial range estimates)
- Table 8: Model comparison (tPGOcc vs stPGOcc): Moran's I, spatial range, MCMC diagnostics

**3.4 Covariate trends 2017-2023** (~200 words)
- Climate: consistent tmmn warming (+1.2 C/yr), weak pr decline, spatially heterogeneous NDVI decline
- Land use: bidirectional cropland change (9-10% gaining, 6% losing), slow urban expansion (+0.004 pp/yr)
- Few individual site trends are significant over the 7-year window
- Figure 8: Trend maps for key covariates (6-panel: NDVI, pr, tmmn, tmmx, LC12, LC13)

**3.5 Demographic asymmetry and extinction debt** (~400 words)
- ε/γ ratio quantifies the fundamental demographic trap: median ratios 201–1,054,932 across species
- Bootstrap 95% CIs confirm: P(ε/γ > 100) = 81–99% for all species
- Delta-γ analysis: even modest recovery (ψ* = 10%) requires 22× (T. tetrax) to 117,215× (P. orientalis) increase in colonisation
- Extinction debt: 5% (T. tetrax) to 100% (P. orientalis) of current sites are transient
- Figure: Isocline plot (γ vs ε on log-log axes) with equilibrium contours + extinction debt bars
  - `figs/pub_fig_isocline_equilibrium.png` (KEY FIGURE)
- Table: ε/γ ratios, equilibrium ψ*, recolonisation timescales, extinction debt fractions

**3.6 Counterfactual attribution** (~400 words)
- Based on colext fitted models: predictions under factorial frozen/observed covariate scenarios
- Effects are small in absolute terms (<0.01 probability units) but informative in relative terms
- Climate dominates for O. tarda, P. alchata, P. orientalis; land use for T. tetrax
- Interaction terms negligible (<5x10^-4 for all species)
- NDVI decomposition sensitivity: decomposition amplifies climate signal for P. orientalis epsilon (3×) but does not change dominant-driver classification
- Spatial patterns: climate effects diffuse, land-use effects concentrated at transition zones
- Table: Cross-species attribution summary (original + decomposed NDVI comparison)

**3.7 Detection bias and model validation** (~300 words)
- Naive colonisation overestimates corrected by 5–22,000×; largest for P. orientalis
- Naive extinction typically lower than corrected (ratios 0.56–2.6)
- Calibration slopes: 0.575–0.623 (all significantly < 1), indicating overdispersed predictions
- Despite miscalibration, discrimination is good (AUC > 0.82)
- Temporal mismatch (ψ₁ 2017 vs multi-year atlas) partly explains calibration bias
- Figure: Calibration curves (4 panels, one per species)
  - `figs/{sp}_calibration_curve.png`
- Table: Calibration slopes with SEs and p-values

**3.8 Cross-framework consistency** (~150 words)
- Weak correlations (|r| < 0.15) between colext-derived attribution maps and stPGOcc-derived turnover hotspots
- Interpretation: spatial random effects capture drivers orthogonal to measured covariates
- This is itself informative: measured climate and land-use covariates explain *temporal* variation (via colext) but not *spatial* variation in turnover
- Table: Correlation matrix

---

### 4. DISCUSSION

**4.1 The colonisation-extinction asymmetry** (~400 words)
- **Key finding**: Extinction rates exceed colonisation by 200–1,000,000× (ε/γ bootstrap medians with 95% CIs)
- Probability that ε/γ > 100 ranges from 81% (P. alchata) to 99% (P. orientalis)
- Delta-γ analysis: achieving even 10% equilibrium occupancy requires 22–117,000× more colonisation
- This asymmetry is the dominant feature of steppe bird dynamics
- Implications: range contractions are effectively irreversible on decadal timescales
- Naive estimates mask the true severity: detection correction reveals colonisation 5–22,000× lower than apparent
- Comparison with other taxa showing similar asymmetries (large-bodied birds, habitat specialists)
- Link to dispersal limitation, Allee effects, conspecific attraction in steppe birds
- This paragraph should be the "memorable finding" of the paper

**4.2 Dynamic occupancy models as a mechanistic framework** (~300 words)
- colext provides explicit demographic parameters (gamma, epsilon) rather than correlative projections
- Key advantage: covariates act on *processes* (colonisation, extinction), not just on static probability
- This enables the counterfactual decomposition -- not possible with SDMs or simple trend models
- Species-specific model selection reveals different ecological mechanisms per species
- Limitations: maximum likelihood estimation sensitive to complete separation with sparse data (ptealc: 16 colonisation events; tettet: 25 transitions from occupied sites)
- No spatial structure in colext -> residual Moran's I in predicted occupancy (motivates stPGOcc)
- Comparison with other dynamic occupancy applications in European bird conservation
- The colext + counterfactual combination is a transferable framework applicable to any multi-season occupancy dataset

**4.3 Climate vs land-use attribution** (~350 words)
- Climate effects are larger for 3/4 species, but driven by different mechanisms per species
- O. tarda: warming + NDVI decline reduces colonisation; urbanisation drives extinction
- P. orientalis: temperature effects create opposing colonisation signals; NDVI decline increases extinction
- T. tetrax: purely land-use driven (model structure constraint, but ecologically valid)
- The small absolute effects over 7 years should not be dismissed -- extrapolation suggests meaningful cumulative impacts over decades
- Comparison with other attribution studies in European birds

**4.4 The value of citizen science for dynamic occupancy modelling** (~300 words)
- eBird coverage enables analysis at scales impossible with conventional surveys
- 3700-4100 sites x 7 years x 10 replicates per species
- Detection modelling accounts for heterogeneous observer effort
- Comparison with atlas-based approaches: continuous monitoring vs snapshots
- Limitations: spatial and temporal coverage gaps, urban/accessible area bias
- Validation against independent atlas data provides confidence in predictions
- Recommendation for integrating citizen science into national monitoring frameworks

**4.5 Spatial structure and its implications** (~300 words)
- Substantial spatial autocorrelation (Moran's I = 0.14-0.24) indicates unmeasured spatial drivers
- stPGOcc spatial random effects reduce this by 59-82%
- Effective spatial ranges (23-334 km) suggest regional-scale processes (landscape connectivity, metapopulation dynamics)
- The weak correlation between attribution maps and stPGOcc hotspots confirms that measured covariates explain *temporal* variation but not *spatial* variation
- Implication: site-level conservation actions must account for spatial context beyond the measured environmental gradients
- Future directions: spatially-varying coefficient models (Doser et al. 2024)

**4.6 Conservation implications** (~250 words)
- **Priority 1**: Prevent extinction at occupied sites (especially from urbanisation for O. tarda, cropland loss for T. tetrax)
- **Priority 2**: Maintain cropland mosaics in steppe regions (protective effect for T. tetrax and P. orientalis)
- **Priority 3**: Do not rely on natural recolonisation -- colonisation rates are too low
- **Priority 4**: Climate adaptation will become increasingly important as trends strengthen over coming decades
- Policy relevance: EU Common Agricultural Policy, Natura 2000 network, IUCN Red List assessments
- Monitoring recommendation: continue eBird-based citizen science monitoring to track attribution effects over longer time horizons

**4.7 Limitations and future directions** (~200 words)
- 7-year window captures only early-stage trends (stronger signals expected over decades)
- NDVI classification as climate-adjacent is debatable
- colext assumes first-order Markov process (no lag effects)
- Convergence issues for rare species (P. alchata, T. tetrax) limit coefficient precision
- No spatial structure in colext likelihood -- stPGOcc addresses this but cannot estimate explicit transitions
- Future: longer time series, spatial dynamic models with explicit transitions (when software matures), integration with population viability analysis

---

### 5. CONCLUSIONS (~150 words)
- Citizen science data enables landscape-scale inference about steppe bird dynamics
- Colonisation-extinction asymmetry is the dominant finding
- Climate effects dominate over land-use effects for most species over 2017-2023
- Spatial structure reveals unmeasured drivers beyond climate and land use
- Conservation priority: prevent extinctions, don't wait for recolonisation

---

### FIGURES (10 main + supplements)

1. **Fig 1**: Study area map with site locations, coloured by survey effort (`figs/pub_map_main_figure.png`)
2. **Fig 2**: Initial occupancy probability maps (4-panel) (`figs/pub_map_occupancy_4species.png`)
3. **Fig 3**: All submodel coefficients (γ + ε + ψ₁) forest plot (`figs/pub_fig_coef_all_submodels.png`)
4. **Fig 4**: Detection probability diagnostics (`figs/pub_fig4_detection_comparison.png`)
5. **Fig 5**: Prevalence trends 2017-2023 (`figs/pub_fig_occupancy_trends_panel.png`)
6. **Fig 6 (KEY)**: Isocline plot (γ vs ε log-log) + extinction debt bars (`figs/pub_fig_isocline_equilibrium.png`)
7. **Fig 7**: Spatial diagnostics (Moran's I reduction, spatial range) (`figs/pub_fig_spatial_moran.png`)
8. **Fig 8**: Covariate trend maps (NDVI, pr, tmmn, tmmx, LC12, LC13)
9. **Fig 9**: Attribution summary (4-species, original + NDVI-decomposed comparison)
10. **Fig 10**: Calibration curves (4 panels) (`figs/{sp}_calibration_curve.png`)

### TABLES (13 main)

1. Study species overview (IUCN status, population estimates, habitat)
2. Covariate descriptions with sources and temporal resolution
3. Model formulas per species (psi, gamma, epsilon, p)
4. Detection model coefficients across species
5. Colonisation model coefficients (beta, SE, z, P)
6. Extinction model coefficients (beta, SE, z, P)
7. **ε/γ ratio with bootstrap CIs** (median, 95% CI, P(>100), P(>1000)) — `results/ratio_bootstrap.csv`
8. **Equilibrium occupancy, extinction debt, recolonisation timescales** — `results/equilibrium_occupancy_table.csv`, `results/extinction_debt_table.csv`
9. **Delta-γ colonisation multiplier** (for ψ*=5% and 10%) — `results/delta_gamma.csv`
10. Attribution summary (cross-species delta-gamma, delta-epsilon) — `results/attribution_table3.csv`
11. **NDVI decomposition comparison** (original vs decomposed attribution) — `results/attribution_comparison_ndvi.csv`
12. **Naive vs corrected transition rates** — `results/naive_vs_corrected_full.csv`
13. **Calibration slopes** (270 km blocks) — `results/calibration_slopes.csv`

### SUPPLEMENTARY

- S1: Full detection model results with interaction effects
- S2: AIC model selection tables for all species (all candidate models)
- S3: Goodness-of-fit diagnostics (bootstrap distributions, chi-sq plots)
- S4: Spatial cross-validation results (AUC, TSS maps by fold)
- S5: Full covariate trend statistics per site
- S6: Individual attribution maps (all 24)
- S7: MCMC diagnostics for stPGOcc (traceplots, Rhat, ESS)
- S8: **NDVI decomposition sensitivity** (decomposed coefficients, model comparison) — `results/attribution_ndvi_decomposed_table.csv`
- S9: Complete separation diagnostics for P. alchata and T. tetrax
- S10: **Naive vs corrected rates** (full table with n events/opportunities) — `results/naive_vs_corrected_full.csv`
- S11: Spatial model comparison (tPGOcc vs stPGOcc: Moran's I, spatial range)
- S12: Cross-framework correlation matrix
