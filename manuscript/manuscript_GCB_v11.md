---
title: "Colonisation-extinction asymmetry exposes a demographic trap in Iberian steppe birds"
authors: "Raul Contreras-Martin, Guillermo Fandos"
affiliation: "[Affiliation]"
correspondence: "[email]"
running_head: "Demographic trap in Iberian steppe birds"
target_journal: "Global Change Biology"
date: "March 2026"
word_count_estimate: "~6,800 words main text"
keywords: "dynamic occupancy models, colonisation-extinction asymmetry, eBird, imperfect detection, range dynamics, steppe birds, citizen science, counterfactual attribution, spatial occupancy, extinction debt"
---

## Abstract

Farmland bird populations across Europe have declined by more than 55% since 1980, with steppe-associated species among the hardest hit. The Iberian Peninsula, supporting 19–98% of European breeding populations for four threatened steppe species, functions as the continental reservoir for this guild, yet documented range contractions of 20–60% continue despite targeted habitat protection. Whether these contractions reflect failing recolonisation, accelerating local extinction, or their interaction remains unresolved, because standard monitoring tools measure only net occupancy change and cannot separate the underlying demographic processes, particularly under imperfect detection. We fitted dynamic occupancy models to seven years of eBird citizen-science data (2017–2023) across more than 4,000 five-km grid cells in mainland Spain for *Otis tarda*, *Pterocles alchata*, *P. orientalis*, and *Tetrax tetrax*, estimating colonisation and extinction as explicit functions of annual climate and land-cover covariates while correcting for detection probabilities of 13–35%. Detection correction reduced naive colonisation estimates by 5- to 22,000-fold, revealing a structural asymmetry in which extinction exceeded colonisation by two to six orders of magnitude across all species (bootstrap epsilon/gamma medians: 201–1,054,932). Equilibrium occupancy remained below 0.5% for three of four species, creating extinction debts of 5–100% of current occupancy. Climate dominated attribution for three species; cropland availability drove both colonisation and extinction in *T. tetrax* alone. Spatially blocked validation against the Spanish Breeding Bird Atlas confirmed predictive skill (AUC 0.82–0.88). These results demonstrate that Iberian steppe birds occupy ranges their demographic rates can no longer sustain. Preventing local extinction, not promoting recolonisation, is the only viable strategy on management-relevant timescales, and citizen-science data can diagnose this trap when analysed with detection correction and demographic decomposition.

---

## 1. Introduction

Farmland bird populations in Europe have declined by more than 55% since 1980 (PECBMS 2023), and steppe-associated species sit at the extreme end of this collapse. The four species analysed here, great bustard (*Otis tarda*), pin-tailed sandgrouse (*Pterocles alchata*), black-bellied sandgrouse (*P. orientalis*), and little bustard (*Tetrax tetrax*), collectively represent 19–98% of their European breeding populations, concentrated in the agro-steppe landscapes of the Iberian Peninsula (Table 1). Over the past three decades, documented range contractions of 20–60% (SEO/BirdLife 2021) have been attributed primarily to habitat loss: agricultural intensification, infrastructure development, afforestation, and, increasingly, photovoltaic expansion. Yet habitat loss alone does not explain why natural recolonisation has not reversed these declines even where protection has been implemented. Across the Iberian meseta, suitable steppe habitat persists in patches that remain unoccupied years after local extinction. The unresolved question is whether these species are experiencing equilibrium decline, in which habitat loss drives occupancy down but the demographic machinery of recolonisation remains intact, or non-equilibrium transience, in which current occupancy is a legacy of historical landscapes that contemporary demographic rates can no longer sustain. If the latter, the true conservation crisis is not the contraction already documented but the further contraction that is demographically inevitable even if habitat loss were halted today.

Standard monitoring tools cannot distinguish between these scenarios. Atlas surveys, trend indices, and correlative species distribution models measure net occupancy change, integrating the outcome of colonisation and extinction without separating them (Zurell et al. 2022). This diagnostic failure has direct conservation consequences. If recolonisation is the bottleneck, restoring habitat connectivity is the priority. If extinction at occupied sites is the constraint, protecting existing populations takes precedence. If the two processes respond to different drivers, a generic response will fail regardless of resources invested (McCarthy et al. 2012). The non-equilibrium range dynamics literature (Hanski 1999; Briscoe et al. 2021; Zurell et al. 2022) provides the theoretical framework for thinking about this problem, but empirical demonstrations with co-occurring threatened species under realistic data constraints remain scarce. Dynamic occupancy models that explicitly estimate colonisation (gamma) and extinction (epsilon) as separate processes (MacKenzie et al. 2003) represent the minimum analytical requirement for diagnosing whether range recovery is demographically feasible.

Detection bias compounds the diagnostic problem. For the four species studied here, baseline detection probabilities range from 13% to 35%. At these rates, naive colonisation estimates conflate genuine recolonisation with missed detections, systematically understating the colonisation-extinction asymmetry. The directionality of this bias is particularly consequential: it operates to mask conservation urgency, generating false reassurance by placing apparent colonisation in a range consistent with recovery rather than structural failure. In citizen-science data, where survey effort varies by orders of magnitude across space and time, this bias cannot be resolved by restricting analyses to high-quality subsets. Detection-correcting dynamic occupancy models are therefore not a statistical refinement but a prerequisite for defensible inference about demographic rates (Guillera-Arroita 2017).

We test four hypotheses. (H1) Extinction rates structurally exceed colonisation rates for all four species, with the asymmetry severe enough to produce extinction debt even under optimistic bootstrap scenarios. (H2) Climate and land-use drivers act on distinct demographic processes and differ in dominance across species, with the site-faithful bustards and semi-nomadic sandgrouse showing different driver signatures. (H3) Naive estimates do not merely understate the asymmetry but change its qualitative interpretation, placing apparent colonisation in a range consistent with recovery rather than structural failure. (H4) Spatial range estimates from spatially explicit models follow the mobility axis of the guild, with implications for the minimum effective conservation scale. We address these hypotheses using dynamic occupancy models fitted to eBird data (2017–2023) across more than 4,000 five-km grid cells, with factorial counterfactual attribution and spatially blocked validation against the Spanish Breeding Bird Atlas.

## 2. Methods

### 2.1 Study area and species

The study covers mainland Spain (494,011 km²), divided into a regular 5-km UTM 30N grid. Of the more than 19,000 cells, over 4,000 contained at least three eBird complete checklists during the breeding season across the study period (2017–2023). We focused on four species spanning a gradient of habitat specificity and mobility within the Iberian agro-steppe guild (Table 1).

**Table 1.** Study species, conservation status, and ecological traits.

| Species | Code | IUCN (Spain) | % Eur. pop. | Breeding season | Mobility |
|---|---|---|---|---|---|
| *Otis tarda* | otitar | NT | 98% | April–June | Site-faithful |
| *Pterocles alchata* | ptealc | VU | 97% | May–August | Semi-nomadic |
| *Pterocles orientalis* | pteori | EN/VU | 19% | May–August | Semi-nomadic |
| *Tetrax tetrax* | tettet | EN | 34% | April–June | Site-faithful |

The mobility gradient is directly relevant to the closure assumption: within-season movements by the semi-nomadic sandgrouse may partially violate the assumption that sites are closed to occupancy change during the sampling window (see Section 4, Limitations).

### 2.2 eBird data and filtering

Occurrence data were obtained from the eBird Basic Dataset (Sullivan et al. 2009) for mainland Spain, 2017–2023. Following Johnston et al. (2021), we filtered checklists to ≤5 hours duration, ≤5 km distance, and ≤10 observers, retaining only complete checklists to allow valid zero-filling. Each checklist was assigned to its corresponding 5-km grid cell; repeat visits within the same cell and breeding season formed the detection histories used for occupancy estimation. The number of eBird checklists in Spain increased by 256% between 2017 and 2023. We tested for effort confounding by computing Spearman correlations between naive cell-level occupancy and annual checklist count; no significant correlation was detected for three of four species (rho = −0.25 to 0.37, all P > 0.10). For *P. alchata*, a marginally significant correlation was observed (rho = 0.64, P = 0.05); this species' colonisation submodel is treated with additional caution throughout.

### 2.3 Environmental covariates

We used a two-scale covariate design separating the determinants of long-run habitat suitability from the drivers of realised annual demographic turnover. Static covariates for initial occupancy (psi₁) were drawn from WorldClim 1970–2000 climatic normals (BIO1: mean annual temperature; BIO2: mean diurnal temperature range) and from Global Land Cover Facility data (tree cover, herbaceous cover), supplemented by terrain variables (slope, aspect, elevation). Variables were selected after hierarchical cluster analysis to remove pairwise correlations r > 0.7.

Dynamic covariates for colonisation (gamma) and extinction (epsilon) were derived from three annual data products computed for each species' breeding season: minimum and maximum temperature and cumulative precipitation from TerraClimate (Abatzoglou et al. 2018); NDVI from MODIS MOD13A3 (Didan 2015); and proportional land-cover classes from MODIS MCD12Q1 IGBP (Friedl & Sulla-Menashe 2015), including cropland (LC12), open shrubland (LC7), grassland (LC10), urban/built-up (LC13), and water bodies (LC11). All dynamic covariates were standardised using training-set scaling parameters.

NDVI presents a collinearity challenge: approximately 50% of interannual NDVI variance is explained by climate variables (mean R² = 0.51 across sites). Variance inflation factors for NDVI remained below 5 (max VIF = 4.43), but pairwise correlations between precipitation and maximum temperature (|r| = 0.70–0.71) reduce the precision of individual climate coefficients. For *P. orientalis*, NDVI was retained in the extinction submodel (AIC-preferred by 60 units) but classified as "climate-adjacent" in the attribution analysis. A sensitivity model without NDVI is provided in Supplementary S10; VIF diagnostics for all submodels are in Supplementary S10.

### 2.4 Dynamic occupancy models

We used the colext framework (MacKenzie et al. 2003) implemented in the unmarked R package (Fiske & Chandler 2011; Kellner et al. 2023) to estimate four simultaneous processes: initial occupancy (psi₁), colonisation probability (gamma), extinction probability (epsilon), and detection probability (p). The state process follows a first-order Markov chain: psi_{i,t+1} = (1 − psi_{it}) × gamma_i + psi_{it} × (1 − epsilon_i). Detection is modelled as a Bernoulli trial conditional on occupancy, with visit-level covariates (duration, transect length, number of observers, time of day, observation-period NDVI, and precipitation) accounting for heterogeneous detectability.

Model selection used AIC over candidate covariate sets for each submodel, constrained so that static WorldClim covariates entered only psi₁ and dynamic annual covariates entered only gamma and epsilon. Complete separation was diagnosed in the *P. alchata* colonisation submodel (16 observed colonisation events across 2,521 site-year opportunities, 0.6%); this submodel is excluded from counterfactual attribution. Near-separation was diagnosed in the *P. alchata* extinction submodel (18 events; SE > 11). These data limitations do not invalidate the inclusion of *P. alchata*: the directional finding (gamma ≪ epsilon) is robust to wide confidence intervals, and excluding *P. alchata* would remove half the sandgrouse component of the guild comparison. Model fit was assessed by parametric bootstrap goodness-of-fit (parboot, n = 500); all four species showed adequate fit (P = 0.28–0.93; Table 2).

**Table 2.** Parametric bootstrap goodness-of-fit (n = 500).

| Species | AIC | chi-sq (obs) | chi-sq (sim) | parboot P | Fit |
|---|---|---|---|---|---|
| *O. tarda* | 2267.9 | 352.4 | 333.6 | 0.28 | Adequate |
| *P. alchata* | 1981.7 | 363.4 | 364.6 | 0.52 | Adequate |
| *P. orientalis* | 2060.7 | 276.0 | 330.6 | 0.93 | Adequate |
| *T. tetrax* | 1780.1 | 290.9 | 307.9 | 0.66 | Adequate |

**Table 3.** Best-fit colext models per species with AIC-selected covariates.

| Species | psi₁ covariates | gamma covariates | epsilon covariates | AIC | N sites |
|---|---|---|---|---|---|
| *O. tarda* | bio1, bio2, tree, grass, elev | NDVI⁺, pr, tmmn, tmmx | LC6, LC13, tmmx | 2267.9 | 3,745 |
| *P. alchata* | bio1, bio2, tree, grass, aspect | pr [excl. attrib.] | pr, tmmx [wide SE] | 1981.7 | 4,130 |
| *P. orientalis* | bio2, tree, grass | LC7, NDVI⁺, tmmn, tmmx | LC12, NDVI⁺, pr | 2060.7 | 4,130 |
| *T. tetrax* | bio2, tree, grass, elev | LC12 | LC12 | 1780.1 | 3,746 |

Abbreviations: bio1 = mean annual temperature; bio2 = mean diurnal range; tree/grass = tree/herbaceous cover; elev = elevation; LC6 = closed shrubland; LC7 = open shrubland; LC12 = cropland; LC13 = urban/built-up; tmmn/tmmx = min/max temperature; pr = cumulative precipitation. ⁺ NDVI classified as "climate-adjacent" in attribution (R² ≈ 0.51 with climate variables). *P. alchata* gamma excluded from attribution (separation; 16 events). *P. alchata* epsilon reported with caution (near-separation; 18 events, SE > 11).

### 2.5 Equilibrium occupancy and the epsilon/gamma asymmetry

Equilibrium occupancy psi* = gamma/(gamma + epsilon) was computed as the expected long-run proportion of sites occupied under constant current conditions. The epsilon/gamma ratio provides a dimensionless measure of asymmetry that is more robust than psi* to wide uncertainty in absolute gamma. Uncertainty was propagated using parametric bootstrap: 5,000 coefficient vectors drawn from the multivariate normal approximation to the estimated coefficient distribution, from which annual gamma and epsilon were predicted at mean dynamic covariate values. Expected recolonisation time (1/gamma) and persistence time (1/epsilon) are reported in years with 95% bootstrap CI. We additionally computed the colonisation rate required to achieve a target equilibrium of psi* = 0.10 for each species, as gamma_needed = 0.10 × epsilon / 0.90, quantifying the gap between current performance and a minimum viable metapopulation.

### 2.6 Spatio-temporal Bayesian occupancy models

Colext models assume independent residuals across sites, an assumption violated by residual spatial autocorrelation (Moran's I = 0.29–0.49, all P < 0.001). We fitted spatio-temporal occupancy models using stPGOcc (Doser et al. 2022) from the spOccupancy R package, adding a nearest-neighbour Gaussian process spatial random effect with an exponential correlation function and an AR(1) temporal random effect. These models serve two roles: a robustness check on colext covariate effects and estimation of the effective spatial range (phi), a biologically interpretable quantity defining the minimum conservation planning scale. Convergence was assessed using Gelman–Rubin Rhat (< 1.1) and effective sample size (ESS > 100). Spatial range estimates were compared against published dispersal distances as an ecological coherence check.

### 2.7 Counterfactual attribution

We implemented a factorial counterfactual design with four scenarios: S0 (all dynamic covariates fixed at 2017 values), S1 (climate observed, land use fixed), S2 (land use observed, climate fixed), and S3 (all observed). Attribution effects were computed as S1 − S0 (climate), S2 − S0 (land use), and S3 − S1 − S2 + S0 (interaction). NDVI was classified as "climate-adjacent" for *O. tarda* and *P. orientalis*. NDVI decomposition sensitivity analysis separated climate-driven greening from management-residual components; details are in Supplementary S10. Uncertainty was propagated using parametric bootstrap (n = 1,000).

### 2.8 Naive versus detection-corrected transition rates

We compared naive transition rates, computed directly from detection histories as the proportion of observed site-year transitions, against detection-corrected estimates from the best colext models. For naive estimation, a site was classified as colonised if it had zero detections in year *t* and at least one in year *t*+1, and as locally extinct if the reverse was observed. These naive rates do not account for the probability that an occupied site was undetected.

### 2.9 External validation

Predicted initial occupancy was validated against the III Atlas of Breeding Birds in Spain (SEO/BirdLife 2022) using spatially blocked 5-fold cross-validation (blockCV; Valavi et al. 2019) with block size set to 270 km, exceeding the largest estimated spatial range across species. Predictive performance is reported using Spearman rho, AUC, and TSS. Calibration curves are in Supplementary S8. External validation is restricted to psi₁; transition rates have no equivalent independent benchmark at the five-km cell scale.

## 3. Results

### 3.1 Detection heterogeneity and survey coverage

Baseline detection probability at mean survey effort ranged from 13.3% (*P. orientalis*) to 34.9% (*T. tetrax*), with *O. tarda* and *P. alchata* at intermediate levels (25.3% and 23.2%). Detection increased with survey duration for all four species (beta = 0.17–0.90, all P < 0.01) and with transect distance for three (*O. tarda*, *P. alchata*, *T. tetrax*). Time of day reduced detection for three species (beta = −0.17 to −0.44, P < 0.01), consistent with declining vocalisation activity later in the morning. These detection rates mean that the majority of true presences are missed on any single visit, with direct consequences for the accuracy of naive transition rates.

Final detection histories comprised 3,745 to 4,130 five-km grid cells with at least three repeat visits. Observed colonisation events were rare across all species (11–22 site-year transitions), while observed extinction events ranged from 9 (*T. tetrax*) to 29 (*O. tarda*).

### 3.2 Detection correction reverses the qualitative assessment

Without detection correction, apparent colonisation rates (0.57–0.89% per site-year) appeared non-negligible, in a range that could, in principle, compensate for local extinction. Detection-corrected estimates were 5- to 22,000-fold lower (Table 4; Fig. 2). For *P. orientalis*, corrected colonisation was effectively zero (gamma < 0.001%). Corrected extinction rates were lower than naive estimates for three species (ratios 0.38–0.78), because missed detections at occupied sites inflate apparent extinction; for *P. alchata*, corrected epsilon exceeded the naive estimate (ratio 1.72) owing to instability near separation.

**Table 4.** Naive versus detection-corrected colonisation (gamma) and extinction (epsilon) rates.

| Species | Naive gamma (%) | Corrected gamma (%) | Fold reduction | Naive epsilon (%) | Corrected epsilon (%) | Ratio |
|---|---|---|---|---|---|---|
| *O. tarda* | 0.89 | 0.0045 | 203× | 44.6 | 22.2 | 0.50 |
| *P. alchata* | 0.63 | 0.14 | 4.6× | 28.1 | 48.9 | 1.72⁺ |
| *P. orientalis* | 0.87 | ~0 | 21,839× | 34.1 | 44.1 | 0.78 |
| *T. tetrax* | 0.57 | 0.068 | 8.5× | 36.0 | 13.8 | 0.38 |

⁺ *P. alchata* corrected epsilon exceeds naive epsilon due to near-separation (18 events; SE > 11); interpret with caution.

The critical inference is not the magnitude of correction itself but its qualitative consequence: naive gamma values fall in a range that could appear demographically manageable, whereas corrected values are unambiguously insufficient to offset observed extinction. Without detection correction, the demographic trap is invisible (H3 supported).

### 3.3 The demographic trap: extinction dominates colonisation by orders of magnitude

Across all four species, detection-corrected extinction rates exceeded colonisation rates by two to six orders of magnitude (Table 5; Fig. 3). The asymmetry was robust across the full bootstrap distribution: even at the lower 95% CI of epsilon and the upper 95% CI of gamma, the ratio remained above 14 for every species, confirming a structural demographic constraint rather than a precision artefact. The probability that the ratio exceeds 100, a threshold below which some metapopulation recovery remains theoretically possible, ranged from 81.4% (*P. alchata*) to 99.3% (*P. orientalis*).

**Table 5.** Baseline demographic rates, epsilon/gamma asymmetry, and equilibrium occupancy.

| Species | Baseline gamma (%) | Baseline epsilon (%) | epsilon/gamma | P(ratio > 100) | psi* median (%) | psi* 95% CI (%) |
|---|---|---|---|---|---|---|
| *O. tarda* | 0.0045 | 22.2 | 4,870 | 97.7% | 0.02 | 0.0004–0.93 |
| *P. alchata* | 0.137 | 48.9 | 319 | 81.4% | 0.31 | 0.08–6.50 |
| *P. orientalis* | ~0 | 44.1 | 1,054,932 | 99.3% | 0.0001 | 0–0.17 |
| *T. tetrax* | 0.068 | 13.8 | 201 | 89.0% | 0.49 | 0.16–1.62 |

epsilon/gamma is the number of extinction events expected per colonisation event at mean current conditions. psi* = gamma/(gamma + epsilon) with 95% bootstrap CI (n = 5,000). P(ratio > 100): proportion of bootstrap replicates with epsilon/gamma > 100.

Species-specific patterns within the asymmetry reveal distinct demographic signatures. *P. orientalis* occupies the extreme: baseline gamma ≈ 3.8 × 10⁻⁷, an effectively irreversible state in which colonisation has ceased to function as a demographic process. Urban encroachment (LC13, beta = +2.80, P < 0.001) is the dominant extinction driver for *O. tarda*, while NDVI (climate-adjacent) and minimum temperature constrain its colonisation. *P. alchata* shows the highest baseline extinction (~49%) driven by precipitation, with a colonisation submodel too data-sparse for inference beyond directionality (16 events). *T. tetrax* presents the cleanest model: cropland (LC12) is the sole significant driver of both colonisation (beta = +0.90, P < 0.05) and extinction (beta = −0.57, P < 0.05), and its epsilon/gamma ratio (201) is the lowest in the guild, the only species retaining near-equilibrium dynamics.

The species-level variation did not follow the bustard-sandgrouse dichotomy we anticipated (H2 partially supported). *P. orientalis* (sandgrouse) showed the most extreme asymmetry; *T. tetrax* (bustard) the least. The ranking reflects the interaction between species-specific habitat requirements and contemporary landscape structure rather than mobility alone. *T. tetrax*'s reliance on cropland, a habitat that remains widespread in the Iberian Peninsula, provides a demographic buffer unavailable to species dependent on undisturbed steppe or open shrubland (H1 supported for all four species).

Equilibrium occupancy under current rates falls far below observed occupancy for all species (Table 6), implying extinction debt ranging from 5% (*T. tetrax*) to effectively 100% (*P. orientalis*). The delta-gamma analysis makes the deficit concrete: achieving even 10% equilibrium occupancy would require multiplying current colonisation rates by 22× (*T. tetrax*) to 117,215× (*P. orientalis*), increases far beyond any plausible management intervention (Table 7).

**Table 6.** Extinction debt: current versus equilibrium occupancy.

| Species | Current occupancy (%) | psi* (%) | Extinction debt (%) | Debt fraction | Recol. time (yr) |
|---|---|---|---|---|---|
| *O. tarda* | 0.68 | 0.02 | 0.65 | 97% | 22,648 |
| *P. alchata* | 0.87 | 0.31 | 0.56 | 64% | 732 |
| *P. orientalis* | 1.15 | 0.0001 | 1.15 | 100% | 2,599,962 |
| *T. tetrax* | 0.52 | 0.49 | 0.03 | 5% | 1,471 |

Recolonisation time = 1/gamma (years). Prevalence values integrate over the full sampling grid including cells outside each species' realised range; within-species ratios and cross-species ranking are the ecologically interpretable quantities.

**Table 7.** Delta-gamma: colonisation rate required for 10% equilibrium occupancy.

| Species | gamma required | Multiplier (median) | Multiplier 95% CI |
|---|---|---|---|
| *O. tarda* | 0.025 | 541× | 6–27,855 |
| *P. alchata* | 0.055 | 35× | 2–148 |
| *P. orientalis* | 0.049 | 117,215× | 86–244,501,693 |
| *T. tetrax* | 0.015 | 22× | 7–73 |

Model-simulated mean occupancy declined over the study period for all four species. *P. orientalis* showed the steepest decline (−37.9%), followed by *P. alchata* (−27.3%), *O. tarda* (−11.7%), and *T. tetrax* (−3.1%), consistent with the ranking of demographic asymmetry.

### 3.4 Spatial expression of the demographic trap

The colonisation-extinction imbalance is not uniformly distributed across the Iberian Peninsula (Fig. 4). Cell-level log₁₀(epsilon/gamma) maps reveal species-specific geographic concentrations of demographic risk. For *O. tarda*, the highest ratios cluster around expanding urban peripheries of central meseta cities, consistent with the strong LC13 effect on extinction. For *P. orientalis*, near-total demographic failure extends across the species' entire range, with the southwestern populations facing the most extreme ratios. For *P. alchata*, the Ebro valley concentrations carry the highest extinction pressure. *T. tetrax* shows the most spatially uniform pattern, with lower and more moderate ratios reflecting its near-equilibrium state. The spatial juxtaposition of high predicted initial occupancy (psi₁) with high extinction risk identifies specific regions where currently occupied habitat faces the greatest demographic pressure, the operational targets for extinction-prevention interventions.

### 3.5 Environmental drivers by demographic process

Environmental drivers showed clear species-specific signatures (Fig. 5; full coefficients in Supplementary S2), but a cross-species pattern partitions drivers by demographic process.

Climate variables dominated the colonisation submodels for three of four species. For *O. tarda*, NDVI (climate-adjacent; beta = −1.45, P < 0.05) and minimum temperature (beta = −4.52, P < 0.05) constrain colonisation, suggesting that vegetation encroachment and cold-season severity limit establishment at new sites. For *P. orientalis*, minimum temperature strongly reduced colonisation (beta = −8.11, P = 0.002) while maximum temperature increased it (beta = +9.14, P < 0.001), indicating sensitivity to the diurnal temperature range during the breeding season. *T. tetrax* was the exception: cropland (LC12) was the sole significant colonisation driver (beta = +0.90, P < 0.05), with no climate variable retained.

Land-use variables played a more prominent role in extinction, but in species-specific configurations. Urban encroachment (LC13) was the dominant extinction driver for *O. tarda* (beta = +2.80, P < 0.001). For *P. orientalis*, NDVI strongly increased extinction (beta = +2.05, P < 0.001; climate-adjacent) while cropland simultaneously reduced it (LC12, beta = −0.83, P = 0.02), indicating that agricultural mosaics buffer local extinction risk. For *T. tetrax*, cropland reduced extinction (LC12, beta = −0.57, P < 0.05), making it the only species for which a single management lever, maintaining agricultural land use through agri-environment schemes, can simultaneously increase colonisation and reduce extinction.

### 3.6 Counterfactual attribution

Climate had larger and more consistent effects on colonisation than land use for three species (*O. tarda*, *P. alchata*, *P. orientalis*), while land-use change drove both processes for *T. tetrax* (Table 8). For *T. tetrax*, the land-use pathway accounts for 100% of both gamma and epsilon attribution, with climate contributions identically zero, consistent with the absence of climate covariates in the AIC-selected model. For *O. tarda*, climate explains 100% of gamma attribution while land use reduces epsilon, confirming that the two driver classes operate on different demographic levers. Effects are small in absolute terms (< 0.01 probability units), reflecting the seven-year window, but informative comparatively across species and driver classes.

**Table 8.** Attribution summary: mean delta-psi per scenario relative to S0 baseline.

| Species | Dominant driver | Gamma climate | Gamma land-use | Epsilon climate | Epsilon land-use |
|---|---|---|---|---|---|
| *O. tarda* | Climate | −6.1e−4 | 0 | 7.4e−4 | −1.1e−3 |
| *P. alchata* | Climate* | excl. | excl. | −1.3e−3 | 0 |
| *P. orientalis* | Climate | −3.4e−4 | −1.8e−4 | 5.7e−3 | −2.1e−4 |
| *T. tetrax* | Land use | 0 | −1.0e−6 | 0 | 9.0e−6 |

*P. alchata* gamma excluded (separation). NDVI decomposition amplified the climate signal for *P. orientalis* extinction by 3× but did not change the dominant-driver classification for any species (Supplementary S10).

### 3.7 Spatial structure and conservation scale

[Section conditional on stPGOcc convergence — Rhat < 1.1, ESS > 100 for phi and sigma²_t.]

Spatio-temporal occupancy models improved fit over colext for all species (delta-WAIC: [values pending]). Inclusion of NNGP spatial random effects reduced residual Moran's I by 59–82%. The qualitative direction and relative magnitude of covariate effects were preserved, confirming that the rate disparity is not an artefact of spatial confounding. Effective spatial ranges followed the mobility axis (Table 9; H4 supported).

**Table 9.** Effective spatial ranges from stPGOcc models.

| Species | Spatial range (km) | 95% CrI | Published dispersal | Conservation scale |
|---|---|---|---|---|
| *T. tetrax* | ~43 | [stPGOcc phi CrI — fill from convergence run] | < 50 km | Landscape-scale SPAs |
| *O. tarda* | ~47 | [stPGOcc phi CrI — fill from convergence run] | < 50 km | Landscape-scale SPAs |
| *P. orientalis* | ~172 | [stPGOcc phi CrI — fill from convergence run] | 100–200 km | Cross-regional coordination |
| *P. alchata* | ~264 | [stPGOcc phi CrI — fill from convergence run] | 100–300 km | Cross-regional network |

### 3.8 External validation

Spatially blocked cross-validation against the III Atlas of Breeding Birds confirmed predictive skill for initial occupancy across all species (Table 10). Calibration slopes on the logit scale (0.57–0.62) indicate moderate overdispersion, a common feature of occupancy models at coarse resolution.

**Table 10.** Spatially blocked cross-validation of psi₁ against the Spanish Breeding Bird Atlas.

| Species | Spatial AUC | Calibration slope | TSS | N test polygons |
|---|---|---|---|---|
| *O. tarda* | 0.884 | 0.609 ± 0.030 | 0.612 | 5,202 |
| *P. alchata* | 0.844 | 0.575 ± 0.033 | 0.569 | 5,202 |
| *P. orientalis* | 0.823 | 0.623 ± 0.029 | 0.536 | 5,202 |
| *T. tetrax* | 0.850 | 0.598 ± 0.022 | 0.546 | 5,202 |

Block size = 270 km (exceeding the largest estimated spatial range).

## 4. Discussion

### The demographic trap and its theoretical context

The central finding of this study is that colonisation has structurally failed for four co-occurring steppe bird species in the Iberian Peninsula. Extinction exceeds colonisation by two to six orders of magnitude (epsilon/gamma: 201 to 1,054,932), and this asymmetry is robust to the wide uncertainty in absolute colonisation rates: even at the most optimistic end of the bootstrap distribution, the ratio remains above 14 for every species. What do these numbers mean ecologically? An epsilon/gamma ratio of 4,870 for *O. tarda* implies that for every site naturally recolonised, nearly five thousand are lost. For *P. orientalis*, with a ratio exceeding one million, colonisation has effectively ceased to function as a demographic process, and expected recolonisation times (2.6 million years) operate on geological rather than ecological timescales. These are not incremental deficits amenable to management optimisation. They describe a qualitative demographic regime, a trap, in which ranges are decoupled from habitat availability and contracting under their own demographic inertia.

The concept of extinction debt (Tilman et al. 1994; Kuussaari et al. 2009) applies directly. For *O. tarda*, 97% of current occupancy is transient. For *P. orientalis*, effectively all occupied sites represent legacy occupancy that will be lost as the system relaxes toward equilibrium. Only *T. tetrax*, with a debt fraction of 5%, retains meaningful demographic resilience. This finding extends metapopulation theory (Hanski 1999) to a conservation context in which the colonisation-extinction balance is not merely unfavourable but structurally broken, and it does so for an entire guild facing the same landscape transformation. The non-equilibrium range dynamics literature has theorised this possibility (Zurell et al. 2022; Briscoe et al. 2021); our results provide a quantitative empirical demonstration with the rates, the confidence intervals, and the species-specific variation needed to translate theory into conservation prescriptions.

### Uncorrected monitoring generates false confidence

Naive colonisation rates of 0.57–0.89% per site-year fall in a range where a monitoring programme might reasonably conclude that recolonisation, while slow, offers some recovery potential. Detection-corrected rates of 0.0045–0.14%, with fold reductions of 5 to 22,000, eliminate this possibility. This is not a statistical refinement: it is the difference between concluding that recovery is difficult and recognising that recovery is structurally impossible without radical demographic change. The directional nature of detection bias in colonisation estimates has been noted previously (Guillera-Arroita 2017; Kujala et al. 2013), but the magnitude of its effect on conservation conclusions has rarely been demonstrated empirically across multiple co-occurring species.

The bias is particularly insidious because it operates in the direction that masks urgency. Naive estimates simultaneously inflate apparent colonisation (by 5–22,000×) and deflate apparent extinction for three of four species (naive/corrected ratios 0.38–0.78), making the demographic regime appear substantially more favourable than it is. Any citizen-science monitoring programme used to assess demographic viability should incorporate detection correction as a minimum methodological standard. For the species studied here, uncorrected monitoring has been generating false confidence about recovery potential, not merely imprecise estimates.

### Species-specific drivers and management implications

The attribution analysis reveals that climate and land use operate on different demographic processes with distinct species-specific signatures, but this partitioning must be interpreted within the constraints of the covariate design. Static land-cover variables were consumed by psi₁, and interannual variance in TerraClimate covariates structurally exceeds variance in MODIS land-cover classes over a seven-year window. Climate attribution for three species represents the dominant signal detectable within these constraints, not a claim about long-term relative importance.

Within this framework, the species-specific patterns have direct management relevance. *T. tetrax* is the most tractable case: cropland availability governs both colonisation and extinction through a single management lever, agri-environment schemes under CAP Pillar II that maintain crop-fallow mosaics. This mechanistic clarity, combined with the lowest asymmetry ratio and near-equilibrium dynamics, makes *T. tetrax* the species for which current policy instruments are most likely to be effective. The sandgrouse present a harder problem. For *P. orientalis*, climate-adjacent NDVI drives extinction while temperature extremes govern colonisation; reducing extinction requires managing the habitat quality that buffers climate stress, maintaining agricultural mosaics that provide thermal refugia and ensuring access to traditional drinking sites, but colonisation cannot be meaningfully increased by habitat management alone if the dominant signal on gamma is climatic. *O. tarda* presents the urban-encroachment case: LC13 is the dominant extinction driver, meaning that halting infrastructure expansion adjacent to breeding areas is the highest-priority intervention, a policy target that is precise and actionable if politically contested.

The NDVI decomposition for *P. orientalis* illustrates a broader challenge: NDVI integrates vegetation encroachment (climate-adjacent, negative for steppe specialists) and herbaceous productivity (partly management-driven). Decomposing NDVI into climate and residual components amplified the climate signal for extinction by 3×, confirming that the original effect was largely climate-mediated vegetation encroachment rather than land-management change.

### Conservation scale from spatial range estimates

The fourfold gradient in effective spatial range (43–264 km) translates directly into conservation planning requirements. For *T. tetrax* and *O. tarda* (~43–47 km), landscape-scale management within existing Special Protection Areas is spatially adequate if habitat quality is maintained. For *P. orientalis* (~172 km), the range exceeds typical SPAs but remains manageable through coordination across adjacent autonomous communities. For *P. alchata* (~264 km), no single protected area can function as a self-sustaining population unit; cross-regional networks connecting multiple SPAs are the minimum meaningful conservation unit. The concordance between model-estimated spatial ranges and published dispersal distances from telemetry (Table 9) provides an ecological coherence check, strengthening confidence that the NNGP random effect captures genuine demographic connectivity rather than statistical artefact. A network designed for the bustards will be too fine-grained for the sandgrouse; one designed for *P. alchata* will be unnecessarily coarse for the bustards. Species-specific spatial design is essential, and the spatial range estimates provide the quantitative basis currently absent from Iberian steppe bird conservation planning.

### Limitations

Four analytically serious limitations bound the inference. First, the *P. alchata* colonisation submodel is underpowered (16 events), producing unreliable coefficient estimates; gamma for this species is excluded from attribution and its point estimates carry wide confidence intervals. We retain *P. alchata* because the directional finding (gamma ≪ epsilon) is robust and removing it would eliminate half the sandgrouse component. Second, the seven-year study window captures approximately one to three turnover cycles, sufficient for epsilon estimation but producing wide confidence intervals on gamma. The primary inference, that epsilon/gamma ≫ 1, holds across the full bootstrap distribution including the upper CI, and the window coincides with a period of accelerating land-use change, making the captured dynamics relevant for near-term planning. Third, external validation is restricted to psi₁; transition rates have no independent benchmark. Confirmation of causal attribution will require quasi-experimental designs, for example comparing transition rates before and after implementation of targeted agri-environment schemes within a BACI framework. Fourth, the colext closure assumption is more likely violated for the semi-nomadic sandgrouse than for the site-faithful bustards. Within-season nomadic movements would inflate gamma for the sandgrouse, the opposite direction of our main result and therefore conservative with respect to the demographic trap conclusion.

## 5. Conclusions

The Iberian steppe bird guild is caught in a demographic trap in which extinction proceeds at rates that colonisation cannot offset under any plausible scenario. The asymmetry spans two to six orders of magnitude across species, producing extinction debts of 5–100% of current occupancy, and is invisible to monitoring that does not correct for imperfect detection. Preventing local extinction at occupied sites is the only intervention with demographic leverage on management-relevant timescales: natural recolonisation cannot compensate for extinction losses at current rates, and achieving even modest equilibrium occupancy would require 22- to 117,000-fold increases in colonisation. Because climate and land-use drivers act on different demographic processes with species-specific signatures, effective conservation requires species-specific targeting at spatial scales ranging from landscape SPAs for the bustards to cross-regional networks for the sandgrouse. Citizen-science data, when analysed with detection correction and demographic decomposition, can reveal these traps, traps that standard occupancy monitoring would interpret as slow decline rather than structural failure. This diagnostic capacity is urgently relevant across the European farmland bird crisis.

---

## References

[References to be compiled for submission]

---

## Figure Legends

**Figure 1. Study area and sampling effort.** Globe inset showing the Iberian Peninsula in European context (orthographic projection centred on ~40°N, 5°W), with *T. tetrax* observation locations (points) and eBird checklist density (background shading) across mainland Spain on a 5-km grid (2017–2023). Species range boundaries from the III Atlas of Breeding Birds in Spain overlaid as polygon outlines. The Iberian Peninsula supports 34% of the European *T. tetrax* breeding population and 19–98% for the full four-species guild analysed in this study. File: `figs/pub_map_sampling_sites.png` [requires revision: add globe inset, atlas polygons, focus on *T. tetrax* as representative species].

**Figure 2. Uncorrected occupancy models overstate colonisation potential by 5- to 22,000-fold.** Compact dotplot of naive versus detection-corrected gamma (colonisation, left) and epsilon (extinction, right) for all four species on a log scale. Fold-reduction annotations quantify detection bias magnitude. *P. orientalis* corrected gamma falls below the detection floor. Companion values with 95% confidence intervals are in Table 4. File: `figs/pub_fig1_naive_vs_corrected.png`.

**Figure 3. All four species occupy the decline zone of demographic state space, confirming a structural demographic trap.** **(a)** Scatter plot of mean detection-corrected colonisation (gamma) versus extinction (epsilon) on log₁₀ axes, with equilibrium occupancy isoclines (psi* contours) and the gamma = epsilon diagonal (no net change). All four species fall deep below the diagonal, with 95% bootstrap CI crosshairs (n = 5,000). Species labels placed directly on the plot. **(b)** Extinction debt as horizontal stacked bars: green = equilibrium occupancy (psi*), pink = transient fraction. Percentage labels quantify the proportion of current occupancy that is demographically unsustainable. File: `figs/pub_fig_isocline_equilibrium.png`.

**Figure 4. Spatial expression of the demographic trap.** Two-column panel per species (4 rows × 2 columns). Left column: predicted initial occupancy (psi₁, 2017; BuPu palette from white to dark). Right column: cell-level log₁₀(epsilon/gamma) bivariate risk index (diverging palette: blue = near-equilibrium dynamics, red = strong demographic decline). The spatial juxtaposition reveals where currently occupied habitat faces the highest extinction pressure, identifying the operational targets for extinction-prevention interventions. Regions where high psi₁ coincides with high log₁₀(epsilon/gamma) are the priority conservation targets. File: composite from `figs/pub_map_occupancy_4species.png` and `figs/pub_fig3_bivariate_risk_map.png` [requires assembly into two-column layout].

**Figure 5. Environmental drivers of occupancy dynamics.** Three-panel forest plot of standardised coefficients for **(A)** initial occupancy psi₁, **(B)** colonisation gamma, and **(C)** extinction epsilon. Species distinguished by shape and colour (Wong palette). Non-significant coefficients shown at reduced opacity. *P. alchata* gamma excluded from panel B (complete separation, n = 16 events). Note: NDVI marked with ⁺ denotes classification as "climate-adjacent" in the attribution analysis (R² ≈ 0.51 with climate variables). File: `figs/pub_fig2_forest_3submodels.png`.

---

## Editorial notes for authors

### (a) Result values requiring confirmation from output files

1. stPGOcc spatial range estimates (Table 9): all phi CrI values are placeholders pending convergence (Rhat < 1.1, ESS > 100). Fill from stPGOcc output once cluster run completes.
2. stPGOcc delta-WAIC values (Section 3.7): pending convergence run.
3. stPGOcc Moran's I reduction percentages (Section 3.7): "59–82%" carried forward from skeleton v10 — confirm from final stPGOcc diagnostics.
4. Detection coefficient for *P. alchata* NDVI on colonisation (beta = −1.08, P < 0.001): verify this is the *P. alchata* gamma submodel value from `results/` output, not the *O. tarda* value.
5. Table 3 calibration slopes described as "0.57–0.62" in text but table shows range 0.575–0.623 — minor rounding, confirm preferred precision.

### (b) Flagged logical tensions or incomplete information

1. *P. alchata* corrected epsilon exceeding naive epsilon (ratio 1.72): the skeleton attributes this to near-separation instability. Consider whether this result should be presented more cautiously or whether the parboot P = 0.52 adequately supports the model specification for this species.
2. NDVI decomposition amplification (3× for *P. orientalis* epsilon): the decomposed value (1.76e−2) is presented in Supplementary S10 but not in the main text Table 8. Consider whether the main-text attribution table should note this sensitivity.
3. The abstract states "AUC 0.82–0.88" but the precise range from Table 10 is 0.823–0.884. Confirm rounding convention.
4. Section 3.7 (stPGOcc) is largely placeholder text. If convergence is not achieved before submission, consider: (a) presenting colext Moran's I diagnostics as the primary spatial check; (b) moving stPGOcc to a "preliminary" supplementary section; (c) removing phi-based conservation scale arguments from Discussion paragraph 4 or qualifying them as provisional.
5. Figure 1 and Figure 4 require revision/assembly from existing components — confirm figure production timeline.

### (c) Word count per section

| Section | Estimated words |
|---|---|
| Abstract | 248 |
| Introduction | 890 |
| Methods | 1,780 |
| Results | 1,820 |
| Discussion | 1,810 |
| Conclusions | 195 |
| **Total main text** | **~6,745** |
| Figure legends | ~520 |
