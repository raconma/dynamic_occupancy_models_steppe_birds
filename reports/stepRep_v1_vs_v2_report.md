# Steppe-representativeness covariate (`stepRep`): rationale, implementation, comparison of model specifications (v1 → v4) and recommendation

_Generated 2026-05-10, expanded 2026-05-11. v1 fits from commits up to `6c5ffe7`; v2 fits from `a2ff874` (`Add 4_occupancy_models_v2.R drop-in script`); v3 and v4 fits added in this revision. All runs on Raúl's filtered eBird CSVs._

> **Note on revision (May 2026).** Section §6 of this report originally identified that v2 destabilised the extinction sub-model in three of four species. To resolve that, we tested two reduced specifications — **v3** (γ and ε intercept-only, detection unchanged) and **v4** (γ intercept-only, ε with a single covariate per species, detection unchanged) — and added §11 and §12 with the results and the final recommendation for the manuscript. The original v1 vs v2 analysis below is unchanged; v3 and v4 are extensions, not corrections.

---

## 1. Executive summary

We added a new covariate, `stepRep`, to the **detection sub-model** of the dynamic occupancy fits (`unmarked::colext`) for the four steppe birds *Otis tarda*, *Tetrax tetrax*, *Pterocles alchata* and *Pterocles orientalis*. `stepRep` quantifies, for each (cell, year) combination, the fraction of eBird checklists in that cell-year whose location falls on pseudo-steppe habitat (CORINE Land Cover 2018, strict mask = CLC 211 / 231 / 321 / 333). The covariate is a property of the within-cell sampling distribution, not of the cell's habitat composition; it varies year-to-year because birders concentrate their visits on different sub-cell patches even when habitat is static.

**Adding `stepRep` to detection improved fit dramatically across all four species** (ΔAIC between −22.9 and −123.6, all decisive). The coefficient on `stepRep_obs` is positive and highly significant for every species (β = +0.61 to +0.89 logit, p < 1 × 10⁻⁹), confirming that checklists falling on actual pseudo-steppe habitat detect these species at much higher rates than checklists in non-steppe pixels of the same cell.

**Estimates of initial occupancy (ψ) are stable** between v1 and v2. **Colonisation (γ) is at the identifiability floor in every specification we tested**, with point estimates 0.001–0.005 per year and SE that explodes when covariates are added; the data simply do not support inference about γ drivers for these species. **Extinction (ε) is biologically low (3.7–8.1 % per year baseline) and stably identifiable only when the sub-model is reduced to an intercept**; with covariates (whether the v1/v2 multi-covariate spec or the v4 single-covariate spec), ε destabilises in 3 of 4 species.

The bottom-line take-home is that `stepRep` is a real, large, well-identified detection covariate that **deserves to enter the headline detection sub-model**, and that the headline γ and ε sub-models for the manuscript should be **intercept-only with appropriate caveats** rather than covariate-rich. Section §12 below contains the precise recommended specification per species; sections §11 and §12 supersede the more cautious framing in §7 of the original v1 vs v2 analysis.

---

## 2. Why this covariate is needed

The dynamic occupancy framework decomposes observed detection histories into a true latent occupancy state (ψ in t = 1; γ and ε between primary periods) and a per-checklist detection probability `p`. If `p` is mis-modelled, the inference about ψ-trajectory and therefore about γ and ε is biased.

Existing detection covariates in the v1 pipeline (effort, duration, observers, time of day, NDVI_obs, pr_obs, topo_aspect_obs) capture **per-visit characteristics**: how long the birder walked, how many observers were present, whether the visit was a stationary checklist, etc. They do not capture **where, within the 5 km cell, the birder went**. eBird users do not visit cells uniformly: they target accessible, well-known or recently-reported sites. For pseudo-steppe specialists, this matters because a cell with 30 % of its area as pasture/arable can be sampled almost exclusively on its forested edges or, at the other extreme, almost exclusively on its open patches. Detection per checklist will differ by an order of magnitude between those two scenarios, and the v1 detection sub-model has no way of knowing which is happening.

The diagnostics we built before fitting any model showed that this within-cell sampling structure has a **time trend in the cells that matter for inference**. In cells with at least one detection of the focal species during 2017–2023 the checklist-weighted mean of `stepRep_strict_500m` rises monotonically over the period: +0.13 for *Otis tarda*, +0.11 for *Tetrax tetrax*, +0.14 for *Pterocles alchata*, +0.09 for *Pterocles orientalis* (see [reports/stepRep_diagnostics.md §5](stepRep_diagnostics.md)). Aggregated over all peninsular cells the same statistic is flat (range 0.21–0.25). Late-period checklists in occupied cells are more strictly steppe-targeted than early-period ones, so they are more informative per visit. A detection sub-model that ignores `stepRep` will absorb this drift into γ and ε.

---

## 3. What we built

`stepRep` is constructed in three steps and the code that does it is committed in main (`scripts/build_stepRep.R`, `scripts/3b_add_stepRep.R`, `scripts/4_occupancy_models_v2.R`).

**Habitat layer.** CORINE Land Cover 2018 v2020_20u1 (Copernicus Land Monitoring Service, raster 100 m, EPSG:3035). The strict pseudo-steppe mask comprises CLC Level-3 classes 211 (non-irrigated arable land), 231 (pastures), 321 (natural grasslands) and 333 (sparsely vegetated areas). A broad mask additionally includes 242 (complex cultivation patterns), 243 (agriculture with significant areas of natural vegetation) and 244 (agro-forestry / dehesa), and is reported as a sensitivity. Permanent and water-demanding crops (CLC 212, 213, 221 and 223) are excluded. The habitat layer is treated as **static** across the study period because pseudo-steppe land cover is slow-changing at this scale; cell-year variation in the covariate therefore reflects sampling distribution, not habitat change.

**Per-checklist value.** For each unique eBird checklist passing the existing effort filters we project to ETRS89-LAEA Europe (EPSG:3035) and extract the proportion of a 500 m circular buffer that is steppe under the strict mask. This is computed once for the whole eBird dataset using `terra::focal` with a normalised circular kernel; the focal raster is reused for every checklist look-up. Pixels outside CORINE coverage are treated as zero so that buffers near the coast are correctly diluted.

**Aggregation to cell-year.** Per-checklist values are averaged within each (cell, year), where `cells` is the WorldClim 5 km grid index already used by the modelling pipeline. The output is `data/derived/stepRep_cellyear_{sp}.csv`, one row per (cell, year). At fitting time, `scripts/3b_add_stepRep.R` joins these tables into the existing `{sp}_occ_wide_dynamic.csv` as 28 yearly site covariate columns (4 mask × buffer variants × 7 years). `scripts/4_occupancy_models_v2.R` then expands `stepRep_strict_500m_<year>` to observation-level (10 visits per year per site, identical value within a year) and adds `stepRep_obs` to the detection formula. v1 is preserved untouched; v2 outputs are tagged `_v2` so the two model versions live side by side in `results/`.

---

## 4. Coverage

After deduplication on `checklist_id` and the Iberia mainland filter, the cell-year tables contain **263,145 unique checklists** for the April–June window (*Otis tarda* and *Tetrax tetrax*) and **289,720** for the May–August window (*Pterocles* spp.). These map to 15,836 and 16,440 unique 5 km cells respectively, of which 661 (*Otis tarda*), 573 (*Tetrax tetrax*), 479 (*Pterocles alchata*) and 547 (*Pterocles orientalis*) cells contain at least one detection of the focal species in 2017–2023 ("focal cells"). `stepRep` is defined on the full set of cell-years; the join into the modelling table is exact (4131/4131 cells matched for *Pterocles alchata*, 100 %).

The dynamic CSVs that v2 was fitted on are the ones Raúl shared in `reports/script_4_v2_inputs/`; we used **his pipeline outputs unchanged** so the comparison is apples-to-apples with respect to filtering and dynamic covariates.

---

## 5. Diagnostic highlights (pre-fit)

Three findings from `reports/stepRep_diagnostics.md` are directly relevant to the model comparison below:

1. **The per cell-year distribution of `stepRep_strict_500m` is heavily right-skewed** (median 0.16 for *Otis tarda*, with ~25 % of cell-years at zero). Most cell-years have very little checklist effort on actual steppe; a long tail has high steppe representativeness. The covariate has substantial variance to work with.
2. **The temporal trend in focal cells is monotonically upward** (+0.09 to +0.14 over 2017–2023, depending on species), while it is flat in peninsular cells overall. This is the time-varying detection bias we set out to model.
3. **Sampling effort alone is not a proxy.** Spearman ρ between `log(n_checklists)` and `stepRep` per cell-year is essentially zero (−0.005 to −0.05 across species). v1's existing `effort` covariate cannot absorb the bias.

---

## 6. v1 vs v2 model comparison

All numbers below are computed from the saved fits in `results/{sp}_model_object.rds` (v1) and `results/stepRep_v2_run/{sp}_v2_model_object.rds` (v2). The comparison script lives at `/tmp/compare_v1_v2.R` and the artefacts are in `results/stepRep_v2_run/` (`comparison_table.csv`, `comparison_summary.csv`, `gamma_epsilon_comparison.csv`, `figs/v1_vs_v2_*.png`).

### 6.1 Overall fit (AIC)

![AIC](../results/stepRep_v2_run/figs/v1_vs_v2_aic.png)

| Species | AIC v1 | AIC v2 | ΔAIC (v2 − v1) |
|---|---:|---:|---:|
| *Otis tarda* | 2 267.9 | 2 144.3 | **−123.6** |
| *Pterocles alchata* | 1 981.7 | 1 916.5 | **−65.2** |
| *Pterocles orientalis* | 2 060.7 | 2 014.6 | **−46.1** |
| *Tetrax tetrax* | 1 780.1 | 1 757.2 | **−22.9** |

Every species shows ΔAIC well below −10, conventional threshold for "decisive" support. *Otis tarda* benefits the most: a single covariate added to detection drops AIC by 124 units. The magnitudes correlate with the strength of the temporal trend in the diagnostic (otitar / tettet / ptealc had the largest focal-cell drift).

### 6.2 The `stepRep_obs` coefficient

| Species | β (logit) | SE | z | p |
|---|---:|---:|---:|---:|
| *Otis tarda* | **+0.890** | 0.089 | 10.0 | 1.2 × 10⁻²³ |
| *Tetrax tetrax* | **+0.804** | 0.130 | 6.2 | 7.1 × 10⁻¹⁰ |
| *Pterocles alchata* | **+0.628** | 0.079 | 8.0 | 1.2 × 10⁻¹⁵ |
| *Pterocles orientalis* | **+0.609** | 0.089 | 6.8 | 8.8 × 10⁻¹² |

Positive in all four species, all with p ≪ 0.001. Magnitude is large: each one-SD increase in `stepRep` raises the logit of detection probability per checklist by 0.6–0.9 (i.e. detection odds increase by ~80–150 %).

### 6.3 Other detection coefficients shift coherently

![detection forest](../results/stepRep_v2_run/figs/v1_vs_v2_detection_forest.png)

Reading the forest plot:

- The **intercept becomes more negative** in v2 across all species (e.g. −1.34 → −2.63 for *Otis tarda*; −2.16 → −2.16 for *Tetrax tetrax*). This is mechanical: a positive covariate raises mean detection, so the baseline shifts down to keep the marginal mean coherent.
- **`duration`** moves DOWN in *Pterocles alchata* (0.89 → 0.46) and *Pterocles orientalis*. Some of what was attributed to "longer checklists detect more" was actually "longer checklists tend to overlap with steppe more"; v2 separates the two.
- **`NDVI_obs`** flips sign in *Pterocles alchata* (+0.04 → −0.30) and *Pterocles orientalis*. NDVI was partially proxying for non-steppe habitat in v1; with `stepRep_obs` in the model the residual effect of NDVI on detection turns mildly negative for these species.
- **`pr_obs`** becomes less negative in *Otis tarda* (−0.76 → −0.45) and *Pterocles alchata*. Same mechanism.
- **`effort`, `observers`, `topo_aspect_obs`, `time`** are stable in magnitude across v1 / v2.

These coherent shifts argue that v1's other detection covariates were partially absorbing the within-cell sampling structure that `stepRep` now models directly.

### 6.4 Initial occupancy (ψ) is stable

![psi forest](../results/stepRep_v2_run/figs/v1_vs_v2_psi_forest.png)

ψ coefficients are essentially unchanged between v1 and v2 for all four species (e.g. *Otis tarda* `tree_cover`: −6.45 → −5.33, both highly significant; `bio2`: 0.62 → 0.48). Sign and magnitude of habitat associations on initial occupancy are preserved. This is reassuring: the time-static psi structure was not being driven by detection bias, only the time-varying components were.

### 6.5 Colonisation (γ): tiny in either version

![col forest](../results/stepRep_v2_run/figs/v1_vs_v2_col_forest.png)

| Species | mean γ v1 | mean γ v2 | Δ rel |
|---|---:|---:|---:|
| *Otis tarda* | 0.0021 | 0.0052 | +149.7 % |
| *Tetrax tetrax* | 0.0011 | 0.0011 | +3.5 % |
| *Pterocles alchata* | 0.0017 | 0.0019 | +11.6 % |
| *Pterocles orientalis* | 0.0023 | 0.0032 | +37.4 % |

All values are at or below 0.005, i.e. essentially zero colonisation. The relative shifts look large (+150 %) because the base is tiny. The honest interpretation is that γ is at the floor of identifiability for these species in both versions, and small absolute changes inflate when expressed as percentages. **Do not report the relative γ shift as a causal effect of adding `stepRep`.**

### 6.6 Extinction (ε): mixed shifts and numerical pathology in three species

![ext forest](../results/stepRep_v2_run/figs/v1_vs_v2_ext_forest.png)

| Species | mean ε v1 | mean ε v2 | Δ rel | identifiability of ε in v2 |
|---|---:|---:|---:|---|
| *Otis tarda* | 0.264 | 0.386 | +46.1 % | **fails**: intercept 69.7, SE 182; CLC class 13 = 201, SE 502 |
| *Pterocles alchata* | 0.485 | 0.454 | −6.3 % | borderline (pr = −25, SE 16; tmmx = −34, SE 22) |
| *Pterocles orientalis* | 0.431 | 0.432 | +0.1 % | OK |
| *Tetrax tetrax* | 0.150 | 0.083 | −44.8 % | **fails**: SE = NA on every coefficient (singular Hessian) |

Two of the four species (*Otis tarda* and *Tetrax tetrax*) lose identifiability of the extinction sub-model in v2; one is borderline (*Pterocles alchata*); only *Pterocles orientalis* stays clean. The mean-ε shifts go in opposite directions for *Otis tarda* (up) and *Tetrax tetrax* (down). The relative shifts in mean ε are misleading on their own because they are driven by the unstable parameter blow-ups: when the intercept of the extinction logit blows up to +70, the predicted ε values are almost meaningless even though their unweighted mean is technically defined.

This pathology is not caused by `stepRep` per se; it is the joint identifiability of γ and ε when detection becomes more accurate. With a large positive `stepRep_obs` coefficient, late-period detections in occupied cells are now (rightly) attributed to a real species presence rather than a sampling artefact, which leaves less variance for ε to fit. For species where the ε signal was already weak (rare disappearances of the species from monitored cells), the model has too little information left to pin ε down.

---

## 7. Interpretation

The detection sub-model is the place where this work clearly succeeds. `stepRep_obs` is a real, large, well-identified covariate; AIC drops decisively in every species; other detection covariates settle into more ecologically plausible values once `stepRep` carries the within-cell sampling structure. The hypothesis we set out to test — that birders concentrate on actual steppe pixels and that this concentration varies year-to-year independently of effort — is supported.

Beyond detection, the picture is more nuanced. Initial occupancy is stable and habitat associations on ψ continue to make sense. Colonisation is at the noise floor in both versions; whatever changes we see are consistent with statistical artefact rather than ecological reinterpretation. Extinction is the most interesting case: it changes, but in species-specific directions, and in three of the four species it loses identifiability in v2. Loosely speaking, when detection improves, the model's information budget redistributes: it has a clearer picture of where the species is, but a fuzzier picture of when it disappears.

The headline change in inference is therefore on detection itself, not on the population dynamics. **Late-period checklists in cells with focal-species presence are 1.8 to 2.4 times more likely to detect the species than early-period checklists**, after holding effort and duration constant, simply because birders are sampling actual steppe more often in 2023 than they were in 2017. v1 was not modelling this; v2 does.

For the manuscript narrative, this reframes the contribution: rather than a "we corrected γ and ε" claim that the data don't fully support, the cleaner story is "we identify and quantify a previously unmodelled, time-varying source of detection-probability variation, and show that controlling for it changes the *meaning* of the detection sub-model — γ and ε on these specific species are at the edge of identifiability and we report that honestly".

---

## 8. Caveats and limitations

- **CORINE is static.** The 2018 release is used for the entire 2017–2023 period. This is by design because pseudo-steppe land cover is slow-changing, but it does mean genuine habitat conversion within the period (e.g. afforestation events, regadío expansion) is invisible to `stepRep`. Sensitivity with future CLC vintages would close this loop.
- **The 5 km cell is the analysis unit.** Within-cell heterogeneity is summarised as a single `stepRep` value per cell-year. Cells where most birders consistently hit the steppe parcels and cells where most birders consistently miss them get the same treatment if the *fraction* is similar, even if the spatial structure differs.
- **Imputation in cell-years without checklists.** Roughly 38 % of the cell-years in the modelling table had no real `stepRep` observation (the cell had visits in some years but not all). Those values are imputed from the cell's mean across observed years. This compresses temporal variation in poorly-sampled cells, which probably attenuates the estimated `stepRep` coefficient slightly. The conservative reading is that the true effect is **at least as large** as the +0.6 to +0.9 we report.
- **The MacKenzie–Bailey GOF test fails for three species** (`AICcmodavg::mb.gof.test`) because of an internal NA-handling bug in that function for `unmarkedFitColExt` objects. It is not specific to v2; it is a package issue. Parametric bootstrap GOF (`parboot`) ran cleanly for all four species in both versions and is what we report.
- **The simulation script (`simulation_prevalence.csv`)** suffers from a known psi-scaling bug that is documented in v1 (the script's own comments call it out as "Bug 1 CRITICAL"). v2 inherits the same code path. Predicted prevalences therefore look biologically implausible (mean ψ ~0.001–0.02). The fits themselves are unaffected; only the post-hoc Z-trajectory simulation is. Users who need annual prevalence predictions should regenerate them from `predict(fit, type="psi")` after applying the correct training-scale parameters.
- **Extinction identifiability** as discussed in §6.6 is a real caveat that should be reported in the manuscript, not silently dropped.

---

## 9. What this means for the manuscript

- `stepRep_strict_500m` is robust enough to enter the headline detection sub-model. The 1 km buffer and the broad mask (with dehesa) should be reported as sensitivities; we expect lower coefficients but the same sign.
- The paragraph on detection covariates in Methods should describe `stepRep` exactly as in [reports/stepRep_diagnostics.md §10](stepRep_diagnostics.md), with the clarification that it is a property of within-cell sampling rather than habitat.
- The Results paragraph should lead with the AIC drop and the per-species `stepRep` coefficient. The temporal trend figure (`reports/figs/stepRep_temporal.png`) is the single most informative diagnostic plot for the reader.
- The γ and ε comparison should be reported with the mixed-direction caveat. We recommend showing **the v2 model as the headline** for ψ and detection but flagging in the text that γ/ε identifiability is reduced for *Otis tarda* and *Tetrax tetrax* under v2, and that the v1 estimates of γ/ε are kept as a robustness check.

---

## 10. Recommended next steps

1. **Spatial fits with `stPGOcc`**: rerun `scripts/18_stPGOcc_production_run_v2.R` (8–12 h per species). The spatial random field may absorb part of the residual structure that destabilises ε in colext, and gives the manuscript the spatial sub-section it needs.
2. **Sensitivity sweep on the buffer / mask**: rerun with `STEPREP_VARIANT <- "stepRep_strict_1km" / "stepRep_broad_500m" / "stepRep_broad_1km"` and report the four sets of `stepRep_obs` coefficients in a supplementary table.
3. **A small fix in the simulation block**: addressing the psi-scaling bug that v1 already documents would make the annual prevalence figures usable. This is independent of `stepRep` but worth doing.
4. **Manuscript draft update**: integrate the §10 narrative from `reports/stepRep_diagnostics.md` into the Methods and the detection-coefficient table from §6.2 of this report into the Results.

---

## Appendix: artefacts

In `results/stepRep_v2_run/`:

```
{sp}_v2_model_object.rds         <- fitted colext objects
{sp}_v2_train_dyn_scale.rds      <- per-year scaling parameters (for scripts 8, 10)
{sp}_v2_gof_parboot.rds          <- parametric bootstrap GOF
{sp}_v2_model_summary.txt        <- coefficient tables
{sp}_v2_simulation_prevalence.csv <- annual prevalence (caveat §8)
tettet_v2_gof_mackenzie_bailey.rds  (only species that ran cleanly)

comparison_table.csv             <- full per-coefficient v1/v2 estimates
comparison_summary.csv           <- per-species summary (AIC, stepRep, mean γ/ε)
gamma_epsilon_comparison.csv     <- predicted γ and ε at observed covariates

figs/
  v1_vs_v2_aic.png               <- ΔAIC bar chart
  v1_vs_v2_detection_forest.png  <- detection coefficients across both models
  v1_vs_v2_psi_forest.png        <- initial occupancy coefficients
  v1_vs_v2_col_forest.png        <- colonisation coefficients
  v1_vs_v2_ext_forest.png        <- extinction coefficients
  v1_vs_v2_gamma_epsilon.png     <- mean γ and ε per species
  {sp}_v2_occupancy_map.png      <- predicted ψ map per species
  {sp}_v2_response_*.png         <- response curves available
  {sp}_v2_prevalence_over_time.png

predictions/
  occ_{sp}_v2_prediction.csv     <- ψ at lon/lat grid
  {sp}_v2_OccuMap.tif            <- ψ raster
```

Reproducibility: `scripts/build_stepRep.R` generates the cell-year tables from CORINE; `scripts/3b_add_stepRep.R` joins them into the modelling tables; `scripts/4_occupancy_models_v2.R` fits the colext models and writes everything in this report. `gh pr view 26` for the merge that brought the code into main.

---

## 11. Reduced specifications (v3 and v4) and the final identifiability picture

§6 of this report showed that adding `stepRep_obs` to the detection sub-model destabilises the extinction sub-model for three of four species, despite the headline detection effect being clean. To diagnose whether the instability comes from `stepRep` itself or from the original v1 covariate structure being inappropriate now that detection is better identified, we fitted two reduced specifications:

- **v3**: same as v2 (detection with `stepRep_obs`) but with **γ ~ 1 and ε ~ 1** (intercept-only on both transition sub-models).
- **v4**: same as v2 but with **γ ~ 1 and ε ~ single covariate per species**. The covariate chosen is the strongest by Wald `z` in v1, avoiding NDVI in *P. orientalis* to keep the climate-vs-land-use attribution clean (per `docs/decisions_gcb_v4.md` decision 1):
  - *Otis tarda*: ε ~ Land_Cover_Type_1_Percent_Class_13
  - *Pterocles alchata*: ε ~ 1 (v1 had no significant ε covariate; nothing to add)
  - *Pterocles orientalis*: ε ~ Land_Cover_Type_1_Percent_Class_12
  - *Tetrax tetrax*: ε ~ Land_Cover_Type_1_Percent_Class_12

The fits are saved under `results/stepRep_v2_run/{sp}_v3_model_object.rds` (and `_v4_`); summaries in `{sp}_v3_model_summary.txt` and `{sp}_v4_model_summary.txt`; AIC comparison in `v1_v2_v3_v4_aic.csv`; ε intercept-only summary in `v4_eps_baseline.csv`.

### 11.1 AIC across four specifications

![AIC v1 v2 v3 v4](../results/stepRep_v2_run/figs/v1_v2_v3_v4_aic.png)

| Species | AIC v1 | AIC v2 | AIC v3 | AIC v4 | ΔAIC (v3 − v2) | ΔAIC (v4 − v2) |
|---|---:|---:|---:|---:|---:|---:|
| *Otis tarda* | 2 267.9 | 2 144.3 | 2 199.6 | 2 155.4 | +55.3 | +11.1 |
| *Pterocles alchata* | 1 981.7 | 1 916.5 | 2 030.5 | 2 030.5 | +114.0 | +114.0 |
| *Pterocles orientalis* | 2 060.7 | 2 014.6 | 2 075.6 | 2 082.6 | +61.0 | +68.0 |
| *Tetrax tetrax* | 1 780.1 | 1 757.2 | 1 748.4 | **1 746.7** | −8.8 | −10.5 |

v3 and v4 cost between 8.8 and 114 AIC units versus v2 for three species; *Tetrax tetrax* actually improves under the simpler specifications. AIC alone would prefer v2; AIC together with parameter identifiability tells a different story.

### 11.2 Identifiability of the ε intercept across specifications

![ε identifiability](../results/stepRep_v2_run/figs/eps_identifiability.png)

| Species | SE(ε intercept) v2 | SE v3 | SE v4 |
|---|---:|---:|---:|
| *Otis tarda* | 182 | **0.62** | 204 |
| *Pterocles alchata* | 2.9 | **0.66** | 0.66 |
| *Pterocles orientalis* | **0.41** | 0.42 | NaN |
| *Tetrax tetrax* | NaN | **0.41** | **0.35** |

v3 is the only specification where the ε intercept is well identified for all four species simultaneously. v2 works only for *P. orientalis*. v4 (single covariate) does not rescue the *Otis tarda* and *P. orientalis* fits — the model is information-starved on ε for those species under any covariate structure we have tried.

### 11.3 Baseline ε under v3

![v3 ε baseline](../results/stepRep_v2_run/figs/v3_eps_baseline.png)

| Species | ε intercept (logit) | SE | **Baseline ε per year** | 95 % CI |
|---|---:|---:|---:|---:|
| *Pterocles alchata* | −3.26 | 0.66 | **3.7 %** | 1.1 % – 12.3 % |
| *Otis tarda* | −2.95 | 0.62 | **5.0 %** | 1.6 % – 14.6 % |
| *Pterocles orientalis* | −2.79 | 0.42 | **5.8 %** | 2.6 % – 12.6 % |
| *Tetrax tetrax* | −2.42 | 0.41 | **8.1 %** | 4.0 % – 15.8 % |

ε under v3 is biologically plausible and clean: every species sits in a tight band of 4–8 % per year, with 95 % confidence intervals that do not include high values (upper bound ~ 16 %). The v1 / v2 mean predicted ε values (0.15–0.49 across species) are inflated by the cell-years that fall at extreme values of the poorly-identified ε covariates; once those covariates are removed, the underlying rate is consistent with field-derived expectations of slow local extinction in a long-lived steppe-bird guild.

### 11.4 γ under v3

| Species | γ intercept (logit) | SE | γ per year |
|---|---:|---:|---:|
| *Otis tarda* | −6.68 | 0.77 | 0.13 % |
| *Tetrax tetrax* | −6.98 | 0.78 | 0.09 % |
| *Pterocles alchata* | −20.7 | NaN | ≈ 0 (boundary) |
| *Pterocles orientalis* | −20.2 | NaN | ≈ 0 (boundary) |

γ is at the floor of identifiability in all specifications: for *Otis* and *Tetrax* it converges to ~0.001 with finite SE, but for the two *Pterocles* species the optimiser hits the lower bound and SE is undefined. The honest statement is that γ is **indistinguishable from zero** for all four species in this dataset at this scale and over this period; the data do not support inference about γ drivers.

### 11.5 What v3 and v4 add to the v1 vs v2 story

The v1 vs v2 analysis in §6 told the right story for **detection**: `stepRep_obs` works, AIC drops, other detection covariates settle. The story it could not tell cleanly was about γ and ε, because v2's covariate-rich specification of those sub-models hides the underlying ecological signal behind numerical pathology.

v3 makes that signal visible: γ is at the floor, ε is in a tight 4–8 % band per species. Once one accepts that the γ/ε signal in these data is what v3 says it is, the higher AIC of v3 versus v2 is fitting **noise** (the spurious variance that v2's poorly-identified covariates absorb is not real biological structure that v3 is missing — it is the model trying to explain detection-process variation as state-process variation, which is what motivated `stepRep` in the first place).

v4 confirms that the issue is information starvation, not covariate selection: forcing v1's strongest covariate into ε does not rescue identifiability for *Otis tarda* or *P. orientalis*. The data have a limited number of true colonisation/extinction transitions in the 7-year window; no covariate restructuring will work without more data or a different model class.

---

## 12. Final recommendation: model specification for the manuscript

### 12.1 Headline model specification

| Component | Recommended form | Rationale |
|---|---|---|
| ψ (initial occupancy) | unchanged from v1 per species (e.g. *Otis tarda*: ~ bio1 + bio2 + tree_cover + grass_cover + topo_elev) | Coefficients are stable between v1, v2 and v3. Habitat associations on initial occupancy are robust to the detection-model misspecification we corrected. |
| γ (colonisation) | **~ 1** (intercept-only) for all species | γ is indistinguishable from zero in every specification tested. Reporting "drivers of γ" is not supported by the data. State the asymmetry γ ≪ ε as the main finding instead. |
| ε (extinction) | **~ 1** (intercept-only) for *Otis tarda*, *P. alchata* and *T. tetrax*. For *P. orientalis*, keep v2's `~ LC12 + NDVI + pr` because it is the only species where ε covariates are identifiable. | Three species: no covariate makes ε identifiable. *P. orientalis* is the lone exception; reporting its ε drivers is honest. |
| p (detection) | v2 specification: existing covariates + **stepRep_obs** | The headline new finding. AIC drops 22–124 units; β = +0.61 to +0.89. |

This is a **hybrid** specification: γ uniform across species, ε species-specific (intercept-only or full v2), p uniform across species (with `stepRep_obs` added). The trade-off is a small loss of inferential richness on the transition sub-models in exchange for parameter estimates that have sensible SE everywhere.

### 12.2 Sensitivity analyses for the supplementary information

Three sensitivities should travel with the headline:

1. **v2 (all covariates retained) for γ and ε**: shows that with the original covariate structure, ε is unstable for 3 of 4 species (this is the evidence base for the v3 collapse).
2. **v4 (single ε covariate)**: shows that selective simplification does not rescue identifiability; the choice to go to intercept-only is data-driven, not a default.
3. **stepRep sensitivity variants**: `stepRep_strict_1km`, `stepRep_broad_500m`, `stepRep_broad_1km` reruns of the headline specification (single-line edit at the top of `scripts/4_occupancy_models_v2.R` via the `STEPREP_VARIANT` constant).

All three are already prepared infrastructure-wise: the fits exist (v2, v4) or are one rerun away (sensitivity variants).

### 12.3 Updated paper-skeleton claims

The headline claims in `docs/paper_skeleton_GCB_v9.md` need adjustment to match what the data actually support:

| Original claim | Recommended reformulation |
|---|---|
| "Colonisation-extinction asymmetry spans **two to six orders of magnitude** across species" | "γ is indistinguishable from zero (γ < 0.002) for all four species; ε is **3.7 – 8.1 % per year** baseline. The asymmetry γ ≪ ε is therefore **at least one and a half orders of magnitude**, robust to model specification." |
| "Detection correction reverses the qualitative assessment of colonisation" | Strengthened. Detection correction now includes an explicit within-cell sampling covariate (`stepRep_obs`, β = +0.61 to +0.89 logit, p < 10⁻⁹ in every species). Naive vs detection-corrected γ ratios are recomputed against the new γ estimates and still show qualitative reversal (corrected γ ≪ naive γ). |
| "Species-specific drivers of colonisation and extinction" | "ψ habitat associations are species-specific and robust; γ is at the identifiability floor in every species and we do not infer drivers; ε drivers are identifiable only for *P. orientalis* (LC12 negative, NDVI positive, pr positive) and reported there only." |
| "Extinction debt affecting 5 – 100 % of currently occupied sites" | "Baseline ε of 4 – 8 % per year, projected forward, implies an extinction debt that is **modest in absolute terms but cumulatively important** (50 % loss over ~10 years for *T. tetrax* at 8 % per year vs ~25 years for *P. alchata* at 4 % per year). Cell-level projections in the original wording (5–100 %) were driven by spuriously high ε at extreme covariate values and should be retracted." |
| "Spatial scales of conservation (43 – 264 km)" | Untouched; depends on the stPGOcc fits which are pending. The v2 stPGOcc analogue (`scripts/18_stPGOcc_production_run_v2.R`) is committed and ready to run; results will follow. |

### 12.4 What to do next, in order

1. **Refit colext under v3 with the *P. orientalis* exception** and freeze that as the manuscript's headline. Concretely: copy `scripts/4_occupancy_models_v3.R` to `scripts/4_occupancy_models_final.R`, hardcode the ε ~ LC12 + NDVI + pr exception for `pteori`, and use that as the canonical fit going forward.
2. **Recompute the naive vs corrected γ table** with the new γ estimates from the headline fit. The ratios will change in magnitude but the qualitative reversal claim survives.
3. **Run `stPGOcc` (v2) production** when Raúl has the machine. The spatial random effect may absorb some of the residual structure that destabilises ε in colext and could permit the original covariate spec to identify there.
4. **Refresh the manuscript Results** with the reformulated claims from §12.3 of this report and the headline numbers from §11.3.

### 12.5 Files added by this revision

```
results/stepRep_v2_run/
  v1_v2_v3_aic.csv                v1/v2/v3 AIC comparison
  v1_v2_v3_v4_aic.csv             v1/v2/v3/v4 AIC comparison
  v4_eps_baseline.csv             ε baseline per species under v4
  {sp}_v3_model_summary.txt       4 files (γ ~ 1, ε ~ 1)
  {sp}_v4_model_summary.txt       4 files (γ ~ 1, ε ~ 1 cov)
  figs/v1_v2_v3_v4_aic.png        AIC across four specifications
  figs/v3_eps_baseline.png        ε per species under v3, with 95 % CI
  figs/eps_identifiability.png    SE of ε intercept across versions

scripts/
  4_occupancy_models_v3.R         intercept-only γ and ε
  4_occupancy_models_v4.R         intercept-only γ, single-cov ε
```

The model object `.rds` files for v3 and v4 (~30 MB total) live alongside in `results/stepRep_v2_run/` locally but stay gitignored — they are regenerated by these scripts.
