# Technical Audit Report — Iberian Steppe Birds Occupancy Dynamics
**Version:** v4 — March 2026
**Target journal:** Global Change Biology
**Author:** Guillermo Fandos
**Repository:** `dynamic_occupancy_models_steppe_birds`
**Branch:** `audit-gcb-v4`

---

## Executive Summary

The analysis pipeline is substantially complete for a GCB submission, with four **blocking gaps** and two **critical methodological issues** that must be resolved:

| Block | Status | Blocking? |
|-------|--------|-----------|
| A. colext models | Complete | No |
| B. Detection & effort | Complete | No |
| C. Equilibrium occupancy psi* | **COMPUTED** (this audit) | Was blocking; now resolved |
| D. Attribution (scripts/8) | Preliminary; scaling bug found | Yes (scripts/10 missing) |
| E. NDVI collinearity diagnostics | Complete | Decision pending (Guillermo) |
| F. stPGOcc spatial models | Not converged | Yes (cluster) |
| G. Validation (blockCV) | Complete but block size too small | Yes (methodological) |
| H. Goodness of fit (parboot) | Failed (formlist slot error) | Yes (cluster) |
| I. Cross-analysis consistency | 1 bug found (scaling), 1 design issue (blockCV) | Yes |

**Key numbers for the Abstract (computed during this audit):**

| Species | psi* median | 95% CI | Recolonisation time (yr) |
|---------|-------------|--------|--------------------------|
| O. tarda | 0.02% | [0.0004%, 0.93%] | 22,648 [427 - 1,168,912] |
| P. alchata | 4.55% | [0.11%, 99.67%] | 680 [303 - 1,515] |
| P. orientalis | 0.0001% | [0%, 0.17%] | 2,599,962 [1,234 - 5.4 billion] |
| T. tetrax | 0.49% | [0.16%, 1.62%] | 1,471 [462 - 4,630] |

**Interpretation:** Three of four species (otitar, pteori, tettet) have psi* < 1% — confirming the demographic trap. ptealc's estimate is unreliable due to near-complete separation in epsilon (SE = 30-41). T. tetrax is the most precisely estimated and the clearest case for the paper's narrative.

---

## 1. Project Map

### 1.1 Directory Structure

```
scripts/
  1_prepare_database_ebird.R       # eBird data preparation
  2_prepare_static_variables.R      # WorldClim, MODIS static
  3_prepare_dynamic_variables.R     # GEE dynamic covariate extraction
  3_download_dynamic_variables_*.js # GEE JavaScript per species
  4_occupancy_models.R              # colext fitting (main pipeline)
  5_validation.R                    # blockCV validation vs atlas
  6_spatial_occupancy_test.R        # spPGOcc single-year + stPGOcc multi-year
  7_spatial_results.R               # stPGOcc results extraction
  8_counterfactual_attribution.R    # Factorial attribution (original models)
  8_effort_confounding.R            # Effort diagnostics
  9_collinearity_ndvi.R             # VIF, stability, NDVI decomposition
  compute_equilibrium.R             # NEW: psi* calculation (this audit)
  publication_*.R                   # Figure/table generation scripts
  test_pipeline_v*.R                # Development/testing scripts

R/
  model_configs.R                   # Centralised model formulas

results/
  {sp}_model_object.rds             # 4 fitted colext models
  {sp}_simulation_prevalence.csv    # Stochastic simulation outputs
  {sp}_model_summary.txt            # Model summaries
  {sp}_validation_summary.txt       # Validation results
  {sp}_validation_cv.csv            # Per-fold CV results
  equilibrium_occupancy_table.csv   # NEW: psi* with CI
  counterfactual_predictions.rds    # Factorial scenarios
  pub_table_attribution_summary.csv # Attribution cross-species table
  pub_table_attribution_revised.csv # 3-way NDVI decomposition
  reviewer_gof_parboot.csv          # GoF results (ALL FAILED)
  results_spatial/                  # stPGOcc + tPGOcc objects and diagnostics
  results2/                         # Intermediate model objects (v4b, v5)
  diagnostics/                      # VIF, stability, NDVI decomposition
  covariate_trends/                 # Linear trend rasters per covariate

data/processed_2023/{sp}/
  {sp}_occ_wide_dynamic.csv         # Site-year data matrices
  {sp}_scaling_params.rds           # Static covariate scaling

docs/
  master_analysis_report.md         # Authoritative synthesis (supersedes others)

manuscript/
  skeleton_GCB_v6.docx              # Paper skeleton
```

### 1.2 Model Objects Inventory

| File | Species | Pipeline | AIC | Status |
|------|---------|----------|-----|--------|
| results/otitar_model_object.rds | O. tarda | v4b | 2267.9 | Verified |
| results/ptealc_model_object.rds | P. alchata | v5 | 1973.9 | Verified (gamma/epsilon flagged) |
| results/pteori_model_object.rds | P. orientalis | v4b | 2060.7 | Verified (epsilon flagged) |
| results/tettet_model_object.rds | T. tetrax | v4b | 1780.1 | Verified (cleanest) |
| results/results_spatial/*_stPGOcc.rds | All 4 spp | stPGOcc | — | Not converged |
| results/results_spatial/*_tPGOcc.rds | All 4 spp | tPGOcc | — | Converged |
| results/pteori_epsilon_revised.rds | — | — | — | MISSING |
| results/otitar_gamma_revised.rds | — | — | — | MISSING |

### 1.3 Figures & Tables: Skeleton v6 Readiness

| Item | Description | Status |
|------|-------------|--------|
| Fig 1A | Study area map | Exists (pub_map_sampling_sites.png) |
| Fig 1B | psi* isoclines gamma vs epsilon | **COMPUTABLE NOW** (equilibrium table computed) |
| Fig 2 | Initial occupancy maps | Exists (pub_map_occupancy_4species.png) |
| Fig 3 | Colonisation coefficient forest plot with bootstrap CI | Exists (pub_fig2_coefficient_forest_plot.png) but **bootstrap CI pending** |
| Fig 4 | Extinction coefficient forest plot | Exists (pub_fig_coef_all_submodels.png) |
| Fig 5 | Predicted gamma/epsilon maps | Exists (pub_fig_spatial_col_ext_rates.png) |
| Fig 6 | Prevalence trends | Exists (pub_fig_occupancy_trends_panel.png) |
| Fig 7 | Spatial diagnostics | Exists (pub_fig_spatial_moran.png) but **stPGOcc not converged** |
| Fig 8 | Covariate trend maps | Exists (covariate_trends/*.tif) |
| Fig 9 | Attribution summary | Exists but **no bootstrap CI; scaling bug** |
| T1 | Species overview | Content available |
| T2 | Covariate descriptions | Content available |
| T3 | Model formulas | Verified in model_configs.R |
| T4 | Detection coefficients | Exists (pub_table4_detection.csv) |
| T5 | Gamma/epsilon coefficients | Exists (pub_table3_coefficients_gamma_epsilon.csv) |
| T7 | Model fit / validation | AUC/TSS exists; **GoF missing** |
| T8 | Spatial model comparison | Exists (pub_spatial_table1_model_comparison.csv) |
| T9 | Attribution summary | Exists but **preliminary (no CI, scaling bug)** |

---

## 2. Audit by Block

### BLOQUE A — COLEXT MODELS

**AIC verification:** All 4 match expectations (otitar=2267.9, ptealc=1973.9, pteori=2060.7, tettet=1780.1).

**Formula verification:** All formulas match R/model_configs.R exactly.

**WorldClim in gamma/epsilon:** NONE found. Correct two-scale design confirmed.

**Coefficient audit — flags by severity:**

**CRITICAL:**
- ptealc epsilon: SE(pr) = 30.67, SE(tmmx) = 40.79 — **quasi-complete separation**. Only ~18 extinction events. Coefficients (-52, -71) are numerically meaningless. The entire epsilon submodel for ptealc is not interpretable.
- pteori gamma(LC7): SE = 6.46, estimate = -9.80 — **probable separation** for open shrubland.

**MODERATE (|z| < 1.5, not interpretable for inference):**
- otitar epsilon(LC6): z = 1.21
- otitar gamma(pr): z = -1.56 (borderline)
- otitar gamma(tmmx): z = 1.54 (borderline)
- ptealc gamma(pr): z = 0.45
- ptealc psi(topo_aspect): z = -0.50
- pteori gamma(NDVI): z = -1.14
- pteori gamma(LC7): z = -1.52
- tettet psi(topo_elev): z = -0.64

**ECOLOGICAL SIGN CHECK:**
- NDVI negative for colonization (otitar z=-2.49, ptealc z=-3.65): counterintuitive but consistent — NDVI in Mediterranean steppe may reflect agricultural intensification. Requires explicit justification in Discussion.
- pteori epsilon(NDVI) = +2.05 (z=3.68): higher NDVI increases extinction. Same framing as above.
- pteori epsilon(pr) = +1.69 (z=2.30): more precipitation increases extinction. Interpretable as wet years driving agricultural intensification.

**vcov finite check:** All 4 species PASS.

**ptealc/gamma collinearity flag status:** NDVI + pr confirmed in current model. The SIGN_CHANGE in pr (from +0.162 to -0.717 when NDVI removed) is verified in coefficient_stability.csv.

**pteori/epsilon collinearity flag status:** LC12 + NDVI + pr confirmed. BETA_SHIFT (LC12: -0.831 to -0.095, shift = 2.8 SE) verified. results/pteori_epsilon_revised.rds does NOT exist.

---

### BLOQUE B — DETECTION AND EFFORT

**Effort confounding verified:**
- O. tarda: no confounding (rho = -0.25, P = 0.59)
- P. alchata: **marginal** (GLM effort P = 0.050). Correctly reported as cautionary.
- P. orientalis: no confounding (rho = 0.64, P = 0.12 — not significant)
- T. tetrax: clearest decoupling (rho = 0.018, P = 0.97)

**Detection covariate signs:** All correct.
- effort/duration: positive (more time = more detection)
- time (hour): negative (dawn/dusk better)
- observers: positive (but weak for ptealc z=1.48, tettet z=1.32)

**Naive gamma vs corrected gamma comparison: DOES NOT EXIST.** This is needed for Discussion 4.4 to demonstrate that detection correction changes conclusions. GAP.

---

### BLOQUE C — EQUILIBRIUM OCCUPANCY psi*

**STATUS: COMPUTED DURING THIS AUDIT.** Output: `results/equilibrium_occupancy_table.csv`

| Species | gamma baseline | epsilon baseline | psi* median | psi* 95% CI | Recol. time |
|---------|---------------|-----------------|-------------|-------------|-------------|
| O. tarda | 0.0045% | 22.24% | 0.02% | [0.0004%, 0.93%] | 22,648 yr |
| P. alchata | 0.15% | 3.05% | 4.55% | [0.11%, 99.67%] | 680 yr |
| P. orientalis | ~0% | 44.06% | 0.0001% | [0%, 0.17%] | 2.6M yr |
| T. tetrax | 0.07% | 13.75% | 0.49% | [0.16%, 1.62%] | 1,471 yr |

**Interpretation notes:**
- ptealc has absurdly wide CI (0.11% to 99.67%) because the epsilon intercept (-3.46) has SE = 4.45 (near-separation). The bootstrap draws include extreme epsilon values near 0, producing psi* near 100%. **This species should be excluded from psi* reporting** or reported with explicit caveat.
- pteori gamma intercept = -14.78 is extremely low (baseline gamma essentially zero). This is the most extreme demographic trap.
- **For the Abstract:** "equilibrium occupancy psi* = 0.02-0.49% for three species (excluding P. alchata due to estimation instability), implying expected recolonisation times of 1,500-2.6 million years."

---

### BLOQUE D — ATTRIBUTION

**Scaling bug found (SEVERITY: HIGH):**
scripts/8 recomputes scaling parameters on drop_na-filtered rows. scripts/4 computed scaling on ALL rows before filtering. The center and scale values differ. This means the counterfactual predictions use slightly different covariate standardization than the model was trained on.

**Fix:** scripts/8 should load the train_dyn_scale from scripts/4, or use a shared scaling parameters file. The stochastic simulation section of scripts/4 (lines 512-545) correctly uses train_dyn_scale — scripts/8 missed this.

**ptealc gamma exclusion:** NOT implemented in scripts/8. Gamma attribution is computed and reported for ptealc even though only 16 colonisation events make it unreliable. The master report says to exclude it.

**NDVI classification:** Consistent — classified as climate-adjacent in scripts/8 line 69.

**scripts/10 status:** DOES NOT EXIST. Missing outputs:
- results/ndvi_decomp_objects.rds — MISSING
- results/attribution_boot_summary.csv — MISSING
- results/pteori_epsilon_revised.rds — MISSING
- results/otitar_gamma_revised.rds — MISSING
- results/attribution_scenarios.rds — MISSING

**Existing but preliminary:** results/pub_table_attribution_revised.csv exists with 3-way NDVI decomposition columns. Provenance unclear — may have been computed ad hoc.

---

### BLOQUE E — NDVI COLLINEARITY (verification only)

All diagnostic outputs exist and are internally consistent:

| File | Exists | Key finding |
|------|--------|-------------|
| results/diagnostics/vif_summary.csv | Yes | All VIF < 5. Max = 4.43 (otitar/gamma/tmmx) |
| results/diagnostics/coefficient_stability.csv | Yes | ptealc/gamma SIGN_CHANGE confirmed; pteori/epsilon BETA_SHIFT confirmed |
| results/diagnostics/ndvi_climate_component.rds | Yes | — |
| results/diagnostics/ndvi_residual_component.rds | Yes | — |
| results/diagnostics/correlation_flags.csv | Yes | 5 SEVERE pairs (pr-tmmx, tmmn-tmmx); 6 MODERATE |
| results/diagnostics/ndvi_model_revision_recommendations.csv | Yes | ptealc/gamma: REMOVE NDVI; pteori/epsilon: REMOVE NDVI |

**NDVI decomposition case study values:** +0.007 (NDVI_climate -> extinction pteori) and -0.011 (NDVI_landuse -> extinction). Bootstrap CI: **NOT computed.** Needed before reporting even as "preliminary."

**Pending decision (Guillermo):** Option A (retain) vs Option B (remove) for ptealc/gamma and pteori/epsilon. This blocks scripts/10.

---

### BLOQUE F — SPATIAL MODELS (stPGOcc)

**WAIC comparison (stPGOcc vs tPGOcc):** Confirmed DELTA_WAIC > 350 for all species.

| Species | tPGOcc WAIC | stPGOcc WAIC | DELTA_WAIC |
|---------|-------------|-------------|------------|
| O. tarda | 2496.9 | 2040.4 | 456.5 |
| P. alchata | 2103.8 | 1655.2 | 448.6 |
| P. orientalis | 2273.7 | 1821.3 | 452.4 |
| T. tetrax | 1890.4 | 1538.8 | 351.6 |

**Moran's I reduction:** Confirmed 118-141% reduction (tPGOcc: 0.14-0.24 -> stPGOcc: -0.04 to -0.06).

**Convergence diagnostics — CATASTROPHIC for 3 species:**

| Species | phi Rhat | sigma2_t Rhat | Intercept Rhat | Classification |
|---------|----------|---------------|----------------|----------------|
| O. tarda | 1.26 | **4.19** | **3.32** | NOT PUBLISHABLE |
| P. alchata | **1.61** | **3.45** | **7.05** | NOT PUBLISHABLE |
| P. orientalis | **5.81** | 1.46 | **5.56** | NOT PUBLISHABLE |
| T. tetrax | 1.31 | **1.84** | 1.16 | USABLE WITH CAUTION |

**Current MCMC settings:** n.iter = 20,000; burn-in = 10,000; thin = 20 -> 500 post-burn samples/chain. Far too few.

**Root cause:** 500 post-thinning samples per chain is insufficient for NNGP + AR(1) models. phi and sigma2_t are slow-mixing parameters that require 5-10x more samples.

**Spatial ranges (PRELIMINARY):** otitar=46.7 km, ptealc=263.9 km, pteori=171.6 km, tettet=42.7 km. These are point estimates only; CIs unreliable due to non-convergence.

**Recommended cluster settings:**
```r
# RAUL: Run in cluster
stPGOcc(
  n.batch      = 4000,       # 4000 x 25 = 100,000 total
  batch.length = 25,
  n.burn       = 50000,
  n.thin       = 20,         # -> 2,500 post-burnin samples
  n.chains     = 3,
  NNGP         = TRUE,
  n.neighbors  = 5,
  ar1          = TRUE
)
# Consider informative phi priors:
# priors.list <- list(phi.unif = c(3/500, 3/10))  # range 10-500 km
# Estimated runtime: 6-24 hours per species
```

---

### BLOQUE G — VALIDATION

**Design:** Spatially-blocked 5-fold using blockCV::cv_spatial(). Correct.

**CRITICAL ISSUE: Block size = 50 km << 264 km (P. alchata spatial range).**

The blockCV block size should be >= the largest spatial autocorrelation range to ensure spatial independence between folds. At 50 km, training and test folds share correlated spatial structure. This makes:
- AUC values potentially optimistic
- Spearman rho P-values anticonservative

**Spearman rho:** Computed both raw (full dataset) and out-of-fold (per-fold mean). Values are very close (difference < 0.01), which itself may reflect the inadequate block size.

**Calibration curves:** Exist for all 4 species (decile observed vs predicted).

**Validation results (potentially optimistic due to block size):**

| Species | AUC | TSS | Spearman rho (CV mean +/- SD) |
|---------|-----|-----|-------------------------------|
| O. tarda | 0.884 | 0.612 | 0.402 +/- 0.073 |
| P. alchata | 0.844 | 0.569 | 0.330 +/- 0.056 |
| P. orientalis | 0.823 | 0.536 | 0.411 +/- 0.051 |
| T. tetrax | 0.850 | 0.546 | 0.523 +/- 0.053 |

---

### BLOQUE H — GOODNESS OF FIT

**STATUS: ALL FAILED.**

parboot failed for all 4 species with error: `"no hay un slot de nombre 'formlist' para ese objeto de clase 'unmarkedFitColExt'"`. This is likely a serialization/version incompatibility — the model objects saved via saveRDS() lose the `formlist` slot when reloaded.

MacKenzie-Bailey GOF failed for 3 of 4 species (otitar, ptealc, pteori) due to row-count mismatch (AICcmodavg bug with 7-year data).

**Current nsim:** 10 (defined at scripts/4_occupancy_models.R line 63). Not reportable regardless.

**c-hat:** Never computed because parboot never succeeded.

**Fix required:**
```r
# RAUL: In scripts/4_occupancy_models.R
# Line 63: change to
NSIM_PARBOOT <- 500
# Line 64: change to
NSIM_MB_GOF  <- 500

# CRITICAL: Models must be re-fitted in the SAME R session where parboot
# runs — do NOT load from RDS. The formlist slot is not preserved during
# serialization. Alternatively, update unmarked to latest version.
# Test with nsim = 5 first.

# Runtime estimate: ~50x current time per species
```

---

### BLOQUE I — CROSS-ANALYSIS CONSISTENCY

**I-1: Scaling parameters (scripts/4 vs scripts/8):**
**BUG CONFIRMED.** scripts/4 computes scaling on ALL rows. scripts/8 computes scaling on drop_na-filtered rows. Different center/scale values. The train_dyn_scale from scripts/4 is NOT passed to scripts/8. Fix: scripts/8 should use saved scaling parameters from scripts/4.

**I-2: blockCV block size vs spatial ranges:**
**ISSUE CONFIRMED.** Block size = 50 km. P. alchata spatial range = 264 km. Folds are not spatially independent. Must increase to >= 264 km (or justify the current choice).

**I-3: Moran's I colext vs stPGOcc:**
**CORRECTLY DIFFERENT.** colext residual Moran's I = 0.358-0.491 (from validation script, on atlas polygons). tPGOcc residual Moran's I = 0.144-0.236 (from stPGOcc comparison, on grid cells). Different measurement units, different models, different spatial support. No error.

**I-4: Attribution formulas vs model_configs.R:**
**NO DIVERGENCE.** scripts/8 correctly sources R/model_configs.R and uses cfg$gamma_vars / cfg$epsilon_vars for predict(). Verified.

---

## 3. Inconsistencies Detected

### Previously identified (confirmed):
- **IC-1:** spPGOcc vs stPGOcc spatial ranges — confirmed, skeleton uses stPGOcc (correct).
- **IC-2:** Attribution uses flagged models — confirmed, ptealc/pteori attribution unreliable.
- **IC-3:** GoF "completed" vs "not reportable" — confirmed, nsim=10, all parboot failed.
- **IC-4:** DELTA_WAIC numbers differ by report — confirmed, skeleton uses stPGOcc numbers (correct).

### Newly identified:
- **IC-5:** Scaling parameter mismatch between scripts/4 and scripts/8 (Bloque I-1).
- **IC-6:** blockCV block size (50 km) inconsistent with largest spatial range (264 km) (Bloque I-2).
- **IC-7:** ptealc psi* unreportable — CI spans 0.1% to 99.7% due to epsilon separation. Not flagged in master report.
- **IC-8:** ptealc gamma not excluded in scripts/8 attribution despite master report instruction.
- **IC-9:** Naive vs corrected gamma comparison missing — needed for Discussion 4.4.
- **IC-10:** parboot failure is a serialization bug, not just a sample size issue — fixing nsim alone won't resolve it.

---

## 4. Action Plan (Prioritised)

### LEVEL 1 — BLOCKING (without this, no submission)

| # | Item | Unblocks | Who | Output |
|---|------|----------|-----|--------|
| 1.1 | **Decision: ptealc/gamma NDVI removal** | scripts/10, attribution | Guillermo | Decision recorded |
| 1.2 | **Fix scaling bug in scripts/8** (use train_dyn_scale from scripts/4) | Attribution correctness | Guillermo Fandos | Corrected scripts/8 or scripts/10 |
| 1.3 | **Refit pteori/epsilon without NDVI** (Task 1 of scripts/10) | Attribution, Table 9 | Raul or Guillermo Fandos | results/pteori_epsilon_revised.rds |
| 1.4 | **Create scripts/10_attribution_revised.R** with Tasks 0-5 | Results 3.7, Fig 9 | Guillermo Fandos + Raul (cluster for bootstrap) | Multiple outputs |
| 1.5 | **stPGOcc convergence** (100,000 iterations) | Results 3.5, Fig 7 | Raul (cluster) | Rhat < 1.3 |
| 1.6 | **parboot GoF** (fix formlist bug + nsim=500) | Methods 2.4 | Raul (cluster) | results/gof/pb_*_nsim500.rds |
| 1.7 | **Increase blockCV block size to >= 264 km** | Validation integrity | Guillermo Fandos | Re-run scripts/5 |

### LEVEL 2 — STRONGLY RECOMMENDED (without this, review will be harsh)

| # | Item | Why |
|---|------|-----|
| 2.1 | Naive vs corrected gamma comparison | Discussion 4.4 — proves detection correction changes conclusions |
| 2.2 | Fig 1B: isocline plot gamma vs epsilon | Core figure — now computable with equilibrium table |
| 2.3 | NDVI decomposition bootstrap CIs for pteori/epsilon | Discussion 4.2 case study reportability |
| 2.4 | Exclude ptealc from psi* reporting (or add explicit caveat) | CI of 0.1-99.7% is not informative |
| 2.5 | Bootstrap CI for gamma/epsilon coefficients (Fig 3, Fig 4) | Currently showing SE-based CI only |

### LEVEL 3 — OPTIONAL BUT VALUABLE

| # | Item | Value |
|---|------|-------|
| 3.1 | Delta-gamma required for psi* = 10% | Discussion 4.1 — recovery threshold |
| 3.2 | Informative phi priors for stPGOcc | May improve convergence |
| 3.3 | Photovoltaic overlap analysis | Policy relevance |
| 3.4 | ptealc epsilon: consider intercept-only model | Currently not interpretable |
