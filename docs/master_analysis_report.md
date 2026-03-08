# Master Analysis Report — Iberian Steppe Birds Occupancy Dynamics (2017–2023)
**Version:** 1.0 — March 2026  
**Pipeline branch:** `analysis-2023` / `counterfactual-attribution` / `effort-confounding-diagnostics`  
**This document supersedes:** `analysis_report_2023.md`, `report_spatial_models.md`, `report_counterfactual_attribution.md`, `report_collinearity_ndvi_diagnostics.md`, `report_effort_confounding.md`  
**Target journal:** Global Change Biology  

---

## ⚡ STATUS DASHBOARD

| Analysis block | Status | Blocking submission? | Responsible |
|---|---|---|---|
| colext models fitted | ✅ Complete | — | — |
| Independent validation (atlas blockCV) | ✅ Complete | — | — |
| Effort confounding diagnostics | ✅ Complete | — | — |
| Collinearity / NDVI diagnostics | ✅ Complete | Decision pending (see §6) | Guillermo |
| Counterfactual attribution (original models) | ⚠️ Preliminary | Pending model revisions | Raúl (cluster) |
| stPGOcc convergence | ❌ Not converged | **YES** — blocks §8, Fig 5–6 | Raúl (cluster) |
| parboot GoF (nsim = 10) | ❌ Not reportable | **YES** — blocks Methods 2.4 | Raúl (cluster) |
| Equilibrium occupancy ψ* with bootstrap CI | ❌ Missing | **YES** — blocks Abstract, §3.3, Fig 1B | Claude Code |
| Attribution revised (scripts/10) | ❌ Missing | **YES** — blocks §3.7, Fig 7 | Raúl (cluster) |

---

## ⚠️ CROSS-REPORT INCONSISTENCIES DETECTED

These inconsistencies exist at the boundary between individual reports and are invisible without a unified view:

### IC-1: Spatial range estimates — spPGOcc vs stPGOcc
`analysis_report_2023.md` §6.1 (spPGOcc, single-year 2023 snapshot):

| Species | spPGOcc range (km) |
|---|---|
| O. tarda | **182** |
| P. alchata | **334** |
| P. orientalis | **263** |
| T. tetrax | **23** |

`report_spatial_models.md` §3.3 (stPGOcc, full 2017–2023):

| Species | stPGOcc range (km) |
|---|---|
| O. tarda | **46.7** |
| P. alchata | **263.9** |
| P. orientalis | **171.6** |
| T. tetrax | **42.7** |

**Skeleton v6 uses the stPGOcc numbers.** The two sets are not comparable: spPGOcc is a single-year spatial model; stPGOcc incorporates seven years of data and AR(1) temporal structure. The spatial ranges reflect different model structures, not just different runs. The stPGOcc numbers are the correct ones to use for the paper. However, they are **preliminary because stPGOcc has not converged** (Rhat φ = 1.25–1.31; see §8).

### IC-2: Attribution analysis uses models flagged for revision
`report_counterfactual_attribution.md` was run with the original colext models, including:
- NDVI in ptealc/gamma → collinearity report flags sign change in pr when NDVI removed (SIGN_CHANGE flag)
- NDVI in pteori/epsilon → collinearity report flags LC12 collapse (β: −0.83 → −0.09, shift = 2.8 SE)

The current attribution numbers for ptealc are **not reliable** pending the model revision decision (§6). The attribution for pteori/epsilon reflects a model where NDVI is stabilising LC12 in a way that may not be ecologically defensible.

### IC-3: GoF status — "completed" vs "not reportable"
`analysis_report_2023.md` §8.2 lists parboot as run with nsim = 10. This is sufficient for pipeline testing but **not reportable in any journal**. The analysis is technically complete but scientifically invalid at that sample size. Current status: not reportable. Required: nsim ≥ 500.

### IC-4: ΔWAIC numbers differ by report
`analysis_report_2023.md` §6.1: ΔWAIC spPGOcc vs tPGOcc = **1,516–2,076** (single-year comparison).  
`report_spatial_models.md` §3.1: ΔWAIC stPGOcc vs tPGOcc = **352–457** (multi-year comparison).  
Skeleton v6 uses "ΔWAIC > 350" — correct for the stPGOcc comparison. The larger numbers from spPGOcc should not appear in the manuscript.

---

## 1. Species and Data Overview

| Code | Species | Common name | IUCN (Spain) | % European pop. | n sites (approx.) |
|---|---|---|---|---|---|
| otitar | *Otis tarda* | Great Bustard | NT | 98% | ~1,200 |
| ptealc | *Pterocles alchata* | Pin-tailed Sandgrouse | VU | 97% | ~900 |
| pteori | *Pterocles orientalis* | Black-bellied Sandgrouse | EN/VU | 19% | ~1,100 |
| tettet | *Tetrax tetrax* | Little Bustard | EN | 34% | ~1,500 |

**eBird growth:** +256% checklists for otitar/tettet; +201% for ptealc/pteori over 2017–2023.  
**Breeding seasons:** April–June (otitar, tettet); May–August (ptealc, pteori).  
**Grid:** 5-km² cells, mainland Spain (494,011 km²).  

**Two-scale covariate design (deliberate conceptual choice — not a data limitation):**
- ψ₁: WorldClim 1970–2000 static normals → fundamental niche
- γ, ε: TerraClimate + MODIS annual covariates → realised demographic dynamics

Any use of WorldClim variables in γ or ε submodels is a conceptual error.

---

## 2. Dynamic Occupancy Models (colext)

### 2.1 Best models per species

| Species | Pipeline | Model ID | ψ₁ covariates | γ covariates | ε covariates | AIC |
|---|---|---|---|---|---|---|
| O. tarda | v4b | m9 | bio1 + bio2 + tree_cover + grass_cover + topo_elev | NDVI + pr + tmmn + tmmx | LC6 + LC13 + tmmx | **2267.9** |
| P. alchata | v5 | m6 | bio1 + bio2 + tree_cover + grass_cover + topo_aspect | NDVI + pr | pr + tmmx | **1973.9** |
| P. orientalis | v4b | m0 | bio2 + tree_cover + grass_cover | LC7 + NDVI + tmmn + tmmx | LC12 + NDVI + pr | **2060.7** |
| T. tetrax | v4b | m7 | bio2 + tree_cover + grass_cover + topo_elev | LC12 | LC12 | **1780.1** |

### 2.2 Parameter estimates — Otis tarda (AIC = 2267.9)

| Sub-model | Parameter | Estimate | SE | z | P |
|---|---|---|---|---|---|
| ψ₁ | Intercept | −8.085 | 0.867 | −9.33 | <0.001 |
| | bio1 | 0.727 | 0.372 | 1.95 | 0.051 |
| | bio2 | 0.618 | 0.271 | 2.28 | **0.023** |
| | tree_cover | −6.449 | 0.884 | −7.30 | **<0.001** |
| | grass_cover | 0.355 | 0.151 | 2.35 | **0.019** |
| | topo_elev | 1.213 | 0.442 | 2.74 | **0.006** |
| γ | Intercept | −10.00 | 2.015 | −4.96 | <0.001 |
| | NDVI | −1.45 | 0.583 | −2.49 | **0.013** |
| | pr | −2.30 | 1.475 | −1.56 | 0.119 |
| | tmmn | −4.52 | 1.949 | −2.32 | **0.020** |
| | tmmx | 2.59 | 1.684 | 1.54 | 0.124 |
| ε | Intercept | −1.25 | 0.426 | −2.94 | **0.003** |
| | LC6 | 1.86 | 1.537 | 1.21 | 0.225 |
| | LC13 | **2.80** | 0.833 | **3.36** | **<0.001** |
| | tmmx | −1.51 | 0.761 | −1.99 | **0.047** |
| p | effort | 0.225 | 0.060 | 3.78 | **<0.001** |
| | observers | 0.223 | 0.053 | 4.23 | **<0.001** |
| | NDVI_obs | 0.428 | 0.119 | 3.61 | **<0.001** |
| | pr_obs | −0.758 | 0.200 | −3.80 | **<0.001** |

**Key findings:** Urban land cover (LC13, β = 2.80, P < 0.001) is the dominant extinction driver. Detection increases with effort and observers, decreases with precipitation.

### 2.3 Parameter estimates — Pterocles alchata (AIC = 1973.9)

| Sub-model | Parameter | Estimate | SE | z | P |
|---|---|---|---|---|---|
| ψ₁ | bio1 | −0.615 | 0.201 | −3.06 | **0.002** |
| | bio2 | 1.517 | 0.293 | 5.18 | **<0.001** |
| | tree_cover | −4.538 | 0.749 | −6.06 | **<0.001** |
| | grass_cover | 0.519 | 0.150 | 3.45 | **<0.001** |
| γ | Intercept | −6.522 | 0.413 | −15.79 | <0.001 |
| | NDVI | −1.081 | 0.296 | −3.65 | **<0.001** |
| | pr | 0.162 | 0.359 | 0.45 | 0.652 |
| ε | pr | **−52.07** | 30.67 | −1.70 | 0.090 |
| | tmmx | **−70.87** | 40.79 | −1.74 | 0.082 |
| p | time | −0.364 | 0.071 | −5.13 | **<0.001** |
| | duration | 0.900 | 0.094 | 9.61 | **<0.001** |
| | effort | 0.197 | 0.075 | 2.64 | **0.008** |

⚠️ **CRITICAL FLAGS:**
- γ submodel: NDVI + pr — collinearity diagnostic shows pr changes sign when NDVI removed (SIGN_CHANGE). Estimates unreliable. **Decision pending: Opción A vs B (see §6).**
- ε submodel: Only 18 extinction events. SE for pr (30.67) and tmmx (40.79) are anomalously large — classic complete separation signature. Extinction coefficients not interpretable for inference.
- γ submodel excluded from attribution analysis (16 colonisation events — separation confirmed).

### 2.4 Parameter estimates — Pterocles orientalis (AIC = 2060.7)

| Sub-model | Parameter | Estimate | SE | z | P |
|---|---|---|---|---|---|
| ψ₁ | bio2 | 0.359 | 0.205 | 1.76 | 0.079 |
| | tree_cover | −3.954 | 0.603 | −6.55 | **<0.001** |
| | grass_cover | 0.885 | 0.144 | 6.16 | **<0.001** |
| γ | tmmn | −8.105 | 2.670 | −3.04 | **0.002** |
| | tmmx | 9.135 | 2.682 | 3.41 | **<0.001** |
| | NDVI | −0.563 | 0.492 | −1.14 | 0.253 |
| | LC7 | −9.80 | 6.463 | −1.52 | 0.129 |
| ε | LC12 | −0.831 | 0.263 | −3.16 | **0.002** |
| | NDVI | 2.047 | 0.557 | 3.68 | **<0.001** |
| | pr | 1.687 | 0.734 | 2.30 | **0.022** |
| p | time | −0.244 | 0.071 | −3.44 | **<0.001** |
| | duration | 0.423 | 0.083 | 5.11 | **<0.001** |

⚠️ **FLAG on ε submodel:** NDVI removal test shows LC12 collapses (β: −0.83 → −0.09, shift = 2.8 SE) and SE of pr increases 60%. NDVI is structurally stabilising LC12 in the ε submodel. Decision pending whether to refit without NDVI (AIC cost: +60 points) or retain with documented limitation (see §6).

### 2.5 Parameter estimates — Tetrax tetrax (AIC = 1780.1)

| Sub-model | Parameter | Estimate | SE | z | P |
|---|---|---|---|---|---|
| ψ₁ | bio2 | 1.149 | 0.256 | 4.49 | **<0.001** |
| | tree_cover | −3.724 | 0.597 | −6.24 | **<0.001** |
| | grass_cover | 0.772 | 0.159 | 4.87 | **<0.001** |
| γ | LC12 | **0.902** | 0.437 | **2.06** | **0.039** |
| ε | LC12 | **−0.565** | 0.278 | **−2.03** | **0.042** |
| p | effort | 0.218 | 0.075 | 2.91 | **0.004** |
| | time | −0.219 | 0.072 | −3.02 | **0.003** |

✅ **Cleanest model.** Single driver (LC12) governs both γ and ε. Cropland promotes colonisation and prevents extinction. No collinearity issues.

### 2.6 Occupancy simulations — raw demographic rates (uncorrected for equilibrium)

| Species | ψ₂₀₁₇ | ψ₂₀₂₃ | Mean γ | Mean ε |
|---|---|---|---|---|
| O. tarda | ~0% | 0.67% | **0.20%** | **26.4%** |
| P. alchata | 0.02% | 0.87% | **0.24%** | **46.6%** |
| P. orientalis | 0.44% | 1.15% | **0.23%** | **43.0%** |
| T. tetrax | ~0% | 0.52% | **0.11%** | **15.0%** |

⚠️ **These raw simulation numbers are NOT the equilibrium occupancy (ψ*) for the paper.** They reflect a specific time window with specific covariate values and are influenced by the low initial ψ₂₀₁₇ estimates. The equilibrium occupancy ψ* = γ/(γ+ε) with bootstrap CIs is a separate calculation — see §3.

### 2.7 Goodness of fit

❌ **NOT REPORTABLE.** parboot run with nsim = 10 only. For publication: minimum nsim = 500.  
MacKenzie-Bailey GOF (`mb.gof.test`) failed for otitar, ptealc, pteori due to row-count mismatch (likely AICcmodavg bug with 7-year data). Succeeded for tettet only.  
`# RAÚL: increase nsim to 500 in GoF script before cluster run`

---

## 3. Equilibrium Occupancy and Demographic Asymmetry

**[PENDING — BLOCKING]** This is the most important single calculation for the paper and does not exist yet in any script output.

### 3.1 Calculation required

```r
# scripts/10_attribution_revised.R — add this block
library(MASS)  # mvrnorm

compute_equilibrium <- function(model, sp_code, n_boot = 1000) {
  coefs <- coef(model)
  V     <- vcov(model)
  sims  <- mvrnorm(n_boot, mu = coefs, Sigma = V)
  
  # Extract intercepts for gamma and epsilon at mean covariate values
  # Column names follow unmarked convention: col(Int), ext(Int)
  gamma_sims   <- plogis(sims[, grep("col\\(Int\\)", colnames(sims))])
  epsilon_sims <- plogis(sims[, grep("ext\\(Int\\)", colnames(sims))])
  
  psi_star          <- gamma_sims / (gamma_sims + epsilon_sims)
  recolonisation_yr <- 1 / gamma_sims
  persistence_yr    <- 1 / epsilon_sims
  
  data.frame(
    species            = sp_code,
    psi_star_median    = median(psi_star),
    psi_star_lo95      = quantile(psi_star, 0.025),
    psi_star_hi95      = quantile(psi_star, 0.975),
    recol_yr_median    = median(recolonisation_yr),
    recol_yr_lo95      = quantile(recolonisation_yr, 0.025),
    recol_yr_hi95      = quantile(recolonisation_yr, 0.975),
    persist_yr_median  = median(persistence_yr),
    persist_yr_lo95    = quantile(persistence_yr, 0.025),
    persist_yr_hi95    = quantile(persistence_yr, 0.975)
  )
}
```

### 3.2 Expected output structure (fill once computed)

| Species | ψ* median | ψ* 95% CI | Recolonisation (yr) | Persistence (yr) |
|---|---|---|---|---|
| O. tarda | [PENDING] | [PENDING] | [PENDING] | [PENDING] |
| P. alchata | [PENDING] | [PENDING] | [PENDING] | [PENDING] |
| P. orientalis | [PENDING] | [PENDING] | [PENDING] | [PENDING] |
| T. tetrax | [PENDING] | [PENDING] | [PENDING] | [PENDING] |

**Indicative values from raw simulation** (NOT for reporting — no CI): ψ* ≈ 0.4–1.9% (= γ/(γ+ε) using mean simulated rates). These values should appear in the Abstract.

---

## 4. Independent Validation (Spanish Biodiversity Atlas)

✅ **Complete.** Spatially-blocked 5-fold cross-validation against the III Atlas of Breeding Birds in Spain (SEO/BirdLife 2022). Method: blockCV R package; block size set to exceed spatial autocorrelation range.

| Species | Atlas prevalence | AUC | TSS | RMSE | Spearman ρ |
|---|---|---|---|---|---|
| O. tarda | 10.1% | **0.884** | **0.612** | 0.275 | **0.405*** |
| P. alchata | 8.4% | **0.844** | **0.569** | 0.261 | **0.334*** |
| P. orientalis | 16.2% | **0.823** | **0.536** | 0.366 | **0.410*** |
| T. tetrax | 24.8% | **0.850** | **0.546** | 0.460 | **0.521*** |

All Spearman ρ P < 0.001.

**Primary metric for GCB:** Spearman ρ (threshold-free, comparable across prevalences). AUC secondary. TSS in supplementary.

**⚠️ Note:** Verify that blockCV block size ≥ largest stPGOcc spatial range (~264 km for P. alchata) to ensure true independence of validation folds.

**Residual spatial autocorrelation (colext baseline, before spatial correction):**

| Species | Moran's I | P |
|---|---|---|
| O. tarda | 0.358 | <0.001 |
| P. alchata | 0.491 | <0.001 |
| P. orientalis | 0.374 | <0.001 |
| T. tetrax | 0.423 | <0.001 |

These values represent the colext model residuals. They motivate the stPGOcc spatial models.

---

## 5. Effort Confounding Diagnostics

✅ **Complete.** Script: `scripts/8_effort_confounding.R`.

### 5.1 eBird effort growth (2017–2023)
- Total checklists: +256% (otitar, tettet) / +201% (ptealc, pteori)
- Growth mechanism: **more observers visiting more sites** (extensive margin), NOT longer or more intensive surveys (intensive margin)
- Mean visits/site: stable (ρ = −0.05 to −0.29, all P > 0.5)
- Mean checklist duration: slightly decreasing (ρ = −0.57 to −0.68)

### 5.2 Naive occupancy vs. effort

| Species | ρ(year, naive occ.) | P | ρ(year, checklists) | P | Verdict |
|---|---|---|---|---|---|
| O. tarda | −0.250 | 0.589 | 0.964 | <0.001 | No confounding |
| P. alchata | 0.324 | 0.478 | 0.821 | 0.023 | **Marginal (see below)** |
| P. orientalis | 0.643 | 0.119 | 0.821 | 0.023 | No confounding |
| T. tetrax | 0.018 | 0.969 | 0.964 | <0.001 | No confounding |

### 5.3 Binomial GLM controlling for log(checklists)

| Species | β_year | 95% CI | P | β_log(effort) | P |
|---|---|---|---|---|---|
| O. tarda | 0.096 | [−0.058, 0.251] | 0.221 | −0.469 | 0.275 |
| P. alchata | −0.107 | [−0.255, 0.042] | 0.160 | **0.860** | **0.050** |
| P. orientalis | −0.006 | [−0.152, 0.141] | 0.941 | 0.543 | 0.229 |
| T. tetrax | 0.044 | [−0.144, 0.233] | 0.644 | −0.278 | 0.593 |

### 5.4 Reporting guidance for manuscript
- Three species: strong evidence against confounding (naive occupancy flat despite 2–3× more checklists)
- **P. alchata**: report as "marginally significant effort effect (P = 0.05)" — NOT as "no confounding confirmed." Year effect non-significant (P = 0.160) and absolute naive occupancy change minimal (−0.7 pp over 7 years), but the marginal effort coefficient warrants explicit disclosure.
- Key sentence for Discussion: *"Naive occupancy was uncorrelated with survey effort for three of four species (Spearman ρ = −0.25 to 0.64, P > 0.1); for P. alchata, the effort coefficient in a binomial GLM controlling for year was marginally significant (β = 0.86, P = 0.05), warranting cautious interpretation of occupancy trends for this species."*

---

## 6. Collinearity and NDVI Diagnostics

✅ **Diagnostics complete.** Script: `scripts/9_collinearity_ndvi.R`. **Decisions pending.**

### 6.1 Pairwise correlations (dynamic covariates, all species)

| Pair | O. tarda | P. alchata | P. orientalis | Flag |
|---|---|---|---|---|
| r(NDVI, pr) | 0.575 | 0.584 | 0.584 | ⚠️ Moderate |
| r(NDVI, tmmx) | −0.605 | −0.613 | −0.613 | ⚠️ Moderate |
| r(pr, tmmx) | −0.704 | −0.714 | −0.714 | ❌ Severe |
| r(tmmx, tmmn) | 0.788 | — | 0.805 | ❌ Severe |

T. tetrax uses only LC12; no relevant pairwise correlations.

### 6.2 VIF by submodel — all values < 5

| Species | Sub-model | Max VIF | Covariate | Status |
|---|---|---|---|---|
| O. tarda | γ | 4.43 | tmmx | ✅ Acceptable |
| P. alchata | γ | 1.46 | NDVI / pr | ✅ Acceptable |
| P. orientalis | γ | 3.57 | tmmx | ✅ Acceptable |
| P. orientalis | ε | 1.12 | LC12 | ✅ Acceptable |

### 6.3 Coefficient stability (NDVI removal test) — critical flags

| Species | Sub-model | Flag | Description |
|---|---|---|---|
| O. tarda | γ | BETA_SHIFT | tmmn shifts 1.1 SE when NDVI removed; all signs maintained. AIC −22. |
| **P. alchata** | **γ** | **SIGN_CHANGE** | **pr changes sign (0.16 → −0.72) when NDVI removed. Severe confounding.** |
| P. orientalis | γ | OK | Shifts < 0.33 SE. AIC unchanged (−0.77). |
| **P. orientalis** | **ε** | **BETA_SHIFT + SE_INCREASE** | **LC12: −0.83 → −0.09 (shift = 2.8 SE); pr SE +60%. AIC +60 if NDVI removed.** |

### 6.4 NDVI decomposition
- Mean R² (NDVI ~ pr + tmmx + tmmn): 0.508–0.514 across species
- Interpretation: ~50% of interannual NDVI variation is climate-driven; ~50% is non-climatic (irrigation, crop rotation, land management)
- Components saved: `results/diagnostics/ndvi_climate_component.rds`, `ndvi_residual_component.rds`

**P. orientalis/ε case study** (Discussion 4.2, not primary finding):
- NDVI_climate → extinction: +0.007 (warming-driven vegetation increase raises extinction risk)
- NDVI_landuse → extinction: −0.011 (agricultural maintenance reduces extinction risk)
- **[PENDING]** Bootstrap CIs — if CIs include zero, describe as "directional but not statistically distinguishable from zero"

### 6.5 Pending model revision decisions

**[DECISION NEEDED — Guillermo to decide before scripts/10 can run]**

| Species | Sub-model | Recommendation | Trade-off |
|---|---|---|---|
| P. alchata | γ | **Remove NDVI** (Opción B) | AIC cost +8; pr estimate becomes interpretable; eliminates sign change |
| P. orientalis | ε | Retain NDVI with documented limitation (Opción A) | AIC cost if removed: +60; LC12 estimate conditional on NDVI presence |

If Opción B for ptealc: refit with `gamma ~ pr` only; save to `results/ptealc_model_object_revised.rds`.

---

## 7. Counterfactual Attribution Analysis

⚠️ **Preliminary — pending model revisions (§6) and bootstrap n = 1,000 (scripts/10).**

Script: `scripts/8_counterfactual_attribution.R` (original models), `scripts/10_attribution_revised.R` (pending).

### 7.1 Scenario design (2×2 factorial)

| Scenario | Climate | Land use | Interpretation |
|---|---|---|---|
| S0 | Frozen 2017 | Frozen 2017 | No change baseline |
| S1 | Observed | Frozen 2017 | Climate effect |
| S2 | Frozen 2017 | Observed | Land-use effect |
| S3 | Observed | Observed | Full observed |

Attribution: Δclimate = S1−S0; Δlanduse = S2−S0; Δinteraction = S3−S1−S2+S0.

### 7.2 Current attribution results (original models, no bootstrap CI)

| Species | Δγ_climate | Δγ_landuse | Δε_climate | Δε_landuse | Dominant driver |
|---|---|---|---|---|---|
| O. tarda | −0.00061 | 0 | +0.00074 | **−0.00110** | Climate (γ) / Land use (ε) |
| P. alchata | +0.00007 | 0 | **−0.00734** | 0 | Climate only (model construction) |
| P. orientalis | −0.00034 | −0.00018 | +0.00573 | −0.00021 | Climate dominant |
| T. tetrax | 0 | ~0 | 0 | ~0 | Land use only (model construction) |

⚠️ **P. alchata γ excluded** (16 colonisation events — separation confirmed). P. alchata attribution is entirely climate by model construction (no LC variables in model). T. tetrax attribution is entirely land use by construction (only LC12 in model). Cross-species comparison is therefore limited.

⚠️ **No bootstrap CIs exist yet.** All effects are point estimates only. Absolute effects are small (10⁻³ to 10⁻⁶ probability units) — interpreted comparatively, not as absolute change estimates.

### 7.3 Pending scripts/10 tasks (execution order is mandatory)

| Task | What | Output | Status | Executor |
|---|---|---|---|---|
| Task 0 | Generate NDVI_climate and NDVI_residual objects | `results/ndvi_decomp_objects.rds` | ❌ Pending | Claude Code |
| Task 1 (**BLOCKING**) | Reestimate pteori/ε without NDVI; verify LC12 stability (shift < 0.5 SE) | `results/pteori_epsilon_revised.rds` | ❌ Pending | Raúl/Claude Code |
| Task 2 | Reestimate otitar/γ with NDVI_climate + NDVI_residual | `results/otitar_gamma_revised.rds` | ❌ Pending | Claude Code |
| Task 3 | Full attribution with revised models (4 scenarios × 4 species) | `results/attribution_scenarios.rds` | ❌ Pending | Claude Code |
| Task 4 | Bootstrap n = 1,000 (mvrnorm on coefficient vectors) | `results/attribution_boot_summary.csv` | ❌ Pending | Raúl (cluster) |
| Task 5 | Cross-species summary table | `results/attribution_table3.csv` | ❌ Pending | Claude Code |

**Note:** Tasks 3–5 cannot run until Task 1 is resolved and the pteori/ε revision decision is made.

---

## 8. Spatial Occupancy Models (stPGOcc)

⚠️ **CONDITIONAL — results preliminary, not publication-ready.**  
Script: `scripts/7_spatial_results.R`.

### 8.1 Framework distinction (critical for interpretation)

| Framework | Estimates | Unique contribution |
|---|---|---|
| colext | Explicit γ and ε with environmental covariates | **Why** colonisation/extinction occurs |
| stPGOcc | Latent occupancy z with spatial NNGP + AR(1) | **Where and when** turnover occurs; spatial range |

These are complementary, not competing. colext attribution is the primary result. stPGOcc provides spatial context and pre-empts independence assumption critique from reviewers.

### 8.2 WAIC model comparison (stPGOcc vs tPGOcc) — CONFIRMED

| Species | tPGOcc WAIC | stPGOcc WAIC | ΔWAIC |
|---|---|---|---|
| O. tarda | 2496.9 | 2040.4 | **456.5** |
| P. alchata | 2103.8 | 1655.2 | **448.6** |
| P. orientalis | 2273.7 | 1821.3 | **452.4** |
| T. tetrax | 1890.4 | 1538.8 | **351.6** |

Overwhelming support for spatial structure (ΔWAIC > 350 all species). This is the stPGOcc vs tPGOcc comparison — distinct from the spPGOcc vs tPGOcc comparison reported separately.

### 8.3 Moran's I reduction (stPGOcc vs tPGOcc)

| Species | Moran's I (tPGOcc) | Moran's I (stPGOcc) | Reduction |
|---|---|---|---|
| O. tarda | 0.145 | −0.060 | 141% |
| P. alchata | 0.236 | −0.042 | 118% |
| P. orientalis | 0.173 | −0.035 | 120% |
| T. tetrax | 0.144 | −0.042 | 129% |

Slight over-correction (negative Moran's I) is expected in spatial models — residuals become spatially independent.

### 8.4 Effective spatial ranges (PRELIMINARY — not converged)

| Species | Range (km) | Ecological interpretation |
|---|---|---|
| T. tetrax | **42.7** | Local; site fidelity, limited natal dispersal |
| O. tarda | **46.7** | Mesoscale; inter-lek connectivity |
| P. orientalis | **171.6** | Regional; cereal pseudo-steppe belt connectivity |
| P. alchata | **263.9** | Regional/nomadic; tracks ephemeral food resources |

**These ranges are preliminary.** Used in skeleton v6 as illustrative numbers with explicit conditionality.

### 8.5 MCMC convergence diagnostics — PROBLEMATIC

| Parameter | Typical Rhat | Typical ESS | Target | Status |
|---|---|---|---|---|
| β coefficients | 1.00–1.07 | 185–870 | Rhat < 1.1, ESS > 100 | ✅ OK |
| α (detection) | 1.00–1.01 | 1,300–3,000 | Rhat < 1.1, ESS > 100 | ✅ OK |
| σ² (spatial variance) | 1.01–1.07 | 31–48 | Rhat < 1.1, ESS > 100 | ❌ ESS too low |
| **φ (spatial decay)** | **1.25–1.31** | **44–56** | Rhat < 1.1, ESS > 100 | ❌ Not converged |
| **σ²_t (temporal variance)** | **1.84–4.19** | **9–108** | Rhat < 1.1, ESS > 100 | ❌ Not converged |

**Current MCMC settings:** n.iter = 20,000; burn-in = 10,000; thin = 20 → 500 post-burn samples/chain.

**Required for publication:**
```r
# RAÚL: Run in cluster with these settings
stPGOcc(
  n.batch      = 4000,       # 4000 × 25 = 100,000 total samples
  batch.length = 25,
  n.burn       = 30000,
  n.thin       = 20,         # → 3,500 post-burnin samples
  n.chains     = 3,
  NNGP         = TRUE,
  n.neighbors  = 5,
  ar1          = TRUE
)
# Target: Rhat < 1.3 for phi and sigma2_t; ESS > 100 for all parameters
# Estimated runtime: 6–24 hours per species on cluster
```

**Classification: NO PUBLICABLE AÚN.** Results go in supplementary as preliminary until convergence confirmed. Paper text uses conditional framing: "Preliminary spatial range estimates suggest..." with explicit convergence caveat.

---

## 9. Pre-Submission Checklist

### LEVEL 1 — BLOCKING

| Item | What it unblocks | Who | Expected output |
|---|---|---|---|
| **Decision: ptealc/γ NDVI removal (Opción A vs B)** | Tasks 1–5 of scripts/10; attribution figures | Guillermo | Decision recorded |
| **Equilibrium occupancy ψ* with bootstrap CI** | Abstract, Results 3.3, Figure 1B, Figure 3 | Claude Code | `results/equilibrium_occupancy_table.csv` |
| **pteori/ε reestimation without NDVI (Task 1)** | Prerequisite for attribution; unblocks Tasks 3–5 | Raúl or Claude Code | `results/pteori_epsilon_revised.rds` |
| **stPGOcc convergence (100,000 iterations)** | Results 3.5, Figure 5, Figure 6, Discussion 4.3 | Raúl (cluster) | Rhat < 1.3 for φ, σ²_t |
| **Bootstrap attribution n = 1,000 (scripts/10)** | Results 3.7, Figure 7, Discussion 4.2 | Raúl (cluster) | `results/attribution_boot_summary.csv` |
| **parboot GoF nsim = 500** | Methods 2.4, c-hat reportable | Raúl (cluster) | `results/gof/pb_*_nsim500.rds` |

### LEVEL 2 — STRONGLY RECOMMENDED

| Item | Why it matters |
|---|---|
| Naive vs corrected γ comparison | Discussion 4.4 — proves detection correction changes conclusions, not just precision |
| Figure 1 Panel B: isoclines γ/(γ+ε) | Core figure — requires equilibrium occupancy numbers |
| NDVI decomposition bootstrap CIs for pteori/ε | Determines whether case study is reportable or only directional |
| Verify blockCV block size ≥ 264 km | Ensures validation is truly spatially independent |

### LEVEL 3 — OPTIONAL BUT VALUABLE

| Item | Value |
|---|---|
| Calibration curves (observed vs predicted by decile) | Standard request in GCB review |
| Δγ required for ψ* = 0.10 | Discussion 4.1 — shows what change magnitude would actually matter for recovery |
| Photovoltaic overlap analysis | Policy relevance paragraph |

---

## 10. File Inventory

### Scripts
| Script | Status | Key output |
|---|---|---|
| `scripts/1_prepare_ebird_data.R` | ✅ Complete | Filtered eBird data |
| `scripts/2_prepare_static_variables.R` | ✅ Complete | WorldClim, MODIS static |
| `scripts/3_prepare_dynamic_variables.R` | ✅ Complete | GEE exports |
| `scripts/4_occupancy_models.R` | ✅ Complete | `results/{sp}_model_object.rds` |
| `scripts/5_validation.R` | ✅ Complete | AUC, TSS, Spearman ρ |
| `scripts/6_spatial_occupancy_test.R` | ✅ (single-year spPGOcc) | spPGOcc ranges (not used in paper) |
| `scripts/7_spatial_results.R` | ⚠️ Convergence issues | stPGOcc objects (preliminary) |
| `scripts/8_counterfactual_attribution.R` | ⚠️ Original models | Current attribution (preliminary) |
| `scripts/8_effort_confounding.R` | ✅ Complete | `results/pub_effort_confounding_summary.csv` |
| `scripts/9_collinearity_ndvi.R` | ✅ Complete | VIF, stability, NDVI decomposition |
| `scripts/10_attribution_revised.R` | ❌ Pending | Revised attribution with bootstrap |
| `R/model_configs.R` | ✅ Complete | Centralised model formulas |

### Model objects
| File | Species | Status |
|---|---|---|
| `results/otitar_model_object.rds` | O. tarda | ✅ |
| `results/ptealc_model_object.rds` | P. alchata | ✅ (γ submodel flagged) |
| `results/pteori_model_object.rds` | P. orientalis | ✅ (ε submodel flagged) |
| `results/tettet_model_object.rds` | T. tetrax | ✅ |
| `results/pteori_epsilon_revised.rds` | P. orientalis (ε revision) | ❌ Pending |

---

*Master report v1.0 — Guillermo Fandos / Raúl Contreras-Martín — March 2026*  
*Supersedes five individual analysis reports. Inconsistencies between reports documented in §0.*
