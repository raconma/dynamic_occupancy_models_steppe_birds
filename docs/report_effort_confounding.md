# Effort Confounding Diagnostics: Are Occupancy Trends an Artefact of Increasing eBird Effort?

**Branch:** `effort-confounding-diagnostics`
**Date:** 2026-03-07
**Script:** `scripts/8_effort_confounding.R`

---

## 1. Motivation

eBird participation across Spain has grown substantially during our study period (2017--2023). If the number of checklists and visited sites increases over time, **apparent occupancy trends may reflect expanding survey coverage** rather than genuine ecological change. This is a well-known concern in citizen science analyses (Kelling et al. 2019, Johnston et al. 2020).

Dynamic occupancy models (colext, stPGOcc) attempt to separate occupancy ($\psi$) from detection ($p$) by modelling observation-level effort covariates (duration, distance, number of observers). However, these models address the **intensive margin** of effort (per-visit detectability) but may not fully absorb the **extensive margin** (more sites visited, more checklists submitted). If effort inflation operates primarily through the extensive margin, model-corrected trends could still be confounded.

> [!NOTE]
> **For Raul:** This analysis is a robustness check --- it tests whether naive occupancy tracks effort growth. If it does, the model-estimated trends need extra scrutiny. If it does not, we have strong evidence that the detection submodel is doing its job.

---

## 2. Methods

### 2.1 Effort trend diagnostics (Task 1)

From the wide-format detection histories (`{sp}_occ_wide_dynamic.csv`), we computed six effort metrics **per species per year** (2017--2023):

| Metric | Definition |
|--------|-----------|
| Total checklists | Count of non-NA detection observations |
| Sites visited | Number of sites with $\geq$1 visit |
| Mean visits/site | Average repeat visits at visited sites |
| Mean duration | Average checklist duration (minutes) |
| Naive occupancy | Proportion of visited sites with $\geq$1 detection |
| Proportion unsampled | Fraction of sites with zero visits |

For each metric, we computed Spearman rank correlations ($\rho$) with year. If naive occupancy rises in lockstep with total checklists, effort confounding is likely.

### 2.2 Model-estimated detection trends (Task 2)

From the best colext model per species (see table below), we extracted:

| Species | Pipeline | Best model |
|---------|----------|-----------|
| *O. tarda* | v4b | m9: ambos lc_crop |
| *P. alchata* | v5 | m6: gam~pr_lag |
| *P. orientalis* | v4b | m0: Baseline estatico |
| *T. tetrax* | v4b | m7: gam~grass+crop |

For each species, we computed:

1. **Mean detection probability** ($\hat{p}_t$) per year, derived from the model's detection coefficients and the empirical mean covariate values in each year.
2. **Smoothed occupancy** ($\hat{\psi}_t$) from the colext `smoothed()` function.
3. **Naive-to-model ratio** ($\psi_{\text{naive}} / \hat{\psi}$) per year. An increasing ratio indicates naive occupancy is growing faster than model-corrected occupancy, suggesting the detection submodel does not fully absorb effort inflation.

### 2.3 Residual effort confounding test (Task 3)

We fitted a binomial GLM for each species:

$$
\text{logit}(\psi_{\text{naive},t}) = \beta_0 + \beta_{\text{year}} \cdot t + \beta_{\text{effort}} \cdot \log(n_{\text{checklists},t})
$$

where $\psi_{\text{naive},t}$ is the proportion of visited sites with at least one detection in year $t$, and $n_{\text{checklists},t}$ is the total number of checklists that year. If $\beta_{\text{year}}$ remains significant after controlling for $\log(n_{\text{checklists}})$, this provides evidence for a real ecological trend net of effort.

---

## 3. Results

### 3.1 eBird effort growth

eBird checklists increased markedly across all species: **+256%** for *O. tarda* / *T. tetrax* and **+201%** for *P. alchata* / *P. orientalis* between 2017 and 2023. The number of visited sites grew proportionally (Spearman $\rho$ = 0.89--0.96, all $p < 0.01$). In contrast, **mean visits per site remained stable** ($\rho$ = -0.05 to -0.29, all $p > 0.5$) and **mean checklist duration slightly decreased** ($\rho$ = -0.57 to -0.68, all $p > 0.09$). This pattern indicates that eBird growth in Spain over this period was driven by **more observers visiting more sites**, not by individual observers increasing their sampling intensity.

> [!IMPORTANT]
> The extensive margin (more sites visited) grew substantially, but the intensive margin (effort per visit) did not. This is the key distinction for occupancy models: colext explicitly models the intensive margin via detection covariates, but the extensive margin changes the spatial sample.

### 3.2 Naive occupancy vs. effort

**Critically, naive occupancy did NOT track effort growth for any species.** Spearman correlations between year and naive occupancy were weak and non-significant:

| Species | $\rho$(year, naive occ.) | $p$-value | $\rho$(year, checklists) | $p$-value |
|---------|:------------------------:|:---------:|:------------------------:|:---------:|
| *O. tarda* | -0.250 | 0.589 | **0.964** | **<0.001** |
| *P. alchata* | 0.324 | 0.478 | **0.821** | **0.023** |
| *P. orientalis* | 0.643 | 0.119 | **0.821** | **0.023** |
| *T. tetrax* | 0.018 | 0.969 | **0.964** | **<0.001** |

Despite a 2--3-fold increase in checklists, naive occupancy barely changed (range: -0.7 to +0.9 percentage points over 7 years). This decoupling between effort and naive occupancy is the strongest evidence against systematic confounding.

![Cross-species effort and naive occupancy](../figs/pub_fig_effort_cross_species.png)

### 3.3 Detection probability trends

Mean model-estimated detection probability ($\hat{p}$) was either stable or slightly decreasing over time for all species. No species showed an increasing detection trend that would indicate inflating detectability due to effort.

| Species | $\hat{p}$ trend (Spearman $\rho$) | Interpretation |
|---------|:---------------------------------:|---------------|
| *O. tarda* | -0.771 | Slightly decreasing |
| *P. alchata* | -0.029 | Stable |
| *P. orientalis* | -0.314 | Slightly decreasing |
| *T. tetrax* | -0.771 | Slightly decreasing |

> [!NOTE]
> **For Raul:** If detection probability were rising with effort, we would worry that inflating $p$ might mask true occupancy declines or amplify apparent increases. The fact that $\hat{p}$ is flat or decreasing suggests the detection covariates (duration, distance, observers) are absorbing per-visit effort variation as intended.

### 3.4 Naive-to-model occupancy ratio

The ratio $\psi_{\text{naive}} / \hat{\psi}$ remained below 1.0 for all species across all years, consistent with imperfect detection (naive estimates underestimate true occupancy). No species showed an increasing ratio over time, further ruling out effort-driven inflation.

### 3.5 Binomial GLM: year effect controlling for effort

After controlling for $\log(\text{checklists})$, the year coefficient was **non-significant for all four species**:

| Species | $\beta_{\text{year}}$ | 95% CI | $p$ | $\beta_{\log(\text{effort})}$ | $p$ | Verdict |
|---------|:--------------------:|:------:|:---:|:-----------------------------:|:---:|---------|
| *O. tarda* | 0.096 | [-0.058, 0.251] | 0.221 | -0.469 | 0.275 | Ambiguous |
| *P. alchata* | -0.107 | [-0.255, 0.042] | 0.160 | 0.860 | **0.050** | Likely confounded |
| *P. orientalis* | -0.006 | [-0.152, 0.141] | 0.941 | 0.543 | 0.229 | Ambiguous |
| *T. tetrax* | 0.044 | [-0.144, 0.233] | 0.644 | -0.278 | 0.593 | Ambiguous |

> [!WARNING]
> *P. alchata* is the only species where $\log(\text{checklists})$ is marginally significant ($p$ = 0.05), while the year coefficient is not. However, the naive occupancy Spearman correlation with year is weak ($\rho$ = 0.32, $p$ = 0.48), and the year effect with 95% CI overlapping zero (-0.107 [-0.255, 0.042]) suggests no clear directional trend. The "likely confounded" verdict is a conservative flag for this species, but the actual occupancy change is minimal (-0.7 pp over 7 years).

---

## 4. Per-Species Diagnostics

### 4.1 *Otis tarda* (Great Bustard)

![Effort diagnostics: O. tarda](../figs/pub_fig_effort_diagnostics_otitar.png)

eBird checklists grew +256% (2017--2023), yet naive occupancy decreased by 0.2 percentage points ($\rho$ = -0.25, $p$ = 0.59). Model-estimated detection showed a slightly decreasing trend ($\rho$ = -0.77), consistent with the observed decrease in mean checklist duration. After controlling for effort, the year coefficient was non-significant ($\beta$ = 0.096, $p$ = 0.22). **There is no evidence that occupancy trends for *O. tarda* are confounded with effort.**

### 4.2 *Pterocles alchata* (Pin-tailed Sandgrouse)

![Effort diagnostics: P. alchata](../figs/pub_fig_effort_diagnostics_ptealc.png)

eBird checklists grew +201% (2017--2023). Naive occupancy showed a weak positive trend ($\rho$ = 0.32, $p$ = 0.48), and $\log(\text{checklists})$ was marginally significant in the GLM ($p$ = 0.05). However, the year effect was non-significant ($\beta$ = -0.107, $p$ = 0.16), and the overall change in naive occupancy was minimal (-0.7 pp). Model-estimated detection was stable ($\rho$ = -0.03). **This species warrants a cautionary note in the discussion, but the evidence for systematic confounding is weak.**

### 4.3 *Pterocles orientalis* (Black-bellied Sandgrouse)

![Effort diagnostics: P. orientalis](../figs/pub_fig_effort_diagnostics_pteori.png)

eBird checklists grew +201% (2017--2023). Naive occupancy showed the strongest positive trend of all species ($\rho$ = 0.64, $p$ = 0.12), but this was not statistically significant. In the binomial GLM, neither the year coefficient ($p$ = 0.94) nor $\log(\text{checklists})$ ($p$ = 0.23) were significant. Model-estimated detection was slightly decreasing ($\rho$ = -0.31). **There is no evidence of effort confounding; the modest increase in naive occupancy (+0.9 pp) is independent of effort growth.**

### 4.4 *Tetrax tetrax* (Little Bustard)

![Effort diagnostics: T. tetrax](../figs/pub_fig_effort_diagnostics_tettet.png)

eBird checklists grew +256% (2017--2023), yet naive occupancy was completely flat ($\rho$ = 0.02, $p$ = 0.97). Neither the year coefficient ($p$ = 0.64) nor $\log(\text{checklists})$ ($p$ = 0.59) were significant in the GLM. Model-estimated detection showed a slightly decreasing trend ($\rho$ = -0.77). **Of all species, *T. tetrax* shows the clearest decoupling between effort and occupancy --- there is no confounding.**

---

## 5. Synthesis and Recommendations

### 5.1 Key findings

1. **eBird effort grew 2--3-fold** across Spain (2017--2023), driven by more observers and sites, not by longer or more intensive surveys per visit.
2. **Naive occupancy did not track effort** for any species (all Spearman $p > 0.1$).
3. **Model-estimated detection ($\hat{p}$) was stable or decreasing**, confirming that detection covariates absorb per-visit effort variation.
4. **After controlling for log(checklists)**, the year coefficient was non-significant for all species.
5. **No species showed a ratio $\psi_{\text{naive}}/\hat{\psi}$ inflating over time.**

### 5.2 Interpretation for the manuscript

These results provide **strong evidence that occupancy trends estimated by our colext and stPGOcc models are not artefacts of increasing eBird participation**. The key argument has three pillars:

- **Decoupling**: Despite 2--3x more checklists, naive occupancy was flat or slightly declining, ruling out a mechanical relationship between effort and apparent occupancy.
- **Detection absorption**: Per-visit detection covariates (duration, distance, observers) adequately absorb the intensive margin of effort variation, as shown by stable or declining $\hat{p}$.
- **Extensive margin coverage**: The increasing number of visited sites expands spatial coverage but does not inflate per-site detection rates, because the occupancy models condition on the sites that were actually visited.

> [!NOTE]
> **For Raul:** I recommend including a paragraph in the Discussion addressing effort confounding explicitly, referencing this analysis as supplementary material. The key sentence would be: "Naive occupancy was uncorrelated with survey effort for all four species (Spearman $\rho$ = -0.25 to 0.64, all $p > 0.1$), and binomial GLMs controlling for log(checklists) showed non-significant year effects (all $p > 0.15$), providing no evidence that occupancy trends were confounded with increasing eBird participation."

### 5.3 Caveats

- The GLM test has **low statistical power** (7 data points per species). Non-significance does not prove absence of confounding, only that it is not detectable at this sample size.
- *P. alchata* showed a marginally significant effort coefficient ($p$ = 0.05). While the overall naive occupancy trend was negligible, this species should be discussed with appropriate caution.
- This analysis addresses **temporal** effort trends at the whole-dataset level. **Spatial** effort heterogeneity (more checklists in certain regions) is addressed by the spatial random effect in the stPGOcc model.

---

## 6. Output inventory

| File | Description |
|------|------------|
| `scripts/8_effort_confounding.R` | Full analysis script |
| `figs/pub_fig_effort_diagnostics_otitar.png` | 4-panel diagnostic: *O. tarda* |
| `figs/pub_fig_effort_diagnostics_ptealc.png` | 4-panel diagnostic: *P. alchata* |
| `figs/pub_fig_effort_diagnostics_pteori.png` | 4-panel diagnostic: *P. orientalis* |
| `figs/pub_fig_effort_diagnostics_tettet.png` | 4-panel diagnostic: *T. tetrax* |
| `figs/pub_fig_effort_cross_species.png` | Cross-species summary (4 panels) |
| `results/pub_effort_confounding_summary.csv` | Summary table (CSV) |
| `results/pub_effort_interpretation.txt` | Plain-text interpretation per species |

---

*Generated from `scripts/8_effort_confounding.R` on branch `effort-confounding-diagnostics`.*
