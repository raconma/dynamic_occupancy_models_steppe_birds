# Figure Plan -- GCB Submission v10

*Last updated: March 2026*

---

## Main Figures

| Figure | Title | Narrative role | Data source | Status | Priority |
|---|---|---|---|---|---|
| **Fig 1** | Study area and sampling effort | Establishes geographic scope, sampling intensity, and species-specific detection coverage before any modelling results are introduced; the reader needs to see where the data come from and that coverage is sufficient | `figs/pub_map_sampling_sites.png` (`scripts/publication_maps.R`) | Complete | Essential |
| **Fig 2** | Predicted initial occupancy (psi_1) maps | Establishes that predicted distributions match known range geography (confirmed by atlas validation, Table 3) before demographic inference is introduced; anchors all subsequent results in plausible spatial predictions | `figs/pub_map_occupancy_4species.png` (`scripts/publication_maps.R`) | Complete | Essential |
| **Fig 3** | Naive vs detection-corrected transition rates | Demonstrates that detection correction changes the qualitative conservation conclusion (Section 4.2), not just precision; the methodological demonstration that justifies the entire DOM approach | `figs/pub_fig1_naive_vs_corrected.png` (`scripts/fig_naive_vs_corrected.R`) | Complete | Essential |
| **Fig 4** | Isocline plot with extinction debt (KEY FIGURE) | The conceptual centrepiece: visually demonstrates that all four species fall deep in the decline zone (epsilon >> gamma) and quantifies the extinction debt as the gap between current and equilibrium occupancy; gains maximum impact after Figure 3 establishes that corrected rates are the credible ones | `figs/pub_fig_isocline_equilibrium.png` (`scripts/fig_isocline_equilibrium.R`) | Complete | Essential |
| **Fig 5** | Spatial distribution of extinction risk (bivariate risk map) | Anchors the demographic asymmetry (Figure 4) in geography: shows WHERE the demographic trap is most severe per species, identifying priority areas for extinction-prevention interventions; essential for translating the aggregate result into spatially explicit conservation targets | `figs/pub_fig3_bivariate_risk_map.png` (`scripts/4_predict_maps.R`) | Complete | Essential |
| **Fig 6** | Spatial hotspots: psi_1 and epsilon | Reveals the spatial mismatch between where species currently persist (psi_1 hotspots) and where they face the highest extinction risk (epsilon hotspots); this mismatch is the spatial expression of the demographic trap and informs site prioritisation | `figs/pub_fig1_hotspots_4species.png` (`scripts/publication_maps.R`) | Needs revision | Essential |
| **Fig 7** | Forest plot of environmental drivers | Identifies the species-specific environmental drivers of each demographic process (psi_1, gamma, epsilon), establishing the empirical basis for the species-specific intervention recommendations in Discussion 4.3 and 4.7 | `figs/pub_fig2_forest_3submodels.png` (`scripts/fig_forest_3submodels.R`) | Complete | Essential |
| **Fig 8** | Response curves for key driver-process combinations | Illustrates the functional form of the most important driver-process relationships, particularly the T. tetrax LC12 result where cropland simultaneously increases colonisation and reduces extinction -- the single most actionable conservation finding | `figs/pub_fig5_response_curves.png` (`scripts/fig_response_curves.R`) | Complete | Important |

## Supplementary Figures

| Figure | Title | Narrative role | Data source | Status | Priority |
|---|---|---|---|---|---|
| **Fig S1** | Forest plot with coefficient heatmap | Provides the comprehensive coefficient overview (forest + heatmap) for readers wanting full detail beyond the streamlined Figure 7; supports reproducibility | `figs/pub_fig4_forest_heatmap_combined.png` (`scripts/fig_forest_heatmap.R`) | Complete | Supplementary |
| **Fig S2** | stPGOcc extinction map | Shows the spatial distribution of predicted extinction rates from the spatially-explicit model, providing visual confirmation that the stPGOcc framework reproduces the colext extinction pattern while accounting for spatial autocorrelation | `figs/pub_fig3b_extinction_only_map.png` (`scripts/6_spatial_occupancy_test.R`) | Blocked | Supplementary |
| **Fig S3** | stPGOcc spatial diagnostics | Documents model convergence (Rhat, ESS) and spatial range estimation (phi), supporting the robustness check described in Section 3.6 | `figs/pub_fig_spatial_moran.png` (`scripts/6_spatial_occupancy_test.R`) | Blocked | Supplementary |
| **Fig S4** | Detection probability diagnostics | Shows how detection probability varies with survey effort covariates, supporting the detection model specification in Section 2.4 and the detection rates reported in Section 3.1 | `figs/pub_fig4_detection_comparison.png` (`scripts/fig_detection_comparison.R`) | Complete | Supplementary |
| **Fig S5** | blockCV design map | Visualises the spatial fold assignment used for cross-validation, allowing readers to verify that training and validation folds are geographically independent at the 270 km scale | pending: `scripts/5_validation.R` | Pending | Supplementary |
| **Fig S6** | Occupancy prevalence trends | Shows the model-simulated temporal trajectory of occupancy for all species with bootstrap CI and equilibrium psi* reference lines; supports the temporal narrative in Section 3.5 | `figs/pub_fig7_prevalence_trends_A.png` (`scripts/fig_prevalence_trends.R`) | Complete | Important |
| **Fig S7** | Full response curves | Complete set of response curves for all species x covariate combinations; serves as reference if Figure 8 is moved to supplementary for space | `figs/pub_fig5_response_curves.png` (`scripts/fig_response_curves.R`) | Complete | Supplementary |
| **Fig S8** | Calibration curves | Observed vs predicted occupancy by decile for each species, showing calibration slope and overdispersion pattern; supports the validation results in Table 3 and Section 3.8 | `figs/{sp}_calibration_curve.png` (`scripts/14_calibration.R`) | Complete | Supplementary |
| **Fig S9** | Model evaluation: atlas validation + temporal trends | Left: predicted psi_1 maps with atlas presence polygons overlaid, demonstrating spatial agreement with independent data. Right: detection-corrected vs naive temporal trends, quantifying detection-induced bias magnitude per species. Supports Section 3.8 and the detection-naive divergence discussion in Section 4.5 | `figs/pub_fig_model_evaluation.png` (`scripts/16_model_evaluation_figure.R`) | Complete | Important |

---

## Critical path to submission

The following Essential figures are not yet Complete, ordered by dependency chain:

### 1. Figure 6 -- Spatial hotspots (Needs revision)

**Current state:** Figure exists (`figs/pub_fig1_hotspots_4species.png`) but in a 4-column x 3-row layout (species x [psi_1, gamma, epsilon]). The v10 specification requires a two-column layout (psi_1 + epsilon only), omitting the gamma column given separation issues in *P. alchata* and near-zero gamma across species.

**What is blocking it:** Script revision needed in `scripts/publication_maps.R` to produce a 4-row x 2-column panel (species as rows, psi_1 and epsilon as columns).

**Effort estimate:** Low (~1 hour). Modify existing plotting code; data is already generated.

**Fallback:** Use existing 3-column figure and note in the caption that gamma hotspots are presented for completeness but should be interpreted with caution due to separation issues.

### 2. Figure S2 -- stPGOcc extinction map (Blocked)

**Current state:** Figure exists from a preliminary stPGOcc run (`figs/pub_fig3b_extinction_only_map.png`), but final values depend on stPGOcc convergence (Rhat < 1.1, ESS > 100 for phi and sigma^2^_t).

**What is blocking it:** stPGOcc MCMC re-run on cluster (Raul). This also blocks Figure S3 (diagnostics), Section 3.6 text, and Discussion 4.4.

**Effort estimate:** Dependent on cluster availability; MCMC run time ~48-72 hours per species.

**Fallback:** If stPGOcc does not converge before submission: (a) use colext-predicted epsilon surface from `scripts/4_predict_maps.R` for Figure S2; (b) present Section 3.6 as preliminary with convergence caveats; (c) add stPGOcc results in revision. The primary demographic inference (epsilon/gamma asymmetry, attribution, extinction debt) rests entirely on colext models and is not affected by stPGOcc convergence.

### 3. Figure S3 -- stPGOcc spatial diagnostics (Blocked)

**Current state:** Partial output exists from preliminary runs. Depends on same stPGOcc convergence as Figure S2.

**What is blocking it:** Same as Figure S2 (stPGOcc MCMC convergence).

**Fallback:** Same as Figure S2 -- present colext Moran's I diagnostics and note that stPGOcc spatial range estimates will be added in revision.

### 4. Figure S5 -- blockCV design map (Pending)

**Current state:** Analysis is complete (`scripts/5_validation.R`), but the map visualising the spatial fold assignment has not been rendered as a standalone figure.

**What is blocking it:** Minor script addition needed to export the fold assignment map from the validation pipeline.

**Effort estimate:** Low (~30 minutes). The blockCV object already contains the fold geometry; needs a simple ggplot/sf map export.

**Fallback:** Describe the block design textually in Supplementary S5 without the map figure. The 270 km block size and 5-fold design are fully documented in Methods 2.9.

---

## Summary

| Category | Count | Status |
|---|---|---|
| Main figures (Essential + Important) | 8 | 7 Complete, 1 Needs revision |
| Supplementary figures | 9 | 5 Complete, 2 Blocked (stPGOcc), 1 Pending, 1 Needs revision |
| **Total blocking items** | **3** | Figure 6 revision (low effort), stPGOcc convergence (cluster), Figure S5 render (low effort) |

**Critical dependency:** stPGOcc convergence is the only external blocker. All other items can be resolved independently. The core demographic inference (Figures 1-5, 7) is complete and does not depend on stPGOcc.
