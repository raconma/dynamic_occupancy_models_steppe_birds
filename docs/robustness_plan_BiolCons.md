# Robustness plan for the Biological Conservation submission

Fixed on 2026-10-04, before any of the analyses below was run. Each analysis
states what will be reported and how the text will change for each possible
outcome, so that the results cannot steer the specification. All analyses use
the primary formulas of script 41; no covariate selection is repeated.

The plan answers two reviews of the 8,000-word version (an internal review and
an external one), which agree on the same points: the static benchmark in the
temporal validation is too weak, the anchored-transition summary may be biased
by selection on observed records, the marginal-habitat result may reflect
seasonal use, the SPA contrasts lack power, and Fig. 3 reads as confirmatory.

## A. Competitive static benchmark and rolling-origin validation (script 51)

- Test years 2021, 2022 and 2023, each predicted from all earlier years
  (training 2017-2020, 2017-2021 and 2017-2022).
- Models: dynamic (posterior state in the last training year propagated with
  the fitted colonisation and extinction); S1 fixed-state over all training
  years (current reference); S2 single-season occupancy model fitted to the
  last training year only (a map updated every year); S3 fixed-state over the
  last three training years (a recent atlas); S4 habitat-only stacked model.
- Formulas: primary formulas, selected on 2017-2023. Repeating the selection
  within each training window is not done; this is stated as a limitation.
- Scores: held-out log-likelihood of the test-year detection histories, with
  95 % intervals from 1,000 resamples of 50 x 50 km blocks, for all surveyed
  cells and for cells with earlier records; AUC and Brier score.
- Reporting rule: all models, years and species are reported in the
  supplement whatever the outcome.
- Text rule: if the dynamic model beats S2 and S3 for the great bustard in the
  pooled test years (interval of the summed difference above zero), the claim
  "a static map, even one updated every year, overstates persistence" stands.
  If S2 or S3 matches the dynamic model, the claim is restricted to maps that
  are not updated, and the title changes accordingly.

## B. Simulation of the anchored-transition summary (new script 52)

- Simulate great bustard detection histories from the fitted primary model
  under two truths: equilibrium (colonisation and extinction rescaled so that
  expected extinctions equal colonisations) and the fitted decline.
- Keep the observed survey pattern: which cell-years were surveyed, the number
  of checklists, and survey continuation that depends on the previous year's
  record at the observed rates (0.37-0.50 after a record, 0.66-0.72 without).
- Refit the primary model to each replicate (100 per truth) and compute the
  anchored ratio of extinctions to colonisations as in script 26.
- Text rule: if the anchored ratio under equilibrium is centred on one (median
  within 0.8-1.25), one sentence in Methods states that the summary does not
  create losses. If it is inflated, the observed ratio is reported against the
  simulated null distribution and the claim of loss is restated accordingly.

## C. Breeding-season window and marginal habitat (script 39 filter, script 26)

- Rebuild great bustard histories using only checklists from 1 April to
  15 May (the lek period, before males leave the leks), refit the primary
  formulas, and recompute the anchored ratio and the core-marginal extinction
  contrast at the 30, 50 and 70 % thresholds.
- Report the number of recorded cells and of cells with anchored transitions
  in each class.
- Text rule: if marginal extinction remains higher (ratio interval above one
  at the 50 % threshold), seasonal use is no longer offered as an explanation.
  If the contrast disappears, the result is reported as seasonal use of
  marginal cells, not as higher extinction.

## D. Power of the SPA contrasts (no new models except one check)

- Minimum detectable effect on extinction (alpha 0.05, power 0.8) from the
  standard errors of the primary SPA model (script 49, psi_spa_t), and the
  range of odds ratios compatible with the data.
- Refit the primary SPA model without SPA cover in detection for the great
  bustard, to show whether the colonisation contrast depends on that term.
- Text rule: the SPA result is described as inconclusive for retention, with
  the detectable effect stated; no sentence implies low effectiveness.

## E. Cells with stable survey effort (post-processing of script 26)

- Restrict the anchored summary to recorded cells whose number of checklists
  varied by less than a factor of two between consecutive surveyed years.
- Text rule: one sentence in Results with the ratio and interval, table in
  the supplement.

## G. Interval check (Tables S6 and S13)

- Verify which intervals include resampling of 50 x 50 km blocks. Correct the
  tables or their captions so that every interval says what it carries.

## Not done

- F (common covariate set for initial occupancy, colonisation and extinction):
  Fig. 3 is labelled as exploratory, model-selected associations instead.
- Covariate reselection within the rolling-origin training windows.
- Year-to-year weather moves to the supplement; no new analysis.
