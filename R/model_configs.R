###############################################################################
# R/model_configs.R
#
# Purpose: Define species-specific model formulas for colext() dynamic
#          occupancy models. Centralises all model specifications.
#
# Each species config contains:
#   - psi_formula:     initial occupancy covariates (static site-level)
#   - gamma_formula:   colonisation covariates (yearly site-level)
#   - epsilon_formula: extinction covariates (yearly site-level)
#   - p_formula:       detection covariates (observation-level)
#   - psi_vars:        variables used in occupancy map prediction
#
# Model formulas were selected via AIC comparison (see guille_feedback.md).
###############################################################################

get_model_config <- function(species_code) {

  configs <- list(

    # ---- Otis tarda ----
    otitar = list(
      psi_formula     = ~ bio1 + bio2 + tree_cover + grass_cover + topo_elev,
      gamma_formula   = ~ NDVI + pr + tmmn + tmmx,
      epsilon_formula = ~ Land_Cover_Type_1_Percent_Class_6 +
                           Land_Cover_Type_1_Percent_Class_13 + tmmx,
      p_formula       = ~ effort + observers + NDVI_obs + pr_obs + topo_aspect_obs,
      psi_vars        = c("bio1", "bio2", "tree_cover", "grass_cover", "topo_elev"),
      gamma_vars      = c("NDVI", "pr", "tmmn", "tmmx"),
      epsilon_vars    = c("Land_Cover_Type_1_Percent_Class_6",
                          "Land_Cover_Type_1_Percent_Class_13", "tmmx")
    ),

    # ---- Pterocles alchata ----
    # FIX: Removed Class_0 from gamma/epsilon (complete separation: 98% zero
    #      at occupied sites, coefficients diverged to ±930).
    #      Removed Class_13 from gamma (only 2/16 colonization events had >0).
    #      Removed tmmn from epsilon (r=0.81 with tmmx, only 18 extinction events).
    # FIX v4: Removed NDVI from gamma. NDVI is ~50% climate-driven; its
    #      inclusion caused pr coefficient sign change (+0.162 → -0.717),
    #      making climate attribution uninterpretable. AIC cost = +8.
    #      Without NDVI, pr is the sole climate covariate in gamma and its
    #      sign is ecologically consistent. Decision record: audit_gcb_v4.md.
    ptealc = list(
      psi_formula     = ~ bio1 + bio2 + tree_cover + grass_cover + topo_aspect,
      gamma_formula   = ~ pr,
      epsilon_formula = ~ pr + tmmx,
      p_formula       = ~ time + duration + effort + observers + NDVI_obs + pr_obs,
      psi_vars        = c("bio1", "bio2", "tree_cover", "grass_cover", "topo_aspect"),
      gamma_vars      = c("pr"),
      epsilon_vars    = c("pr", "tmmx")
    ),

    # ---- Pterocles orientalis ----
    # FIX: Removed Class_7 and Class_13 from epsilon (quasi-complete separation:
    #      Class_13=0 in 43/44 transitions, estimate=+109; Class_7 near-zero
    #      at most occupied sites, estimate=-5.4 with NaN SE).
    # FIX v4: Removed NDVI from epsilon. NDVI is ~50% climate-driven,
    #      creating attribution ambiguity. Without NDVI, LC12 coefficient
    #      should recover to a more negative value (currently -0.09),
    #      strengthening the land-use signal. Decision: audit_gcb_v4.md.
    pteori = list(
      psi_formula     = ~ bio2 + tree_cover + grass_cover,
      gamma_formula   = ~ Land_Cover_Type_1_Percent_Class_7 + NDVI + tmmn + tmmx,
      epsilon_formula = ~ Land_Cover_Type_1_Percent_Class_12 + pr,
      p_formula       = ~ time + duration + observers + pr_obs + topo_aspect_obs,
      psi_vars        = c("bio2", "tree_cover", "grass_cover"),
      gamma_vars      = c("Land_Cover_Type_1_Percent_Class_7", "NDVI", "tmmn", "tmmx"),
      epsilon_vars    = c("Land_Cover_Type_1_Percent_Class_12", "pr")
    ),

    # ---- Tetrax tetrax ----
    # FIX v3: Only 25 transitions from occupied sites (9 ext, 16 persist).
    #   Original formula had Class_0 + Class_12 + Class_13 + NDVI + tmmx
    #   (6 params). Class_0 and Class_13 had zero variance → -Inf log-lik.
    #   Progressive testing revealed `pr` causes complete separation in
    #   gamma (ANY model with pr in gamma/epsilon fails to converge).
    #   `tmmx` also fails in epsilon. Only Class_12 works in both submodels.
    #   Best model: gam=~LC12, eps=~LC12 (AIC=1780.1).
    #   Simplified p formula to reduce overall parameter count.
    tettet = list(
      psi_formula     = ~ bio2 + tree_cover + grass_cover + topo_elev,
      gamma_formula   = ~ Land_Cover_Type_1_Percent_Class_12,
      epsilon_formula = ~ Land_Cover_Type_1_Percent_Class_12,
      p_formula       = ~ effort + observers + time,
      psi_vars        = c("bio2", "tree_cover", "grass_cover", "topo_elev"),
      gamma_vars      = c("Land_Cover_Type_1_Percent_Class_12"),
      epsilon_vars    = c("Land_Cover_Type_1_Percent_Class_12")
    )
  )

  if (!species_code %in% names(configs)) {
    stop("Unknown species code: ", species_code,
         ". Available: ", paste(names(configs), collapse = ", "))
  }

  configs[[species_code]]
}
