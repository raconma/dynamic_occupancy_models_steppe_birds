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
    ptealc = list(
      psi_formula     = ~ bio1 + bio2 + tree_cover + grass_cover + topo_aspect,
      gamma_formula   = ~ Land_Cover_Type_1_Percent_Class_0 +
                           Land_Cover_Type_1_Percent_Class_13 + NDVI + pr,
      epsilon_formula = ~ Land_Cover_Type_1_Percent_Class_0 + pr + tmmn + tmmx,
      p_formula       = ~ time + duration + effort + observers + NDVI_obs + pr_obs,
      psi_vars        = c("bio1", "bio2", "tree_cover", "grass_cover", "topo_aspect"),
      gamma_vars      = c("Land_Cover_Type_1_Percent_Class_0",
                          "Land_Cover_Type_1_Percent_Class_13", "NDVI", "pr"),
      epsilon_vars    = c("Land_Cover_Type_1_Percent_Class_0", "pr", "tmmn", "tmmx")
    ),

    # ---- Pterocles orientalis ----
    pteori = list(
      psi_formula     = ~ bio2 + tree_cover + grass_cover,
      gamma_formula   = ~ Land_Cover_Type_1_Percent_Class_7 + NDVI + tmmn + tmmx,
      epsilon_formula = ~ Land_Cover_Type_1_Percent_Class_7 +
                           Land_Cover_Type_1_Percent_Class_12 +
                           Land_Cover_Type_1_Percent_Class_13 + NDVI + pr,
      p_formula       = ~ time + duration + observers + pr_obs + topo_aspect_obs,
      psi_vars        = c("bio2", "tree_cover", "grass_cover"),
      gamma_vars      = c("Land_Cover_Type_1_Percent_Class_7", "NDVI", "tmmn", "tmmx"),
      epsilon_vars    = c("Land_Cover_Type_1_Percent_Class_7",
                          "Land_Cover_Type_1_Percent_Class_12",
                          "Land_Cover_Type_1_Percent_Class_13", "NDVI", "pr")
    ),

    # ---- Tetrax tetrax ----
    tettet = list(
      psi_formula     = ~ bio2 + tree_cover + grass_cover + topo_elev,
      gamma_formula   = ~ Land_Cover_Type_1_Percent_Class_12 + pr,
      epsilon_formula = ~ Land_Cover_Type_1_Percent_Class_0 +
                           Land_Cover_Type_1_Percent_Class_12 +
                           Land_Cover_Type_1_Percent_Class_13 + NDVI + tmmx,
      p_formula       = ~ effort + observers + time + NDVI_obs + topo_elev_obs,
      psi_vars        = c("bio2", "tree_cover", "grass_cover", "topo_elev"),
      gamma_vars      = c("Land_Cover_Type_1_Percent_Class_12", "pr"),
      epsilon_vars    = c("Land_Cover_Type_1_Percent_Class_0",
                          "Land_Cover_Type_1_Percent_Class_12",
                          "Land_Cover_Type_1_Percent_Class_13", "NDVI", "tmmx")
    )
  )

  if (!species_code %in% names(configs)) {
    stop("Unknown species code: ", species_code,
         ". Available: ", paste(names(configs), collapse = ", "))
  }

  configs[[species_code]]
}
