###############################################################################
# fig_index.R
#
# Registry mapping each publication figure to its data source, generating
# script, and GCB manuscript figure number.
###############################################################################

fig_index <- tibble::tribble(
  ~fig_number, ~filename,                              ~data_source,                          ~script,
  "Fig1",      "pub_map_main_figure.png",              "eBird checklists",                    "scripts/publication_maps.R",
  "Fig2",      "pub_fig_forest_gamma_epsilon.png",     "colext model coefficients",           "scripts/fig_forest_coefs.R",
  "Fig3",      "pub_fig_spatial_moran.png",            "stPGOcc residuals",                   "scripts/7_spatial_results.R",
  "Fig4",      "pub_fig4_detection_comparison.png",    "colext detection submodel",           "scripts/publication_figures_models.R",
  "Fig5",      "pub_fig_isocline_equilibrium.png",     "isocline_plot_data.csv, extinction_debt_table.csv", "scripts/fig_isocline_equilibrium.R",
  "Fig6",      "pub_fig_maps_3process_4species.png",   "OccuMap rasters, colext predictions", "scripts/fig_maps_3process.R",
  "Fig7",      "pub_fig_prevalence_trends.png",        "simulation_prevalence.csv, bootstrap_draws_5000.rds", "scripts/fig_prevalence_trends.R",
  "Fig8",      "pub_fig_spatial_range.png",            "stPGOcc phi estimates",               "scripts/7_spatial_results.R",
  "Fig9",      "pub_fig_response_curves.png",          "colext model objects",                "scripts/fig_response_curves.R",
  "FigS1",     "otitar_calibration_curve.png",         "calibration_slopes.csv",              "scripts/14_calibration.R"
)

print(fig_index)
