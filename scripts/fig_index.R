###############################################################################
# fig_index.R
#
# Registry mapping each publication figure to its data source, generating
# script, and GCB manuscript figure number.
###############################################################################

fig_index <- tibble::tribble(
  ~fig_number, ~filename,                              ~data_source,                          ~script,
  # Main figures
  "Fig1",      "pub_map_sampling_sites.png",           "eBird checklists",                    "scripts/publication_maps.R",
  "Fig2",      "pub_map_occupancy_4species.png",       "OccuMap rasters",                     "scripts/publication_maps.R",
  "Fig3",      "pub_fig3b_extinction_only_map.png",    "stPGOcc ext rasters",                 "scripts/fig_bivariate_risk_map.R",
  "Fig4",      "pub_fig1_naive_vs_corrected.png",      "naive_vs_corrected_full.csv",         "scripts/fig_naive_vs_corrected.R",
  "Fig5",      "pub_fig_isocline_equilibrium.png",     "isocline_plot_data.csv, extinction_debt_table.csv", "scripts/fig_isocline_equilibrium.R",
  "Fig6",      "pub_fig2_forest_3submodels.png",       "colext model objects",                "scripts/fig_forest_heatmap.R",
  "Fig7",      "pub_fig5_response_curves.png",         "colext model objects",                "scripts/fig_response_curves_v2.R",
  "Fig8",      "pub_fig7_prevalence_trends_A.png",     "simulation_prevalence.csv, bootstrap_draws_5000.rds", "scripts/fig_prevalence_trends_v2.R",
  # Supplementary figures
  "FigS1",     "pub_fig3_bivariate_risk_map.png",      "stPGOcc col/ext rasters",             "scripts/fig_bivariate_risk_map.R",
  "FigS2",     "pub_fig4_forest_heatmap_combined.png", "colext model objects",                "scripts/fig_forest_heatmap.R",
  "FigS3",     "pub_fig1_hotspots_4species.png",       "OccuMap + stPGOcc rasters",           "scripts/fig_hotspots.R",
  "FigS4",     "pub_fig_spatial_moran.png",            "stPGOcc residuals",                   "scripts/7_spatial_results.R",
  "FigS5",     "pub_fig4_detection_comparison.png",    "colext detection submodel",           "scripts/publication_figures_models.R",
  "FigS6",     "otitar_calibration_curve.png",         "calibration_slopes.csv",              "scripts/14_calibration.R"
)

print(fig_index)
