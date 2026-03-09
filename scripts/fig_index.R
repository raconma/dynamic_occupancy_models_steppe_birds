###############################################################################
# fig_index.R
#
# Registry mapping each publication figure to its data source, generating
# script, and GCB manuscript figure number.
###############################################################################

fig_index <- tibble::tribble(
  ~fig_number, ~filename,                              ~data_source,                          ~script,
  "Fig1",      "pub_fig1_naive_vs_corrected.png",      "naive_vs_corrected_full.csv",         "scripts/fig_naive_vs_corrected.R",
  "Fig2",      "pub_fig2_forest_3submodels.png",       "colext model objects",                "scripts/fig_forest_heatmap.R",
  "Fig3",      "pub_fig3_bivariate_risk_map.png",      "stPGOcc col/ext hotspot rasters",     "scripts/fig_bivariate_risk_map.R",
  "Fig4",      "pub_fig4_forest_heatmap_combined.png", "colext model objects",                "scripts/fig_forest_heatmap.R",
  "Fig5",      "pub_fig_isocline_equilibrium.png",     "isocline_plot_data.csv, extinction_debt_table.csv", "scripts/fig_isocline_equilibrium.R",
  "Fig6",      "pub_fig1_hotspots_4species.png",       "OccuMap + stPGOcc rasters",           "scripts/fig_hotspots.R",
  "Fig7",      "pub_fig_maps_3process_4species.png",   "OccuMap rasters, colext predictions", "scripts/fig_maps_3process.R",
  "Fig8",      "pub_fig7_prevalence_trends_A.png",     "simulation_prevalence.csv, bootstrap_draws_5000.rds", "scripts/fig_prevalence_trends_v2.R",
  "Fig9",      "pub_fig5_response_curves.png",         "colext model objects",                "scripts/fig_response_curves_v2.R",
  "Fig10",     "pub_fig_spatial_moran.png",            "stPGOcc residuals",                   "scripts/7_spatial_results.R",
  "Fig11",     "otitar_calibration_curve.png",         "calibration_slopes.csv",              "scripts/14_calibration.R"
)

print(fig_index)
