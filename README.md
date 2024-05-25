# Dynamic occupancy models for steppe birds

Just a quick first draft of the README.

## How to

Run preprare_sp_ebird_steppe.R to filter the raw data and get the information of our species, it will create two files. We are creating two files instead of just one so Google Earth Engine (GEE) doesn't time out while processing them.

Now we will append the variables to our data. In this case minimum and maximum temperature and precipitation, NDVI, EVI and land cover.
Go to [GEE code editor](https://code.earthengine.google.com/) and create four files:
- 1_GEE_script_landcover_part1
- 1_GEE_script_landcover_part2
- 2_GEE_script_climate_ndvi_evi_part1
- 2_GEE_script_climate_ndvi_evi_part2

Paste the code of the scripts with the same names.
Then, go to "Assets" -> "NEW" -> CSV file (.csv) and select both files (part1 and part2) that we created with our R code.

Run all four scripts in order. Download from Google Drive the two files that we have just created.

Then, run bind_csv.R to merge both files. Our data is prepared to run the dynamic models.
