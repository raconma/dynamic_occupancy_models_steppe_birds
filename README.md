# Dynamic occupancy models for steppe birds with citizen science

This projects aims to.........................

## How to

### 1. Prepare the database
1_prepare_database_ebird.R:

This is the first script you will have to run. It takes the raw data downloaded from eBird, filters it and creates a new database for each species.

### 2. Prepare the static variables
2_prepare_static_variables_XXX.R:

For each one of the species: first it formats the data and transforms it from long to wide format, then, it appends the bioclimatic variables. The new database in wide format woth the bioclimatic variables is stored. Lastly, a new database is created just containing the sites, this file will be used to dowload the dynamic variables later on.

### 3. Prepare the dynamic variables
Dynamic variables need an special treatment, it is required for them to have the database in wide format.
First of all, we upload "3_dowload_dynamic_variables.js" in Google Earth Engines's (GEE) API https://code.earthengine.google.com/. The easiest way is to copy-paste the code in the code editor. Then, in the Assets' tab, we upload the output of the previous scripts, the ones named "XXX_occ_wide_latlong.csv". We can do this just by clicking in New - CSV file (.csv) - SELECT - Upload.
We update the script by changing the name of the file at the beginning and at the end of the script.
This script uses GEE to create a file with the sites and the asociated value for each dynamic variable each year.

3_prepare_dynamic_variables_XXX.R:
Now, we are ready to add the dynamic variables to our databases and that's what these scripts do.

### 4. Models
4_occupancy_XXX.R:

For each speciesirst, we format the variables to use them in the unmarkedMultFrame. When the unmarkedMultFrame is created we can start designing the model. After obtaining the best model (with the AIC), we can analize the Goodness of Fit or how well the results fit the observations.

Then we analize the response variables for the occupancy, colonization, extinction and detection, the initial occupancy maps, and the predicted mean occupancy, colonization and extinction.

### 5. Richness

5_richness.R:

This script transforms the initial occupancy map of each species in a binary map (raster values just 1 or 0, 1 meaning occupied and 0 absence) and merges every map in a richness map, with each pixel having a value from 0 to 6, 0 meaning complete absence and 6 being the presence of the 6 steppe birds.

### 6. Threats

6_threats.R:

With this script we are able to calculate how much steppe birds' space has also presence of photovoltaic plants.

### 7. Natura 2000

7_overlap_richness_RN2000.R:

Finnally, with this script we overlap the richness' map with the Red Natura 2000's map. We can analize then how much of the space occupied by the 6 steppe birds is protected by this network. 
