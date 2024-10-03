# Dynamic Occupancy Models for Steppe Birds with Citizen Science Data

This project aims .....TODO.....

## Prerequisites
- eBird account
- R
- Google Earth Engine (GEE) account

---

## Step-by-Step Guide

### 0. Download eBird data

.....TODO.....

### 1. Prepare the Database
**Script:** `1_prepare_database_ebird.R`

This script processes the raw eBird data. It filters the data and creates a new database for each steppe bird species.

### 2. Prepare the Static Variables
**Script:** `2_prepare_static_variables_XXX.R`

This script performs the following tasks for each species:
- Formats the data and transforms it from long to wide format.
- Appends bioclimatic variables to the data.
- Saves the new wide-format databaseS.
- Creates a separate file containing only the sites, which will be used to download dynamic variables later on.

### 3. Prepare the Dynamic Variables

1. **Upload to Google Earth Engine:**  
   Upload the script `3_download_dynamic_variables.js` to Google Earth Engine's (GEE) code editor (https://code.earthengine.google.com/).  
   - In the GEE editor, upload the file generated in the previous step (e.g., `XXX_occ_wide_latlong.csv`) by going to the "Assets" tab and clicking "New" -> "CSV File (.csv)" -> "SELECT" -> "Upload".
   This script merges the site data with dynamic variable values for each year.
   - Update the script manually by changing the file names at the beginning and end of the script.

2. **Add Dynamic Variables to the Database:**  
   **Script:** `3_prepare_dynamic_variables_XXX.R`  
    This script adds the dynamic variables to the previous databases.

### 4. Fit the Models
**Script:** `4_occupancy_XXX.R`

This script formats the variables for use in the `unmarkedMultFrame` object and designs the occupancy model for each species. The workflow includes:
- Analyzing the goodness of fit of the model.
- Evaluating response variables for occupancy, colonization, extinction, and detection.
- Generating initial occupancy maps and predicted mean values for occupancy, colonization, and extinction.

### 5. Create a Richness Map
**Script:** `5_richness.R`

This script converts each species' occupancy map into a binary map, where 1 indicates occupancy and 0 indicates absence. It, then, merges all species maps into a single richness map. The resulting map shows a value between 0 (no species present) and 6 (all six species present).

### 6. Assess Threats to Steppe Birds
**Script:** `6_threats.R`

This script calculates how much of the steppe birds' occupied space overlaps with photovoltaic plants. This helps identify potential threats from renewable energy installations.

### 7. Overlap with Natura 2000
**Script:** `7_overlap_richness_RN2000.R`

This script overlaps the steppe bird species richness map with the Red Natura 2000 map. It calculates the proportion of the steppe bird habitat that is protected under the Natura 2000 network.

---