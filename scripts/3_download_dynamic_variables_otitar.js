// Import the first part of the data
var coords = ee.FeatureCollection('projects/ee-rcontr03/assets/otitar_occ_wide_latlong');

// Visualize the data on the map
Map.addLayer(coords, {color: 'green'}, 'Presence - Absence');
Map.centerObject(coords, 11);

// Define the time range for filtering the MODIS image collection (NDVI/EVI, Land Cover, and Temperature)
var time_start_ndviEvi = '2017-04-01';
var time_end_ndviEvi = '2022-06-30';

var time_start_landCover = '2017-01-01';
var time_end_landCover = '2022-12-31';

var time_start_temperature = '2017-04-01';
var time_end_temperature = '2022-06-30';

// Import the MODIS NDVI and EVI image collection
var ndviEviCollection = ee.ImageCollection("MODIS/061/MOD13A3").filterDate(time_start_ndviEvi, time_end_ndviEvi);
var ndviEviSubset = ndviEviCollection.select(['NDVI', 'EVI']);

// Import the MODIS land cover image collection and filter by date
var landCoverCollection = ee.ImageCollection("MODIS/061/MCD12C1").filterDate(time_start_landCover, time_end_landCover);
var landCoverSubset = landCoverCollection.select([
  'Majority_Land_Cover_Type_1',
  'Land_Cover_Type_1_Percent_Class_0',
  'Land_Cover_Type_1_Percent_Class_6',
  'Land_Cover_Type_1_Percent_Class_7',
  'Land_Cover_Type_1_Percent_Class_10',
  'Land_Cover_Type_1_Percent_Class_12',
  'Land_Cover_Type_1_Percent_Class_13',
  'Land_Cover_Type_1_Percent_Class_14'
]);

// Import the TerraClimate image collection
var temperatureCollection = ee.ImageCollection("IDAHO_EPSCOR/TERRACLIMATE").filterDate(time_start_temperature, time_end_temperature);
var temperatureSubset = temperatureCollection.select(['tmmn', 'tmmx', 'pr']);

var years = [2017, 2018, 2019, 2020, 2021, 2022];

// Function to process each point and extract statistics for NDVI/EVI, Land Cover, and Temperature
var processPoint = function(point) {
  var properties = {
    'cells': point.get('cells'), // Add point ID property
    'longitude': point.geometry().coordinates().get(0),
    'latitude': point.geometry().coordinates().get(1)
  };

  years.forEach(function(year) {
    // NDVI/EVI Processing
    var dateRangeStart_ndviEvi = ee.Date.fromYMD(year, 4, 1);
    var dateRangeEnd_ndviEvi = ee.Date.fromYMD(year, 6, 30);
    var ndviEviImage = ndviEviSubset.filterDate(dateRangeStart_ndviEvi, dateRangeEnd_ndviEvi).mean();
    var ndviEviData = ndviEviImage.reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: point.geometry(),
      scale: 1,
      bestEffort: true
    });

    // Land Cover Processing
    var dateRangeStart_landCover = ee.Date.fromYMD(year, 1, 1);
    var dateRangeEnd_landCover = ee.Date.fromYMD(year, 12, 31);
    var landCoverImage = landCoverSubset.filterDate(dateRangeStart_landCover, dateRangeEnd_landCover).first();
    var landCoverStats = landCoverImage.reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: point.geometry(),
      scale: 1,
      bestEffort: true
    });

    // Temperature Processing
    var dateRangeStart_temperature = ee.Date.fromYMD(year, 4, 1);
    var dateRangeEnd_temperature = ee.Date.fromYMD(year, 6, 30);
    var temperatureImage = temperatureSubset.filterDate(dateRangeStart_temperature, dateRangeEnd_temperature).mean();
    var temperatureData = temperatureImage.reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: point.geometry(),
      scale: 1,
      bestEffort: true
    });

    // Add extracted statistics to properties
    properties['NDVI_' + year] = ndviEviData.get('NDVI');
    properties['EVI_' + year] = ndviEviData.get('EVI');
    properties['Majority_Land_Cover_Type_1_' + year] = landCoverStats.get('Majority_Land_Cover_Type_1');
    properties['Land_Cover_Type_1_Percent_Class_0_' + year] = landCoverStats.get('Land_Cover_Type_1_Percent_Class_0');
    properties['Land_Cover_Type_1_Percent_Class_6_' + year] = landCoverStats.get('Land_Cover_Type_1_Percent_Class_6');
    properties['Land_Cover_Type_1_Percent_Class_7_' + year] = landCoverStats.get('Land_Cover_Type_1_Percent_Class_7');
    properties['Land_Cover_Type_1_Percent_Class_10_' + year] = landCoverStats.get('Land_Cover_Type_1_Percent_Class_10');
    properties['Land_Cover_Type_1_Percent_Class_12_' + year] = landCoverStats.get('Land_Cover_Type_1_Percent_Class_12');
    properties['Land_Cover_Type_1_Percent_Class_13_' + year] = landCoverStats.get('Land_Cover_Type_1_Percent_Class_13');
    properties['Land_Cover_Type_1_Percent_Class_14_' + year] = landCoverStats.get('Land_Cover_Type_1_Percent_Class_14');
    properties['tmmn_' + year] = temperatureData.get('tmmn');
    properties['tmmx_' + year] = temperatureData.get('tmmx');
    properties['pr_' + year] = temperatureData.get('pr');
  });

  // Create a new feature with the accumulated properties and original point geometry
  return ee.Feature(point.geometry(), properties);
};

// Apply the processing function to each point
var pointstats = coords.map(processPoint);

// Export the feature collection to CSV
Export.table.toDrive({
  collection: pointstats,
  description: 'otitar_dynamic_variables',
  fileFormat: 'CSV'
});