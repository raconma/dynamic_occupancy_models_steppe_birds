// input
var coords = ee.FeatureCollection('projects/ee-rcontr03/assets/ebd_otitar_breeding_spain_zf_part2_landcover');

Map.addLayer(coords, {color: 'green'}, 'Otis tarda. Presence - Absence');
Map.centerObject(coords, 11);

// Define the input parameters
var time_start = '2010-04-01';
var time_end = '2023-07-30';

// Import the MODIS NDVI and EVI image collection
var ndviEviCollection = ee.ImageCollection("MODIS/061/MOD13A3").filterDate(time_start, time_end);
var ndviEviSubset = ndviEviCollection.select(['NDVI', 'EVI']);

// Import the TerraClimate image collection
var temperatureCollection = ee.ImageCollection("IDAHO_EPSCOR/TERRACLIMATE").filterDate(time_start, time_end);
var temperatureSubset = temperatureCollection.select(['tmmn', 'tmmx', 'pr']);

// Map over the points, filter the collection for the date range of every point, and get the mean data for that point
var pointstats = coords.map(function(point) {
  // Get the year from the point feature
  var year = ee.Number(point.get('year')).toInt();
  
  // Define the date range, our study season
  var dateRangeStart = ee.Date.fromYMD(year, 4, 1);
  var dateRangeEnd = ee.Date.fromYMD(year, 7, 30);
  
  // Filter the NDVI/EVI collection for the date range and calculate the mean image
  var ndviEviImage = ndviEviSubset.filterDate(dateRangeStart, dateRangeEnd).mean();
  
  // Reduce the NDVI/EVI image to get the mean NDVI and EVI values for the point
  var ndviEviData = ndviEviImage.reduceRegion({
    reducer: ee.Reducer.mean(), // Use mean reducer to get the mean value at the point
    geometry: point.geometry(),
    scale: 1, // Point scale
    bestEffort: true
  });
  
  // Filter the temperature collection for the date range and calculate the mean image
  var temperatureImage = temperatureSubset.filterDate(dateRangeStart, dateRangeEnd).mean();
  
  // Reduce the temperature image to get the mean temperature values for the point
  var temperatureData = temperatureImage.reduceRegion({
    reducer: ee.Reducer.mean(), // Use mean reducer to get the mean value at the point
    geometry: point.geometry(),
    scale: 1, // Point scale
    bestEffort: true
  });
  
  // Extract coordinates from geometry
  var coords = point.geometry().coordinates();
  
  // Merge the NDVI/EVI and temperature data into a single feature
  var data = ee.Feature(point.setMulti(ndviEviData.combine(temperatureData)))
  // Export lat long as separate columns not as a geojson object
                 .set('longitude', coords.get(0)) 
                 .set('latitude', coords.get(1));
  
  return data;
});

// Export the results to CSV
Export.table.toDrive({
  collection: ee.FeatureCollection(pointstats),
  description: 'ebd_otitar_breeding_spain_zf_part2_variables',
  fileFormat: 'CSV'
});
