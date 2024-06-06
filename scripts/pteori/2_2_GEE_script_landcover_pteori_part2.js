// Import the second part of the data
var coords = ee.FeatureCollection('projects/ee-rcontr03/assets/ebd_pteori_breeding_spain_zf_part2');

// Visualize the data on the map
Map.addLayer(coords, {color: 'green'}, 'Pterocles orientalis. Presence - Absence');
Map.centerObject(coords, 11);

// Define the time range for filtering the MODIS image collection
var time_start = '2010-01-01';
var time_end = '2023-12-31';

// Import the MODIS land cover image collection and filter by date
var collection = ee.ImageCollection("MODIS/061/MCD12C1").filterDate(time_start, time_end);

// Select the band from the collection
var collectionSubset = collection.select(['Majority_Land_Cover_Type_1','Land_Cover_Type_1_Percent_Class_0','Land_Cover_Type_1_Percent_Class_13']);

// Function to process each point
var pointstats = coords.map(function(point){
  // Get the actual year for that point
  var year = ee.Number(point.get('year')).toInt();
  
  // MODIS/061/MCD12Q1 is a yearly collection, the date range must be the whole year
  var dateRangeStart = ee.Date.fromYMD(year, 1, 1);
  var dateRangeEnd = ee.Date.fromYMD(year, 12, 31);
  
  // Filter the collection for the date range and select the first image
  var image = collectionSubset.filterDate(dateRangeStart, dateRangeEnd).first();
  
  // Reduce the image to a single value at the point's location
  var data = image.reduceRegion({
      reducer: ee.Reducer.first(), // Use first to get the land cover type at the point
      geometry: point.geometry(),
      scale: 1 // 1 so it just uses the exact point
  });
  
  // Return the point with the new data properties added
  return point.setMulti(data);
});

// .toAsset to store it locally
Export.table.toAsset({
  collection: pointstats,
  assetId: 'projects/ee-rcontr03/assets/ebd_pteori_breeding_spain_zf_part2_landcover',
  description: 'ebd_pteori_breeding_spain_zf_part1_landcover'
});
