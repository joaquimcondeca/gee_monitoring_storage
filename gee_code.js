/////////// JAVASCRIPT CODE DEVELOPED ////////////////////////////////
////////// Monitoring the storage volume of water reservoirs using Google Earth Engine /////////////////////////////////////
////////// Joaquim Condeça, João Nascimento and Nuno Barreiras /////////////////////////////////////////////////////////
///////////WATER AREA(m2) FOR ALVITO, CAIA, MARANHÃO AND ROXO RESERVOIRES - LANDSAT 4,5 e 8/////////////////////////////////
///////////EXEMPLE FOR MNDWI, XU,2006///////////
////////// NDWI = (green - mir) / (green + mir)/////////////////////////////////////////////////////////
////////// green: Banda B2   //////////////////////////////////////////////////////////////////////////
////////// mir:   Banda B5  //////////////////////////////////////////////////////////////////////////
///////// THE var albufeira must be include as Assets (shapefile geometry of the reservoirs)


/// BEGIN Functions cloudMaskL457 and maskL8sr ///
/// Functions to mask clouds based on the pixel_qa band of Landsat SR data
/// SOURCE: https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LT05_C01_T1_SR#bands
/// The recommendation is the use of two functions (cloudMaskL457, for Landsat 4 and 5 images, and maskL8sr, for Landsat 8 images), 
/// which use the “pixel_qa” band of Landsat SR images to mask the clouds (GEE, 2020).
/// pixel_qa -> Pixel quality attributes generated from the CFMASK algorithm.

var cloudMaskL457 = function(image) {
  var qa = image.select('pixel_qa');
  // If the cloud bit (5) is set and the cloud confidence (7) is high
  // or the cloud shadow bit is set (3), then it's a bad pixel.
  var cloud = qa.bitwiseAnd(1 << 5)
                  .and(qa.bitwiseAnd(1 << 7))
                  .or(qa.bitwiseAnd(1 << 3));
  // Remove edge pixels that don't occur in all bands
  var mask2 = image.mask().reduce(ee.Reducer.min());
  return image.updateMask(cloud.not()).updateMask(mask2);
};

/**
 * Function to mask clouds based on the pixel_qa band of Landsat 8 SR data.
 * @param {ee.Image} image input Landsat 8 SR image
 * @return {ee.Image} cloudmasked Landsat 8 image
 */
function maskL8sr(image) {
  // Bits 3 and 5 are cloud shadow and cloud, respectively.
  var cloudShadowBitMask = (1 << 3);
  var cloudsBitMask = (1 << 5);
  // Get the pixel QA band.
  var qa = image.select('pixel_qa');
  // Both flags should be set to zero, indicating clear conditions.
  var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)
                 .and(qa.bitwiseAnd(cloudsBitMask).eq(0));
  return image.updateMask(mask);
}

/// END Functions cloudMaskL457 and maskL8sr ///



/// BEGIN Function rejectImagesWithValidAreaUnderLimit ///
/// Function to removes images from the collection (var img) 
/// that contain pixels affected by clouds (already masked) 
/// that remain within the reservoir area
/// The application of this function allows obtaining a new collection
/// of images (var img_00) with adequate quality for the calculation of water indices
/// It retrieves the number of valid pixels of each image over the reservoir and corresponding area. 
/// Then, dividing it by the reservoir area, it obtains an index that varies between 0 and 1,
/// where 0 means that clouds completely cover the reservoir, and 1 that there are no clouds.
///////// THE var albufeira must be include as Assets (shapefile geometry of the reservoirs)

function rejectImagesWithValidAreaUnderLimit(image) { 
  var limite = 1; 
  var scaleforfunction = 30;
  var pixelscount = image.select('B3').reduceRegion({ //Count pixel number over revervoir geometry
  reducer: ee.Reducer.count(),
  geometry: albufeira,
  crs: 'EPSG:3763', // PORTUGAL
  scale: 30,
  maxPixels: 1e13
  });
  var numpixelscount = ee.Number(pixelscount.get('B3'));
  var areanumpixelscount = ee.Number(numpixelscount).multiply(900);//set valid pixel area over the geometry reservoir
  var areascore = areanumpixelscount.divide(ee.Number(areaalb)); //set relation between valid pixel area and the reservoir area
  var usar = areascore.gte(limite);//set "0" or "1" to var usar
  // 0 means that clouds completely cover the reservoir
  // 1 means that there are no clouds
  var usarnumero = ee.Number(usar);
  var valor = (image.select('B3').multiply(0).toFloat()).add(usar)
  .rename('B99') //set "0" or "1" to B99 band
  var time0 = image.get('system:time_start');
  var ntime0 = ee.Number(image.get('system:time_start'));
  return(image.addBands(valor).set('system:time_start',ntime0.multiply(usarnumero))); 
  // Images with the null 'system:time_start' property will be discarded 
  // When applying the function to an image collection with a time filter, 
  //all images with 'system:time_start' null will be ignored and removed from the resulting collection

}

/// END rejectImagesWithValidAreaUnderLimit ///


/// BEGIN Function rejectImagesThatDontFullyOverlapROI ///
/// This function returns images with the null 'system:time_start' 
/// property in case of images to be discarded. 
/// When applying the function to an image collection with a time filter, 
/// all images with 'system:time_start' null will be ignored and removed from the resulting collection.

function rejectImagesThatDontFullyOverlapROI(image) {
  var geometriaimagem = image.geometry();
  var geometriaalbufeira = albufeira.geometry();
  var intersection = geometriaalbufeira.intersection(geometriaimagem, ee.ErrorMargin(1));
  var areaintersection = intersection.area();
  var areaalbufeira = albufeira.geometry().area();
  // usar = 1 , good image
  // usar = 0 , reject image
  var usar = areaintersection.gte(areaalbufeira);
  var usarnumero = ee.Number(usar);
  var valor = (image.select('B3').multiply(0).toFloat()).add(usar)
  .rename('B98')
  var time0 = image.get('system:time_start');
  var ntime0 = ee.Number(image.get('system:time_start'));
  return(image.addBands(valor).set('system:time_start',ntime0.multiply(usarnumero)));
}

/// END Function rejectImagesThatDontFullyOverlapROI ///


//////////////////////////BEGIN WATER FUNCTION///////////////////////////////////////////////////////

//EXEMPLE FOR MNDWI INDEX
///////////MNDWI, SEGUNDO XU,2006///////////
////////// MNDWI = (green - mir) / (green + mir)/////////////////////////////////////////////////////////
////////// green: Banda B2   //////////////////////////////////////////////////////////////////////////
////////// mir:   Banda B5  //////////////////////////////////////////////////////////////////////////


var waterThreshold = 0; // Threshold

var waterfunction = function(image){

var ndwi = image.expression(
    '(green - mir) / (green + mir)',
    {
        green: image.select('B2'),   
        mir: image.select('B5')
    }
  )
.rename('NDWI');

// To use other indices replace in var ndwi the respective formulas:
///////////NDWI, MCFEETERS,1996///////////
////////// NDWI = (green - nir) / (green + nir)/////////////////////////////////////////////////////////
////////// green: Banda B2   //////////////////////////////////////////////////////////////////////////
////////// nir:   Banda B4  //////////////////////////////////////////////////////////////////////////

///////////NDWI, GAO,1996///////////
////////// NDWI = (nir - mir) / (nir + mir)/////////////////////////////////////////////////////////
////////// nir: Banda B4   //////////////////////////////////////////////////////////////////////////
////////// mir:   Banda B5  //////////////////////////////////////////////////////////////////////////


  //Indices greater than Threshold
  var water01 = ndwi.gt(waterThreshold);

  //Pixel mask
  image = image.updateMask(water01).addBands(ndwi).clip(albufeira);

  //Area for the water pixels, in m2
  var area = ee.Image.pixelArea();
  var waterArea = water01.multiply(area).rename('waterArea');
  image = image.addBands(waterArea);

  var stats = waterArea.reduceRegion({
    reducer: ee.Reducer.sum(), 
    geometry: albufeira, 
    scale: 30,
  });

  return image.set(stats);
};

//////////////////////////END WATER FUNCTION///////////////////////////////////////////////////////

/*
/**************************************************************************************************
*
*       END OF FUNCTIONS
*
***************************************************************************************************
*/


var visParams = { bands: ['B3','B2','B1'],min: 0, max: 3000, gamma: 1.4, };
Map.centerObject(albufeira, 8);

var startdate = ee.Date('1982-01-01');
var enddate = ee.Date('2019-12-31');


var L8 = ee.ImageCollection('LANDSAT/LC08/C01/T1_SR')
  .filterBounds(albufeira)// geometry filter
  .map(function(image){
  return image.select(
      ['B1', 'B2', 'B3', 'B4', 'B5','B6', 'pixel_qa'],
      ['B00', 'B1', 'B2', 'B3', 'B4','B5', 'pixel_qa']);
  })
  .select(['B1', 'B2', 'B3', 'B4', 'B5', 'pixel_qa']);

var L5 = ee.ImageCollection('LANDSAT/LT05/C01/T1_SR')
.filterBounds(albufeira) // geometry filter
.select(['B1', 'B2', 'B3', 'B4', 'B5', 'pixel_qa']);

var L4 = ee.ImageCollection('LANDSAT/LT04/C01/T1_SR')
.filterBounds(albufeira) // geometry filter
.select(['B1', 'B2', 'B3', 'B4', 'B5', 'pixel_qa']);

var img = ee.ImageCollection(L4.merge(L5.merge(L8)));


var imgL45 = ee.ImageCollection(L4.merge(L5));


///// BEGIN Reservoir Area
//////////////////////////
var AreaAlb0 = ee.Image(1);
var scaleforAlb = 30;
var reducer = AreaAlb0.reduceRegion({
  reducer: ee.Reducer.sum(),
  geometry: albufeira,
  crs: 'EPSG:3763', // PORTUGAL
  scale: scaleforAlb,
  maxPixels: 1E13
});
// value in m2
var areaalb = ee.Number(reducer.get('constant')).multiply(scaleforAlb).multiply(scaleforAlb);
print('Area da albufeira calculated using pixel count method: ', areaalb.getInfo() + ' m2');

///// END Reservoir Area


// Remove images without data 
// They does not overpal the entire reservoir
//
// 1a - Apply filter area to var L45
var L45_00 = imgL45.map(rejectImagesThatDontFullyOverlapROI)
.filterDate(startdate, enddate);
print('L45_00 Collection: Number of images: ', L45_00.size());
// 1b - Apply filter area to var L8
var L8_00 = L8.map(rejectImagesThatDontFullyOverlapROI)
.filterDate(startdate, enddate);
print('L8_00 Collection: Number of images: ', L8_00.size());

// 2a - Apply filter mask to Landsat 4 and 5 collection
var L45_0 = L45_00.map(cloudMaskL457);
print('L45_0 Collection: Number of images: ', L45_0.size());
// 2a - Apply filter mask to Landsat 8 collection
var L8_0 = L8_00.map(maskL8sr);
print('L8_0 Collection: Number of images: ', L8_0.size());

// Merge var L45_0 and L8_0
var img = ee.ImageCollection(L45_0.merge(L8_0));
print('img Collection: Number of images: ', img.size());

/// Apply function rejectImagesWithValidAreaUnderLimit
/// To the collection (var img) 
/// And remove images that contain pixels affected by clouds (already masked)
/// that remain within the reservoir area

var img_00 = img.map(rejectImagesWithValidAreaUnderLimit)
.filterDate(startdate, enddate);
print('Fltered Collection: Number of images: ', img_00.size());

/////////////////////////////////////////////////////////////////////////
/// APPLY WATER FUNCTION TO IMG_00 COLLECTION ///
/////////////////////////////////////////////////////////////////////////

var collection = img_00.map(waterfunction);
print(collection);

// Define chart customization options.
var options = {
  lineWidth: 1,
  pointSize: 2,
  hAxis: {title: 'Date'},
  vAxis: {title: 'Area (m2)'},
  title: 'Roxo Reservoir'
};

var chart = ui.Chart.image.series({
  imageCollection: collection.select('waterArea'), 
  region: albufeira, 
  reducer: ee.Reducer.sum(), 
  scale: 30})
  .setOptions(options);
  
print(chart);


// Exportar para csv
Export.table.toDrive({
  collection: collection,
  description: 'description',
  folder: 'folder',
  fileFormat: 'CSV'
}); 


//////////////////////////////////////////////////////////////////

Map.addLayer(albufeira);

print('  --  The END  --  ');
