/************
 * SCRIPT FINAL COMPLET GEE — ORAN 2025 (CORRIGÉ)
 * Surveillance côtière : hydrocarbures, algues, turbidité
 * Sentinel-1 / Sentinel-2 / MODIS / ERA5
 ************/

/****
 * 0. PARAMÈTRES
 ****/
var startDate = '2025-01-01';
var endDate   = '2025-06-30';

/****
 * 1. WILAYA D’ORAN
 ****/
var gaul = ee.FeatureCollection('FAO/GAUL/2015/level1');
var oran = gaul.filter(
  ee.Filter.and(
    ee.Filter.eq('ADM0_NAME', 'Algeria'),
    ee.Filter.eq('ADM1_NAME', 'Oran')
  )
);
Map.centerObject(oran, 9);
Map.addLayer(oran, {color: 'red'}, 'Wilaya Oran');

/****
 * 2. BUFFER CÔTIER 20 km
 ****/
var buffer20km = oran.geometry().buffer(20000);

/****
 * 3. MASQUE MARIN
 ****/
var water = ee.Image('JRC/GSW1_4/GlobalSurfaceWater')
  .select('occurrence')
  .gt(0);

var marineBuffer = water.updateMask(water).clip(buffer20km);
Map.addLayer(marineBuffer, {palette:['blue']}, 'Zone marine 20 km');

/****
 * 4. SENTINEL-2 — PRÉTRAITEMENT
 ****/
function maskS2clouds(image) {
  var scl = image.select('SCL');
  var mask = scl.neq(3)
    .and(scl.neq(7))
    .and(scl.neq(8))
    .and(scl.neq(9))
    .and(scl.neq(10));
  return image.updateMask(mask)
    .divide(10000)
    .copyProperties(image, ['system:time_start']);
}

var s2 = ee.ImageCollection('COPERNICUS/S2_SR')
  .filterDate(startDate, endDate)
  .filterBounds(marineBuffer.geometry())
  .map(maskS2clouds);

var s2_median = s2.median().clip(marineBuffer.geometry());

Map.addLayer(
  s2_median.select(['B4','B3','B2']),
  {min:0, max:0.3},
  'Sentinel-2 RGB'
);

/****
 * 5. INDICES SPECTRAUX
 ****/
// NDWI
var ndwi = s2_median.normalizedDifference(['B3','B8']).rename('NDWI');

// NDCI
var ndci = s2_median.normalizedDifference(['B5','B4']).rename('NDCI');

// FAI (formule scientifique correcte)
var fai = s2_median.expression(
  'B8 - (B4 + (B11 - B4) * ((833 - 665) / (1610 - 665)))',
  {
    B4: s2_median.select('B4'),
    B8: s2_median.select('B8'),
    B11: s2_median.select('B11')
  }
).rename('FAI');

Map.addLayer(ndwi, {min:-1, max:1, palette:['brown','blue']}, 'NDWI');

/****
 * 6. MASQUE EAU
 ****/
var waterMask = ndwi.gt(0);
var ndci_w = ndci.updateMask(waterMask);
var fai_w  = fai.updateMask(waterMask);

/****
 * 7. ALGUES
 ****/
var algaeMask = fai_w.gt(0).or(ndci_w.gt(0.05));
Map.addLayer(algaeMask.updateMask(algaeMask),
  {palette:['green']}, 'Algues');

/****
 * 8. TURBIDITÉ
 ****/
var turbidityMask = ndwi.lt(0.1).and(waterMask);
Map.addLayer(turbidityMask.updateMask(turbidityMask),
  {palette:['yellow']}, 'Turbidité');

/****
 * 9. SENTINEL-1 — PRÉTRAITEMENT
 ****/
function removeBorderNoise(img) {
  var edge = img.lt(0.0001);
  var mask = edge.focal_min(1).focal_max(1).not();
  return img.updateMask(mask);
}

// 9.2 Lee Filter (sur valeurs linéaires)
function leeFilterLinear(img) {
  var kernel = ee.Kernel.square(3);
  var mean = img.reduceNeighborhood(ee.Reducer.mean(), kernel);
  var variance = img.reduceNeighborhood(ee.Reducer.variance(), kernel);
  var noiseVar = ee.Number(0.25);
  var weight = variance.subtract(noiseVar)
    .divide(variance)
    .clamp(0, 1);
  return mean.add(weight.multiply(img.subtract(mean)));
}

// 9.3 Conversion en dB
function toDB(img) {
  return img.log10().multiply(10);
}

// 9.4 Collection Sentinel-1 propre
var s1 = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterDate(startDate, endDate)
  .filterBounds(marineBuffer.geometry())
  .filter(ee.Filter.eq('instrumentMode', 'IW'))
  .filter(ee.Filter.eq('orbitProperties_pass', 'ASCENDING'))
  .filter(ee.Filter.listContains(
    'transmitterReceiverPolarisation', 'VV'))
  .select('VV')
  .map(removeBorderNoise)
  .map(leeFilterLinear)
  .map(toDB);

// Image moyenne (référence SAR)
var s1_mean = s1.mean().clip(marineBuffer.geometry());

// 
Map.addLayer(
  s1_mean,
  {min:-30, max:-5},
  'Sentinel-1 VV (référence)'
);
/****
 * 10. HYDROCARBURES
 ****/
/** 10. HYDROCARBURES (S1) **/

// 1. Sélection des composantes du vent et moyenne
var windCollection = ee.ImageCollection('ECMWF/ERA5/HOURLY')
  .filterDate('2025-03-01', '2025-03-31')
  .filterBounds(marineBuffer.geometry())
  .select(['u_component_of_wind_10m', 'v_component_of_wind_10m']);

// 2. Moyenne pour obtenir une seule image
var windImage = windCollection.mean().clip(marineBuffer.geometry());

// 3. Calcul de la vitesse du vent
var windSpeed = windImage.select('u_component_of_wind_10m')
  .hypot(windImage.select('v_component_of_wind_10m'));

// 4. Seuil SAR
var oilThreshold = -17;
var oilMask_raw = s1_mean.lt(oilThreshold);

// 5. Filtrage par vent faible et eau
var oilMask = oilMask_raw
  .updateMask(windSpeed.lt(6))
  .updateMask(waterMask);

// 6. Affichage final
Map.addLayer(oilMask, {palette: ['black']}, 'Hydrocarbures détectés (S1)');

/** FLÈCHES DE VENT — SAFE VERSION GEE **/

var wind = ee.ImageCollection('ECMWF/ERA5/HOURLY')
  .filterDate('2025-03-01', '2025-03-31')
  .filterBounds(marineBuffer.geometry())
  .select([
    'u_component_of_wind_10m',
    'v_component_of_wind_10m'
  ])
  .mean()
  .clip(marineBuffer.geometry());

//
var points = wind.sample({
  region: marineBuffer.geometry(),
  scale: 30000,   // 30 km
  geometries: true
});

// 
var arrows = points.map(function(f) {

  var u = ee.Number(f.get('u_component_of_wind_10m'));
  var v = ee.Number(f.get('v_component_of_wind_10m'));

  var speed = u.hypot(v);
  var angle = v.atan2(u);

  var length = speed.multiply(3000); 
  var headLength = 0.2;           
  var headAngle = 0.3;              

  //
  var start = f.geometry();

  // 
  var start3857 = start.transform('EPSG:3857', 1);
  var coords = ee.List(start3857.coordinates());
  var x = ee.Number(coords.get(0));
  var y = ee.Number(coords.get(1));

  // 
  var dx = length.multiply(angle.cos());
  var dy = length.multiply(angle.sin());
  var end3857 = ee.Geometry.Point([x.add(dx), y.add(dy)], 'EPSG:3857');

  // 
  var leftX = x.add(length.multiply(headLength).multiply(angle.cos().multiply(Math.cos(headAngle)).subtract(angle.sin().multiply(Math.sin(headAngle)))));
  var leftY = y.add(length.multiply(headLength).multiply(angle.sin().multiply(Math.cos(headAngle)).add(angle.cos().multiply(Math.sin(headAngle)))));
  var rightX = x.add(length.multiply(headLength).multiply(angle.cos().multiply(Math.cos(headAngle)).add(angle.sin().multiply(Math.sin(headAngle)))));
  var rightY = y.add(length.multiply(headLength).multiply(angle.sin().multiply(Math.cos(headAngle)).subtract(angle.cos().multiply(Math.sin(headAngle)))));

  var end4326 = end3857.transform('EPSG:4326', 1);
  var left = ee.Geometry.Point([leftX, leftY], 'EPSG:3857').transform('EPSG:4326', 1);
  var right = ee.Geometry.Point([rightX, rightY], 'EPSG:3857').transform('EPSG:4326', 1);

  //
  var mainLine = ee.Geometry.LineString([start.coordinates(), end4326.coordinates()]);
  var leftLine = ee.Geometry.LineString([end4326.coordinates(), left.coordinates()]);
  var rightLine = ee.Geometry.LineString([end4326.coordinates(), right.coordinates()]);

  var arrowGeom = mainLine.union(leftLine).union(rightLine);

  return ee.Feature(arrowGeom, {speed: speed});
});

// 
Map.addLayer(
  arrows,
  {color: 'black', width: 2},
  'Flèches du vent ERA5'
);

/****
 * 12. BATHYMÉTRIE
 ****/
var bathy = ee.Image('NOAA/NGDC/ETOPO1')
  .select('bedrock')
  .clip(marineBuffer.geometry());

Map.addLayer(bathy,
  {min:-4000, max:0, palette:['darkblue','blue','cyan']},
  'Bathymétrie');

/****
 * 13. MODIS — SST
 ****/
var modis_sst = ee.ImageCollection(
  'NASA/OCEANDATA/MODIS-Aqua/L3SMI')
  .filterDate(startDate, endDate)
  .select('sst')
  .mean()
  .clip(marineBuffer.geometry());

Map.addLayer(modis_sst,
  {min:12, max:30, palette:['blue','cyan','yellow','red']},
  'MODIS SST');

/************
 * PORTS DE LA WILAYA D’ORAN (MANUEL – COMPLET)
 * Types : Commercial, Industriel, Militaire, Pêche
 ************/
var ports = ee.FeatureCollection([

  // --- Ports commerciaux & industriels ---
  ee.Feature(
    ee.Geometry.Point([-0.6510, 35.7040]),
    {name: 'Port d’Oran', type: 'Commercial'}
  ),

  ee.Feature(
    ee.Geometry.Point([-0.3180, 35.8470]),
    {name: 'Port d’Arzew', type: 'Industriel / Hydrocarbures'}
  ),

  ee.Feature(
    ee.Geometry.Point([-0.2460, 35.7980]),
    {name: 'Port de Bethioua', type: 'Pétrochimique'}
  ),

  // --- Port militaire ---
  ee.Feature(
    ee.Geometry.Point([-0.7140, 35.7340]),
    {name: 'Port de Mers El Kébir', type: 'Militaire'}
  ),

  // --- Ports de pêche ---
  ee.Feature(
    ee.Geometry.Point([-0.7700, 35.7430]),
    {name: 'Port d’Aïn El Turk', type: 'Pêche'}
  ),

  ee.Feature(
    ee.Geometry.Point([-0.8340, 35.7210]),
    {name: 'Port de Bousfer', type: 'Pêche'}
  ),

  ee.Feature(
    ee.Geometry.Point([-0.5665, 35.7745]),
    {name: 'Port de Kristel', type: 'Pêche'}
  )
]);

Map.addLayer(
  ports,
  {color: 'black', pointSize: 6},
  'Ports – Wilaya d’Oran'
);

Map.addLayer(ports, {color: 'black'}, 'Ports');

/****
 * 15. EXPORTS
 ****/
Export.image.toDrive({
  image: oilMask,
  description: 'OilSpill_S1_Oran',
  scale: 10,
  region: marineBuffer.geometry(),
  maxPixels: 1e13
});

Export.image.toDrive({
  image: algaeMask,
  description: 'Algae_S2_Oran',
  scale: 10,
  region: marineBuffer.geometry(),
  maxPixels: 1e13
});

Export.image.toDrive({
  image: turbidityMask,
  description: 'Turbidity_S2_Oran',
  scale: 10,
  region: marineBuffer.geometry(),
  maxPixels: 1e13
});