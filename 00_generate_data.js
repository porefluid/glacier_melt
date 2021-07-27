var s1 = ee.ImageCollection("COPERNICUS/S1_GRD"),
    s2 = ee.ImageCollection("COPERNICUS/S2"),
    srtm = ee.Image("USGS/SRTMGL1_003"),
    hkh = ee.FeatureCollection("users/coreyscher/HKHoutline"),
    geometry = 
    /* color: #d63000 */
    /* shown: false */
    ee.Geometry.MultiPolygon(
        [[[[70.22190717002036, 40.71381359186029],
           [70.22190717002036, 25.7205660764414],
           [103.88401654502036, 25.7205660764414],
           [103.88401654502036, 40.71381359186029]]],
         [[[86.25697167791463, 28.125581230691978],
           [86.25697167791463, 28.00561256288194],
           [86.54398950018026, 28.00561256288194],
           [86.54398950018026, 28.125581230691978]]]], null, false),
    geometry2 = 
    /* color: #98ff00 */
    /* shown: false */
    ee.Geometry.Polygon(
        [[[86.14593025401427, 28.85806553283277],
          [86.14593025401427, 27.307275851833403],
          [87.20061775401427, 27.307275851833403],
          [87.20061775401427, 28.85806553283277]]], null, false),
    gamdam = ee.FeatureCollection("users/coreyscher/GAMDAM2018"),
    dry_area = 
    /* color: #d63000 */
    /* shown: false */
    ee.Geometry.MultiPolygon(
        [[[[76.6541304425802, 35.63827172939908],
           [76.67575977607629, 35.626830990892714],
           [76.67816303535363, 35.629202985862335],
           [76.68279789253137, 35.62459846075183],
           [76.69464252754113, 35.62529613311733],
           [76.70168064399621, 35.63715563187438],
           [76.70992039008996, 35.64092239855066],
           [76.70322559638879, 35.64538648469737],
           [76.68726108833215, 35.64524698578049],
           [76.68005131050012, 35.643433477696625],
           [76.67421482368371, 35.64566548180045],
           [76.68211124702356, 35.653895458030796],
           [76.69069431587121, 35.66031146140305],
           [76.68399952217004, 35.66407713635255],
           [76.6628851728048, 35.65542976659926],
           [76.6515555219259, 35.654592874670335],
           [76.65687702461145, 35.644828487568425]]],
         [[[76.75695560737512, 35.68555203495484],
           [76.7478575543966, 35.69154724082967],
           [76.73275135322473, 35.68569146346471],
           [76.74579761787317, 35.67802253356749],
           [76.74390934272668, 35.67105013950879],
           [76.73755787177942, 35.66337980253849],
           [76.74871586128137, 35.65501132173098],
           [76.75661228462121, 35.648455399246295],
           [76.77137516303918, 35.65333752033456],
           [76.77875660224817, 35.649850321433306],
           [76.78150318427942, 35.66017198854992],
           [76.77635334297082, 35.6665874876537],
           [76.7698302106466, 35.66379820355779],
           [76.76296375556848, 35.674954755225706]]]]),
    himap = ee.FeatureCollection("users/coreyscher/boundary_mountain_regions_hma_v3_clean-polygon"),
    high_karakorm = 
    /* color: #98ff00 */
    /* shown: false */
    /* displayProperties: [
      {
        "type": "rectangle"
      }
    ] */
    ee.Geometry.Polygon(
        [[[76.45630175552125, 35.78447099314777],
          [76.45630175552125, 35.50183912226047],
          [77.0025282569861, 35.50183912226047],
          [77.0025282569861, 35.78447099314777]]], null, false),
    geometry3 = /* color: #d63000 */ee.Geometry.Point([86.69954899741792, 28.024359267611104]);

//###################################################################
//#####  Snowmelt onset from Sentinel-1 in Hindu Kush Himalayas ####### 
//#####                     A final product                     #######  
//###################################################################

// --------------Runtime parameters ---------------- //

var dryStart = '2018-01-01';
var dryEnd = '2018-02-28';

var meltStart = '2017-07-01';
var meltEnd = '2017-08-30';

var studyStart = '2017-01-01';
var studyEnd = '2018-01-01';

var zLim = 2

var band = ['VH']

var threshold = 3 //dB

// -------------------------------------------------------------------
// Define region 

himap = himap.filterBounds(hkh)

var himap_glac = gamdam.filterBounds(himap);

var gMask = ee.Image().paint(himap_glac, 1).eq(0)

Map.addLayer(gMask)
// -------------------------------------------------------------------------------
// Bring in Radar data -----------------------------------------------------------
// ---- Import Sentinel-1 image collection -------
s1 = s1.filterBounds(hkh)
// Filter to get images with VV and VH dual polarization.
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
// Filter to get images collected in interferometric wide swath mode.
    .filter(ee.Filter.eq('instrumentMode', 'IW'))//.map(maskAngGT30).map(maskAngLT45);

//----Separate S1 by orbit type ---------------------
var s1A = s1.filter(ee.Filter.eq('orbitProperties_pass', 'ASCENDING'))
var s1D = s1.filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'))

// --------------------------------------------------------------------
// Make list of relative orbit numbers
//
/*
var orbs = function(image){
  image = ee.Image(image);
  var start = image.get('relativeOrbitNumber_start');
  return start;
};

var orbitsD = s1D.toList(s1D.size()).map(orbs).distinct();
print(orbitsD);
*/


var orbitsD = ee.List([107,121,122,135,136,150,151,164,165,4,5,19,20,33,34,48,49,62,63,77,78,92,106]);
var orbitsA = ee.List([12,13,26,27,41,42,55,56,70,71,85,86,99,100,114,115,128,129,143,144,158,172,173])


// ----------------------------------------------------------------------
// Functions

// Descending track algorithm

var algDesc = function(orbitNumber){
  
  // Calculate relevant statistics for threshold detection 
  
  var winterMean = s1D.filterDate(dryStart, dryEnd).select(band)
                  .filterMetadata('relativeOrbitNumber_start', 'equals', orbitNumber)
                  .mean();
  var summerMean = s1D.filterDate(meltStart, meltEnd).select(band)
                  .filterMetadata('relativeOrbitNumber_start', 'equals', orbitNumber)
                  .mean();
  var winterStd = s1D.filterDate(dryStart, dryEnd).select(band)
                  .filterMetadata('relativeOrbitNumber_start', 'equals', orbitNumber)
                  .reduce(ee.Reducer.stdDev());
                  
  //Filter image collection by orbit track
  
  var imColl = s1D.filterDate(studyStart, studyEnd).select(band)
                .filterMetadata('relativeOrbitNumber_start', 'equals', orbitNumber);
                
  var zScore = winterMean.subtract(summerMean).divide(winterStd);
  var zMask = zScore.gt(zLim);
  
  //Classify melt
  var meltColl = imColl.map(function(i){
    i = ee.Image(i);
    var date = ee.Image(ee.Number(i.date().getRelative('day', 'year')));

    i = i.lt(winterMean.subtract(threshold))
            .set('system:time_start', i.get('system:time_start'));
    var finalIm = i.addBands(date).updateMask(zMask)
                    .addBands(zScore).addBands(srtm.select(['elevation']))
                    .rename(['meltObs', 'DOY', 'zScore', 'elevation'])
            .toInt()
            .set('layer', 'melt');
    return finalIm
            .updateMask(
              ee.Image(
                finalIm.select(['meltObs']))
                .eq(1));
                
  }); 

  return meltColl;
};

// Ascending track algorithm

var algAsc = function(orbitNumber){
  
  // Calculate relevant statistics for threshold detection 
  
  var winterMean = s1A.filterDate(dryStart, dryEnd).select(band)
                  .filterMetadata('relativeOrbitNumber_start', 'equals', orbitNumber)
                  .mean();
  var summerMean = s1A.filterDate(meltStart, meltEnd).select(band)
                  .filterMetadata('relativeOrbitNumber_start', 'equals', orbitNumber)
                  .mean();
  var winterStd = s1A.filterDate(dryStart, dryEnd).select(band)
                  .filterMetadata('relativeOrbitNumber_start', 'equals', orbitNumber)
                  .reduce(ee.Reducer.stdDev());
                  
  //Filter image collection by orbit track
  
  var imColl = s1A.filterDate(studyStart, studyEnd).select(band)
                .filterMetadata('relativeOrbitNumber_start', 'equals', orbitNumber);
                
  var zScore = winterMean.subtract(summerMean).divide(winterStd);
  var zMask = zScore.gt(zLim);
  
  //Classify melt
  var meltColl = imColl.map(function(i){
    i = ee.Image(i);
    var date = ee.Image(ee.Number(i.date().getRelative('day', 'year')));

    i = i.lt(winterMean.subtract(threshold))
            .set('system:time_start', i.get('system:time_start'));
    var finalIm = i.addBands(date).updateMask(zMask)
                    .addBands(zScore).addBands(srtm)
                    .rename(['meltObs', 'DOY', 'zScore', 'elevation'])
                    .toInt()
            .set('layer', 'melt');
    return finalIm
            
            .updateMask(
              ee.Image(
                finalIm.select(['meltObs']))
                .eq(1));
  
  }); 

  return meltColl;
};

// -------------- Final product generation ------------ //

var productD = orbitsD.map(algDesc);
var productA = orbitsA.map(algAsc);

var product = productD.cat(productA);


var melt = product.map(function(imageCol){
  imageCol = ee.ImageCollection(imageCol);
  
  // melt images
  var zScore = imageCol.filterMetadata('layer', 'equals', 'melt').select(['zScore']).mean();
  var meltCol = imageCol.filterMetadata('layer', 'equals', 'melt').select(['meltObs', 'DOY']);
  var onset = meltCol.select(['DOY']).reduce(ee.Reducer.firstNonNull()).toInt();
  var offset = meltCol.select(['DOY']).reduce(ee.Reducer.lastNonNull()).toInt();
  var elev = imageCol.select(['elevation']);
  var count = ee.Image(imageCol.select(['meltObs']).map(function(im){return im.unmask()}).count());
  var daysMeltObs = meltCol.select(['meltObs']).sum();//.divide(imageCol.select(['meltObs']).count());
  var meltFrac = daysMeltObs.divide(count);
  
  return ee.Image(onset.addBands(offset)
              .addBands(daysMeltObs)
              .addBands(zScore).addBands(elev.mean()).addBands(count).addBands(ee.Image.pixelArea())
        .rename(['meltOn', 'meltOff', 'meltObs', 'zScore', 'elev','count','area'])
        .updateMask(gMask.not()));
});

// Merge collections
var finalProduct = ee.ImageCollection(melt).qualityMosaic('zScore');



var meltFraction = finalProduct.select(['meltObs']).divide(finalProduct.select(['count'])).multiply(100).rename(['meltFrac']).toInt()
finalProduct = finalProduct.addBands(meltFraction)
              .focal_median(9, "square", "pixels").updateMask(gMask.not())
              //.updateMask(finalProduct.select('zScore').gt(threshold))
              .addBands(gMask).toInt();
var meltVis = {"opacity":1,"bands":['meltOff'],  "min":250,"max":350,"palette":["00ffdb","005aff","002387"]};
var fracVis = {"opacity":1,"bands":['meltFrac'],  "min":0,"max":100,"palette":["00ffdb","005aff","002387"]};


Map.addLayer(finalProduct, meltVis)


Export.image.toAsset({image: finalProduct, 
                      description: '2017_glacMelt_90m_072721_zmask', 
                      scale: 90, 
                      region: hkh.geometry().bounds().buffer(1000),
                      maxPixels: 1e13
});








//Map.addLayer(finalProduct, meltVis, 'final melt product');

//Map.addLayer(meltFraction, fracVis)

//Reduce regions by the grid and count melt area, onset, offset, zscore, etc

//grid = grid.filterBounds(hhk)
//Map.addLayer(grid)
/*
var gridStats = finalProduct.select(['area']).rename(['areaDetected'])
            .reduceRegions({
              collection: grid,
              reducer: ee.Reducer.sum(),
              scale: 90,
              tileScale: 2
            });
            
gridStats = ee.Image.pixelArea().rename('glacierArea')
        .updateMask(gMask)
        .reduceRegions({
          collection: grid, 
          reducer: ee.Reducer.sum(), 
          scale: 90,
          tileScale: 2
        });

gridStats = finalProduct.select(['meltOn', 'meltOff', 'meltFrac', 'zScore']).rename(['onSt', 'offSt', 'fracSt', 'zSt'])
            .reduceRegions({
              collection: grid,
              reducer: ee.Reducer.stdDev(),
              scale: 90,
              tileScale: 2
            });
            
gridStats = finalProduct.select(['meltOn', 'meltOff', 'meltFrac', 'zScore'])
            .reduceRegions({
              collection: grid,
              reducer: ee.Reducer.mean(),
              scale: 90,
              tileScale: 2
            });
            

Export.table.toDrive({collection: gridStats, description: 'gridded_stats_2019_smaller', folder: 'gridded_melt_stats'});

 Map.addLayer(hhk)
*/

/*
// Export elevation statistics

var bins = ee.List.sequence(2000, 9000, 100);
var i = 60

var bin = ee.Number(bins.get(i));
var n = ee.String('meltElevationBin2019Agg_').cat(bin.getInfo().toString());
var elMask = srtm.lt(bin).or(srtm.gt(bin.add(100))).not()

var im = finalProduct;
//var im = finalProduct.addBands(slope);

var area = im.select(['elev', 'area']).rename(['elev', 'areaDetected'])
            .reduceRegions({
              collection: glims,
              reducer: ee.Reducer.sum(),
              scale: 90,
              tileScale: 2
            });

area = im.select(['meltOn', 'meltOff', 'meltFrac']).rename(['onSt', 'offSt', 'fracSt'])
            .reduceRegions({
              collection: area,
              reducer: ee.Reducer.stdDev(),
              scale: 90,
              tileScale: 2
            });

area = im.select(['meltOn', 'meltOff', 'meltFrac'])
            .reduceRegions({
              collection: area,
              reducer: ee.Reducer.mean(),
              scale: 90,
              tileScale: 2
            });


var area = im.select(['zScore']).rename(['zScore_count'])
            .reduceRegions({
              collection: glims,
              reducer: ee.Reducer.count(),
              scale: 90,
              tileScale: 2
            });
var area = im.select(['zScore']).rename(['zScore_sum'])
            .reduceRegions({
              collection: area,
              reducer: ee.Reducer.sum(),
              scale: 90,
              tileScale: 2
            });
            
Map.addLayer(im)

print(area.first())
Export.table.toAsset({collection: area, description: 'z_score_sum_count_2017'});
*/
