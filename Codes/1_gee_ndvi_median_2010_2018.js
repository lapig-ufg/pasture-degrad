var L8SR = ee.ImageCollection("LANDSAT/LC08/C01/T1_SR"),
    L5SR = ee.ImageCollection("LANDSAT/LT05/C01/T1_SR"),
    LANDSAT_GRID = ee.FeatureCollection("users/vieiramesquita/LAPIG-PASTURE/VECTORS/LANDSAT_GRID_V2_PASTURE"),
    brasil = ee.FeatureCollection("users/claudineisan/limites_geopoliticos/pa_br_brasil"),
    Biomas = ee.FeatureCollection("users/vieiramesquita/lm_bioma_250");


//------------------------------------------------------------------------//
//-------------------Pasture Quality Analalyse Approach-------------------//
//--------------------------NDVIq and CVP---------------------------------//
//------------------------------------------------------------------------//

// Downlaod da série temporal de dados de NDVI (Mediana Bianual corrigida e equalizada mensalmente)
// Anos de 2010-2018, utilizando a respectiva máscaras de pastagens;

var biomList = ['Pampa','Pantanal','Cerrado','Amazônia','Caatinga','Mata Atlântica']

for (var biom = 0; biom < 6; biom++){

var monthBand = function(img){
  
  var monthImg = ee.Image(ee.Number(ee.Date(img.get('SENSING_TIME')).get('month')))
  
  return img.addBands(monthImg.clip(img.geometry()).rename('month_band'))
}

var factorA = 0.99082355;
var factorB = 0.02726984;

var LANDSAT_GRID_List = LANDSAT_GRID.toList(LANDSAT_GRID.size());

var getBRNDVI = LANDSAT_GRID_List.map(function(feature){
  
  var wrsProp = ee.String(ee.Feature(feature).get('SPRNOME'));
  var WRS_PATH = ee.Number.parse(ee.String(wrsProp.split('/').get(0)));
  var WRS_ROW = ee.Number.parse(ee.String(wrsProp.split('/').get(1)));
  
  var col2010 = L5SR.filterDate('2009-07-01','2011-06-30')
                                .filterBounds(Biomas.filter(ee.Filter.eq('Bioma',biomList[biom])))
                                .filterMetadata('WRS_PATH', 'equals',WRS_PATH)
                                .filterMetadata('WRS_ROW', 'equals',WRS_ROW)
                                .filterMetadata('CLOUD_COVER','less_than',80)
                                .map(monthBand);
  
  var col2018 = L8SR.filterDate('2017-07-01','2019-06-30')
                                .filterBounds(Biomas.filter(ee.Filter.eq('Bioma',biomList[biom])))
                                .filterMetadata('WRS_PATH', 'equals',WRS_PATH )
                                .filterMetadata('WRS_ROW', 'equals',WRS_ROW )
                                .filterMetadata('CLOUD_COVER','less_than',80)
                                .map(monthBand);
  
  var monthList = ee.List.sequence(1, 12);
  
  var fillMonths = monthList.map(function(month){
    
    var fillMonth_2010 = col2010.filter(ee.Filter.calendarRange({start:month, field: 'month'}))
      .map(function(img){
        var mask = img.select(['pixel_qa']).lte(68)
        var fillNegative = img.normalizedDifference(['B4','B3']).gte(0)
        return img.normalizedDifference(['B4','B3'])
                  .mask(fillNegative)
                  .addBands(img.select('month_band'))
                  .mask(mask)
                  .copyProperties(img);
      })
      .reduce(ee.Reducer.mean())
      .rename([ee.String('L5_NDVI_2010'),ee.String('month_band_2010')])//.cat(ee.Number(month).toInt()))
      .set({Month: month,WRS_PATH: WRS_PATH,WRS_ROW: WRS_ROW});
      
    var fillMonth_2018 = col2018.filter(ee.Filter.calendarRange({start:month, field: 'month'}))
      .map(function(img){
        var mask = img.select(['pixel_qa']).lte(324);
        var fillNegative = img.normalizedDifference(['B4','B3']).gte(0)
        return img.normalizedDifference(['B5','B4'])
                  .mask(fillNegative)
                  .addBands(img.select('month_band'))
                  .mask(mask)
                  .copyProperties(img);
      })
      .reduce(ee.Reducer.mean())
      .rename([ee.String('L8_NDVI_2018'),ee.String('month_band_2018')])//.cat(ee.Number(month).toInt()))
      .set({Month: month,WRS_PATH: WRS_PATH,WRS_ROW: WRS_ROW});
    return ee.Image([fillMonth_2010,fillMonth_2018]);
  });
  
  fillMonths = ee.ImageCollection(fillMonths).map(function(img){
    var bandSize = ee.Number(img.bandNames().size());
    return img.clip(feature).set({BandNumber: bandSize});
  });
  
  fillMonths = fillMonths.filter(ee.Filter.equals('BandNumber',4));
  
  return fillMonths;
});

var mergeColBRNDVI = getBRNDVI.iterate(function(ImgCols,mainCol){
  return ee.ImageCollection(mainCol).merge(ee.ImageCollection(ImgCols));
},ee.ImageCollection([]));

print(mergeColBRNDVI)

var mosaicColBRNDVI = ee.List.sequence(1, 12).map(function(month){
  
  var bandName2010 = ee.String('L5_NDVI_2010');//.cat(ee.Number(month).toInt());
  var mosaic2010 = ee.ImageCollection(mergeColBRNDVI)
      .filter(ee.Filter.eq('Month',month));
  
  var bandName2018 = ee.String('L8_NDVI_2018');//.cat(ee.Number(month).toInt())
  var mosaic2018 = ee.ImageCollection(mergeColBRNDVI)
      .filter(ee.Filter.eq('Month',month));
  
  var month_pixel_mask_2010 = mosaic2010.select('month_band_2010').mosaic()
  var month_pixel_mask_2018 = mosaic2018.select('month_band_2018').mosaic()
  
  var month_pixel_mask_intersect = month_pixel_mask_2010.and(month_pixel_mask_2018)
  
  mosaic2010 = mosaic2010.select(bandName2010).mosaic()
  mosaic2018 = mosaic2018.select(bandName2018).mosaic()
  
  var mfilled_mosaic2010 = mosaic2010.updateMask(month_pixel_mask_intersect)
  var mfilled_mosaic2018 = mosaic2018.updateMask(month_pixel_mask_intersect)
  
  return mfilled_mosaic2010.addBands(mfilled_mosaic2018);
  //return mosaic2010.addBands(mosaic2018);
});

mosaicColBRNDVI = ee.ImageCollection(mosaicColBRNDVI);

//var teste_masks = mosaicColBRNDVI.select('month_band_2010')

//------------------------------------------------------------------------//
// mediana bianual corrigida equalizada mensalmente 2010
var Median2010 = ee.Image(mosaicColBRNDVI.select('L5_NDVI_2010')
                          .reduce(ee.Reducer.median())
                          .multiply(factorA)
                          .add(factorB))
                          .updateMask(ee.Image('users/vieiramesquita/PASTAGEM3/pasture_lapig_2010'))
                          .clip(brasil);

Export.image.toDrive({image: Median2010,
description:'pa_br_ndvi_med_binaual_cor_eqlmes_l8sr_2010'+biomList[biom],
folder:'00_pa_br_ndvi_bianual_eqlmes_cor',
fileNamePrefix: 'pa_br_ndvi_med_binaual_cor_eqlmes_l8sr_2010'+biomList[biom],
region: LANDSAT_GRID.geometry().bounds(),
scale: 30, 
maxPixels: 1e13
});

//------------------------------------------------------------------------//
// mediana bianual corrigida equalizada mensalmente 2018
var Median2018 = mosaicColBRNDVI.select('L8_NDVI_2018')
                                .reduce(ee.Reducer.median())
                                .updateMask(ee.Image('users/vieiramesquita/PASTAGEM3_1/pasture_lapig_2018'))
                                .clip(brasil);

Export.image.toDrive({image: Median2018,
description:'pa_br_ndvi_med_binaual_cor_eqlmes_l8sr_2018'+biomList[biom],
folder:'00_pa_br_ndvi_bianual_eqlmes_cor',
fileNamePrefix: 'pa_br_ndvi_med_binaual_cor_eqlmes_l8sr_2018'+biomList[biom],
region: LANDSAT_GRID.geometry().bounds(),
scale: 30, 
maxPixels: 1e13
});

var visParam = {min: 0.25, max: 1, palette: [
    'FFFFFF', 'CE7E45', 'DF923D', 'F1B555', 'FCD163', '99B718', '74A901',
    '66A000', '529400', '3E8601', '207401', '056201', '004C00', '023B01',
    '012E01', '011D01', '011301'
  ],};

Map.addLayer(Median2010,visParam,'2010-'+biomList[biom]);
Map.addLayer(Median2018,visParam,'2018-'+biomList[biom]);

}
