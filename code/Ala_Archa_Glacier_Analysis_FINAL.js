var MAIN_THRESHOLD = 0.40;
var SENSITIVITY_THRESHOLD_LOW = 0.30;
var SENSITIVITY_THRESHOLD_MID = 0.35;
var SENSITIVITY_THRESHOLD_HIGH = 0.45;
var MIN_VALID_OBSERVATIONS = 2;
var CLOUD_LIMIT = 30;
var SCALE = 30;

var STUDY_START = 1984;
var STUDY_END = 2025;


// ============================================================
// 1. STUDY AREA
// ============================================================

var aoi = ee.Geometry.Rectangle([
  74.40, 42.40,
  74.68, 42.65
]);

var basin10 = ee.FeatureCollection(
  'WWF/HydroSHEDS/v1/Basins/hybas_10'
);

var catchment = basin10
  .filter(ee.Filter.eq(
    'HYBAS_ID',
    4101277960
  ));

var catchmentGeometry = catchment.geometry();

var catchmentAreaKm2 = catchmentGeometry
  .area(1)
  .divide(1e6);

print('==================================================');
print('STUDY AREA CHECK');
print('Selected HYBAS_ID:', 4101277960);
print('HydroSHEDS level:', 'Level 10');
print('Ala-Archa catchment area (km2):', catchmentAreaKm2);
print('==================================================');

Map.centerObject(catchment, 10);

Map.addLayer(
  catchment.style({
    color: 'red',
    fillColor: '00000000',
    width: 3
  }),
  {},
  'Ala-Archa catchment'
);


// ============================================================
// 2. LANDSAT COLLECTIONS
// ============================================================

function getCollection(year) {

  year = ee.Number(year);

  var startDate = ee.Date.fromYMD(
    year,
    8,
    15
  );

  var endDate = ee.Date.fromYMD(
    year,
    9,
    16
  );

  var l5 = ee.ImageCollection(
    'LANDSAT/LT05/C02/T1_L2'
  )
  .filterBounds(catchmentGeometry)
  .filterDate(startDate, endDate)
  .filter(
    ee.Filter.lt(
      'CLOUD_COVER',
      CLOUD_LIMIT
    )
  );

  var l7 = ee.ImageCollection(
    'LANDSAT/LE07/C02/T1_L2'
  )
  .filterBounds(catchmentGeometry)
  .filterDate(startDate, endDate)
  .filter(
    ee.Filter.lt(
      'CLOUD_COVER',
      CLOUD_LIMIT
    )
  );

  var l8 = ee.ImageCollection(
    'LANDSAT/LC08/C02/T1_L2'
  )
  .filterBounds(catchmentGeometry)
  .filterDate(startDate, endDate)
  .filter(
    ee.Filter.lt(
      'CLOUD_COVER',
      CLOUD_LIMIT
    )
  );

  var l9 = ee.ImageCollection(
    'LANDSAT/LC09/C02/T1_L2'
  )
  .filterBounds(catchmentGeometry)
  .filterDate(startDate, endDate)
  .filter(
    ee.Filter.lt(
      'CLOUD_COVER',
      CLOUD_LIMIT
    )
  );

  var collection;

  collection = ee.ImageCollection(
    ee.Algorithms.If(
      year.lte(1998),
      l5,
      ee.Algorithms.If(
        year.lte(2012),
        l5.merge(l7),
        ee.Algorithms.If(
          year.lte(2020),
          l8,
          l8.merge(l9)
        )
      )
    )
  );

  return collection;
}


// ============================================================
// 3. LANDSAT SURFACE REFLECTANCE SCALING
// ============================================================

function applyScaleFactors(image) {

  var optical = image
    .select('SR_B.')
    .multiply(0.0000275)
    .add(-0.2);

  return image.addBands(
    optical,
    null,
    true
  );
}


// ============================================================
// 4. QA MASK + NDSI
// ============================================================

function imageToNDSI(image) {

  var scaled = applyScaleFactors(image);

  var qa = scaled.select('QA_PIXEL');

  var qaRadsat = scaled.select('QA_RADSAT');

  var fill = qa
    .bitwiseAnd(1 << 0)
    .eq(0);

  var dilatedCloud = qa
    .bitwiseAnd(1 << 1)
    .eq(0);

  var cloud = qa
    .bitwiseAnd(1 << 3)
    .eq(0);

  var cloudShadow = qa
    .bitwiseAnd(1 << 4)
    .eq(0);

  var oldSensor = ee.String(
    image.get('SPACECRAFT_ID')
  ).match(
    'LANDSAT_5|LANDSAT_7'
  );

  var oldMask = qaRadsat
    .bitwiseAnd(1 << 9)
    .eq(0);

  var newMask = qa
    .bitwiseAnd(1 << 2)
    .eq(0)
    .and(
      qaRadsat
        .bitwiseAnd(1 << 11)
        .eq(0)
    );

  var sensorMask = ee.Image(
    ee.Algorithms.If(
      oldSensor,
      oldMask,
      newMask
    )
  );

  var clean = scaled.updateMask(
    fill
      .and(dilatedCloud)
      .and(cloud)
      .and(cloudShadow)
      .and(sensorMask)
  );

  var spacecraft = ee.String(
    image.get('SPACECRAFT_ID')
  );

  var ndsi = ee.Image(
    ee.Algorithms.If(
      spacecraft.match('LANDSAT_5|LANDSAT_7'),
      clean.normalizedDifference([
        'SR_B2',
        'SR_B5'
      ]),
      clean.normalizedDifference([
        'SR_B3',
        'SR_B6'
      ])
    )
  );

  return ndsi
    .rename('NDSI')
    .copyProperties(
      image,
      [
        'system:time_start',
        'system:index',
        'CLOUD_COVER',
        'SPACECRAFT_ID'
      ]
    );
}


// ============================================================
// 5. AREA CALCULATION
// ============================================================

function calculateAreaKm2(mask) {

  var areaImage = mask
    .selfMask()
    .multiply(
      ee.Image.pixelArea()
    )
    .rename('area');

  var stats = areaImage.reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: catchmentGeometry,
    scale: SCALE,
    maxPixels: 1e13
  });

  var areaM2 = ee.Number(
    ee.Algorithms.If(
      stats.contains('area'),
      stats.get('area'),
      0
    )
  );

  return areaM2.divide(1e6);
}


// ============================================================
// 6. ANALYZE ONE YEAR
// ============================================================

function analyzeYear(year) {

  year = ee.Number(year);

  var collection = getCollection(year);

  var imageCount = collection.size();

  var emptyResult = ee.Feature(
    null,
    {
      year: year,
      area_km2: null,
      area_status: 'NO_DATA',
      image_count: 0,
      selected_date: null,
      scene_dates: null,
      selected_cloud_cover: null,
      sensor: null,
      quality: 'NO_DATA'
    }
  );

  var availableResult = ee.Feature(
    ee.Algorithms.If(
      imageCount.gt(0),
      analyzeAvailableYear(
        year,
        collection
      ),
      emptyResult
    )
  );

  return availableResult;
}


// ============================================================
// 7. ANALYZE AVAILABLE YEAR
// ============================================================

function analyzeAvailableYear(
  year,
  collection
) {

  var ndsiCollection = collection.map(
    imageToNDSI
  );

  var seasonalNDSI =
    ndsiCollection.median();

  var observationCount =
    ndsiCollection.count();

  var glacierMask =
    seasonalNDSI.gte(
      MAIN_THRESHOLD
    );

  var sufficientObservations =
    observationCount.gte(
      MIN_VALID_OBSERVATIONS
    );

  var strictGlacier =
    glacierMask.and(
      sufficientObservations
    );

  var imageCount =
    collection.size();

  var areaKm2 = ee.Algorithms.If(
    imageCount.gte(
      MIN_VALID_OBSERVATIONS
    ),
    calculateAreaKm2(
      strictGlacier
    ),
    null
  );

  var dates = ee.List(
    collection.aggregate_array(
      'system:time_start'
    )
  );

  var formattedDates = dates.map(
    function(time) {
      return ee.Date(time)
        .format('YYYY-MM-dd');
    }
  );

  var uniqueDates = ee.List(
    formattedDates
  )
  .distinct()
  .sort();

  var sceneDates =
    uniqueDates.join(', ');

  var leastCloudy =
    collection.sort(
      'CLOUD_COVER'
    );

  var firstImage = ee.Image(
    leastCloudy.first()
  );

  var selectedDate =
    ee.Date(
      firstImage.get(
        'system:time_start'
      )
    ).format('YYYY-MM-dd');

  var selectedCloud =
    firstImage.get(
      'CLOUD_COVER'
    );

  var sensors = ee.List(
    collection.aggregate_array(
      'SPACECRAFT_ID'
    )
  )
  .distinct()
  .sort()
  .join(', ');

  var quality = ee.String(
    ee.Algorithms.If(
      imageCount.eq(1),
      'SINGLE_SCENE',
      'MULTI_SCENE'
    )
  );

  var status = ee.String(
    ee.Algorithms.If(
      imageCount.gte(
        MIN_VALID_OBSERVATIONS
      ),
      'VALID',
      'INSUFFICIENT_SCENES'
    )
  );

  return ee.Feature(
    null,
    {
      year: year,
      area_km2: areaKm2,
      area_status: status,
      image_count: imageCount,
      selected_date: selectedDate,
      scene_dates: sceneDates,
      selected_cloud_cover: selectedCloud,
      sensor: sensors,
      quality: quality
    }
  );
}


// ============================================================
// 8. RUN 1984–2025
// ============================================================

var years = ee.List.sequence(
  STUDY_START,
  STUDY_END
);

var results = ee.FeatureCollection(
  years.map(analyzeYear)
);

print(
  'ALL YEARS — MAIN SEASONAL MEDIAN NDSI RESULTS:',
  results
);


// ============================================================
// 9. CLEAN RESULTS
// ============================================================

var cleanResults = results.map(
  function(feature) {

    return ee.Feature(
      null,
      {
        year: feature.get('year'),
        area_km2: feature.get('area_km2'),
        area_status: feature.get('area_status'),
        image_count: feature.get('image_count'),
        selected_date: feature.get('selected_date'),
        scene_dates: feature.get('scene_dates'),
        selected_cloud_cover:
          feature.get(
            'selected_cloud_cover'
          ),
        sensor: feature.get('sensor'),
        quality: feature.get('quality')
      }
    );

  }
);

print(
  'YEAR / AREA / IMAGE COUNT / DATES / SENSOR / QUALITY / STATUS:',
  cleanResults
);


// ============================================================
// 10. VALID AND MISSING YEARS
// ============================================================

var validAreaResults =
  cleanResults.filter(
    ee.Filter.notNull([
      'area_km2'
    ])
  );

var validYears =
  ee.List(
    validAreaResults.aggregate_array(
      'year'
    )
  );

var missingResults =
  cleanResults.filter(
    ee.Filter.eq(
      'area_status',
      'INSUFFICIENT_SCENES'
    )
  );

var missingYears =
  ee.List(
    missingResults.aggregate_array(
      'year'
    )
  );

print(
  'VALID GLACIER AREA YEARS:',
  validYears
);

print(
  'NUMBER OF VALID GLACIER AREA YEARS:',
  validYears.size()
);

print(
  'MISSING / UNAVAILABLE AREA YEARS:',
  missingYears
);

print(
  'MISSING / UNAVAILABLE YEAR DIAGNOSTICS:',
  missingResults
);


// ============================================================
// 11. TIME SERIES
// ============================================================

var chartResults =
  validAreaResults.sort(
    'year'
  );

var chart =
  ui.Chart.feature.byFeature(
    chartResults,
    'year',
    'area_km2'
  )
  .setChartType('LineChart')
  .setOptions({
    title:
      'Ala-Archa Glacier Area Over Time',
    hAxis: {
      title: 'Year'
    },
    vAxis: {
      title: 'Glacier Area (km2)'
    },
    pointSize: 4,
    lineWidth: 2,
    legend: {
      position: 'none'
    }
  });

print(chart);


// ============================================================
// 12. ANNUAL AREA CHANGE
// ============================================================

var sortedResults =
  validAreaResults.sort(
    'year'
  );

var sortedList =
  sortedResults.toList(
    sortedResults.size()
  );

var changeFeatures =
  ee.List.sequence(
    0,
    ee.Number(
      sortedResults.size()
    ).subtract(1)
  )
  .map(
    function(i) {

      i = ee.Number(i);

      var current =
        ee.Feature(
          sortedList.get(i)
        );

      var previous =
        ee.Feature(
          sortedList.get(
            i.subtract(1)
          )
        );

      var currentYear =
        ee.Number(
          current.get('year')
        );

      var previousYear =
        ee.Number(
          previous.get('year')
        );

      var consecutive =
        currentYear.subtract(
          previousYear
        ).eq(1);

      return ee.Feature(
        current
      ).set(
        'previous_year',
        previousYear
      )
      .set(
        'previous_year_area_km2',
        previous.get('area_km2')
      )
      .set(
        'area_change_km2',
        ee.Algorithms.If(
          consecutive,
          ee.Number(
            current.get(
              'area_km2'
            )
          ).subtract(
            ee.Number(
              previous.get(
                'area_km2'
              )
            )
          ),
          null
        )
      )
      .set(
        'area_change_status',
        ee.Algorithms.If(
          consecutive,
          'CONSECUTIVE_YEAR',
          'GAP_IN_SERIES'
        )
      );

    }
  );

var annualChange =
  ee.FeatureCollection(
    changeFeatures
  )
  .filter(
    ee.Filter.notNull([
      'area_change_km2'
    ])
  );

print(
  'ANNUAL GLACIER AREA CHANGE DATA:',
  annualChange
);

print(
  'NUMBER OF CONSECUTIVE-YEAR AREA-CHANGE OBSERVATIONS:',
  annualChange.size()
);

print(
  'YEARS WITH ANNUAL AREA CHANGE:',
  annualChange.aggregate_array(
    'year'
  )
);


// ============================================================
// 13. RGI v7
// ============================================================

var rgi =
  ee.FeatureCollection(
    'projects/intense-reactor-455007-a9/assets/RGI2000-v7-0-G-13_central_asia'
  );

var rgiStudyArea =
  rgi.filterBounds(
    catchmentGeometry
  );

print(
  'Number of RGI glaciers in study area:',
  rgiStudyArea.size()
);

var rgiWithArea =
  rgiStudyArea.map(
    function(feature) {

      var clipped =
        feature.geometry()
          .intersection(
            catchmentGeometry,
            1
          );

      return ee.Feature(
        clipped
      ).set(
        'rgi_area_km2',
        clipped.area(1)
          .divide(1e6)
      );

    }
  );

var totalRGIArea =
  ee.Number(
    rgiWithArea.aggregate_sum(
      'rgi_area_km2'
    )
  );

var meanRGIArea =
  totalRGIArea.divide(
    rgiStudyArea.size()
  );

print(
  'Total RGI glacier area inside catchment (km2):',
  totalRGIArea
);

print(
  'Mean RGI glacier area (km2):',
  meanRGIArea
);

// RGI v7 stores the outline source date in src_date.
var rgiDated =
  rgiStudyArea
    .filter(
      ee.Filter.notNull([
        'src_date'
      ])
    )
    .map(function(feature) {
      return feature.set(
        'src_year',
        ee.Date(
          feature.get('src_date')
        ).get('year')
      );
    });

print(
  'RGI glaciers with source dates:',
  rgiDated.size()
);

print(
  'RGI glaciers without source dates:',
  rgiStudyArea.size()
    .subtract(rgiDated.size())
);

print(
  'RGI source date — earliest:',
  rgiDated.aggregate_min('src_date')
);

print(
  'RGI source date — latest:',
  rgiDated.aggregate_max('src_date')
);

print(
  'RGI source year — mean:',
  rgiDated.aggregate_mean('src_year')
);

print(
  'RGI source year — minimum:',
  rgiDated.aggregate_min('src_year')
);

print(
  'RGI source year — maximum:',
  rgiDated.aggregate_max('src_year')
);

Map.addLayer(
  rgiStudyArea.style({
    color: 'red',
    fillColor: '00000000',
    width: 2
  }),
  {},
  'RGI glacier outlines'
);


// ============================================================
// 14. 2020 MAIN ANALYSIS
// ============================================================

var collection2020 =
  getCollection(2020);

var ndsi2020Collection =
  collection2020.map(
    imageToNDSI
  );

var ndsi2020 =
  ndsi2020Collection.median();

var observationCount2020 =
  ndsi2020Collection.count();

var sufficient2020 =
  observationCount2020.gte(
    MIN_VALID_OBSERVATIONS
  );

var glacier030 =
  ndsi2020
    .gte(SENSITIVITY_THRESHOLD_LOW)
    .and(sufficient2020);

var glacier035 =
  ndsi2020
    .gte(SENSITIVITY_THRESHOLD_MID)
    .and(sufficient2020);

var glacier040 =
  ndsi2020
    .gte(0.40)
    .and(sufficient2020);

var glacier045 =
  ndsi2020
    .gte(0.45)
    .and(sufficient2020);

var area030 =
  calculateAreaKm2(
    glacier030
  );

var area035 =
  calculateAreaKm2(
    glacier035
  );

var area040 =
  calculateAreaKm2(
    glacier040
  );

var area045 =
  calculateAreaKm2(
    glacier045
  );

var thresholdSpread =
  area030.subtract(
    area045
  );

var validPixelStats2020 =
  sufficient2020.selfMask().reduceRegion({
    reducer: ee.Reducer.count(),
    geometry: catchmentGeometry,
    scale: SCALE,
    maxPixels: 1e13
  });

var validPixelCount2020 =
  ee.Number(
    ee.Algorithms.If(
      validPixelStats2020.contains('NDSI'),
      validPixelStats2020.get('NDSI'),
      0
    )
  );

var observationStats2020 =
  observationCount2020
    .updateMask(sufficient2020)
    .reduceRegion({
      reducer: ee.Reducer.mean()
        .combine({
          reducer2: ee.Reducer.max(),
          sharedInputs: true
        }),
      geometry: catchmentGeometry,
      scale: SCALE,
      maxPixels: 1e13
    });

var differenceRGI =
  area040.subtract(
    totalRGIArea
  );

var percentDifferenceRGI =
  differenceRGI
    .abs()
    .divide(
      totalRGIArea
    )
    .multiply(100);

print('==================================================');
print('2020 MAIN WINDOW');
print(
  '2020 suitable scenes:',
  collection2020.size()
);

print(
  '2020 glacier area NDSI 0.30 (km2):',
  area030
);

print(
  '2020 glacier area NDSI 0.35 (km2):',
  area035
);

print(
  '2020 glacier area NDSI 0.40 (km2):',
  area040
);

print(
  '2020 glacier area NDSI 0.45 (km2):',
  area045
);

print(
  '2020 threshold sensitivity range (km2):',
  thresholdSpread
);

print(
  '2020 pixels with >= 2 valid observations:',
  validPixelCount2020
);

print(
  '2020 mean valid observations per pixel (among >=2):',
  observationStats2020.get('NDSI_mean')
);

print(
  '2020 maximum valid observations per pixel:',
  observationStats2020.get('NDSI_max')
);

print(
  '2020 Landsat seasonal-composite glacier area (km2):',
  area040
);

print(
  '2020 RGI glacier area inside catchment (km2):',
  totalRGIArea
);

print(
  '2020 Landsat vs RGI difference (km2):',
  differenceRGI
);

print(
  '2020 Landsat vs RGI difference (%):',
  percentDifferenceRGI
);

Map.addLayer(
  ndsi2020.clip(
    catchmentGeometry
  ),
  {
    min: -1,
    max: 1,
    palette: [
      'brown',
      'white',
      'blue'
    ]
  },
  '2020 Seasonal Median NDSI'
);

Map.addLayer(
  glacier040.selfMask()
    .clip(catchmentGeometry),
  {
    palette: ['cyan']
  },
  '2020 Glacier Mask NDSI 0.40'
);

// 2020 true-colour diagnostic for visual assessment of omitted area.
var trueColour2020 =
  collection2020
    .map(applyScaleFactors)
    .median()
    .clip(catchmentGeometry);

Map.addLayer(
  trueColour2020,
  {
    bands: ['SR_B4', 'SR_B3', 'SR_B2'],
    min: 0.02,
    max: 0.35
  },
  '2020 True Colour'
);

Map.addLayer(
  rgiStudyArea.style({
    color: 'yellow',
    fillColor: '00000000',
    width: 2
  }),
  {},
  '2020 RGI Outlines'
);

// Combined visual diagnostic: Landsat mask over true-colour imagery.
Map.addLayer(
  glacier040.selfMask()
    .clip(catchmentGeometry),
  {
    palette: ['00FFFF']
  },
  '2020 Landsat Mask Overlay'
);


// ============================================================
// 15. 2020 LATE-SEASON DIAGNOSTIC
// ============================================================

var lateStart =
  ee.Date('2020-09-01');

var lateEnd =
  ee.Date('2020-09-16');

var lateCollection =
  ee.ImageCollection(
    'LANDSAT/LC08/C02/T1_L2'
  )
  .filterBounds(
    catchmentGeometry
  )
  .filterDate(
    lateStart,
    lateEnd
  )
  .filter(
    ee.Filter.lt(
      'CLOUD_COVER',
      CLOUD_LIMIT
    )
  );

var lateNDSICollection =
  lateCollection.map(
    imageToNDSI
  );

var lateCount =
  lateNDSICollection.size();

var lateArea = ee.Number(
  ee.Algorithms.If(
    lateCount.gte(
      MIN_VALID_OBSERVATIONS
    ),
    calculateAreaKm2(
      lateNDSICollection
        .median()
        .gte(MAIN_THRESHOLD)
        .and(
          lateNDSICollection.count()
            .gte(
              MIN_VALID_OBSERVATIONS
            )
        )
    ),
    0
  )
);

var lateSufficient =
  lateNDSICollection.count()
    .gte(MIN_VALID_OBSERVATIONS);

var lateValidPixelStats =
  lateSufficient.selfMask()
    .reduceRegion({
      reducer: ee.Reducer.count(),
      geometry: catchmentGeometry,
      scale: SCALE,
      maxPixels: 1e13
    });

var lateValidPixelCount =
  ee.Number(
    ee.Algorithms.If(
      lateValidPixelStats.contains('NDSI'),
      lateValidPixelStats.get('NDSI'),
      0
    )
  );

var lateObservationStats =
  lateNDSICollection.count()
    .updateMask(lateSufficient)
    .reduceRegion({
      reducer: ee.Reducer.mean()
        .combine({
          reducer2: ee.Reducer.max(),
          sharedInputs: true
        }),
      geometry: catchmentGeometry,
      scale: SCALE,
      maxPixels: 1e13
    });

print('==================================================');
print('2020 LATE-SEASON DIAGNOSTIC');
print(
  'Late-season window:',
  '1 September – 15 September'
);

print(
  'Late-season suitable scenes:',
  lateCount
);

print(
  'Late-season pixels with >= 2 valid observations:',
  lateValidPixelCount
);

print(
  'Late-season mean valid observations per pixel (among >=2):',
  lateObservationStats.get('NDSI_mean')
);

print(
  'Late-season maximum valid observations per pixel:',
  lateObservationStats.get('NDSI_max')
);

print(
  'Late-season 2020 glacier area NDSI 0.40 (km2):',
  lateArea
);

print(
  'Main-window 2020 glacier area NDSI 0.40 (km2):',
  area040
);

print(
  'Late-season minus main-window area (km2):',
  lateArea.subtract(area040)
);

print(
  'Late-season vs RGI difference (km2):',
  lateArea.subtract(totalRGIArea)
);

print(
  'Late-season vs RGI difference (%):',
  lateArea.subtract(totalRGIArea)
    .abs()
    .divide(totalRGIArea)
    .multiply(100)
);


// ============================================================
// 16. 1993 INDEPENDENT JULY VALIDATION
// ============================================================

var validationCol =
  ee.ImageCollection(
    'LANDSAT/LT05/C02/T1_L2'
  )
  .filterBounds(
    catchmentGeometry
  )
  .filterDate(
    '1993-07-01',
    '1993-07-11'
  )
  .filter(
    ee.Filter.lt(
      'CLOUD_COVER',
      CLOUD_LIMIT
    )
  );

var validationCount =
  validationCol.size();

var validationNDSICollection =
  validationCol.map(
    imageToNDSI
  );

var validationArea =
  ee.Number(
    ee.Algorithms.If(
      validationCount.gt(0),
      calculateAreaKm2(
        validationNDSICollection
          .median()
          .gte(MAIN_THRESHOLD)
          .and(
            validationNDSICollection
              .count()
              .gte(
                MIN_VALID_OBSERVATIONS
              )
          )
      ),
      0
    )
  );

print('==================================================');
print('1993 VALIDATION DIAGNOSTIC');

print(
  '1993 validation scenes:',
  validationCount
);

print(
  '1993 validation area NDSI 0.40 (km2):',
  validationArea
);

print(
  '1993 validation is an independent July diagnostic and is not treated as a direct comparison with the late-season annual series.'
);


// ============================================================
// 17. 2009 QUALITY CONTROL
// ============================================================

var collection2009 =
  getCollection(2009);

var ndsi2009Collection =
  collection2009.map(
    imageToNDSI
  );

var glacier2009 =
  ndsi2009Collection
    .median()
    .gte(MAIN_THRESHOLD)
    .and(
      ndsi2009Collection
        .count()
        .gte(
          MIN_VALID_OBSERVATIONS
        )
    );

var area2009 =
  calculateAreaKm2(
    glacier2009
  );

var image2009 =
  ee.Image(
    collection2009
      .sort('CLOUD_COVER')
      .first()
  );

print('==================================================');
print('2009 QUALITY CONTROL');

print(
  '2009 suitable scenes:',
  collection2009.size()
);

print(
  '2009 strict glacier area NDSI 0.40 (km2):',
  area2009
);

print(
  '2009 least-cloudy scene date:',
  ee.Date(
    image2009.get(
      'system:time_start'
    )
  ).format('YYYY-MM-dd')
);


// ============================================================
// 18. TOPOGRAPHY
// ============================================================

var dem =
  ee.Image(
    'USGS/SRTMGL1_003'
  );

var elevationStats =
  dem.reduceRegion({
    reducer:
      ee.Reducer.minMax(),
    geometry:
      catchmentGeometry,
    scale:
      30,
    maxPixels:
      1e13
  });

print('==================================================');
print('STUDY AREA TOPOGRAPHY');

print(
  'Minimum elevation (m):',
  elevationStats.get(
    'elevation_min'
  )
);

print(
  'Maximum elevation (m):',
  elevationStats.get(
    'elevation_max'
  )
);


// ============================================================
// 19. ERA5 SUMMER TEMPERATURE — FINAL CORRELATION DATA
// ============================================================

var era5 =
  ee.ImageCollection(
    'ECMWF/ERA5/HOURLY'
  )
  .select(
    'temperature_2m'
  );

var correlationData =
  annualChange
    .filter(
      ee.Filter.gte(
        'year',
        2002
      )
    )
    .filter(
      ee.Filter.lte(
        'year',
        2025
      )
    )
    .map(function(feature) {

      var year =
        ee.Number(
          feature.get('year')
        );

      var start =
        ee.Date.fromYMD(
          year,
          6,
          1
        );

      var end =
        ee.Date.fromYMD(
          year,
          9,
          1
        );

      var summerTemperature =
        era5
          .filterDate(
            start,
            end
          )
          .mean()
          .subtract(
            273.15
          );

      var temperatureStats =
        summerTemperature.reduceRegion({
          reducer:
            ee.Reducer.mean(),
          geometry:
            catchmentGeometry,
          scale:
            27830,
          maxPixels:
            1e13
        });

      return feature.set(
        'summer_temp_c',
        temperatureStats.get(
          'temperature_2m'
        )
      );

    })
    .filter(
      ee.Filter.notNull([
        'area_change_km2',
        'summer_temp_c'
      ])
    );

print(
  '=================================================='
);

print(
  'FINAL CORRELATION DATA — 2002–2025:',
  correlationData
);

print(
  'FINAL CORRELATION N:',
  correlationData.size()
);

print(
  'FINAL CORRELATION YEARS:',
  correlationData.aggregate_array(
    'year'
  )
);
// ============================================================
// 20. PEARSON
// ============================================================

var pearson =
  correlationData.reduceColumns({
    reducer:
      ee.Reducer.pearsonsCorrelation(),
    selectors: [
      'area_change_km2',
      'summer_temp_c'
    ]
  });

print(
  '=================================================='
);

print(
  'PEARSON CORRELATION — AREA CHANGE VS SUMMER TEMPERATURE:',
  pearson
);


// ============================================================
// 21. SPEARMAN
// ============================================================

var spearman =
  correlationData.reduceColumns({
    reducer:
      ee.Reducer.spearmansCorrelation(),
    selectors: [
      'area_change_km2',
      'summer_temp_c'
    ]
  });

print(
  'SPEARMAN CORRELATION — AREA CHANGE VS SUMMER TEMPERATURE:',
  spearman
);


// ============================================================
// 22. ABSOLUTE AREA DIAGNOSTIC
// ============================================================

var absoluteAreaData =
  validAreaResults
    .filter(
      ee.Filter.gte(
        'year',
        2002
      )
    )
    .filter(
      ee.Filter.lte(
        'year',
        2025
      )
    )
    .map(
    function(feature) {

      var year =
        ee.Number(
          feature.get('year')
        );

      var start =
        ee.Date.fromYMD(
          year,
          6,
          1
        );

      var end =
        ee.Date.fromYMD(
          year,
          9,
          1
        );

      var summerTemperature =
        era5
          .filterDate(
            start,
            end
          )
          .mean()
          .subtract(
            273.15
          );

      var temperatureStats =
        summerTemperature.reduceRegion({
          reducer:
            ee.Reducer.mean(),
          geometry:
            catchmentGeometry,
          scale:
            27830,
          maxPixels:
            1e13
        });

      return feature.set(
        'summer_temp_c',
        temperatureStats.get(
          'temperature_2m'
        )
      );

    }
  )
  .filter(
    ee.Filter.notNull([
      'area_km2',
      'summer_temp_c'
    ])
  );

var oldPearson =
  absoluteAreaData.reduceColumns({
    reducer:
      ee.Reducer.pearsonsCorrelation(),
    selectors: [
      'area_km2',
      'summer_temp_c'
    ]
  });

var oldSpearman =
  absoluteAreaData.reduceColumns({
    reducer:
      ee.Reducer.spearmansCorrelation(),
    selectors: [
      'area_km2',
      'summer_temp_c'
    ]
  });

print('==================================================');
print(
  'FINAL DIAGNOSTIC — ABSOLUTE AREA VS SUMMER TEMPERATURE (2002–2025)'
);

print(
  'Absolute-area Pearson:',
  oldPearson
);

print(
  'Absolute-area Spearman:',
  oldSpearman
);

print(
  'ABSOLUTE-AREA DIAGNOSTIC N:',
  absoluteAreaData.size()
);

print(
  'ABSOLUTE-AREA DIAGNOSTIC YEARS:',
  absoluteAreaData.aggregate_array('year')
);


// ============================================================
// 23. FINAL SUSPICIOUS-YEAR DIAGNOSTIC
// ============================================================

var suspiciousYears =
  ee.List([
    1998,
    2004,
    2006,
    2013,
    2021
  ]);

var diagnosticFeatures =
  suspiciousYears.map(
    function(year) {

      year = ee.Number(year);

      var collection =
        getCollection(year);

      var count =
        collection.size();

      var ndsiCollection =
        collection.map(
          imageToNDSI
        );

      var ndsiCount =
        ndsiCollection.size();

      var hasData =
        ndsiCount.gt(0);

      var ndsi =
        ee.Image(
          ee.Algorithms.If(
            hasData,
            ndsiCollection.median(),
            ee.Image.constant(
              -999
            ).rename('NDSI')
          )
        );

      var observations =
        ee.Image(
          ee.Algorithms.If(
            hasData,
            ndsiCollection.count(),
            ee.Image.constant(
              0
            ).rename('NDSI')
          )
        );

      var sufficient =
        observations.gte(
          MIN_VALID_OBSERVATIONS
        );

      var area030 =
        ee.Number(
          ee.Algorithms.If(
            hasData,
            calculateAreaKm2(
              ndsi.gte(0.30)
                .and(sufficient)
            ),
            0
          )
        );

      var area035 =
        ee.Number(
          ee.Algorithms.If(
            hasData,
            calculateAreaKm2(
              ndsi.gte(0.35)
                .and(sufficient)
            ),
            0
          )
        );

      var area040 =
        ee.Number(
          ee.Algorithms.If(
            hasData,
            calculateAreaKm2(
              ndsi.gte(0.40)
                .and(sufficient)
            ),
            0
          )
        );

      var area045 =
        ee.Number(
          ee.Algorithms.If(
            hasData,
            calculateAreaKm2(
              ndsi.gte(0.45)
                .and(sufficient)
            ),
            0
          )
        );

      var dates =
        ee.List(
          collection.aggregate_array(
            'system:time_start'
          )
        )
        .map(
          function(time) {
            return ee.Date(time)
              .format(
                'YYYY-MM-dd'
              );
          }
        )
        .join(', ');

      var clouds =
        ee.List(
          collection.aggregate_array(
            'CLOUD_COVER'
          )
        )
        .join(', ');

      var sensors =
        ee.List(
          collection.aggregate_array(
            'SPACECRAFT_ID'
          )
        )
        .distinct()
        .join(', ');

      return ee.Feature(
        null,
        {
          year: year,
          image_count: count,
          scene_dates: dates,
          cloud_cover: clouds,
          sensors: sensors,
          area_030_km2: area030,
          area_035_km2: area035,
          area_040_km2: area040,
          area_045_km2: area045
        }
      );

    }
  );

var suspiciousDiagnostic =
  ee.FeatureCollection(
    diagnosticFeatures
  );

print(
  '=================================================='
);

print(
  'FINAL DIAGNOSTIC — SUSPICIOUS YEARS:',
  suspiciousDiagnostic
);


// ============================================================
// 24. FINAL SUMMARY
// ============================================================

print('==================================================');
print('FINAL SUMMARY');

print(
  'Study period:',
  '1984–2025'
);

print(
  'Selected catchment HYBAS_ID:',
  4101277960
);

print(
  'HydroSHEDS level:',
  'Level 10'
);

print(
  'Catchment area (km2):',
  catchmentAreaKm2
);

print(
  'RGI glacier count:',
  rgiStudyArea.size()
);

print(
  'RGI glacier area inside catchment (km2):',
  totalRGIArea
);

print(
  'Mean RGI glacier area (km2):',
  meanRGIArea
);

print(
  'Valid glacier-area years:',
  validYears
);

print(
  'Missing glacier-area years:',
  missingYears
);

print(
  'Consecutive annual-change observations:',
  annualChange.size()
);

print(
  'Correlation observations:',
  correlationData.size()
);

print(
  'Primary Pearson:',
  pearson
);

print(
  'Primary Spearman:',
  spearman
);

print(
  'Absolute-area Pearson (2002–2025):',
  oldPearson
);

print(
  'Absolute-area Spearman (2002–2025):',
  oldSpearman
);

print(
  '2020 Landsat area NDSI 0.30 (km2):',
  area030
);

print(
  '2020 Landsat area NDSI 0.40 (km2):',
  area040
);

print(
  '2020 RGI area (km2):',
  totalRGIArea
);

print(
  '2020 Landsat minus RGI (km2):',
  differenceRGI
);

print(
  '2020 threshold sensitivity range (km2):',
  thresholdSpread
);

print(
  '2020 pixels with >= 2 valid observations:',
  validPixelCount2020
);

print(
  'RGI source-date earliest:',
  rgiDated.aggregate_min('src_date')
);

print(
  'RGI source-date latest:',
  rgiDated.aggregate_max('src_date')
);

print(
  'RGI source-year mean:',
  rgiDated.aggregate_mean('src_year')
);

print(
  '2020 late-season area NDSI 0.40 (km2):',
  lateArea
);

print(
  '1993 validation area (km2):',
  validationArea
);

print(
  '2009 QC area NDSI 0.40 (km2):',
  area2009
);

print(
  'Main annual analysis window:',
  '15 August – 15 September'
);

print(
  'Late-season diagnostic window:',
  '1 September – 15 September'
);

print(
  'Composite method:',
  'Pixel-wise median of masked NDSI scenes'
);

print(
  'Minimum valid observations per pixel:',
  MIN_VALID_OBSERVATIONS
);

print(
  'Main glacier threshold:',
  'NDSI >= 0.40'
);

print(
  'Sensitivity thresholds:',
  'NDSI 0.30 / 0.35 / 0.40 / 0.45'
);

print(
  'Spatial resolution:',
  '30 m'
);

print(
  'RGI dataset:',
  'RGI v7 Central Asia'
);

print(
  'ERA5 variable:',
  'temperature_2m'
);

print(
  'ERA5 temporal window:',
  'June–August'
);

print(
  'ERA5 analysis units:',
  'Celsius'
);

print(
  'Primary correlation:',
  'Annual glacier area change vs June–August temperature'
);

print('==================================================');


// ============================================================
// 25. EXPORT ANNUAL RESULTS
// ============================================================

Export.table.toDrive({

  collection:
    cleanResults,

  description:
    'Ala_Archa_Glacier_Area_1984_2025_FINAL',

  fileNamePrefix:
    'Ala_Archa_Glacier_Area_1984_2025_FINAL',

  fileFormat:
    'CSV',

  selectors: [
    'year',
    'area_km2',
    'area_status',
    'image_count',
    'selected_date',
    'scene_dates',
    'selected_cloud_cover',
    'sensor',
    'quality'
  ]

});


// ============================================================
// 26. EXPORT CORRELATION DATA
// ============================================================

Export.table.toDrive({

  collection:
    correlationData,

  description:
    'Ala_Archa_Glacier_Area_ERA5_FINAL',

  fileNamePrefix:
    'Ala_Archa_Glacier_Area_ERA5_FINAL',

  fileFormat:
    'CSV',

  selectors: [
    'year',
    'previous_year',
    'previous_year_area_km2',
    'area_change_km2',
    'area_change_status',
    'summer_temp_c'
  ]

});


// ============================================================
// 26B. EXPORT ABSOLUTE-AREA DIAGNOSTIC
// ============================================================

Export.table.toDrive({

  collection:
    absoluteAreaData,

  description:
    'Ala_Archa_Absolute_Area_ERA5_DIAGNOSTIC',

  fileNamePrefix:
    'Ala_Archa_Absolute_Area_ERA5_DIAGNOSTIC',

  fileFormat:
    'CSV',

  selectors: [
    'year',
    'area_km2',
    'summer_temp_c'
  ]

});


// ============================================================
// 27. CARTOGRAPHIC ELEMENTS
// ============================================================

var northArrow = ui.Panel({
  widgets: [
    ui.Label(
      'N',
      {
        fontSize: '18px',
        fontWeight: 'bold',
        textAlign: 'center',
        margin: '0 0 2px 0'
      }
    ),
    ui.Label(
      '↑',
      {
        fontSize: '28px',
        fontWeight: 'bold',
        textAlign: 'center',
        margin: '0'
      }
    )
  ],
  style: {
    position: 'top-left',
    padding: '5px',
    backgroundColor:
      'rgba(255,255,255,0.85)'
  }
});

Map.add(northArrow);


// ============================================================
// 28. KYRGYZSTAN INSET
// ============================================================

var countries =
  ee.FeatureCollection(
    'USDOS/LSIB_SIMPLE/2017'
  );

var kyrgyzstan =
  countries.filter(
    ee.Filter.eq(
      'country_na',
      'Kyrgyzstan'
    )
  );

var alaArchaPoint =
  ee.Feature(
    catchmentGeometry.centroid()
  );

var insetImage =
  ee.Image(0)
    .paint(
      kyrgyzstan,
      1,
      2
    )
    .paint(
      alaArchaPoint,
      2,
      8
    );

var inset =
  ui.Thumbnail({
    image: insetImage,
    params: {
      region:
        kyrgyzstan.geometry(),
      dimensions: 180,
      format: 'png',
      min: 0,
      max: 2,
      palette: [
        'ffffff',
        '000000',
        'ff0000'
      ]
    },
    style: {
      width: '180px',
      height: '140px',
      margin: '0'
    }
  });

var insetPanel =
  ui.Panel({
    widgets: [
      ui.Label(
        'Ala-Archa study area in Kyrgyzstan',
        {
          fontWeight: 'bold',
          fontSize: '11px',
          margin: '0 0 4px 0'
        }
      ),
      inset
    ],
    style: {
      position: 'bottom-left',
      padding: '5px',
      backgroundColor:
        'rgba(255,255,255,0.85)'
    }
  });

Map.add(insetPanel);
