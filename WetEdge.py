"""
This is a conversion of a Google Earth Engine JavaScript to a Python script (convert by me and ChatGPT), the original preamble is below:

This code was used in the manuscript "Hidden becomes clear: optical remote sensing of vegetation reveals water table dynamics in northern peatlands."
submitted to Remote Sensing of Environment
Here the intercept and slop of the wet edge are calculated for Männikjärve (EE_MAN) peatland located in Estonia. 

The code authors: Iuliia Burdun, Viacheslav Komisarenko
"""


import ee
ee.Authenticate()
ee.Initialize()

# CLOUD MASKING AND FILTERING
s2Sr = ee.ImageCollection('COPERNICUS/S2_SR')
s2Clouds = ee.ImageCollection('COPERNICUS/S2_CLOUD_PROBABILITY')

# Set start and end dates, cloud probability, and area of interest
START_DATE = ee.Date('2018-03-01')
END_DATE = ee.Date('2021-12-30')
MAX_CLOUD_PROBABILITY = 5
table = ee.FeatureCollection('users/iuliiaburdun/OPTRAM_FI/EE_MAN_4326')
region = ee.FeatureCollection(table)  # Assuming 'table' is a variable or constant containing the region

# Define maskClouds function to get the images with the lowest cloud probability
def maskClouds(img):
    clouds = ee.Image(img.get('cloud_mask')).select('probability')
    isNotCloud = clouds.lt(MAX_CLOUD_PROBABILITY)
    return img.updateMask(isNotCloud)

# The masks for the 10m bands sometimes do not exclude bad data at
# scene edges, so we apply masks from the 20m (B8A) and 60m  (B9) bands as well.
# More information: https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_S2_CLOUD_PROBABILITY
def maskEdges(s2_img):
    return s2_img.updateMask(
        s2_img.select('B8A').mask().updateMask(s2_img.select('B9').mask()))

# Filter s2Sr and s2Clouds collections by START_DATE, END_DATE, and region.
criteria = ee.Filter.And(
    ee.Filter.bounds(region),
    ee.Filter.date(START_DATE, END_DATE))

s2Sr = (s2Sr
        .filter(criteria)
        .map(maskEdges))
s2Clouds = s2Clouds.filter(criteria)

# Join s2Sr with s2Clouds dataset to add a cloud mask.
s2SrWithCloudMask = (ee.Join.saveFirst('cloud_mask')
                     .apply(primary=s2Sr,
                            secondary=s2Clouds,
                            condition=ee.Filter.equals(leftField='system:index', rightField='system:index')))

# Apply the cloud mask and filter the collection based on the DOY
s2CloudMasked2 = (ee.ImageCollection(s2SrWithCloudMask)
                  .map(maskClouds)
                  .filter(ee.Filter.dayOfYear(90, 270)))


# CLIP DATA AND CALCULATE VEGETATION INDICES (VIs)

# Clip data to polygon
def clip_image(image):
    clip = image.clip(table)
    return clip

S2_boa = s2CloudMasked2.map(clip_image).select(['B2', 'B3', 'B5', 'B4', 'B6', 'B8', 'B12'])  # Bands needed for VIs

# Apply scale factor to the bands
def scale_bands(image):
    Blue = image.select(['B2']).rename('Blue').divide(10000)
    Green = image.select(['B3']).rename('Green').divide(10000)
    Red = image.select(['B4']).rename('Red').divide(10000)
    RE = image.select(['B5']).rename('RE').divide(10000)
    RE2 = image.select(['B6']).rename('RE2').divide(10000)
    NIR = image.select(['B8']).rename('NIR').divide(10000)
    SWIR = image.select(['B12']).rename('SWIR').divide(10000)
    new_bands = image.addBands([Blue, Red, RE, RE2, NIR, SWIR, Green])
    return new_bands.select(['Blue', 'Red', 'RE', 'RE2', 'NIR', 'SWIR', 'Green'])

S2_boa = S2_boa.map(scale_bands)

# Function to calculate NDVI
def calculate_ndvi(image):
    ndvi = image.normalizedDifference(['NIR', 'Red']).rename('NDVI')
    return image.addBands(ndvi)

# Function to calculate RENDVI
def calculate_rendvi(image):
    rendvi = image.normalizedDifference(['RE2', 'RE']).rename('RENDVI')
    return image.addBands(rendvi)

# Function to calculate kNDVI
def calculate_kndvi(image):
    kndvi2 = image.expression(
        '(tanh((NIR-Red)/(NIR+Red))**2)', {
            'NIR': image.select('NIR'),
            'Red': image.select('Red')
        })
    kndvi = kndvi2.rename('kNDVI')
    return image.addBands(kndvi)

# Function to calculate EVI
def calculate_evi(image):
    evi2 = image.expression(
        '(2.5*(B5-B4)/(B5+6*B4-7.5*B2+1))', {
            'B2': image.select('Blue'),
            'B4': image.select('Red'),
            'B5': image.select('NIR')
        })
    evi = evi2.rename('EVI')
    return image.addBands(evi)

# Function to calculate MNDWI
def calculate_mndwi(image):
    mndwi = image.normalizedDifference(['Green', 'NIR']).rename('MNDWI')
    return image.addBands(mndwi)

# Function to calculate STR
def calculate_str(image):
    str2 = image.expression(
        '((1-SWIR)*(1-SWIR))/(2*SWIR)', {
            'SWIR': image.select('SWIR')
        })
    str_band = str2.rename('STR')
    return image.addBands(str_band)

S2_boa = (S2_boa
          .map(calculate_ndvi)
          .map(calculate_rendvi)
          .map(calculate_kndvi)
          .map(calculate_evi)
          .map(calculate_mndwi)
          .map(calculate_str))

# Function to calculate the number of NDVI pixels for each image
def calculate_ndvi_pixels(image):
    count = image.gt(-1).reduceRegion(
        reducer=ee.Reducer.sum(),
        geometry=table,
        maxPixels=1e13
    ).values().get(0)
    return image.set({'NDVI_pixel_count': count})

S2_boa = S2_boa.map(calculate_ndvi_pixels)

# Keep images that have at least 1 valid NDVI pixel
S2_boa = S2_boa.filter(ee.Filter.gt('NDVI_pixel_count', 1))

"""
# Based on VIs and STR we apply additional filters:
1) NDVI >= 0 - to exclude water bodies and bare land
2) MNDWI <= -0.2 - to exclude water bodies
3) STR [0,20] and EVI [-1, 1] - to exclude the pixel affected by shadows
"""

def apply_additional_filters(image):
    full_cover = image.select('NDVI').gte(0)
    no_water = image.select('MNDWI').lte(-0.2)
    low_str = image.select('STR').gte(0)
    high_str = image.select('STR').lte(20)
    high_evi = image.select('EVI').lte(1)
    low_evi = image.select('EVI').gte(-1)
    return (image
            .updateMask(full_cover)
            .updateMask(no_water)
            .updateMask(low_str)
            .updateMask(high_str)
            .updateMask(high_evi)
            .updateMask(low_evi))

S2_boa = S2_boa.map(apply_additional_filters)

# KERNAL SMOOTH

"""
S2_boa still may contain some outliers (water bodies, shadows, clouds).
To minimize the impact of those pixels - we apply Kernal smooth
"""

def boxcar_smooth(image):
    boxcar = ee.Kernel.square(radius=10, units='pixels', normalize=True)
    image_smoothed = image.convolve(boxcar)
    return image_smoothed

S2_boa = S2_boa.map(boxcar_smooth)

# Convert STR and NDVI values to array
array = S2_boa.select(['STR', 'NDVI']).toArray()
axes = {'image': 0, 'band': 1}
str_array = array.arraySlice(axes['band'], 0, 1)  # creates a subarray of STR
ndvi_array = array.arraySlice(axes['band'], 1, 2)  # Creates a subarray of NDVI

# Create the empty variables for the further loop
ndvi_interval = ee.List([])
ndvi_interval2 = ee.List([])
str_interval = ee.List([])
str_interval2 = ee.List([])

# Set the initial parameters for Wet edge estimation
j = 0  # the min number of subinterval
minimal_ndvi = 0.1  # min NDVI value for which the Wet edge will be calculated
maximal_ndvi = 0.7  # max NDVI value for which the Wet edge will be calculated
sub_number = 10  # number of subintervals in each interval
step = 0.001  # step for intervals
reduced_ndvi = ndvi_array.arrayReduce(ee.Reducer.median(), ee.List([0]))  # one median NDVI image to list

"""
Within the NDVI-STR space we derive the max STR value for each NDVI subinterval
This max STR value is arranged with the median NDVI value of each NDVI subinterval
"""

for i in range(int(minimal_ndvi * 1000), int(maximal_ndvi * 1000), int(step * 1000)):
    i = i / 1000.0
    # assign the working subinterval
    temp_mask1 = reduced_ndvi.gt(ee.Number(i))
    temp_mask2 = reduced_ndvi.lt(ee.Number(i + step))

    # derive the max STR within the working subinterval
    temp_mask = temp_mask1.And(temp_mask2)
    temp_mask = temp_mask.arrayProject(ee.List([0]))
    temp_mask = temp_mask.arrayFlatten([['array']])
    temp_masked_array = str_array.updateMask(temp_mask)
    temp_masked_array2 = temp_masked_array.arrayProject(ee.List([0]))
    temp_masked_array2 = temp_masked_array2.arrayReduce(ee.Reducer.max(), ee.List([0]))
    temp_masked_array2 = temp_masked_array2.arrayFlatten([['array']])

    # derive the median NDVI value of the working subinterval
    temp_median_ndvi_per_str = temp_masked_array2.reduceRegion(
        reducer=ee.Reducer.median(),
        geometry=table,
        scale=10
    )
    added_list = ndvi_interval.add(i + step / 2)
    ndvi_interval = ee.List(ee.Algorithms.If(temp_median_ndvi_per_str.get('array'), added_list, ndvi_interval))
    value_to_add = ee.Number(ee.Algorithms.If(
        temp_median_ndvi_per_str.get('array'),
        ee.Number(temp_median_ndvi_per_str.get('array')),
        ee.Number(0)
    ))
    value_to_add2 = ee.Number(ee.Algorithms.If(temp_median_ndvi_per_str.get('array'), i + step / 2, ee.Number(0)))
    str_interval = str_interval.add(value_to_add)
    ndvi_interval2 = ndvi_interval2.add(value_to_add2)
    j += 1

# min STR values within each NDVI subinterval
# print(str_interval, "max STR within each NDVI subinterval")


# NDVI AND STR VALUES FOR WET EDGE ESTIMATION

# Create the empty variables for the further loop
points = ee.List([])
xValues = ee.List([])
yValues = ee.List([])
xValues2 = ee.List([])
yValues2 = ee.List([])

total_sub_number = j  # number of subintervals

# Each NDVI interval has subintervals within which we derived max STR values.
# Now, within each interval we calculate the median and std of max STR values.
# Within each interval, we filter out max STR values that are bigger than median max STR + std max STR.
# The remained max STR values are averaged (median) and associated with the median NDVI value within each interval.
for i in range(0, total_sub_number, sub_number):
    interval_str = ee.List([])
    interval_str2 = ee.List([])
    ndvi_range = ee.List([])

    # Within each interval find the max STR values that are lower than median max STR + std max STR
    for j in range(i, i + sub_number):
        added_val = interval_str.add(str_interval.get(j))
        ndvi_range2 = ndvi_range.add(ndvi_interval2.get(j))
        ndvi_range = ee.List(ee.Algorithms.If(str_interval.get(j), ndvi_range2, ndvi_range))
        interval_str = ee.List(ee.Algorithms.If(str_interval.get(j), added_val, interval_str))
        interval_str2 = interval_str2.add(ee.Number(ee.Algorithms.If(str_interval.get(j), str_interval.get(j), 0)))

    median_value = interval_str.reduce(ee.Reducer.median())
    std_value = interval_str.reduce(ee.Reducer.stdDev())
    threshold1 = ee.Number(median_value)
    threshold2 = ee.Number(median_value)
    added_th = threshold2.add(ee.Number(std_value))
    threshold2 = ee.Number(ee.Algorithms.If(threshold1, added_th, -1000))
    remaining_interval = ee.List([])

    for k in range(sub_number):
        temp_val = ee.Number(interval_str2.get(k))
        added_interval = remaining_interval.add(temp_val)
        remaining_interval = ee.List(ee.Algorithms.If(
            temp_val.lt(threshold2).And(temp_val.gt(ee.Number(0))),
            added_interval,
            remaining_interval
        ))

    # Average remained max STR values for each interval and calculate median NDVI within each interval
    mean_of_interest = remaining_interval.reduce(ee.Reducer.median())
    mean_ndvi = ndvi_range.reduce(ee.Reducer.median())
    extended_points = points.add([mean_ndvi, mean_of_interest])
    points = ee.List(ee.Algorithms.If(mean_of_interest, extended_points, points))
    ex_xValues = xValues.add(mean_ndvi)
    ex_yValues = yValues.add(mean_of_interest)
    yValues = ee.List(ee.Algorithms.If(mean_of_interest, ex_yValues, yValues))
    xValues = ee.List(ee.Algorithms.If(mean_of_interest, ex_xValues, xValues))

    xval_to_add = ee.Algorithms.If(mean_of_interest, mean_ndvi, 0)
    yval_to_add = ee.Algorithms.If(mean_of_interest, mean_of_interest, 0)
    yValues2 = yValues2.add(yval_to_add)
    xValues2 = xValues2.add(xval_to_add)

# NDVI and STR values that are used for further Wet edge estimation
# print("NDVI and STR for Wet edge", points)


# PARAMETERS FOR WET EDGE

# Linear model with all the max STR and NDVI values
linearFit = ee.Dictionary(points.reduce(ee.Reducer.linearFit()))
# print(linearFit);
# print('y-intercept (preliminary):', linearFit.get('offset'))
# print('Slope (preliminary):', linearFit.get('scale'))

# To ensure that there are not STR outliers for Wet edge estimation -
# we calculate RMSE of the resulted points. We filter out the max STR
# values that do not lie within the range modelled max STR ± 2 * RMSE

# Model the max STR values using parameters from linearFit
def prediction(x):
    temp = ee.Number(x)
    temp = temp.multiply(ee.Number(linearFit.get('scale')))
    temp = temp.add(ee.Number(linearFit.get('offset')))
    return temp


preds = xValues2.map(prediction)
rmse = ee.List([])

for k in range(int(total_sub_number / sub_number)):
    res = ee.Algorithms.If(yValues2.get(k), yValues2.get(k), ee.Number(0))
    res1 = ee.Number(res)
    y_pred = ee.Algorithms.If(preds.get(k), preds.get(k), ee.Number(0))
    res1 = res1.subtract(y_pred)
    res1 = res1.multiply(res1)
    added_rmse = rmse.add(res1)
    rmse = ee.List(ee.Algorithms.If(res, added_rmse, rmse))

# Calculate double RMSE
rmse_reduced = ee.Number(rmse.reduce(ee.Reducer.mean()))
rmse_reduced = rmse_reduced.sqrt()
doubled_rmse = rmse_reduced.multiply(2)
# print("Double RMSE", doubled_rmse)

# Create the empty variables for the further loop
reduced_points = ee.List([])
reduced_X = ee.List([])
reduced_Y = ee.List([])

# Filter out the max STR values that do not lie within the range: modelled max STR + 2 * RMSE
for k in range(int(total_sub_number / sub_number)):
    y_true = yValues2.get(k)
    y_true = ee.Number(y_true)
    y_pred = ee.Number(preds.get(k))
    diff = ee.Number(y_true.subtract(y_pred))
    diff2 = ee.Number(y_pred.subtract(y_true))
    added_reduced_points = reduced_points.add([xValues2.get(k), yValues2.get(k)])
    condition = y_pred.lt(y_true)
    condition2 = y_pred.gt(y_true)
    condition = condition.And(diff.lt(doubled_rmse))
    condition2 = condition2.And(diff2.lt(doubled_rmse))
    condition = condition.And(y_true)
    condition2 = condition2.And(y_true)
    condition = condition.Or(condition2)

    reduced_points = ee.List(ee.Algorithms.If(condition, added_reduced_points, reduced_points))

# Final points that can be used for the Wet edge estimation
# print("Reduced NDVI and STR for Wet edge:", reduced_points)

# Wet edge's parameters
linearFit2 = ee.Dictionary(reduced_points.reduce(ee.Reducer.linearFit()))
# print('y-intercept:', linearFit2.get('offset'))
# print('Slope:', linearFit2.get('scale'))

# Export lm parameters
par1 = ee.Number(linearFit2.get('offset')).getInfo()
print("y-intercept", par1)
par2 = ee.Number(linearFit2.get('scale')).getInfo()
print("Slope", par2)

par_list = ee.List([[par1, par2]])
# print(par_list)
par_fc = ee.FeatureCollection(par_list.map(lambda point: ee.Feature(None, {'coef': point})))
table_name = ee.String(table.get("system:id")).split('/').get(3)

# # Export the table to Google Drive
# task = ee.batch.Export.table.toDrive(
#     collection=par_fc,
#     description=table_name.getInfo(),
#     folder='LM_parameters_NDVI_wetedge_5cloudprov',
#     fileFormat='CSV'
# )

# # Start the export task
# task.start()
