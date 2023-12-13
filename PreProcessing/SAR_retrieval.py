import requests
import pystac_client
import ee
ee.Authenticate()
ee.Initialize()

def SAR_retrieval(bbox, time_range, sitename):

    min_lon, min_lat, max_lon, max_lat = bbox

    aoi = {
        "type": "Polygon",
        "coordinates": [
            [
            [min_lon, min_lat],
            [max_lon, min_lat],
            [max_lon, max_lat],
            [min_lon, max_lat],
            [min_lon, min_lat],
        ]
        ],
    }

    sentinel1_viewingangle = get_sentinel1_viewingAngles_GEE(bbox, aoi, time_range)
    sentinel1_rtc = get_sentinel1_collection_MPC(bbox, time_range)

    print("Number of Sentinel-1 images from GEE:", sentinel1_viewingangle.size().getInfo())
    print("Number of Sentinel-1 images from MPC:", len(sentinel1_rtc))

    sentinel1 = MergeLayers(sentinel1_rtc, sentinel1_viewingangle)

    return sentinel1


def clip_image(image, aoi):
    return image.clipToBoundsAndScale(geometry=aoi, scale=30)

# Google Earth Engine S1 Extraction
def get_sentinel1_viewingAngles_GEE(bbox, aoi, time_range):
    region = ee.Geometry.Rectangle(bbox)

    sentinel1 = (ee.ImageCollection('COPERNICUS/S1_GRD')
                 .filterBounds(region)
                 .filterDate(time_range[0], time_range[1])
                 .select('angle')
                 )

    sentinel1 = sentinel1.map(add_orbit_number)

    return sentinel1.map(lambda image: clip_image(image, aoi))

def get_sentinel1_collection_MPC(bbox, time_range):

    catalog = pystac_client.Client.open(
    "https://planetarycomputer.microsoft.com/api/stac/v1")

    collection = ["sentinel-1-rtc"]
    search = catalog.search(collections=collection, bbox=bbox, datetime=time_range)
    sentinel1 = search.item_collection()

    return sentinel1

def add_orbit_number(image):
    # Get the platform number, relative orbit number start, and orbit properties pass
    platform_number = image.get('platform_number')
    relative_orbit_number_start = image.get('relativeOrbitNumber_start')
    orbit_properties_pass = image.get('orbitProperties_pass')
    
    # Concatenate the strings to form the platform_relorbit property
    platform_relorbit = (
        ee.String(platform_number)
        .cat('_')
        .cat(ee.Number(relative_orbit_number_start).format('%.0f'))
        .cat('_')
        .cat(orbit_properties_pass)
    )
    
    # Return the image with the platform_relorbit property added
    return image.set('platform_relorbit', platform_relorbit)


def MergeLayers(sentinel1_rtc, sentinel1_viewingangle):

    for image in sentinel1_rtc:

        orbit_no = image.properties['sat:relative_orbit']

        Selected_AngleImage = sentinel1_viewingangle.filter(ee.Filter.eq('orbit_number', orbit_no)).first()

        image.assets['view_angle'] = Selected_AngleImage

    return sentinel1_rtc