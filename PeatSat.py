import requests
import pystac
import pystac_client
import folium
import shapely.geometry
from IPython.display import Image
import ee
import geemap
ee.Authenticate()
ee.Initialize()

bbox = [-4.4375, 56.0702778, -3.7752778, 58.3775]
time_range = ["2020-12-01","2022-12-31"]
min_lon, min_lat, max_lon, max_lat = bbox
sitename = 'FlowCountry'

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

def clip_image(image):
    return image.clipToBoundsAndScale(geometry=aoi, scale=30)

# Google Earth Engine S1 Extraction
def get_sentinel1_collection_GEE(bbox, time_range):
    region = ee.Geometry.Rectangle(bbox)

    sentinel1 = (ee.ImageCollection('COPERNICUS/S1_GRD')
                 .filterBounds(region)
                 .filterDate(time_range[0], time_range[1])
                 .select('angle')
                 )

    return sentinel1.map(clip_image)

def get_sentinel1_collection_MPC(bbox, time_range):

    catalog = pystac_client.Client.open(
    "https://planetarycomputer.microsoft.com/api/stac/v1")

    collection = ["sentinel-1-rtc"]
    search = catalog.search(collections=collection, bbox=bbox, datetime=time_range)
    sentinel1 = search.item_collection()

    return sentinel1


sentinel1_viewingangle = get_sentinel1_collection_GEE(bbox, time_range)

print("Number of Sentinel-1 images from GEE:", sentinel1_viewingangle.size().getInfo())

print("Number of Sentinel-1 images from MPC:", len(get_sentinel1_collection_MPC(bbox, time_range)))

