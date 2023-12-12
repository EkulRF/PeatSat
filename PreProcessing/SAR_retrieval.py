import requests
import pystac
import pystac_client
import folium
import shapely.geometry
from IPython.display import Image
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

    sentinel1_viewingangle = get_sentinel1_collection_GEE(bbox, aoi, time_range)
    sentinel1_rtc = get_sentinel1_collection_MPC(bbox, time_range)

    print("Number of Sentinel-1 images from GEE:", sentinel1_viewingangle.size().getInfo())
    print("Number of Sentinel-1 images from MPC:", len(sentinel1_rtc))

    return


def clip_image(image, aoi):
    return image.clipToBoundsAndScale(geometry=aoi, scale=30)

# Google Earth Engine S1 Extraction
def get_sentinel1_collection_GEE(bbox, aoi, time_range):
    region = ee.Geometry.Rectangle(bbox)

    sentinel1 = (ee.ImageCollection('COPERNICUS/S1_GRD')
                 .filterBounds(region)
                 .filterDate(time_range[0], time_range[1])
                 .select('angle')
                 )

    return sentinel1.map(lambda image: clip_image(image, aoi))

def get_sentinel1_collection_MPC(bbox, time_range):

    catalog = pystac_client.Client.open(
    "https://planetarycomputer.microsoft.com/api/stac/v1")

    collection = ["sentinel-1-rtc"]
    search = catalog.search(collections=collection, bbox=bbox, datetime=time_range)
    sentinel1 = search.item_collection()

    return sentinel1