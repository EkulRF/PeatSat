import requests
import pystac
import folium
import shapely.geometry
from IPython.display import Image
import ee
ee.Authenticate()
ee.Initialize()

bbox = [-4.4375, 56.0702778, -3.7752778, 58.3775]
sitename = 'FlowCountry'

import pystac_client

catalog = pystac_client.Client.open(
    "https://planetarycomputer.microsoft.com/api/stac/v1"
)

aoi = {
    "type": "Polygon",
    "coordinates": [
        [
            [29.036865234375, 7.857940257224196],
            [31.4813232421875, 7.857940257224196],
            [31.4813232421875, 10.055402736564236],
            [29.036865234375, 10.055402736564236],
            [29.036865234375, 7.857940257224196],
        ]
    ],
}


collection = "sentinel-2-l2a"
query = {"eo:cloud_cover": {"lt": "10"}}
search = catalog.search(intersects=aoi, collections="sentinel-1-rtc", query=query)