import os
import time
import shutil
import subprocess
from pystac_client import Client
import requests
from typing import List
from typing import Tuple
import shapely
from shapely.geometry import box, shape, mapping
import numpy as np
from osgeo import gdal
import osgeo_utils.gdal_merge
import geopandas as gpd
import urllib.parse
from planetary_computer import sign_url

from phase_0_fetch.src.download_dem import snakemake_type_exists, get_shapefile_extent

gdal.UseExceptions()


def query_imagery_stac(
    intersects: dict = None,
    datetime: str = None,
    collections: List[str] = None,
    query: dict = None,
    max_items: int = None,
) -> List[dict]:
    """Queries the STAC server for data with the following input parameters and returns a list of scenes

    Parameters
    ----------
    intersects: GeoJSON object
        A GeoJSON object representing the geographical area of interest
    daterange: string
        A string specifying the date range for the query in the format 'YYYY-MM-DD/YYYY-MM-DD'
    collections: string of list of strings
        A string or a list of strings specifying the collections to search. Data found here: https://planetarycomputer.microsoft.com/catalog
    query: dictionary
        dictionary for additional query parameters
    max_items: integer
        An integer specifying the maximum number of items to return

    Returns
    -------
    query: list
        list of scenes from the stac search

    """

    stac = Client.open("https://planetarycomputer.microsoft.com/api/stac/v1")

    if intersects or datetime or collections or query is not None:
        query = stac.search(
            collections=collections,
            max_items=max_items,
            intersects=intersects,
            datetime=datetime,
            query=query,
        )
        return query
    else:
        return print(
            "Must set at least one of the following parameters: collections, intersects, datetime or query before continuing."
        )


def calculate_coverage_ratio(scene_geom, aoi_geom):
    """Compute what area of the bounding box geometry is covered by the scene geometry.

    Parameters
    ----------
    scene_geom: geometry
        geometry of the scene
    aoi_geom: geometry
        geometry of the remaining area that isn't covered by a scene

    Returns
    -------
    proportion: float
        proportion of the area in the aoi_geometry that is covered by the scene

    """

    intersection = scene_geom.intersection(aoi_geom)

    if aoi_geom.area == 0.0:
        return 1.0
    else:
        return intersection.area / aoi_geom.area


def get_items_list(chosen_item, query_items, bounding_box_geom):
    """get list of items that covers entire bounding box geometry

    Parameters
    ----------
    chosen item: scene item
        scene item that has the most coverage over the bounding box geometry
    query_items: list of scene items
        list of scene items from the query
    bounding_box_geom: geometry
        geometry of the remaining bounding box that isn't covered by a scene

    Returns
    -------
    items_list: list of scene items
        list of scenes from the stac search that will be used to cover the bounding box geometry

    """

    # make a list with the first item
    items_list = [chosen_item]

    # get the geometry of the first item
    downloaded_area = shape(chosen_item.geometry)

    # determine how much area in the bounding box is required to download
    remaining_area = bounding_box_geom.difference(downloaded_area)

    # while loop that adds image to the items_list until the remaining area is either zero or there's not enough images
    not_enough_scenes = False
    while remaining_area.area > 0.0 and not_enough_scenes == False:

        # calculate the overlapped area that of each scene and the reamaining area
        remaining_overlaps = np.array(
            [
                calculate_coverage_ratio(shape(item.geometry), remaining_area)
                for item in query_items
            ]
        )

        # if there is no overlapped area, end the loop
        if max(remaining_overlaps) == 0.0:
            not_enough_scenes = True
            print(
                "couldn't find enough landsat scenes to cover area with both platforms, consider changing date range and platforms to search"
            )

        else:
            # chose the item that has the most overlap
            next_chosen_item = query_items[np.argmax(remaining_overlaps)]

            # add it to the list
            items_list += [next_chosen_item]

            # update how much area is covered by a scene
            downloaded_area = downloaded_area.union(shape(next_chosen_item.geometry))

            # find remaining area that has no scene coverage
            remaining_area = bounding_box_geom.difference(downloaded_area)

    return items_list


def download_scenes(extent, items, target_asset_key, output_file):
    """Download all the scenes in the extent with gdalwrap

    Parameters
    ----------
    extent: list
        list of minx, miny, maxx, maxy values of the bounding box
    items: list of scene items
        list of scene items to download
    target_asset_key: string
        name of target asset to download, e.g., band to download. You'll need to look up the names of the assets
        for the particular collection: https://planetarycomputer.microsoft.com/catalog
    output_file: string
        path of the image output

    Returns
    -------
    none

    """

    # initialize a list of the stac urls to be downloaded
    stac_urls = []
    for item in items:
        for asset_key, asset in item.assets.items():
            # Check if the asset is one of the desired bands
            if asset_key == target_asset_key:
                stac_urls += [sign_url(asset.href)]

    # GDAL command
    command = (
        [
            "gdalwarp",
            "-te",
            str(extent[0]),
            str(extent[1]),
            str(extent[2]),
            str(extent[3]),
            "-t_srs",
            "EPSG:4326",
        ]
        + ["/vsicurl/" + stac_url for stac_url in stac_urls]
        + [output_file]
    )

    # Run the command
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # Check for errors
    if result.returncode != 0:
        print("Error running gdalwarp:", result.stderr.decode())
    else:
        print("gdalwarp completed successfully:", result.stdout.decode())


def get_imagery(
    extent,
    out_dir,
    collection_id,
    query,
    target_asset_keys,
    start_query,
    end_query,
    suffix,
):
    """Gets imagery for a from PlanetaryComputer

    Parameters
    ----------
    extent: list
        list of minx, miny, maxx, maxy values of the bounding box
    out_dir: string
        path where files will be downloaded
    collection_id: string
        name of the collection, see here: https://planetarycomputer.microsoft.com/catalog
    target_asset_keys: list of strings
        name of target asset to download, e.g., band to download. You'll need to look up the names of the assets
        for the particular collection: https://planetarycomputer.microsoft.com/catalog
    query: dictionary
        dictionary used to filter the query search
    start_query: datetime
        YYYY-MM-DD date that brackets the start the search. You can also add time in this format: YYYY-MM-DDTHH:MM:SSZ"
    start_query: datetime
        YYYY-MM-DD date that brackets the end the search. You can also add time in this format: YYYY-MM-DDTHH:MM:SSZ"
    suffix: string
        string suffix to add to the end of the file name

    Returns
    -------
    none

    """

    # make out_dir
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Load extent shapefile
    extent_shp = gpd.read_file(extent_shpfile).to_crs("EPSG:4326")
    west, east, length, south, north, height = get_shapefile_extent(extent_shp)

    # make geo box for geojson
    geo_bbox = (west, south, east, north)

    # Get the geometries of the items
    bounding_box_geom = box(*geo_bbox)

    # Set intersects using polygon geojson
    geojson_obj = mapping(bounding_box_geom)

    # query the stac api to scene matches
    query_return = query_imagery_stac(
        collections=[collection_id],
        query=query,
        intersects=geojson_obj,
        datetime=start_query + "/" + end_query,
    )

    # Store the items in a list
    items = list(query_return.item_collection())

    # Filter out items with banding issues
    # Assuming 'banding_issues' is a flag in the properties (this will depend on the dataset)
    QA_items = [
        item
        for item in items
        if "banding_issues" not in item.properties
        or item.properties["banding_issues"] != True
    ]

    # calculate how much the site overlaps
    overlaps = np.array(
        [
            calculate_coverage_ratio(shape(item.geometry), bounding_box_geom)
            for item in QA_items
        ]
    )

    # choose item with most site overlap
    chosen_item = QA_items[np.argmax(overlaps)]

    # if the AOI is covered by a single scene
    if max(overlaps) == 1.0:
        chosen_items = [chosen_item]
    # if multiple scenes are needed
    else:
        chosen_items = get_items_list(chosen_item, QA_items, bounding_box_geom)

    # download the scenes for each target asset in the collection
    for target_asset_key in target_asset_keys:
        download_scenes(
            [west, south, east, north],
            chosen_items,
            target_asset_key,
            out_dir + "/" + target_asset_key + suffix + ".tif",
        )


if __name__ == "__main__":
    extent_shpfile = snakemake.input["extent_shpfile"]
    out_bands = snakemake.output["out_bands"]
    out_dir = os.path.dirname(out_bands[0])
    collection_id = snakemake_type_exists(
        snakemake.params, "collection_id", "landsat-c2-l2"
    )
    query = snakemake_type_exists(snakemake.params, "query", {})
    target_asset_keys = snakemake_type_exists(
        snakemake.params, "target_asset_keys", ["red", "green", "blue"]
    )
    start_query = snakemake_type_exists(snakemake.params, "start_query", "2013-01-01")
    end_query = snakemake_type_exists(snakemake.params, "end_query", "2025-06-01")
    suffix = snakemake_type_exists(snakemake.params, "suffix", "")

    get_imagery(
        extent_shpfile,
        out_dir,
        collection_id,
        query,
        target_asset_keys,
        start_query,
        end_query,
        suffix,
    )
