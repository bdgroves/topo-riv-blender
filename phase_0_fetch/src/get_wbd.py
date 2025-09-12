import os
import pandas as pd
import geopandas as gpd
from pygeohydro import WBD


def wbd_download(huc_list, extent_shpfile):
    """Downloads wbd data (.shp) according to a list of huc ids

    Parameters
    ----------
    huc_list: list of strings
        strings of the huc ids to download
    extent_shpfile: string
        path of the extent shapefile to be outputted

    Returns
    -------
    none

    """

    # Make out_dir if it doesn't exist
    out_dir = os.path.dirname(extent_shpfile)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Create list of geopandas shapefiles holding HUC extents
    wbd_list = []

    # Run through HUC or HUCs
    for huc in huc_list:

        # Determine HUC level
        huc_level = str(len(huc))

        # Set up WBD
        wbd = WBD("huc" + huc_level)

        # Add the HUC to the list
        wbd_list += [wbd.byids("huc" + huc_level, [huc]).dissolve()]

    # Merge the HUCs in to a single geometry
    extent = gpd.GeoDataFrame(pd.concat(wbd_list)).dissolve()

    # Save WBD data as a shapefile
    extent.to_file(extent_shpfile)


if __name__ == "__main__":
    huc_list = snakemake.params["huc_list"]
    extent_shpfile = snakemake.output["extent_shpfile"]
    wbd_download(huc_list, extent_shpfile)
