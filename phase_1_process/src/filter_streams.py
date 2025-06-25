import os
import geopandas as gpd

from phase_0_fetch.src.download_dem import snakemake_type_exists


def filter_rivers(nhd_shpfile, min_stream_order, filtered_nhd_shpfile):
    """Filters out stream orders smaller than a user-specified stream order

    Parameters
    ----------
    nhd_shpfile: string
        path of of the river shapefile to be outputted
    min_stream_order: int
        minimum stream order to keep
    filtered_nhd_shpfile: string
        path of of the river shapefile to be outputted

    Returns
    -------
    none

    """

    # Make out_dir if it doesn't exist
    out_dir = os.path.dirname(filtered_nhd_shpfile)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Load extent shapefile, ensure it's in geographic coordinates
    nhd_shp = gpd.read_file(nhd_shpfile)

    # Filter out unwanted stream orders
    filtered_nhd_shp = nhd_shp[nhd_shp["StreamOrde"] >= min_stream_order]

    # Save NHD data as a shapefile
    filtered_nhd_shp.to_file(filtered_nhd_shpfile)


if __name__ == "__main__":
    nhd_shpfile = snakemake.input["nhd_shpfile"]
    min_stream_order = snakemake_type_exists(snakemake.params, "min_stream_order", 3)
    filtered_nhd_shpfile = snakemake.output["filtered_nhd_shpfile"]
    filter_rivers(nhd_shpfile, min_stream_order, filtered_nhd_shpfile)
