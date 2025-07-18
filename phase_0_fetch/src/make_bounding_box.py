import os
import geopandas as gpd
from shapely.geometry import Polygon

from phase_0_fetch.src.download_dem import snakemake_type_exists


def make_bbox(UL_corner, LR_corner, extent_shpfile, crs):
    """Creates a shape file based on coordinates of a bounding box

    Parameters
    ----------
    UL_corner : 2-element tuple
        lat, long coordinate of the upper-left corner
    LR_corner: 2-element tuple
        lat, long coordinate of the lower-right corner
    extent_shpfile: string
        path of of the extent shapefile of the bbox to be outputted
    crs: string
        coordinate reference system, defaults to EPSG:4326

    Returns
    -------
    none

    """

    # Create a polygon from the bounding box
    domain_geom = Polygon(
        zip(
            [UL_corner[1], UL_corner[1], LR_corner[1], LR_corner[1], UL_corner[1]],
            [UL_corner[0], LR_corner[0], LR_corner[0], UL_corner[0], UL_corner[0]],
        )
    )

    # Create a geopandas dataframe from the polygon and set the CRS to WGS84
    domain_polygon = gpd.GeoDataFrame(index=[0], crs=crs, geometry=[domain_geom])

    # Save the geopandas dataframe as a shapefile
    domain_polygon.to_file(extent_shpfile)


if __name__ == "__main__":
    UL_corner = snakemake.params["UL_corner"]
    LR_corner = snakemake.params["LR_corner"]
    extent_shpfile = snakemake.output["extent_shpfile"]
    crs = snakemake_type_exists(snakemake.params, "crs", "EPSG:4326")

    make_bbox(UL_corner, LR_corner, extent_shpfile, crs)
