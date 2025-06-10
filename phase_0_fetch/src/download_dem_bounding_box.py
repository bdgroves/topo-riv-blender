import os
import geopandas as gpd
from bmi_topography import Topography
from shapely.geometry import Polygon

from phase_0_fetch.src.download_dem import (
    snakemake_type_exists,
    buffer_shapefile,
    opentopography_api_download,
    determine_dem_product,
    estimate_area,
)


def bmi_download_dem_bbox(
    UL_corner,
    LR_corner,
    extent_shpfile,
    demfile,
    data_product,
    dem_product,
    buffer,
    render_pixels,
):
    """Downloads a DEM (.tif) according to a bounding box and makes a shapefile of the extent of the bbox.

    Parameters
    ----------
    UL_corner : 2-element tuple
        lat, long coordinate of the upper-left corner
    LR_corner: 2-element tuple
        lat, long coordinate of the lower-right corner
    extent_shpfile: string
        path of of the extent shapefile of the bbox to be outputted
    demfile: string
        path of the dem file to be downloaded
    data_product: string
        name of the data product used by the opentopography api
    dem_product: string
        name of the dem product used by the opentopography api
    buffer: float
        a percent buffer around the extent; maintains original aspect
    render_pixels: int
        number of pixels in the final render image, used to determine the dem resolution

    Returns
    -------
    none

    """

    # Buffer those boundaries (expand) using a user-defined buffering percentage
    buff_west, buff_east, buff_south, buff_north = buffer_shapefile(
        float(UL_corner[1]),
        float(LR_corner[1]),
        float(LR_corner[1]) - float(UL_corner[1]),
        float(LR_corner[0]),
        float(UL_corner[0]),
        float(UL_corner[0]) - float(LR_corner[0]),
        buffer,
    )

    # Determine the correct data and dem product if set to Auto-Mode
    temp_data_product, temp_dem_product = determine_dem_product(
        data_product,
        dem_product,
        buff_west,
        buff_east,
        buff_south,
        buff_north,
        render_pixels,
    )
    print("Using " + temp_data_product + " with " + temp_dem_product)

    # Get Default parameter set from BMI-Topography
    params = Topography.DEFAULT.copy()

    # Set parameter values to user-defined settings
    params["data_type"] = temp_data_product
    params["dem_type"] = temp_dem_product
    params["output_format"] = "GTiff"
    params["west"] = buff_west
    params["south"] = buff_south
    params["east"] = buff_east
    params["north"] = buff_north
    params["cache_dir"] = os.path.dirname(demfile)

    # Create a polygon from the bounding box
    domain_geom = Polygon(
        zip(
            [UL_corner[1], UL_corner[1], LR_corner[1], LR_corner[1], UL_corner[1]],
            [UL_corner[0], LR_corner[0], LR_corner[0], UL_corner[0], UL_corner[0]],
        )
    )

    # Create a geopandas dataframe from the polygon and set the CRS to WGS84
    domain_polygon = gpd.GeoDataFrame(
        index=[0], crs="EPSG:4326", geometry=[domain_geom]
    )

    # Save the geopandas dataframe as a shapefile
    domain_polygon.to_file(extent_shpfile)

    # Send a request to the OpenTopography API with the user-defined parameters
    opentopography_api_download(params, demfile)


if __name__ == "__main__":
    data_product = snakemake_type_exists(
        snakemake.params, "data_product", "/API/globaldem"
    )
    dem_product = snakemake_type_exists(snakemake.params, "dem_product", "NASADEM")
    render_pixels = snakemake_type_exists(snakemake.params, "render_pixels", 1)
    buffer = snakemake_type_exists(snakemake.params, "buffer", 1.0)
    UL_corner = snakemake.params["UL_corner"]
    LR_corner = snakemake.params["LR_corner"]
    extent_shpfile = snakemake.output["extent_shpfile"]
    demfile = snakemake.output["demfile"]
    download_mode = snakemake_type_exists(snakemake.params, "download_mode", "bmi")
    bmi_download_dem_bbox(
        UL_corner,
        LR_corner,
        extent_shpfile,
        demfile,
        data_product,
        dem_product,
        buffer,
        render_pixels,
    )
