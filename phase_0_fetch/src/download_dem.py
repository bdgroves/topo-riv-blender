import os
import geopandas as gpd
import numpy as np
from bmi_topography import Topography
import requests
from osgeo import gdal
import osgeo_utils.gdal_merge


def snakemake_type_exists(snakemake_type, string, default_input):
    """Checks if string exists in a snakemake type and applies a default value if it does not exist.

    Parameters
    ----------
    snakemake_type: snakemake type
        e.g. input, output, params, wildcard
    string: string
        name of the string
    default_input:
        default value if the string does not exist

    Returns
    -------
    the original value if exists and the default if not

    """

    if hasattr(snakemake_type, string):
        return snakemake_type[string]
    else:
        return default_input


def get_shapefile_extent(shapefile):
    """Gets a extent coodinates and the height + length according to a shapefile

    Parameters
    ----------
    shapefile: geopandas data frame
        geopandas data frame holding the extent geometry

    Returns
    -------
    west: float
        west long coordinate
    east: float
        east long coordinate
    length: float
        length of the total extent
    south: float
        south lat coordinate
    north: float
        north lat coordinate
    height: float
        height of the total extent

    """

    # Get the total bounds of the shapefile
    west, south, east, north = shapefile.geometry.total_bounds

    # Calculate the x-length of the bounding box
    length = east - west

    # Calculate the y-length of the bounding box
    height = north - south

    return west, east, length, south, north, height


def buffer_shapefile(west, east, length, south, north, height, buffer):
    """Gets a new extent coordinates according to a buffer

    Parameters
    ----------
    west: float
        west long coordinate
    east: float
        east long coordinate
    length: float
        length of the total extent
    south: float
        south lat coordinate
    north: float
        north lat coordinate
    height: float
        height of the total extent
    buffer: float
        a percent buffer around the extent; maintains original aspect

    Returns
    -------
    west: float
        west long coordinate with buffer
    east: float
        east long coordinate with buffer
    south: float
        south lat coordinate with buffer
    north: float
        north lat coordinate with buffer

    """

    # Add a buffer to the total bounds according to the product of a buffer percentage and the corresponding length or height of the bounding box
    west -= buffer / 100.0 * length
    south -= buffer / 100.0 * height
    east += buffer / 100.0 * length
    north += buffer / 100.0 * height

    return west, east, south, north


def estimate_area(west, east, south, north, render_pixels):
    """Determines the approximate land area within a render pixel

    Parameters
    ----------
    west: float
        west long coordinate
    east: float
        east long coordinate
    south: float
        south lat coordinate
    north: float
        north lat coordinate
    render_pixels: int
        number of pixels in the final render image, used to determine the dem resolution

    Returns
    -------
    area_m_per_pixel: float
        approximate land area within a render pixel

    """

    # Radius of the Earth in kilometers
    R = 6371.0  

    # Convert degrees to kilometers, y-distance
    height_km = (north - south) * (R * (np.pi / 180.0))

    # Convert degrees to kilometers, x-distance
    width_km = (
        (east - west) * (R * (np.pi / 180.0)) * np.cos((north + south) * np.pi / 360.0)
    )

    # Calculate area per render pixel in square meters
    return width_km * height_km * 1000000.0 / float(render_pixels)


def determine_dem_product(
    data_product, dem_product, west, east, south, north, render_pixels
):
    """Automatically determine the data and dem project based on the size of the extent and the number of pixels in the final render

    Parameters
    ----------
    data_product: string
        initial name of the data product used by the opentopography api
    dem_product: string
        initial name of the dem product used by the opentopography api
    west: float
        west long coordinate
    east: float
        east long coordinate
    south: float
        south lat coordinate
    north: float
        north lat coordinate
    render_pixels: int
        number of pixels in the final render image, used to determine the dem resolution

    Returns
    -------
    data_product: string
        automated name of the data product used by the opentopography api
    dem_product: string
        automated name of the dem product used by the opentopography api

    """

    # Determine dem product automatically
    if dem_product == "auto":

        # Calculate area per render pixel in square meters
        area_m_per_pixel = estimate_area(west, east, south, north, render_pixels)

        # Initialize switch from usgs to global as false
        switch_to_global = False

        # Require either usgsdem or global dem
        if data_product != "/API/usgsdem" and data_product != "/API/globaldem":
            raise Exception(
                "Incorrect data product. Please use '/API/usgsdem' or '/API/globaldem'."
            )

        # Determine USGS dem product if usgsdem is chosen
        if data_product == "/API/usgsdem":

            # 3DEP 10m
            if area_m_per_pixel < 900.0:
                return "/API/usgsdem", "USGS10m"

            # 3DEP 30m    
            elif 900.0 <= area_m_per_pixel and area_m_per_pixel < 8100.0:
                return "/API/usgsdem", "USGS30m"

            # Switch to global data, area per render pixel is too large
            else:
                print(
                    "This is a pretty large area, consider using data_product = '/API/globaldem' instead of '/API/usgsdem' next time. Switching to '/API/globaldem' for now."
                )
                switch_to_global = True 

        # Determine global dem product if globaldem is chosen or switched to
        if data_product == "/API/globaldem" or switch_to_global == True:

            # STRM 30m
            if area_m_per_pixel < 8100.0:
                return "/API/globaldem", "SRTMGL1"

            # STRM 90m
            elif 8100.0 <= area_m_per_pixel and area_m_per_pixel < 72900.0:
                return "/API/globaldem", "SRTMGL3"

            # STRM 500m + Bathymetry
            elif 72900.0 <= area_m_per_pixel and area_m_per_pixel < 2250000.0:
                return "/API/globaldem", "SRTM15Plus"

            # GEDI DTM 1000m
            else:
                print("You're using the coarsest DEM, this might be a lot of data.")
                return "/API/globaldem", "GEDI_L3"
    
    # Return user-specified parameters
    else:
        return data_product, dem_product


def opentopography_api_download(params, demfile):
    """Creates an API request to OpenTopography to download the dem

    Parameters
    ----------
    params: dict
        params for the api request
    demfile: string
        path of the dem file to be downloaded
    data_product: string
        name of the data product used by the opentopography api
    dem_product: string
        name of the dem product used by the opentopography api

    Returns
    -------
    none

    """

    # Make out_dir if it doesn't exist
    out_dir = os.path.dirname(demfile)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Set api key if exists (you need this when you run out of free requests)
    if os.path.isfile("opentopography_api_key.txt"):
        with open("opentopography_api_key.txt", "r") as f:
            api_key = f.read()
        params["api_key"] = api_key

    # Setup request and fetch
    topo_request = Topography(**params)

    # First try to download the DEM with the basic parameters
    try:
        topo_request.fetch()
        
        # Rename the file
        os.rename(
            out_dir
            + "/"
            + str(params["dem_type"])
            + "_"
            + str(params["south"])
            + "_"
            + str(params["west"])
            + "_"
            + str(params["north"])
            + "_"
            + str(params["east"])
            + ".tif",
            demfile,
        )

    # If you request a DEM that is too large, you'll need to download the dems in strips and then merge them.
    except requests.exceptions.HTTPError as e:
        
        # If the error is status code 400, figure out how large the initial request was and the maximum requested area
        if e.response.status_code == 400:
            print ("DEM is too large for one request, we have to download it in strips.")

            err = e.response.text

            # Get maximum area that can be requested
            max_area = float(
                err[
                    err.find(params["dem_type"])
                    + len(params["dem_type"])
                    + 4 : err.find(" km2. The")
                ].replace(",", "")
            )

            # Get the total area that was initial requested
            total_area = float(
                err[
                    err.find("selected area is") + 17 : err.find(" km2.</error>")
                ].replace(",", "")
            )

            # Make a temp directory to hold the dems strips
            tmp_dir = (
                out_dir.replace("out", "tmp") + "/" + os.path.basename(demfile)[:-4]
            )
            if not os.path.exists(tmp_dir):
                os.makedirs(tmp_dir)

            # Get initial request parameters
            data_product = params["data_type"]
            dem_product = params["dem_type"]
            buff_west =  params["west"]
            buff_south =  params["south"]
            buff_east =  params["east"] 
            buff_north =  params["north"]       

            # Calculate the number of strips needed to keep the area request under the max area
            strips = int(total_area / max_area) + 1

            # Calculate the total width of the request
            total_width = buff_east - buff_west

            # Calculate the width of each strip
            strip_width = total_width / float(strips)

            # Determine the left boundaries of the strip
            strip_left = buff_west + np.linspace(0.0, total_width - strip_width, strips)
            strip_right = buff_west + np.linspace(strip_width, total_width, strips)

            # Add a tiny addition to the right boundary of every strip except the last one, to ensure there are no gaps between the strips
            strip_right[:-1] += strip_width * 0.001

            # Make a list container to hold the dem strip paths
            dem_strips = []

            # Create a OpenTopography API request for each dem strip
            for i in range(0, strips):

                # Only redownload if the strip doesn't exist. When dealing with large datasets this prevents the code from
                # redownloading already retrieved data. However, this may cause issues if you download a huge dem and only
                # change the extent slightly. 
                if os.path.isfile(tmp_dir + "/dem_part_" + str(i) + ".tif") == False:

                    # Get Default parameter set from BMI-Topography
                    params = Topography.DEFAULT.copy()

                    # Set parameter values to user-defined settings
                    params["data_type"] = data_product
                    params["dem_type"] = dem_product
                    params["output_format"] = "GTiff"
                    params["west"] = strip_left[i]
                    params["south"] = buff_south
                    params["east"] = strip_right[i]
                    params["north"] = buff_north
                    params["cache_dir"] = tmp_dir
                    params["api_key"] = api_key

                    # Print dem strip downloading status
                    print(
                        "downloading... "
                        + "dem_part_"
                        + str(i)
                        + ".tif: "
                        + str(i + 1)
                        + " out of "
                        + str(strips)
                    )

                    # Setup request and fetch
                    topo_request = Topography(**params)
                    topo_request.fetch()

                    # Rename the file
                    os.rename(
                        tmp_dir
                        + "/"
                        + str(params["dem_type"])
                        + "_"
                        + str(params["south"])
                        + "_"
                        + str(params["west"])
                        + "_"
                        + str(params["north"])
                        + "_"
                        + str(params["east"])
                        + ".tif",
                        tmp_dir + "/dem_part_" + str(i) + ".tif",
                    )

                # Add path to dem strip path list    
                dem_strips += [tmp_dir + "/dem_part_" + str(i) + ".tif"]

            # Use GDAL to open the first dem strip
            srs = gdal.Open(tmp_dir + "/dem_part_0.tif")
            
            # Get the NaN value
            no_data_val = srs.GetRasterBand(1).GetNoDataValue()

            # Set GDAL parameters for a merge
            parameters = (
                ["", "-o", demfile]
                + dem_strips
                + [
                    "-co",
                    "COMPRESS=LZW",
                    "-co",
                    "BIGTIFF=YES",
                    "-n",
                    str(no_data_val),
                    "-a_nodata",
                    str(no_data_val),
                ]
            )

            # Merge the strips according to the parameters
            osgeo_utils.gdal_merge.main(parameters)


def bmi_download_dem_extent(
    extent_shpfile, demfile, data_product, dem_product, buffer, render_pixels
):
    """Downloads a DEM (.tif) according to the extent of a shapefile.

    Parameters
    ----------
    extent_shpfile: string
        path of of the extent shapefile
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

    # Load extent shapefile
    shp = gpd.read_file(extent_shpfile).to_crs("EPSG:4326")

    # Get the boundaries of the extent files assuming a WGS84 CRS
    west, east, length, south, north, height = get_shapefile_extent(shp)

    # Buffer those boundaries (expand) using a user-defined buffering percentage
    buff_west, buff_east, buff_south, buff_north = buffer_shapefile(
        west, east, length, south, north, height, buffer
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

    # Send a request to the OpenTopography API with the user-defined parameters
    opentopography_api_download(params, demfile)


if __name__ == "__main__":
    data_product = snakemake_type_exists(
        snakemake.params, "data_product", "/API/globaldem"
    )
    dem_product = snakemake_type_exists(snakemake.params, "dem_product", "NASADEM")
    render_pixels = snakemake_type_exists(snakemake.params, "render_pixels", 1)
    buffer =snakemake_type_exists(snakemake.params, "buffer", 1.0)
    extent_shpfile = snakemake.input["extent_shpfile"]
    demfile = snakemake.output["demfile"]

    bmi_download_dem_extent(
        extent_shpfile, demfile, data_product, dem_product, buffer, render_pixels
    )
