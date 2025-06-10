import os, csv
import numpy as np
from osgeo import gdal
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from shapely.geometry import Polygon
from scipy.interpolate import RegularGridInterpolator
from astropy.convolution import Gaussian2DKernel, convolve

from phase_0_fetch.src.download_dem import get_shapefile_extent


def get_raster_extent(dataset):
    """Gets extent coodinates of the raster image

    Parameters
    ----------
    dataset: gdal data set
        gdal data set holding the dem

    Returns
    -------
    ulx: float
        upper left x-coordinate
    lrx: float
        lower right x-coordinate
    lry: float
        lower right y-coordinate
    uly: float
        upper left y-coordinate

    """

    # Get the raster extent from gdal
    ulx, xres, xskew, uly, yskew, yres = dataset.GetGeoTransform()

    # Calculate the lower right x-coordinate
    lrx = ulx + (dataset.RasterXSize * xres)

    # Calculate the lower right y-coordinate
    lry = uly + (dataset.RasterYSize * yres)

    return ulx, lrx, lry, uly


def clean_border(ax, west, east, south, north, edgecolor, linewidth, crs):
    """Adds a black line around the extent to remove anomalous values from edge.
    Because matplotlib saves the texture and height maps as images, the data
    shape and image pixel shape may not fully line up. This causes there to be a
    white line on the border of the images. This function cleans that up.

    Parameters
    ----------
    ax: matplotlib axis
        axis handle to modify
    west: float
        west long coordinate
    east: float
        east long coordinate
    south: float
        south lat coordinate
    north: float
        north lat coordinate
    edgecolor: rgba list
        color of the edge
    linewidth: float
        width of the line to clean up the edges
    crs: string
        coordinate reference system

    Returns
    -------
    none

    """

    # Create a geopandas dataframe based on the extent
    clean_polygon = gpd.GeoDataFrame(
        index=[0],
        crs=crs,
        geometry=[
            Polygon(
                zip(
                    [west, west, east, east, west],
                    [north, south, south, north, north],
                )
            )
        ],
    )

    # Plot this as a line around the extent
    clean_polygon.plot(
        ax=ax,
        facecolor="none",
        edgecolor=edgecolor,
        linewidth=linewidth,
        capstyle="round",
        zorder=1000,
    )


def fig_setup(fignum, mask, min_res, facecolor):
    """Sets up the matplotlib figure

    Parameters
    ----------
    fignum: int
        figure number
    mask: geopandas dataframe
        geopandas dataframe that contains the geometry of the height + texture map
    min_res: int
        minimum pixel resolution of the short side of the image
    facecolor: rbga list
        color of the background of the image

    Returns
    -------
    fig: matplotlib figure
        matplotlib figure of the height and texture map
    ax: matplotlib axis
        axis handle to modify
    dpi: float
        dots per inch

    """

    # Default fig size in inches, 10 in
    fig_width = 10.0

    # Get the boundaries of the extent files assuming a map_crs
    west, east, mask_length, south, north, mask_height = get_shapefile_extent(mask)

    # Calculate the dpi based on the short side and the minimum resolution specified
    if mask_height / mask_length >= 1.0:
        dpi = min_res / fig_width
    else:
        dpi = min_res * mask_length / mask_height / fig_width

    # Make figure
    fig = plt.figure(
        fignum,
        figsize=(fig_width, fig_width * mask_height / mask_length),
        facecolor=facecolor,
    )

    # Make the Axis
    ax = fig.add_axes([0, 0, 1, 1])

    # Set bounds
    ax.set_aspect("equal")
    ax.set_xlim(west, east)
    ax.set_ylim(south, north)

    # Remove axis labels and ticks
    ax.axis("off")

    return fig, ax, dpi


def project(file, map_crs, tmp_dir):
    """Projects a raster file to another crs

    Parameters
    ----------
    file: string
        path of geotiff
    map_crs: string
        coordinate reference system to use for the maps
    tmp_dir: string
        path of the tmp directory

    Returns
    -------
    proj_array: numpy array
        projected numpy array
    proj_ds: projected gdal dataset
        axis handle to modify
    min_array: float
        minimum value of proj_array
    max_array: float
        maximum value of proj_array

    """

    # Open the geotiff
    ds = gdal.Open(file)

    # Read in the data as a numpy array
    array = np.array(ds.GetRasterBand(1).ReadAsArray())

    # Reproject dem
    proj_ds = gdal.Warp(
        tmp_dir + "/proj_" + os.path.basename(file),
        ds,
        dstSRS=map_crs,
        resampleAlg="bilinear",
    )

    # Read in the projected data as a numpy array
    proj_array = np.array(proj_ds.GetRasterBand(1).ReadAsArray())

    # Get min and max of array, excluding NaNs
    NaN = proj_ds.GetRasterBand(1).GetNoDataValue()
    if NaN != None:
        if np.isnan(NaN):
            min_array = np.nanmin(proj_array[~np.isnan(proj_array)])
            max_array = np.nanmax(proj_array[~np.isnan(proj_array)])
        else:
            min_array = np.nanmin(proj_array[proj_array != NaN])
            max_array = np.nanmax(proj_array[proj_array != NaN])
    else:
        min_array = np.nanmin(proj_array)
        max_array = np.nanmax(proj_array)

    return proj_array, proj_ds, min_array, max_array


def contour(ds, dem, coarsen, stddev):
    """Creates a smoothed dem that will look better as a contour

    Parameters
    ----------
    ds: gdal dataset
        gdal dataset holding the dem
    dem: numpy array
        dem values in numpy array
    coarsen: float
        divisor that coarsens the data for contouring
    stdddev: float
        standard deviation in grid-lengths for applying the convolution filter

    Returns
    -------
    astropy_conv: numpy array
        convoluted dem values in numpy array

    """

    # Get extent from raster
    rast_extent = list(get_raster_extent(ds))

    # Make a spatial grid for the raster
    rast_x = np.linspace(rast_extent[0], rast_extent[1], dem.shape[1])
    rast_y = np.linspace(rast_extent[3], rast_extent[2], dem.shape[0])

    # Initialize the interpolator
    interpol = RegularGridInterpolator((rast_y, rast_x), dem, method="linear")

    # Make coarser spatial grid for the raster
    xx = np.linspace(rast_extent[0], rast_extent[1], int(dem.shape[1] / coarsen))
    yy = np.linspace(rast_extent[3], rast_extent[2], int(dem.shape[0] / coarsen))
    Y, X = np.meshgrid(yy, xx, indexing="ij")

    # Initilize the gaussian kernal
    kernel = Gaussian2DKernel(x_stddev=stddev)

    # Use the guassian convolution on the coarsened data to create a smooth contoured topography
    astropy_conv = convolve(interpol((Y, X)), kernel)

    return astropy_conv


def make_labels_file(labels_shpfile, labels_file, extent, map_crs):
    """Sets up labels file for rendering.

    Parameters
    ----------
    map_crs: string
        coordinate reference system to use for the maps
    extent: list
        list of west, east, south, north boundaries
    labels_shpfile: string
        path of the labels shapefile
    labels_file: string
        path to the labels file with a list of the coordinates and labels

    Returns
    -------
    none

    """

    # Load labels shapefile
    labels_shp = gpd.read_file(labels_shpfile).to_crs(map_crs)

    labels = [[-9999.0, -9999.0, ""] for i in range(0, len(labels_shp))]

    for i in range(0, len(labels_shp)):
        labels[i] = [
            (labels_shp.geometry[i].x - 0.5 * (extent[0] + extent[1]))
            / (extent[1] - extent[0]),
            (labels_shp.geometry[i].y - 0.5 * (extent[2] + extent[3]))
            / (extent[3] - extent[2]),
            labels_shp["labels"][i],
            labels_shp["p_colors"][i],
            labels_shp["l_colors"][i],
        ]

    # Write the 2D list to a CSV file
    with open(labels_file, mode="w", newline="") as file:
        writer = csv.writer(file)
        writer.writerows(labels)


def make_custom_cmap(cstops, nstops):
    """Creates a custom matplotlib cmap based on color stops. Defaults to evenly spaced color stops unless specified otherwise.

    Parameters
    ----------
    map_crs: string
        coordinate reference system to use for the maps
    extent: list
        list of west, east, south, north boundaries
    labels_shpfile: string
        path of the labels shapefile
    labels_file: string
        path to the labels file with a list of the coordinates and labels

    Returns
    -------
    cmap: matplotlib cmap
        the custom matplotlib color map

    """
    # make evently spaced n stops if not specified
    if nstops == []:
        number_of_stops = len(cstops)
        nstops = [
            int(256.0 * float(i) / float(number_of_stops - 1))
            for i in range(0, number_of_stops)
        ]

    # create empty color ramp
    vals = np.ones((256, 4))

    # linearly interpolate between color stops
    for i in range(0, 3):
        for j in range(1, len(cstops)):
            vals[nstops[j - 1] : nstops[j], i] = np.linspace(
                cstops[j - 1][i] / 255.0,
                cstops[j][i] / 255.0,
                int(nstops[j] - nstops[j - 1]),
            )

    # return the matplotlib colormap
    return mcolors.ListedColormap(vals)
