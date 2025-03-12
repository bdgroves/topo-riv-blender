import pandas as pd
import geopandas as gpd

def make_shpfile(
        coords,
        labels,
        pin_colors,
        label_colors,
        crs,
        shpfile):
    """Creates a ESRI shapefile from a list of coordinates and labels.

    Parameters
    ----------
    coords: 2D list of floats
        list of coordinates of the shapefile
    labels: 1D list of strings
        list of labels for corresponding labels
    crs: string
        coordinate reference system
    
    Returns
    -------
    none

    """
        
    x_coords = []
    y_coords = []

    for i in range(0,len(labels)):
        x_coords += [coords[i][1]]
        y_coords += [coords[i][0]]

    df = pd.DataFrame(
        {
            "labels": labels,
            "p_colors": pin_colors,
            "l_colors": label_colors,
            "x_coords": x_coords,
            "y_coords": y_coords,
        }
    )

    gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.x_coords, df.y_coords), crs=crs)

    gdf.to_file(shpfile, driver='ESRI Shapefile')

if __name__ == "__main__":
    coords = snakemake.params["coords"]
    labels = snakemake.params["labels"]
    pin_colors = snakemake.params["pin_colors"]
    label_colors = snakemake.params["label_colors"]
    crs = snakemake.params["crs"]
    shpfile = snakemake.output["shpfile"]
    make_shpfile(
        coords,
        labels,
        pin_colors,
        label_colors,
        crs,
        shpfile,
    )