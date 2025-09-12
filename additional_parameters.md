## Additional Parameters
The process function contains more parameters than are specified in the `blender_params_auto.py` file. This was done to avoid overwhelming you with all the options. The following list contains parameters that I have personally used for other visualizations. If there are other things that you want to visualize, you will need to modify the `process.py` or `render.py` script. This workflow is meant to get you started on that.

#### How to add additional parameters

##### Snakefile
The default Snakefile currently contains a rule named: `create_heightmap_texturemap`. Below are the parameters currently being used by the `process.py` script.
```python
rule create_heightmap_texturemap:
    params:
        buffer = parameters.buffer,
        topo_cmap = parameters.topo_cmap,
        background_color = parameters.background_color,
        wall_color = parameters.wall_color,
        river_color = parameters.river_color
```
However, within `process.py`, you will see that there are many `arguments` in the `setup_blender_data` function. Many of these use a default value unless specified by the user. 

##### What to modify
If you would like to add one of these parameters, add the parameter to the `blender_params_auto.py` file (the default workflow) and `param_Y = parameters.param_Y` to the Snakefile. The parameters file (e.g., `blender_params_auto.py`) should have the following line added:
```python
param_Y = VALUE
```
and the Snakefile should now look like this:
```python
rule create_heightmap_texturemap:
    params:
        buffer = parameters.buffer,
        topo_cmap = parameters.topo_cmap,
        background_color = parameters.background_color,
        wall_color = parameters.wall_color,
        river_color = parameters.river_color,
        param_Y = parameters.param_Y
```
##### List of parameters
###### map
- `map_crs`: [ESPG](https://epsg.io/) coordinate reference system to use for the maps
    - Type: *string*
    - Default: `NULL`, gets closest UTM projection
- `background_color`: rgba color of the background
    - Type: *rgba list*
    - Default: `[0.5, 0.5, 0.5, 1.0]`
- `min_res`: the minimum pixel resolution for the height map and texture map
    - Type: *integer*
    - Default: `2000`
###### topography
- `mask_boolean`: option that masks out all topography outside of the extent shapefile
    - Type: *boolean*
    - Default: `True`
- `custom_min_dem`: lowest value for DEM color map
    - Type: *float*
    - Default: `NULL`, finds minimum based on the DEM
- `wall_color`: rgba color of the wall of the topography
    - Type: *rgba list*
    - Default: `[0.2, 0.133, 0.0667, 1.0]`
- `wall_thickness`: thickness of outline around topography
    - Type: *float*
    - Default: `1.0`
    - Note: *This is in units of points (1/72 in). The width of the texture map is set to 10 inches.*
- `contour_boolean`: option that changes the topography into contours
    - Type: *boolean*
    - Default: `False`
- `contour_levels`: number of contour levels
    - Type: *integer*
    - Default: `20`
- `convolution_coarsen`: coarsens data by this factor
    - Type: *float*
    - Default: `5.0`
- `convolution_stdddev`: standard deviation in grid-lengths for applying the convolution filter
    - Type: *float*
    - Default: `2.0`
###### river
- `river_color`: rgba color of the rivers
    - Type: *rgba list*
    - Default: `[0.1294, 0.2275, 0.3608, 1.0]`
- `river_width`: thickness of river lines
    - Type: *float*
    - Default: `auto`, determines this value automatically
    - Note: *This is in units of points (1/72 in). The width of the texture map is set to 10 inches.*
###### ocean
- `ocean_boolean`: option that adds for ocean layer at ocean elevation
    - Type: *boolean*
    - Default: `False`
    - Note: *In order to use this option, you need to set the `number_of_layers` parameter to `1` or greater in the blender-parameters (e.g., `blender_params_lake.py`). See the `Snakefile_lake` workflow and `blender_params_lake.py` as an example. You will also need to set `layers_alpha` in the blender-parameters to set the alpha value of the ocean in blender.*
- `ocean_elevation`: elevation of the ocean layer
    - Type: *float*
    - Default: `0.0`
- `ocean_color`: rgba color of the ocean
    - Type: *list of floats*
    - Default: `[0.1294, 0.2275, 0.3608, 1.0]`
- `oceanfloor_cmap`: matplotlib cmap of the ocean floor
    - Type: *string*
    - Default: `"Greys"`
###### additional layers - *Experimental, more documentation to come*
- `layers_cmap`: list of matplotlib cmaps of each additional layer
    - Type: *list of string*
    - Example: `["Greys", "viridis"]`, uses the `Greys` color map for the first additional layer and the `viridis` color map for the second additional layer
- `layers_vlim`: list of two element long lists
    list of the min and max values for the color maps for each layer
    - Example: `["NULL", [0.0,1.0]]`, uses the minimum and maximum of the first additional layer and user-specified minimum value of 0 to maximum value of 1 for the second additional later

##### Inputs
When adjusting the parameters, you may need to adjust the inputs.
- `extent_shpfile`: path of the extent shapefile
    - Type: *string*
    - Note: *This shapefile can be a HUC boundary, the extent of a bounding box, or can be user-specified.*
- `demfile`: path of the DEM raster file
    - Type: *string*
    - Note: *The workflows will download this file for you, but you can supply custom DEMs here as well.*
- `waterbody_shpfile`: path of the waterbody shapefile
    - Type: *string*
    - Note: *The example workflow in the default Snakefile uses the NHD waterbodies dataset. You can also create your own.*
- `flowlines_shpfile`: path of the flowline shapefile
    - Type: *string*
    - Note: *The example workflow in the default Snakefile uses the NHD flowlines dataset. You can also create your own.*
    - Type: *string*
- `labels_shpfile`: path of the labels shapefile
    - Type: *string*
    - Note: *The `Snakefile_global` workflow shows an example of how to create this shape file.*    
###### *Experimental, more documentation to come*
- `layerfiles` list of paths of the labels shapefile
    - Type: *list of strings*
    - Note: *This is a list of additional layers (elevations) to be renders. For example, if you wanted a snow layer, you would include the elevation of the snow layers surface and set its color to white by using the `Greys` color map and setting `layer_vlim` to `[-10000,-9999]` to make any positive value is white. This is a little wonky, so you now see why it's experimental.*
- `aerialfiles`: list of the paths of the aerial imagery r, g, b bands
    - Type: *list of strings*
    - Default: `NULL`
    - Note: *Manual downloading of aerial imagery required for now.*
##### Outputs
When adjusting the parameters, you may need to adjust the outputs.
- `dimensions_file`: path to the file containing the length, width, and height of the landscape
    - Type: *string*
- `heightmap_file`: path to the image file showing the height map
    - Type: *string*
- `texturemap_file`: path to the image file showing the texture map
    - Type: *string*
- `apronmap_file`:  path to the image file showing the texture map of the background
    - Type: *string*
- `labels_file`: path to the labels file with a list of the coordinates and labels   
    - Type: *string*
###### *Experimental, more documentation to come*
- `heightmap_layerfiles`: list of paths to the height maps for the additional layers
    - Type: *list of strings*
- `texturemap_layerfiles`:  list of paths to the texture maps for the additional layers
    - Type: *list of strings*