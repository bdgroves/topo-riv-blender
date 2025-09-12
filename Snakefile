# Render Parameters File
module_name = "blender_parameters.blender_params_auto"

# Import parameters file
import os, sys, importlib
render_params = module_name.replace(".", os.sep) + ".py"
parameters = importlib.import_module(module_name)

# Define the final outputs
rule all:
    input:
        "phase_2_visualize/out/render_140100.png",
        "phase_2_visualize/out/render_02050302.png",
        "phase_2_visualize/out/render_070400030309.png"

# Downloads a shapefile of the HUC boundary
rule get_wbd:
    params:
        huc_list = lambda wildcards: f"{wildcards.huc_id}".split("and")
    output:
        extent_shpfile = "phase_0_fetch/out/{huc_id}/extent.shp"
    script:
        "phase_0_fetch/src/get_wbd.py"

# Downloads NHD flowline data within extent
rule get_nhd_fl:
    params:
        nhd_type = parameters.nhd_flowline,
        render_pixels = parameters.res_x * parameters.res_y
    input:
        extent_shpfile = "phase_0_fetch/out/{huc_id}/extent.shp"
    output:
        nhd_shpfile = "phase_0_fetch/out/{huc_id}/rivers.shp"
    script:
        "phase_0_fetch/src/get_nhd.py"

# Downloads NHD waterbody data within extent
rule get_nhd_wb:
    params:
        nhd_type = parameters.nhd_waterbody,
        render_pixels = parameters.res_x * parameters.res_y 
    input:
        extent_shpfile = "phase_0_fetch/out/{huc_id}/extent.shp"
    output:
        nhd_shpfile = "phase_0_fetch/out/{huc_id}/waterbody.shp"
    script:
        "phase_0_fetch/src/get_nhd.py"

# Downloads a digital elevation model within extent 
rule download_dem:
    params:
        dem_product = parameters.dem_product,
        buffer = parameters.buffer,
        render_pixels = parameters.res_x * parameters.res_y
    input:
        extent_shpfile = "phase_0_fetch/out/{huc_id}/extent.shp"
    output:
        demfile = "phase_0_fetch/out/{huc_id}/dem.tif"
    script:
        "phase_0_fetch/src/download_dem.py"

# Creates a grayscale height map and texture map for Blender to render.
# The height map is used to determine how high the landscape should be 
# in the render, and the texture map determines the color of the landscape.
rule create_heightmap_texturemap:
    params:
        buffer = parameters.buffer,
        topo_cmap = parameters.topo_cmap,
        background_color = parameters.background_color,
        wall_color = parameters.wall_color,
        river_color = parameters.river_color
    input:
        extent_shpfile = "phase_0_fetch/out/{huc_id}/extent.shp",
        demfile =  "phase_0_fetch/out/{huc_id}/dem.tif",
        flowlines_shpfile = "phase_0_fetch/out/{huc_id}/rivers.shp",
        waterbody_shpfile = "phase_0_fetch/out/{huc_id}/waterbody.shp"
    output:
        dimensions_file = "phase_1_process/out/{huc_id}/dimensions.npy",
        heightmap_file = "phase_1_process/out/{huc_id}/heightmap.png",
        texturemap_file = "phase_1_process/out/{huc_id}/texturemap.png",
        apronmap_file = "phase_1_process/out/{huc_id}/apronmap.png"
    script:
        "phase_1_process/src/process.py"

# Using the height map and texture map, Blender sets up the scene 
# (topography, lighting, and camera) and renders a photorealistic image.
rule render:
    input:
        dimensions_file = "phase_1_process/out/{huc_id}/dimensions.npy",
        heightmap_file = "phase_1_process/out/{huc_id}/heightmap.png",
        texturemap_file = "phase_1_process/out/{huc_id}/texturemap.png",
        apronmap_file = "phase_1_process/out/{huc_id}/apronmap.png",
        blender_params = render_params
    output:
        output_file = "phase_2_visualize/out/render_{huc_id}.png"
    shell:
        "Blender -b -P phase_2_visualize/src/render.py -- "
        "{input.blender_params} "
        "{input.dimensions_file} {input.heightmap_file} {input.texturemap_file} {input.apronmap_file} NULL "
        "{output.output_file} "