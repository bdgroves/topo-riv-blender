# For more details about the parameters, please refer to the README in the Github repository:
# https://github.com/DOI-USGS/topo-riv-blender/blob/main/README.md

########## data management ##########
# OpenTopography Products: https://portal.opentopography.org/apidocs/
dem_product = "SRTM15Plus"  # "GEBCOIceTopo"  # determine automatically
buffer = 0.2  # buffer space around watershed

########## map visualization ##########
background_color = [0.1, 0.1, 0.1, 1.0]
wall_color = [0.2, 0.133, 0.0667, 1.0]

# topography colormap
topo_cmap = "custom"
topo_cstops = [[42, 64, 42], [88, 70, 56], [200, 200, 200]]
topo_nstops = [0, 200, 256]

# "ocean" parameters
oceanfloor_cmap = "bone"
ocean_elevation = 455.5  # m

########## render hardware ##########
GPU_boolean = False  # if you have a GPU, set this to 1 and this will run faster

########## blender scene ##########
# max dimension on plane
plane_size = 1.0  # meters

########## additional data layers ##########
number_of_layers = 1  # number of additional data layers
layers_alpha = [0.75]

########## camera settings ##########
# color management
view_transform = "Filmic"  # 'Standard'
exposure = 0.0  # default is 0
gamma = 0.85  # default is 1

# camera type
camera_type = "orthogonal"  # orthogonal or perspective

# orthogonal
ortho_scale = 1.3  # when using orthogonal scale, increase to "zoom" out

# perspective
focal_length = 50.0  # mm when using perspective camera, increase to zoom in
shift_x = 0.0  # you may need to shift the camera to center the topo in the frame
shift_y = 0.0  # you may need to shift the camera to center the topo in the frame

# camera location
camera_distance = (
    1.0  # meters from the center (maximum horizontal axis is assumed to be 1 meter)
)
camera_tilt = (
    45.0  # degrees from horizontal: 0.0 is profile view, 90.0 is planform view
)
camera_rotation = 150.0  # camera location degrees clockwise from North: 0.0 is North, 90.0 is East, 180.0, is South, and 270.0 is West.

# depth of field
use_depth_of_field = False
dof_distance = camera_distance  # where the focal plane is
f_stop = 100.0  # affects depths of field, lower for a shallow dof, higher for wide dof

########## sun properties ##########
sun_tilt = 20.0  # degrees from horizontal
sun_rotation = 315.0  # degrees clockwise from North
sun_intensity = 0.5  # sun intensity
sun_strength = 1.0  # sun strength

########## landscape representation ###########
min_res = 4000  # minimum resolution of the heightmap, larger value increases detail but takes longer to render
number_of_subdivisions = 2000  # number of subdivisions, larger value increases detail but takes longer to render
exaggeration = 5.0  # vertical exaggeration
displacement_method = "DISPLACEMENT"  # "BOTH" #for more exaggerated shadows

########## render settings ##########
res_x = 2000  # x resolution of the render, larger value increases detail but takes longer to render
res_y = 1500  # y resolution of the render, larger value increases detail but takes longer to render
samples = 10  # number of samples that decides how "good" the render looks, larger value increases detail but takes longer to render
