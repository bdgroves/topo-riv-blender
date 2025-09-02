# For more details about the parameters, please refer to the README in the Github repository:
# https://github.com/DOI-USGS/topo-riv-blender/blob/main/README.md

########## data management ##########
# OpenTopography Products: https://portal.opentopography.org/apidocs/
dem_product = "USGS10m"  # usgs 10m
buffer = 0.2  # buffer space around watershed

########## map visualization ##########

# topography colormap
topo_cmap = "custom"
topo_cstops = [
    [70, 140, 131],
    [134, 166, 144],
    [185, 202, 110],
    [237, 224, 132],
    [242, 140, 119],
    [242, 177, 71],
]

########## render hardware ##########
GPU_boolean = False  # if you have a GPU, set this to 1 and this will run faster

########## blender scene ##########
# max dimension of plane
plane_size = 1.0  # meters

########## additional data layers ##########
number_of_layers = 0  # number of additional data layers

########## camera settings ##########
# color management
view_transform = "Filmic"  # 'Standard'
exposure = -2.0  #'default is 0
gamma = 0.5  # default is 1

# camera type
camera_type = "orthogonal"  # orthogonal or perspective

# orthogonal
ortho_scale = 0.8  # when using orthogonal scale, increase to "zoom" out

# perspective
focal_length = 50.0  # mm when using perspective camera, increase to zoom in
shift_x = 0.0  # you may need to shift the camera to center the topo in the frame
shift_y = 0.05  # you may need to shift the camera to center the topo in the frame

# camera location
camera_distance = (
    1.0  # meters from the center (maximum horizontal axis is assumed to be 1 meter)
)
camera_tilt = (
    30.0  # degrees from horizontal: 0.0 is profile view, 90.0 is planform view
)
camera_rotation = 180.0  # camera location degrees clockwise from North: 0.0 is North, 90.0 is East, 180.0, is South, and 270.0 is West.

# depth of field
use_depth_of_field = True
dof_distance = camera_distance  # where the focal plane is
f_stop = 100.0  # affects depths of field, lower for a shallow dof, higher for wide dof

########## sun properties ##########
sun_tilt = 45.0  # degrees from horizontal
sun_rotation = 0.0  # degrees clockwise from North
sun_intensity = 0.5  # sun intensity
sun_strength = 1.0  # sun strength

########## landscape representation ###########
min_res = 2000  # minimum resolution of the height map, larger value increases detail but takes longer to render
number_of_subdivisions = 2000  # number of subdivisions, larger value increases detail but takes longer to render
exaggeration = 5.0  # vertical exaggeration
displacement_method = "DISPLACEMENT"  # "BOTH" #for more exaggerated shadows

########## render settings ##########
res_x = 2000  # x resolution of the render, larger value increases detail but takes longer to render
res_y = 940  # y resolution of the render, larger value increases detail but takes longer to render
samples = 10  # number of samples that decides how "good" the render looks, larger value increases detail but takes longer to render
