########## data management ##########
# OpenTopography Products: https://portal.opentopography.org/apidocs/
data_product = (
    "/API/usgsdem"  # Within US use'/API/usgsdem', Global use '/API/globaldem'
)
dem_product = "auto"  # determine automatically (see options here: https://portal.opentopography.org/apidocs/#/)
# NHD products
nhd_flowline = "flowline_auto"  # use auto, mr, or hr
nhd_waterbody = "waterbody_auto"  # use auto, mr, or hr
# Data Buffer
buffer = 1.0  # buffer space around watershed

########## map visualization ##########
background_color = [0.5, 0.5, 0.5, 1.0]  # rgba color of the background
wall_color = [0.2, 0.133, 0.0667, 1.0]  # rgba color of the wall edge
river_color = [0.1294, 0.2275, 0.3608, 1.0]  # rgba color of the rivers
topo_cmap = "copper"  # color map of the topography

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
exposure = 0.0  # 'default is 0
gamma = 0.65  # default is 1

# camera type
camera_type = "orthogonal"  # orthogonal or perspective

# orthogonal
ortho_scale = 1.0  # when using orthogonal scale, increase to "zoom" out

# perspective
focal_length = 50.0  # mm when using perspective camera, increase to zoom in

# camera location
camera_distance = (
    1.0  # meters from the center (maximum horizontal axis is assumed to be 1 meter)
)
camera_tilt = (
    45.0  # degrees from horizontal: 0.0 is profile view, 90.0 is planform view
)
camera_rotation = 180.0  # camera location degrees CW from North: 0.0 is North, 90.0 is East, 180.0, is South, and 270.0 is West.
shift_x = 0.0  # you may need to shift the camera to center the topo in the frame
shift_y = 0.0  # you may need to shift the camera to center the topo in the frame

# depth of field
use_depth_of_field = False
dof_distance = camera_distance  # where the focal plane is
f_stop = 100.0  # affects depths of field, lower for a shallow dof, higher for wide dof

########## sun properties ##########
sun_tilt = 20.0  # degrees from horizontal
sun_rotation = 315.0  # degrees CW from North
sun_intensity = 0.5  # sun intensity
sun_strength = 1.0  # sun strength

########## landscape representation ###########
min_res = 2000  # minimum resolution of the heightmap
number_of_subdivisions = 2000  # number of subdivisions, more increases the detail
exaggeration = "auto"  # vertical exaggeration
displacement_method = "DISPLACEMENT"  # "BOTH" #for more exaggerated shadows

########## render settings ##########
res_x = 2000  # x resolution of the render
res_y = 2000  # y resolution of the render
samples = 10  # number of samples that decides how "good" the render looks. more is better but takes longer
