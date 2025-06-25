# For more details about the parameters, please refer to the README in the Github repository:
# https://github.com/DOI-USGS/topo-riv-blender/blob/main/README.md

########## data management ##########
# OpenTopography Products: https://portal.opentopography.org/apidocs/
data_product = (
    "/API/globaldem"  # Within US use'/API/usgsdem', Global use '/API/globaldem'
)
dem_product = "auto"  # determine automatically
buffer = 0.2  # buffer space around watershed

########## map visualization ##########
background_color = [0.1, 0.1, 0.1, 1.0]
wall_color = [0.1, 0.0667, 0.0333, 1.0]

########## Label Parameters ##########
pin_sides = 4  # sides of pin
pin_height = 0.015  # height of pin
pin_radius1 = 0.0  # bottom of pin
pin_radius2 = 0.00667  # top of pin
pin_label_buffer = 0.05  # space between pin and label
label_scale = 0.035  # size of label
label_font = (
    "Barlow Condensed"  # Specify Google Fonts, Browse here: https://fonts.google.com/
)
shadow_boolean = False  # alllow shadows for the label

########## render hardware ##########
GPU_boolean = False  # if you have a GPU, set this to 1 and this will run faster

########## blender scene ##########
# max dimension on plane
plane_size = 1.0  # meters

########## additional data layers ##########
number_of_layers = 0  # number of additional data layers

########## camera settings ##########
# color management
view_transform = "Filmic"  # 'Standard'
exposure = 0.0  # default is 0
gamma = 0.85  # default is 1

# camera type
camera_type = "orthogonal"  # orthogonal or perspective

# orthogonal
ortho_scale = 1.2  # when using orthogonal scale, increase to "zoom" out

# perspective
focal_length = 50.0  # mm when using perspective camera, increase to zoom in
shift_x = 0.0  # you may need to shift the camera to center the topo in the frame
shift_y = 0.0  # you may need to shift the camera to center the topo in the frame

# camera location
camera_distance = (
    1.0  # meters from the center (maximum horizontal axis is assumed to be 1 meter)
)
camera_tilt = (
    30.0  # degrees from horizontal: 0.0 is profile view, 90.0 is planform view
)
camera_rotation = 180.0  # camera location degrees clockwise from North: 0.0 is North, 90.0 is East, 180.0, is South, and 270.0 is West.

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
min_res = 4000  # minimum resolution of the heightmap
number_of_subdivisions = 2000  # number of subdivisions, more increases the detail
exaggeration = "auto"  # vertical exaggeration
displacement_method = "DISPLACEMENT"  # "BOTH" #for more exaggerated shadows

########## render settings ##########
res_x = 3000  # x resolution of the render
res_y = 1500  # y resolution of the render
samples = 10  # number of samples that decides how "good" the render looks. more is better but takes longer
