from RKJ_utils import *
from yade import *
from yade import plot
from math import *

import json
import sys, os

sys.path.append(".")

from brush_vs_target_lib import *
# ---------------------------------------------------------------------------- input parameter
simulation_mode = "batch"
input_data_start_index = 0

if simulation_mode == 'single':
    ID  = 1
    file_name = "Exp_20_04_C001H001S0001"
    x_pos = 13.9e-3

elif simulation_mode == 'batch':
    readParamsFromTable(
        ID  = 1,
        file_name = "Exp_20_04_C001H001S0001",
        x_pos = 13.9e-3
        )
    from yade.params.table import *

# ----------------------------------------------------- Experiment data
bristle_tip_pos_y  =  x_pos
bristle_tip_pos_z  =  63.5e-3


with open(
    'conditioned_data/' + file_name + '.json', 
    'r') as jsonfile: 
    # data.json, target_rotation_angle
    exp_data = json.load(jsonfile)

with open('damping_model.json', 'r') as jsonfile: 
    # data.json, target_rotation_angle
    damping_data = json.load(jsonfile)

fps = exp_data["fps"]
dt  = 1/fps 
# ----------------------------------------------------- brush
bristle_radius_root    = 1e-3 
bristle_radius_diff    = 0.5e-3
bristle_length         = 35e-3
bristle_no_of_segments = 12

bristle_young    = 3.2e6
bristle_density  = 1250
bristle_poisson  = 0.48
bristle_friction = radians(44)
bristle_damping  = 0

bristles_dx      = 2.4375e-3
bristles_dy      = 2.4375e-3
bristle_x_no     = 17
bristle_y_no     = 5

brush_height = 12.0e-3
covar        = 0e-7
clearance    = 1e-5
# ----------------------------------------------------- target
target_mass = 69.0e-3
target_side = 0.0602
p_radius    = 1e-3

target_poisson  = 0.3
target_friction = radians(24)
target_young    = bristle_young

target_pos         = [0.0, 0.0, 70e-3]
target_pfacet_mass = 0.0

target_shape      = '60mm.gts'
geometry_location = './geometry/'
target_shape      = geometry_location + target_shape 

target_axle_radius  = 7.9e-3 
target_axle_mass    = 55.8e-3
# ---------------------------------------------------------------------------- Materials
brush_int_mat  = 'cyl_int_mat'
brush_ext_mat  = 'cyl_ext_mat'

# ----------------------------------------------------- brush
add_official_brush_material(
    young   = bristle_young, 
    poisson = bristle_poisson, 
    density = bristle_density, 
    friction_Angle   = bristle_friction,
    friction_Angle_2 = radians(30),
    ext_mat_name   = brush_ext_mat, 
    int_mat_name   = brush_int_mat,
    )

# ---------------------------------------------------------------------------- Engines
O.engines = [
                ForceResetter(),

                InsertionSortCollider([
                    Bo1_GridConnection_Aabb(),
                    Bo1_PFacet_Aabb(),
                    Bo1_Sphere_Aabb(),
                ]),

                InteractionLoop(
                    [
                        Ig2_PFacet_PFacet_ScGeom(),
                        Ig2_GridConnection_GridConnection_GridCoGridCoGeom(),
                        Ig2_GridNode_GridNode_GridNodeGeom6D(),
                        Ig2_GridConnection_PFacet_ScGeom(),
                        Ig2_Sphere_PFacet_ScGridCoGeom(),
                    ],
                    [
                        Ip2_FrictMat_FrictMat_FrictPhys(),
                        Ip2_CohFrictMat_CohFrictMat_CohFrictPhys(
                            setCohesionNow = True, 
                            setCohesionOnNewContacts = False
                            ),
                    ],
                    [
                        Law2_GridCoGridCoGeom_FrictPhys_CundallStrack(),
                        Law2_ScGeom_FrictPhys_CundallStrack(),
                        Law2_ScGridCoGeom_FrictPhys_CundallStrack(),
                        Law2_ScGeom6D_CohFrictPhys_CohesionMoment(),
                    ],
                ),
                PyRunner(command = "damper()", iterPeriod = 1),
                NewtonIntegrator(gravity = (0,0,-10), damping = 0.05)
            ]
# ---------------------------------------------------------------------------- objects 
# # ----------------------------------------------------- target
target = Target(
        target_shape,
        p_radius,
        target_pos,
        target_young,
        target_poisson,
        target_friction,
        target_mass,
        target_pfacet_mass,
        target_side,
        target_axle_mass,
        target_axle_radius
        )

# ----------------------------------- getting the clump orientations
# target initial orientation
target_init_orientation_quarternion = O.bodies[target.clump_id].state.ori 
target_init_orientation_matrix      = target_init_orientation_quarternion.toRotationMatrix


# ----------------------------------------------------- brush
bristle_radius_tip = bristle_radius_root - bristle_radius_diff
brush_x_length     = bristles_dx * bristle_x_no
brush_y_length     = bristles_dy * bristle_y_no


brush_position = [
    0.0, 
    bristle_tip_pos_y + bristle_length, 
    bristle_tip_pos_z + brush_height/2.0
    ]
orientation    = Quaternion((0,1,0), 0) * Quaternion((0,0,1), -pi/2) * Quaternion((0,1,0), pi/2)


brush_x_no_density = 1 / bristles_dx
brush_y_no_density = 1 / bristles_dy
bristle_tip_spread_covariance = [[covar, 0.0],[0.0, covar]]

brush = ccp_brush(
        bristle_radius   = bristle_radius_tip,
        bristle_radius_2 = bristle_radius_root, 
        bristle_length   = bristle_length,
        no_of_segment    = bristle_no_of_segments, 

        bristle_external_material = brush_ext_mat,
        bristle_internal_material = brush_int_mat,        

        x_lenght = brush_x_length, 
        y_length = brush_y_length, 
        x_no_density = brush_x_no_density, 
        y_no_density = brush_y_no_density, 

        brush_base_pos = brush_position,
        clearance      = clearance,        
        orientation    = orientation,

        root_fixed_Q = True, # If true the roots will be fixed to the ground
        bristle_tip_spread_covariance = bristle_tip_spread_covariance,# 'default',
    )
node_ids_array, cyl_ids_array, root_node_ids_array = brush.generate()


# ----------------------------------------------------------------------------- rotation data
# getting the orientation, relative to be first angle
exp_observed_ori = -(np.array(exp_data["filetered_orientation"])) * pi / 180
exp_observed_ori =   exp_observed_ori[input_data_start_index:]




# generating the time stamps. JSON data has problems with time. 
# At larger values of time the dt is small so it was rounded.
# This means at large values some consecutive time stamps are equal 
len_exp_data  = len(exp_observed_ori)
exp_time_data = np.arange(len_exp_data)
exp_time_data = exp_time_data * dt




# the target starts at an angle which could be in contact with the brush
# If this is the start then there will be a collision before the start of the sim
# So we need to slowly set it to that location
offset_time      = 1e-2
exp_time_data    = [0.0] + list(offset_time + exp_time_data)
exp_observed_ori = np.array([0.0] + list(exp_observed_ori))


stop_mover_index               = int(exp_data["start_idx"])
exp_observed_ori_lerper        = DataInterpolator(exp_time_data, list(exp_observed_ori))
index_where_target_driven_ends = stop_mover_index - input_data_start_index + 1
time_when_target_driven_ends   = exp_time_data[index_where_target_driven_ends]




# Creating a data linear interpolator for the raw data. This will be compared
# against the simulations 
exp_observed_raw_ori = -(np.array(exp_data["target_orientation"])) * pi / 180
exp_observed_raw_ori =   exp_observed_raw_ori[input_data_start_index:]

exp_observed_raw_ori =   np.array([0.0] + list(exp_observed_raw_ori))
exp_observed_raw_ori_lerper = DataInterpolator(exp_time_data, list(exp_observed_raw_ori))
# ----------------------------------------------------------------------------- Additional engines
nextStep         = 'move'
motionStartTime  = 0.0

O.engines += [
    PyRunner(
        command    = 'mover()',
        iterPeriod = 1, 
        label      = 'moverEngine'),

    PyRunner(
        command    = 'plotter()',
        iterPeriod = 1000,
        label      = 'PlotterEngine'
    )
]

# ----------------------------------------------------- plot
plot.plots = {'t': 'theta', 't2':('wx_og', 'wx')}


# ----------------------------------------------------- damper
def damper():
    """
    Adds the damping force caluclated from the damping model 
    """
    # getting orientation **********************
    state   = O.bodies[target.clump_id].state
    ori     = state.ori
    ori_mat = ori.toRotationMatrix()

    # getting the inertia in global space ******
    mmoi_vec          = state.inertia
    mmoi_local_space  = getInertiaMatrixFromVec(mmoi_vec)
    mmoi_global_space = getGlobalInertiaMatrix(ori_mat, mmoi_local_space)
    
    # getting the velocity *********************
    omega = state.angVel
    wx    = omega[0]

    # acceleration and torque calculation ******
    alpha_x = dampingModelFunction(wx)
    tau_x   =   alpha_x * mmoi_global_space[0][0]
    tau     =   Vector3([tau_x, 0, 0])

    # appliying torque *************************
    O.forces.addT(id = target.hole_sphere_id, t = tau)

# ***************************************************** dampingModelFunction
def dampingModelFunction(wx):
    """
    Returns the damping acceleration using the model created
    from the experiment

    Args:
        wx (float): x component of angular velocity

    Returns:
        float: x component of angular Acceleration
    """
    return (
        damping_data["damping_slope"]*wx 
      - damping_data["friction_const"]*signum(wx)
      )

# ***************************************************** getInertiaMatrixFromVec
def getInertiaMatrixFromVec(mmoi_vec):
    """
    Returns the inertia matrix from the MMOI vector. 
    This vector represents the diagonal of the principal inertia matrix.
    This 

    Args:
        mmoi_vec (Vector3): Inertia vector

    Returns:
        Matrix3: Inertia matrix
    """
    mmoi_local_space  = Matrix3(
        mmoi_vec[0], 0, 0,
        0, mmoi_vec[1], 0,
        0, 0, mmoi_vec[2]
        )
    return mmoi_local_space

# ***************************************************** getMMOIGlobal
def getGlobalInertiaMatrix(ori_mat, mmoi_local_space):
    """
    Converts the inertia matrix from local to global space

    Args:
        ori_mat (Matrix3): Orientation matrix
        mmoi_local_space (Matrix3): Inertia matrix in local frame

    Returns:
        Matrix3: Inertia matrix in global frames
    """
    mmoi_global_space = (
        ori_mat.transpose()*mmoi_local_space*ori_mat
        )
    return mmoi_global_space

# ----------------------------------------------------- mover()
track  = Tracker(target_init_orientation_quarternion, Vector3(0,1,0), Vector3(1,0,0))
switch = True

def mover():
    """
    Moves the target according to the experimental observation.
    This simulates the part of the experiment where the target is
    manually moved.
    """
    angle, ang_vel = getAngleAngVel()


    # driving the target
    if O.time < time_when_target_driven_ends:
        # Applying the initial motion to the target
        applyMotionToTarget(angle, ang_vel)


    # target driving by dynamics
    else:
        print("\aSwitch iter:",  O.iter)
        print("Switch time:"  ,  O.time)
        moverEngine.dead =  True
        # The target was not affected by forces before.
        # Now change physics settings to allow forces to
        # affect it 
        assignCorrectPhysics()        

# ***************************************************** assignCorrectPhysics
def assignCorrectPhysics():
    """
    The initial target motion must have the dynamics turned off for it to work,
    once the manual forcing is done, the dynamics is turned on using this function.
    It also assignes correct constrains
    """
    for i in target.moving_id_list:
        O.bodies[i].dynamic = True
    del i

    for i in target.moving_id_list:
        O.bodies[i].state.blockedDOFs = "xyzYZ"
    del i

# ***************************************************** getAngleAngVel
def getAngleAngVel():
    """Calculates the angular position and angular velocity from the lerper

    Returns:
        (float, float): angular position, angular velocity
    """
    angle   =  exp_observed_ori_lerper.out(O.time)
    ang_vel = (exp_observed_ori_lerper.out(O.time + O.dt) - exp_observed_ori_lerper.out(O.time)) / O.dt
    return angle, ang_vel

# ***************************************************** plotSimExpOmega
def plotSimExpOmega():
    """
    Plots the angular position and the simulate angular velocity and
    the simulate experimentally observed angular velocity
    """
    ori               = O.bodies[target.moving_id_list[0]].state.ori
    track.addOri(ori)
    angle             = track.getLatestAngle()            
    w                 = O.bodies[target.moving_id_list[0]].state.angVel

    plot.addData(
        t  = O.time, theta = angle,
        t2 = O.time, wx    = w[0], 
        wx_og = (
            exp_observed_raw_ori_lerper.out(O.time + O.dt) - 
            exp_observed_raw_ori_lerper.out(O.time)
            ) / O.dt
    )

# ***************************************************** applyMotionToTarget
def applyMotionToTarget(angle, ang_vel):
    """a
    Assigns the angular position and angular velocity to all the bodies in the target. 
    This is based on the experimental data.

    Args:
        angle (float): Angular position
        ang_vel (float): Angular velocity
    """
    for id in target.moving_id_list:
        O.bodies[id].state.ori    = Quaternion((1,0,0), angle ) * target_init_orientation_quarternion
        O.bodies[id].state.angVel = [ang_vel, 0.0, 0.0]

# ----------------------------------------------------- plotter
def plotter():
    plotSimExpOmega()

# ----------------------------------------------------- calcNumberIter
def calcNumberIter():
    number_of_data = len(exp_observed_ori)
    time_sim_end   = number_of_data*dt
    no_iter        = int((time_sim_end + offset_time)/ O.dt)
    return no_iter
    
# ----------------------------------------------------- stopper
O.dt = 1e-6
end_iter = calcNumberIter()
O.engines += [
    PyRunner(
        command    = 'stopper()',
        iterPeriod = end_iter,
        label      = 'StopperEngine' 
    )
]

# ----------------------------------------------------- stopper
O.engines += [
    PyRunner(
        command    = 'stopper()',
        iterPeriod = 1600000,
        label      = 'StopperEngine' 
    )
]

def stopper():
    result_dict               = {}
    result_dict['parameters'] = getParameters(globals())
    result_dict['data']       = plot.data
    
    file_json = './res/Sim_result_ID_' + str(ID) + '.json'
    with open(file_json, 'w') as f:
        json.dump(result_dict, f, indent=4, cls=VectorEncoder)

    O.pause()
# ----------------------------------------------------------------------------- simulation controls
O.dt = 1e-6
if simulation_mode == 'batch':
    O.run()
    waitIfBatch()

else:
    O.saveTmp()

plot.plot()

