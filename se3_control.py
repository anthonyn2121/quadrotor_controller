import numpy as np
from scipy.spatial.transform import Rotation
from math import cos, sin


class SE3Control(object):
    """

    """
    def __init__(self, quad_params):
        """
        This is the constructor for the SE3Control object. You may instead
        initialize any parameters, control gain values, or private state here.

        For grading purposes the controller is always initialized with one input
        argument: the quadrotor's physical parameters. If you add any additional
        input arguments for testing purposes, you must provide good default
        values!

        Parameters:
            quad_params, dict with keys specified by crazyflie_params.py

        """

        # Quadrotor physical parameters.
        self.mass            = quad_params['mass'] # kg
        self.Ixx             = quad_params['Ixx']  # kg*m^2
        self.Iyy             = quad_params['Iyy']  # kg*m^2
        self.Izz             = quad_params['Izz']  # kg*m^2
        self.arm_length      = quad_params['arm_length'] # meters
        self.rotor_speed_min = quad_params['rotor_speed_min'] # rad/s
        self.rotor_speed_max = quad_params['rotor_speed_max'] # rad/s
        self.k_thrust        = quad_params['k_thrust'] # N/(rad/s)**2
        self.k_drag          = quad_params['k_drag']   # Nm/(rad/s)**2

        # You may define any additional constants you like including control gains.
        self.inertia = np.diag(np.array([self.Ixx, self.Iyy, self.Izz])) # kg*m^2 body inertia 
        self.g = 9.81 # m/s^2

        # STUDENT CODE HERE
        self.Kp = np.diag(np.array([1,1,1]))
        self.Kd = np.diag(np.array([1,1,1]))
        self.Kr = np.diag(np.array([250, 250, 250]))
        self.Kw = np.diag(np.array([20, 20, 20]))
        self.grav = np.array([0,0,-self.mass * self.g])
        self.d = self.k_drag/self.k_thrust

        
    def update(self, t, state, flat_output):
        """
        This function receives the current time, true state, and desired flat
        outputs. It returns the command inputs.

        Inputs:
            t, present time in seconds
            state, a dict describing the present state with keys
                x, position, m
                v, linear velocity, m/s
                q, quaternion [i,j,k,w]
                w, angular velocity, rad/s
            flat_output, a dict describing the present desired flat outputs with keys
                x,        position, m
                x_dot,    velocity, m/s
                x_ddot,   acceleration, m/s**2
                x_dddot,  jerk, m/s**3
                x_ddddot, snap, m/s**4
                yaw,      yaw angle, rad
                yaw_dot,  yaw rate, rad/s

        Outputs:
            control_input, a dict describing the present computed control inputs with keys
                cmd_motor_speeds, rad/s
                cmd_thrust, N (for debugging and laboratory; not used by simulator)
                cmd_moment, N*m (for debugging; not used by simulator)
                cmd_q, quaternion [i,j,k,w] (for laboratory; not used by simulator)
        """
        cmd_motor_speeds = np.zeros((4,))
        cmd_thrust = 0
        cmd_moment = np.zeros((3,))
        cmd_q = np.zeros((4,))

        # STUDENT CODE HERE
        ''' 
        inner controller STATE = [roll, pitch, yaw, roll rate, pitch rate, yaw rate]
            roll, pitch yaw given from: state['q'] converted to euler angles 
            roll, pitch, yaw RATES given from: state['w']
        ''' 

        (f1, f2, f3, f4) = self.get_forces(state, flat_output)
        print(f1)

        control_input = {'cmd_motor_speeds':cmd_motor_speeds,
                         'cmd_thrust':cmd_thrust,
                         'cmd_moment':cmd_moment,
                         'cmd_q':cmd_q}
        return control_input



    def get_forces(self, state, flat_output): 
        quat = Rotation.from_quat(state['q'])
        euler = quat.as_euler('zyx', degrees=False)
        R = quat.as_matrix()
        
        r_des = flat_output['x_ddot'] + -self.Kd @ (state['v'] - flat_output['x_dot']) + -self.Kp @ (state['x']-flat_output['x'])
        F_des = self.mass*r_des + self.grav 

        b3_des = F_des/np.linalg.norm(F_des) # align b3 body axis so that it is oriented along desried thrust
        a_psi = np.asarray([np.cos(flat_output['yaw']), np.sin(flat_output['yaw']), 0])
        b2_des = np.cross(b3_des, a_psi)/np.linalg.norm(np.cross(b3_des, a_psi))
        b1_des = np.cross(b2_des, b3_des)        
        R_des = np.concatenate((b1_des.reshape((3,1)), b2_des.reshape((3,1)), b3_des.reshape((3,1))), axis=1)
        
        eR = .5* (R_des.T @ R - R.T @ R_des)
        eR = np.array([eR[2, 1], eR[0, 2], eR[1, 0]])
        eR = eR.reshape((3,1))
        w_des = np.zeros((3))
        eW = state['w'] - w_des
        eW = eW.reshape((3,1))
        u2 = np.matmul(self.inertia,(np.matmul(-self.Kr,eR) - np.matmul(self.Kw,eW)))
        u1 = (R @ np.array([0,0,1]).reshape((3,1))).T @ F_des

        u = (np.vstack((u1,u2)))
        
        F = np.linalg.inv(np.array([[1,1,1,1],
                    [0, self.arm_length, 0, -self.arm_length],
                    [-self.arm_length, 0, self.arm_length, 0],
                    [self.d, -self.d, self.d, -self.d]]))   @ u
        print(F)
        # print(F)
        return F
