import numpy as np
from scipy.spatial.transform import Rotation

def S(x): 
    return np.sin(x)

def C(x): 
    return np.cos(x)

def T(x):
    return np.tan(x)



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
        self.ct = self.k_thrust
        self.cq = self.k_drag
        self.m = self.mass
        self.l = self.arm_length
        self.gamma = self.k_drag/self.k_thrust
        
        # STUDENT CODE HERE
        self.Kp = np.diag(np.array([2, 2.3, 2]))
        self.Kd = np.diag(np.array([50, 65, 10]))
        self.Kr = np.diag(np.array([1, 1, 10]))
        self.Kw = np.diag(np.array([40, 40, 10]))
        
        
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

        # Compute the commanded acceleraton r_ddot_des 
        r_ddot_des = flat_output['x_ddot'] - self.Kd @ (state['v'] - flat_output['x_dot']) - self.Kp @ (state['x'] - flat_output['x'])
        
        # Compute the desired angles phi_des, theta_des from the desired acceleration 
        phi_des = (S(flat_output['yaw'])*r_ddot_des[0] - r_ddot_des[1]*C(flat_output['yaw']))/self.g  # roll 
        theta_des = (S(flat_output['yaw'])*r_ddot_des[0] + r_ddot_des[1]*C(flat_output['yaw']))/self.g # pitch
        psi_des = flat_output['yaw']
        euler_des = np.array([phi_des, theta_des, psi_des])

        # Find the total upward thrust from equation of motion 
        u1 = (r_ddot_des[2] + self.g) * self.mass
        
        # Find the moment vector 
        euler_ang = Rotation.from_quat(state['q']).as_euler('xyz', degrees=False)
        omega_des = np.array([0, 0, flat_output['yaw_dot']])
        u2  = self.inertia @ (-self.Kr @ (euler_ang - euler_des) - self.Kw @ (state['w'] - omega_des))
        
        # From sys inputs, find real inputs 
        gamma = self.k_drag/self.k_thrust
        l = self.arm_length 
        A = np.array([
            [1, 1, 1, 1], 
            [0, l, 0, -l], 
            [-l, 0, l, 0], 
            [gamma, -gamma, gamma, -gamma]
        ])
        u = np.append(u1, u2)
        
        F = np.linalg.inv(A) @ u 
        print(F)

        
        cmd_motor_speeds = np.sqrt(F/self.k_thrust)

        for i in range(len(cmd_motor_speeds)): 
            if cmd_motor_speeds[i] < self.rotor_speed_min: 
                cmd_motor_speeds[i] = self.rotor_speed_min
            if cmd_motor_speeds[i] > self.rotor_speed_max: 
                cmd_motor_speeds[i] = self.rotor_speed_max


        # print(np.sqrt(np.sum(np.square(flat_output['x'] - state['x']))))



        control_input = {'cmd_motor_speeds':cmd_motor_speeds,
                         'cmd_thrust':cmd_thrust,
                         'cmd_moment':cmd_moment,
                         'cmd_q':cmd_q}
        return control_input

