import numpy as np 
from math import * 
from scipy.spatial.transform import Rotation as R

quad_params = {
    'mass': 0.030,   # kg
    'Ixx':  1.43e-5, # kg*m^2
    'Iyy':  1.43e-5, # kg*m^2
    'Izz':  2.89e-5, # kg*m^2
    'arm_length': 0.046, # meters
    'rotor_speed_min': 0,    # rad/s
    'rotor_speed_max': 2500, # rad/s
    'k_thrust': 2.3e-08, # N/(rad/s)**2
    'k_drag':   7.8e-11, # Nm/(rad/s)**2
}

def C(x): 
    return np.cos(x)

def S(x): 
    return np.sin(x)

def T(x): 
    return np.tan(x)


class System(object): 
    def __init__(self, quad_params): 
        (m, Ixx, Iyy, Izz, l, Kt, Kd, d, g) = self.get_params()
        
        self.m = m # mass 
        self.Ixx = Ixx 
        self.Iyy = Iyy 
        self.Izz = Izz 
        self.l = l # arm length 
        self.Kt = Kt 
        self.Kd = Kd 
        self.d = Kt/Kd
        self.g = g # gravity 
        ############ 
        self.Kp = np.array([1, 1, 1]) # proportional gains ~ position 
        self.Kd = np.array([1, 1, 1]) # derivative gains ~ position 
        self.Kr = np.array([1, 1, 1]) # proportional gains ~ attitude 
        self.Kw = np.array([1, 1, 1,]) # derivative gains ~ attitude 

    def get_params(self):
        m = quad_params['mass']
        Ixx = quad_params['Ixx']
        Iyy = quad_params['Iyy']
        Izz = quad_params['Izz']
        l = quad_params['arm_length']
        Kt = quad_params['k_thrust']
        Kd = quad_params['k_drag']
        d = Kd/Kt
        g = 9.81
        return m, Ixx, Iyy, Izz, l, Kt, Kd, d, g 

    def PD_control(self, current_state, goal_state): 
        '''
        state = {
            x, position, m 
            v, linear velocity, m/s 
            q, quaternion, [i,j,k,w]
            w, angular velocity, rad/s
        }
        
        Given desired outputs: 
            flat_output = { 
             x, position, m 
             x_dot,  velocity, m/s 
             x_ddot, acceleration, m/s^2
             x_dddot, jerk, m/s^3
             x_ddddot, snap, m/s^4
             yaw 
             yaw_dot   
            }
        '''

        # Find errors in position and velocity 
        error_pos = flat_output['x'] - state['x']
        error_vel = flat_output['x_dot'] - state['v']

        # Altitude control with constant feedforward compensation for gravity 
        # Total thrust 
        T = self.m*self.g - (self.Kp[2]*error_pos[2] + self.Kd[2]*error_vel)

        # Find desired roll and pitch 
        roll_des = self.Kp[0]*error_pos + self.Kd[0]*error_vel # x-axis 
        pitch_des = self.Kp[1]*error_pos + self.Kd[1]*error_vel # y-axis 

        # Find error in attitude 
        euler_angles = R.from_quat(state['quaternion']).as_euler('xyz', degrees=False)
        error_att = np.array([roll_des, pitch_des, flat_output['yaw']]) - euler_angles
        
        
        
        
        # Attitude control, given by torque matrix 
        tau1 = self.Kr[0]*error_att + self.Kw[0]*

        return T, tau 
if __name__ == "__main__": 
    quad = System(quad_params)
