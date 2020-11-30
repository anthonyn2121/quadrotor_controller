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
        
        self.ct = 1 
        self.cq = 1 
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

    def PD_control(self, state, flat_output): 
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
        T_tot = self.m*self.g - (self.Kp[2]*error_pos[2] + self.Kd[2]*error_vel[2])
        
        # Find desired roll and pitch 
        roll_des = self.Kp[0]*error_pos[0] + self.Kd[0]*error_vel[0] # x-axis 
        pitch_des = self.Kp[1]*error_pos[1] + self.Kd[1]*error_vel[1] # y-axis 

        # Find error in attitude 
        euler_angles = Rotation.from_quat(state['q']).as_euler('xyz', degrees=False)
        error_att = np.array([roll_des, pitch_des, flat_output['yaw']]) - euler_angles        

        # Attitude control, given by torque matrix
        # desired values 
        psi_dot_des = (state['w'][2]-flat_output['yaw_dot'])/C(roll_des)    
        theta_dot_des = (state['w'][0] - psi_dot_des*S(roll_des)*S(pitch_des))/C(pitch_des)

        euler_rate_des = np.array([theta_dot_des, psi_dot_des, flat_output['yaw_dot']]) # into array 

        
        # current values 
        psi_dot_curr = (state['w'][1] + (state['w'][0]/C(euler_angles[2]))) / (S(euler_angles[0])*C(euler_angles[2]) + T(euler_angles[0]/euler_angles[2])*T(euler_angles[2]))
        theta_dot_curr = (state['w'][0] - psi_dot_curr*S(euler_angles[0])*S(euler_angles[2]))/C(euler_angles[2])
        phi_dot_curr = (state['w'][2] - psi_dot_curr*C(euler_angles[0])) 

        euler_rate_curr = np.array([theta_dot_curr, psi_dot_curr, phi_dot_curr]) # into array 

        # error matrix 
        euler_rate_err = euler_rate_des - euler_rate_curr
        
        # tau matrix
        tau1 = self.Kr[0]*error_att[0] + self.Kw[0]*euler_rate_err[0]
        tau2 = self.Kr[1]*error_att[1] + self.Kw[1]*euler_rate_err[1]
        tau3 = self.Kr[2]*error_att[2] + self.Kw[2]*euler_rate_err[2]
        # tau = np.array([tau1, tau2, tau3])  


        # Solve for omega 
        u = np.array([tau1, tau2, tau3, T_tot])
        A = np.array([
            [0, -self.ct*self.l, 0, self.ct*self.l], 
            [-self.ct*self.l, 0, self.ct*self.l, 0], 
            [-self.cq, self.cq, -self.cq, self.cq], 
            [self.ct, self.ct, self.ct, self.ct]
        ])
        w_squared =  np.linalg.inv(A) @ u
        print(w_squared)
        return w 


if __name__ == "__main__": 
    quad = System(quad_params)
    quad.PD_control()
