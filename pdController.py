import numpy as np 
from scipy.spatial.transform import Rotation as R 
import matplotlib.pyplot as plt 
import control

def C(x): 
    return np.cos(x)

def S(x): 
    return np.sin(x)

def T(x):
    return np.tan(x)

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

m = quad_params['mass']
Ixx = quad_params['Ixx']
Iyy = quad_params['Iyy']
Izz = quad_params['Izz']
l = quad_params['arm_length']
Kt = quad_params['k_thrust']
Kd = quad_params['k_drag']
d = Kd/Kt
g = 9.81

initial_state = {'x': np.array([0, 0, 0]),
                 'v': np.zeros(3,),
                 'q': np.array([0, 0, 0, 1]), # [i,j,k,w]
                 'w': np.zeros(3,)}

quat = R.from_quat(initial_state['q'])
euler = quat.as_euler('xyz', degrees=False)

x_init = np.asarray([euler, initial_state['w'], initial_state['x'], initial_state['v']])

#Mapping Matrix from u to F

A = np.array([[1,1,1],[-l,l,0],[d,d,-d]])

#Desired Positions
x_d = 0
y_d = 0
z_d = 1

x1_d = 0
x2_d = 0

xdot_d = 0
xdot1_d = 0
xdot2_d = 0

ydot_d = 0
zdot_d = 0

phi_d = 0
theta_d = 0 
psi_d = m*g*d / Kd

#initialize total force

uf = 0
u = np.array([0,0,0])

#Initial States Seperated into changing variables

w = initial_state['w']
x = initial_state['x']
v = initial_state['v']

#Gain values

k_px = k_py = k_pz = k_dx = k_dy = k_dz = k_pphi = k_dphi = k_ptheta = k_dtheta = 1

#Time Parameters

t = 1 #Time step of one second
T = 60 #Time period of sixty seconds

#Dynamics

for i in T:
  
  phidot = w[0] + w[1]*S(euler[0])*T(euler[1]) +w[2]*C(euler[0])*T(euler[1])
  thetadot = w[1]*C(euler[0]) - w[2]*S(euler[0])

  xdot = (C(euler[0])*S(euler[1])*C(euler[2]) + S(euler[0])*S(euler[2]))*uf - (Kt/m)*x[0]
  ydot = (C(euler[0])*S(euler[1])*S(euler[2]) - S(euler[0])*C(euler[2]))*uf - (Kt/m)*x[1]
  zdot = (uf*(euler[0])*C(euler[1]) - Kt*x[2] - m*g)/m

  # #outer loop (position) nu 1
  
  ux = k_px*(x_d - x[0]) + k_dx*(xdot_d - xdot)

  uy = k_py*(y_d - x[1]) + k_dy*(ydot_d - ydot)


  # #inner loop (orientation) nu 2

  uphi = k_pphi*(ux - euler[0]) + k_dphi*(xdot1_d - phidot)

  utheta =  k_ptheta*(uy - euler[1]) + k_dtheta*(xdot2_d - thetadot)

  #Altitude Controller

  uf = k_pz*(z_d - x[2]) + k_dz *(zdot_d - zdot)

  #update Quadrotor

  u = np.array([[uf , uphi, utheta]])
  print(np.linalg.inv(A),u.T)

  f = np.linalg.inv(A) @ u.T

  #Send forces to individual rotors (If using a simulation)

  #Update Positions
  x[0] += xdot*t
  x[1] += ydot*t
  x[2] += zdot*t