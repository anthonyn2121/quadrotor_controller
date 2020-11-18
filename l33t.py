from scipy.io import loadmat 
import numpy as np 
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

x_init = np.hstack((euler.reshape(3), initial_state['w'], initial_state['x'], initial_state['v']))
print(x_init)
x_init = np.ones((5))*.1
print(x_init)

def C(x): 
    return np.cos(x)

def S(x): 
    return np.sin(x)

def T(x): 
    return np.tan(x)

x_hat = x_init[:5]
x = x_hat 

K = loadmat('/home/anthony/Downloads/Cmat_new.mat')
Cmat = K['KnewC']

alpha = np.zeros(3)
alpha[0] = -2*l*Izz*Kd**2*(Ixx*(-.25*C(x[1])**(-2)*Izz*(2*d**2*g**2*m**2*(S(2*x[0])+x[0]+C(2*x[0])*(-4*S(x[0])+x[0])-S(2*x[1])*x[1]*T(x[0]))+ \
    2*d*g*m*Kd*(-3*x[3]-4*C(2*x[1])*S(x[0])*x[4]+C(x[0])**(-1)*(2*x[3]-C(3*x[0])*x[3]+2*C(2*x[0])*C(2*x[1])*x[3]+S(3*x[0])*x[4])+(2*S(2*x[1])*x[2]+x[4])*T(x[0]))+ \
        Kd**2*(-6*x[3]*x[4]-C(3*x[0])**(-1)*(2*S(2*x[1])*x[2]*x[3]+2*C(3*x[0])*x[3]*x[4]+S(3*x[0])*(x[3]-x[4])*(x[3]+x[4]))+ \
            (x[4]**2-3*x[3]**2+2*C(2*x[1])*x[3]**2)*T(x[0]))-C(x[0])**(-1)*Kd**3*x[4]*T(x[1]))+Izz*(Kd**2*(C(x[0])**(-1)*Kd*x[4]*T(x[1])-d*g*m*(x[1]+C(x[0])**(-1)*\
                (-1+S(x[0])*x[0])*T(x[0])))) + Izz*(-Kd**2*x[3]*x[4]+d*g*m*Kd*(-x[3]+C(x[0])*x[3]-S(x[0])*x[4]+x[2]*T(x[0])*T(x[1])+S(x[0])*(x[4]+x[3]*T(x[0]))*T(x[1])**2))+\
                    d**2*g**2*m**2*(x[0]-x[1]*T(x[0])*T(x[1]) + S(x[0])*(-1+T(x[1])**2))))

alpha[1] = Izz*Kd**2*(Ixx*(-Kd**3*x[4]*T(x[0])+ Izz*(d**2*g**2*m**2*(C(x[0])**(-1)*x[1]+(-2+C(x[0]))*T(x[1]))+Kd**2*(2*x[4]+x[3]*T(x[0])+\
    C(x[0])**(-1)*(S(x[0])*x[3]+C(x[0])*x[4])**2*T(x[1]))+2*d*g*m*Kd*(x[2]-C(x[0])**(-1)*x[2]+(S(x[0])*x[3]-x[4]+C(x[0])*x[4]-x[3]*T(x[0]))*T(x[1])))+\
        Izz*(C(x[0])**(-1)*Kd**2*(d*g*m*(S(x[0])*x[0])+S(x[0])*Kd*x[4])+Izz*(-Kd**2*x[2]*x[4]+d**2*g**2*m**2*(-C(x[0])**(-1)*x[1]+T(x[1]))+\
            d*g*m*Kd*((-1+C(x[0])**(-1)*x[2]+(x[4]+x[3]*T(x[0]))*T(x[1])))))))

beta = np.array([
    [1, -(2*T(x[0])*T(x[1]))/l, -(2*C(x[0])**(-1)*Ixx*T(x[1]))/(l*Izz)],
    [0, S(C(x[0])), (Ixx*T(x[0]))/Izz],
    [0, 0, 1]
])

v = Cmat @ x
u = alpha + beta @ v
print(u)
# print(Cmat.shape)