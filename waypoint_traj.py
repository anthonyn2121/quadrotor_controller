import numpy as np
import scipy.interpolate as sp

class WaypointTraj(object):
    """

    """
    def __init__(self, points):
        """
        This is the constructor for the Trajectory object. A fresh trajectory
        object will be constructed before each mission. For a waypoint
        trajectory, the input argument is an array of 3D destination
        coordinates. You are free to choose the times of arrival and the path
        taken between the points in any way you like.

        You should initialize parameters and pre-compute values such as
        polynomial coefficients here.

        Inputs:
            points, (N, 3) array of N waypoint coordinates in 3D
        """

        # self.v = 1 #m/s
        # self.points = points
        # self.t = np.zeros(len(points),)

        # if np.shape(self.points) == (3,) or np.shape(self.points) == (1,3):
        #     pass
        # elif np.shape(self.points) != (3,) or np.shape(self.points) != (1,3):
        #     for i in range(len(self.t)-1):
        #         self.t[(i+1)] = np.linalg.norm((points[(i+1)]-points[i]))/self.v

        #     self.point_t = np.zeros(len(points),)
        #     for i in range(int(len(self.t)-1)):
        #         self.point_t[(i+1)] = self.point_t[i] + self.t[i+1]

        #     self.f = sp.CubicSpline(self.point_t,self.points,axis = 0)

        self.points = points
        self.dist = np.zeros((len(self.points)-1,3)) 
        for i in range(len(self.points)-1):
            self.dist[i] = points[i+1] - points[i]

        self.times = np.linalg.norm(self.dist,axis=1)
        self.times = np.cumsum(self.times)
        self.times = np.insert(self.times,0,0,axis = 0)
        self.f = sp.CubicSpline(self.times,self.points,axis=0)

        self.flag = False
        self.vnew = np.array([0,0,0])
        self.xP = np.array([0,0,0])
        self.tP = 0
        self.xi = 0
        
        self.ctr = 0
        
        # STUDENT CODE HERE

    def update(self, t):
        """
        Given the present time, return the desired flat output and derivatives.

        Inputs
            t, time, s
        Outputs
            flat_output, a dict describing the present desired flat outputs with keys
                x,        position, m
                x_dot,    velocity, m/s
                x_ddot,   acceleration, m/s**2
                x_dddot,  jerk, m/s**3
                x_ddddot, snap, m/s**4
                yaw,      yaw angle, rad
                yaw_dot,  yaw rate, rad/s
        """
        x        = np.zeros((3,))
        x_dot    = np.zeros((3,))
        x_ddot   = np.zeros((3,))
        x_dddot  = np.zeros((3,))
        x_ddddot = np.zeros((3,))
        yaw = 0
        yaw_dot = 0



        # STUDENT CODE HERE

        # if np.shape(self.points) == (3,) or np.shape(self.points) == (1,3):
        #     x = np.reshape(self.points,(3,))
        # elif np.shape(self.points) != (3,) or np.shape(self.points) != (1,3):
        #     if t > self.point_t[-1]:
        #         x = self.points[-1]
        #     else:
        #         x = self.f(t)

        # a = np.array([1/.03, 1/.03, 0])
        # vi = np.array([0,0,0])
        # for i in range():
        #     distance = vi*t + 1/2 *a* self.times**2
        #     v = vi + 2*a*distance
        
        
        const_a = np.linalg.norm(self.dist)*2/(self.times[1]**2)
        #a = np.array([const_a, const_a, 0])
        a = np.array([.1,.1,0])
        disturbance = np.array([.05,0,0])
        #disturbance = np.array([0,0,0])
        a -= disturbance
        #a = np.array([0,0,0])


        v = a*t #self.times[self.ctr]

        if t == np.inf:
            x = np.array([0,0,0])
            v = np.array([0,0,0])
        
        elif np.linalg.norm(self.points[-1]-self.xP) < .5 : 
            
            x = self.points[-1]
            v = np.array([0,0,0])
            a = np.array([0,0,0])

        elif (np.linalg.norm((self.points[3] + self.points[2])/2 -self.xP) < .01) & (self.flag ==False):
            
            self.vnew = a*t
            x = (1/2)*a*t**2
            self.flag = True
            self.tP = t
            self.xi = x

        elif self.flag == True:
            a = -1*a
            x = self.xi + self.vnew*(t-self.tP) + (1/2)*a*(t-self.tP)**2
            v = self.vnew + a*(t-self.tP)

        else:
            x = (1/2)*a*t**2
   
        self.xP = x

        x_dot = v 
        x_ddot = a

        if self.ctr < 5:
            self.ctr += 1

        print(x, x_dot, x_ddot)

        flat_output = { 'x':x, 'x_dot':x_dot, 'x_ddot':x_ddot, 'x_dddot':x_dddot, 'x_ddddot':x_ddddot,
                        'yaw':yaw, 'yaw_dot':yaw_dot}
        return flat_output
