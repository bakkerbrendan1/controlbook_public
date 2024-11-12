import numpy as np 
import VTOLParam as P


class VTOLDynamics:
    def __init__(self, alpha=0.0):
        # Initial state conditions
        self.state = np.array([
            [P.z0],  # initial x
            [P.h0],  # initial height
            [P.theta0],  # initial angle
            [P.zdot0],  # initial velocity in x
            [P.hdot0], # initial velocity in y
            [P.thetadot0] # initial angular velocity
        ])
        # simulation time step
        self.Ts = P.Ts
        # define three masses
        self.Jc = P.Jc * (1.+alpha*(2.*np.random.rand()-1.))
        self.mr = P.mr * (1.+alpha*(2.*np.random.rand()-1.))
        self.mc = P.mc * (1.+alpha*(2.*np.random.rand()-1.))
        # length of arms
        self.d = P.d * (1.+alpha*(2.*np.random.rand()-1.))  
        # damping coefficient
        self.mu = P.mu * (1.+alpha*(2.*np.random.rand()-1.))
        # define gravity, force limit
        self.g = P.g
        self.force_limit = P.fmax

    def update(self, u):
        # This is the external method that takes the input u at time
        # t and returns the output y at time t.
        # saturate the input torque
        u[0][0] = saturate(u[0][0], self.force_limit)
        u[1][0] = saturate(u[1][0], self.force_limit)
        self.rk4_step(u)  # propagate the state by one time sample
        y = self.h()  # return the corresponding output
        return y

    def f(self, state, u):
        # Return xdot = f(x,u)
        zv = state[0][0]
        h = state[1][0]
        theta = state[2][0]
        zvdot = state[3][0]
        hdot = state[4][0]
        thetadot = state[5][0]
        fl = u[0][0]
        fr = u[1][0]
        # The equations of motion.
        # M = np.array([[self.Js, 0],
        #               [0, self.Jp]])
        # C = np.array([[tau - self.b*(thetadot-phidot)-self.k*(theta-phi)],
        #               [-self.b*(phidot-thetadot)-self.k*(phi-theta)]])
        # tmp = np.linalg.inv(M) @ C
        # thetaddot = tmp[0][0]
        # phiddot = tmp[1][0]
        zvddot = -(self.mu*zvdot + (fl + fr)*np.sin(theta) + P.F_wind) / (self.mc + 2.0*self.mr)
        hddot = (self.g*(-self.mc - 2*self.mr) + (fl+fr)*np.cos(theta)) / (self.mc + 2.0*self.mr)
        thetaddot = (self.d*(fl-fr)) / (self.Jc + 2.0*self.d**2 * self.mr)
        # print(zvdot)
        # print(hdot)
        # print(thetadot)
        # print(zvddot)
        # print(hddot)
        # print(thetaddot)
        # build xdot and return
        xdot = np.array([[zvdot], [hdot], [thetadot], [zvddot], [hddot], [thetaddot]])
        return xdot

    def h(self):
        # return y = h(x)
        zv = self.state[0][0]
        h = self.state[1][0]
        theta = self.state[2][0]
        y = np.array([[zv], [h], [theta]])
        return y

    def rk4_step(self, u):
        # Integrate ODE using Runge-Kutta RK4 algorithm
        F1 = self.f(self.state, u)
        F2 = self.f(self.state + self.Ts / 2 * F1, u)
        F3 = self.f(self.state + self.Ts / 2 * F2, u)
        F4 = self.f(self.state + self.Ts * F3, u)
        self.state += self.Ts / 6 * (F1 + 2 * F2 + 2 * F3 + F4)

        
def saturate(u, limit):
    if abs(u) > limit:
        u = limit*np.sign(u)
    return u
