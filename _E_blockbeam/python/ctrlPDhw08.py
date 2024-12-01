import numpy as np
import blockbeamParam as P

class ctrlPD:
    def __init__(self):
        ####################################################
        #       PD Control: Time Design Strategy
        ####################################################
        # tuning parameters
        tr_z = 1             # Rise time for inner loop (z)
        zeta_th = 0.707       # inner loop Damping Coefficient
        zeta_z = 0.707        # outer loop Damping Coefficient

        #---------------------------------------------------
        #                 Inner Loop: theta
        #---------------------------------------------------
        # parameters of the open loop transfer function
        ze = P.length/2.0 # equilibrium position - center of beam
        b0 = P.length/(P.m2*P.length**2/3.0+P.m1*ze**2)
        M = 10.0              # Time scale separation 
        tr_th = tr_z/M # rise time for inner loop
        # coefficients for desired inner loop
        wn_th = 2.2 / tr_th     # Natural frequency
        # compute gains
        self.kp_th = wn_th**2/b0 # 1.83
        self.kd_th = 2.0*zeta_th*wn_th/b0 # 1.17
        DC_gain = 1.0
        #---------------------------------------------------
        #                   Outer Loop: z
        #---------------------------------------------------
        # coefficients for desired outer loop
        # desired rise time, s, defined in "tuning parameters"
        wn_z = 2.2 / tr_z  # desired natural frequency
        # compute gains
        self.kd_z = -wn_z**2/P.g # -0.0317
        self.kp_z = -2.0*zeta_z*wn_z/P.g # -0.00494
        # print control gains to terminal
        print('DC_gain', DC_gain)
        print('kp_th: ', self.kp_th)
        print('kd_th: ', self.kd_th)
        print('kp_z: ', self.kp_z)
        print('kd_z: ', self.kd_z)
        #---------------------------------------------------
        #                    zero canceling filter
        #---------------------------------------------------
        self.filter = zeroCancelingFilter(DC_gain)

    def update(self, z_r, state):
        z = state[0][0]
        theta = state[1][0]
        zdot = state[2][0]
        thetadot = state[3][0]
        # the reference angle for theta comes from the
        # outer loop PD control

        theta_r = self.kp_z * (z_r - z) - self.kd_z * zdot
        Fe = P.m2*P.g/2.0 + P.m1*P.g*z / P.length # define equilibrium force
        F = self.kp_th * (theta_r - theta) - self.kd_th * thetadot + Fe
        return saturate(F, P.Fmax)

def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u

class zeroCancelingFilter:
    def __init__(self, DC_gain):
        self.a = -3.0 / (2.0 * P.length * DC_gain)
        self.b = np.sqrt(3.0 * P.g / (2.0 * P.length))
        self.state = 0.0

    def update(self, input):
        # integrate using RK1
        self.state += P.Ts * (-self.b * self.state + self.a * input)
        return self.state
