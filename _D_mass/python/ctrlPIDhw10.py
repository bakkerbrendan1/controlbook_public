import numpy as np
import massParam as P

class ctrlPID:
    def __init__(self):
        # tuning parameters
        tr = 2.0
        zeta = 0.7
        self.ki = 0.3 # integrator gain

        # PD gains
        # open loop char polynomial and poles
        a1 = P.b/P.m
        a0 = P.k/P.m
        wn = 2.2/tr
        alpha1 = 2.0 * zeta * wn
        alpha0 = wn**2
        
        self.kp = 5.05 # P.m * (alpha0 - a0) # 5.05
        self.kd = 7.2 # P.m * (alpha1 - a1) # 7.2
        print('kp: ', self.kp)
        print('kd: ', self.kd)
        print('ki: ', self.ki)

        ####################################################
        # dirty derivative gains
        self.sigma = 0.05 # from homework requirements
        self.beta = (2.0 * self.sigma - P.Ts) \
                / (2.0 * self.sigma + P.Ts)

        # initialize variables for integrator and differentiator
        self.zdot = 0.0 # estimated derivative of z
        self.z_d1 = 0.0 # z delayed by one sample
        self.error_dot = 0.0 # estimated derivative of error
        self.error_d1 = 0.0 # error delayed by one sample
        self.integrator = 0.0 # integrator

    def update(self, z_r, x):
        # update state variables, estimate zdot with dirty derivative
        z = x[0][0]
        
        # define integrator for PID control
        error = z_r - z
        self.integrator = self.integrator + (P.Ts/2) * (error + self.error_d1)
        self.zdot = self.beta * self.zdot \
                + (1 - self.beta) * ((z - self.z_d1) / P.Ts)

        # feedback linearized force
        F_fl = P.k * z
        # # compute the linearized force using PD control
        # F_tilde = self.kp * (z_r - z) - self.kd * self.zdot
        # PID control
        F_tilde = self.kp * error \
            + self.ki * self.integrator \
            - self.kd * self.zdot

        # compute the final force and saturate
        F_unsat = F_fl + F_tilde
        F = saturate(F_unsat, P.F_max)

        # integrator anti-windup
        if self.ki != 0.0:
            self.integrator = self.integrator + P.Ts / self.ki * (F - F_unsat)

        # update delayed variables
        self.error_d1 = error
        self.z_d1 = z
        return F


def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u







