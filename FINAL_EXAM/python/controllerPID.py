import numpy as np
import massParam as P

class controllerPID:

    def __init__(self):
        self.kp = 2.5
        self.ki = 1.75
        self.kd = 1.497
        self.limit = P.F_max
        self.sigma = 0.05 # from test document and email
        self.beta = (2.0 * self.sigma - P.Ts) \
                / (2.0 * self.sigma + P.Ts)
        self.Ts = P.Ts

        self.z_d1 = 0.
        self.z_dot = 0.
        self.error_d1 = 0.
        self.integrator = 0.
        self.F_e = P.F_e

    def update(self, z_r, y):
        z = y[0][0]

        # define integrator for PID control:
        error = z_r - z
        self.integrator = self.integrator + (P.Ts/2) * (error + self.error_d1)
        self.z_dot = self.beta * self.z_dot \
                + (1 - self.beta) * ((z - self.z_d1) / P.Ts)
        
        # feedback linearized force
        F_fl = P.k1 * z
        # # compute the linearized force using PD control
        # F_tilde = self.kp * (z_r - z) - self.kd * self.zdot
        # PID control
        F_tilde = self.kp * error \
            + self.ki * self.integrator \
            - self.kd * self.z_dot

        # compute the final force and saturate
        F_unsat = F_fl + F_tilde + P.F_e
        F = self.saturate(F_unsat)

        # integrator anti-windup
        if self.ki != 0.0:
            self.integrator = self.integrator + P.Ts / self.ki * (F - F_unsat)

        # update delayed variables
        self.error_d1 = error
        self.z_d1 = z

        return F

    def saturate(self, u):
        if abs(u) > self.limit:
            u = self.limit*np.sign(u)
        return u







