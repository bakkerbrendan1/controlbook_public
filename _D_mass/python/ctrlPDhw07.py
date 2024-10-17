import numpy as np
import massParam as P

class ctrlPD:
    def __init__(self):
        # PD gains
        self.kp = 4.5
        self.kd = 12
        print('kp: ', self.kp)
        print('kd: ', self.kd)

    def update(self, z_r, x):
        z = x[0][0]
        zdot = x[1][0]
        # feedback linearized force
        F_fl = P.k * z
        # no equilibrium force
        z_e = 0.0
        F_e = P.k * z_e
        # compute the linearized force using PD control
        F_tilde = self.kp * (z_r - z) - self.kd * zdot
        # compute total force
        F = F_fl + F_tilde
        #tau = tau_e + tau_tilde
        # always saturate to protect hardware
        F = saturate(F, P.F_max)
        return F


def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u







