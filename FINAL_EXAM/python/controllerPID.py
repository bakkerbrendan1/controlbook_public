import numpy as np
import massParam as P

class controllerPID:

    def __init__(self):
        self.kp =
        self.ki =
        self.kd =
        self.limit =
        self.beta =
        self.Ts =
        self.z_d1 =
        self.z_dot =
        self.error_d1 =
        self.integrator =
        self.F_e =

    def update(self, z_r, y):
        return F

    def saturate(self, u):
        if abs(u) > self.limit:
            u = self.limit*np.sign(u)
        return u







