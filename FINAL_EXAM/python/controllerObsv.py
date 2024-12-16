import numpy as np
import massParam as P

class controllerObsv:
    def __init__(self):
        self.observer_state =
        self.F_d1 =
        self.integrator =
        self.error_d1 =
        self.K =
        self.ki =
        self.L =
        self.A =
        self.B =
        self.C =
        self.limit =
        self.F_e =
        self.Ts =

    def update(self, z_r, y):
        return F, x_hat

    def update_observer(self, y):
        return x_hat

    def observer_f(self, x_hat, y):
        return xhat_dot

    def integrateError(self, error):
        pass

    def saturate(self,u):
        if abs(u) > self.limit:
            u = self.limit*np.sign(u)
        return u

