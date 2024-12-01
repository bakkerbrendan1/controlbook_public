import numpy as np
import massParam as P
import control as cnt

class ctrlObserver:
    def __init__(self):
        # tuning parameters
        tr = 0.25
        zeta = 0.707
        integrator_pole = 1.2
        tr_obs = tr/10  # rise time frequency for observer
        zeta_obs = 0.707  # damping ratio for observer

        # State Space Equations
        # xdot = A*x + B*u
        # y = C*x
        self.A = P.A
        self.B = P.B
        self.C = P.C
        Cr = np.array([[1.0, 0.0]])

        # form augmented system
        A1 = np.vstack((np.hstack((self.A, np.zeros((np.size(self.A,1),1)))), 
                        np.hstack((-Cr, np.array([[0.0]]))) ))
        B1 = np.vstack( (self.B, 0.0) )

        # gain calculation
        wn = 2.2 / tr
        des_char_poly = np.convolve([1, 2*zeta*wn, wn**2],
                                    [1, -integrator_pole])
        des_poles = np.roots(des_char_poly)

        # compute the gains if the system is controllable
        if np.linalg.matrix_rank(cnt.ctrb(A1, B1)) != 3:
            print('The system is not controllable!')
        else:
            K1 = cnt.place(A1, B1, des_poles)
            self.K = K1[0][0:2]
            self.ki = K1[0][2]

        # observer design
        wn_obs = 2.2 / tr_obs
        des_obsv_char_poly = [1, 2*zeta_obs*wn_obs, wn_obs**2]
        des_obsv_poles = np.roots(des_obsv_char_poly)
        # Compute the gains if the system is controllable
        if np.linalg.matrix_rank(cnt.ctrb(self.A.T, self.C.T)) != 2:
            print('The system is not observable')
        else:
            self.L = cnt.acker(self.A.T, self.C.T, des_obsv_poles).T

        print('K: ', self.K)
        print('kr: ', self.ki)
        print('L^T: ', self.L)
        #--------------------------------------------------
        # variables to implement the integrator
        self.integrator = 0.0  # integrator
        self.error_d1 = 0.0  # error delayed by 1 sample
        self.x_hat = np.array([
            [0.0], # z_hat_0
            [0.0], # zdot_hat_0
        ])
        self.F_d1 = 0.0 # control force delayed by 1 sample

    def update(self, z_r, y):
        x_hat = self.update_observer(y)
        z_hat = x_hat[0][0]

        # integrate error
        error = z_r - z_hat
        self.integrator = self.integrator \
                          + (P.Ts / 2.0) * (error + self.error_d1)
        self.error_d1 = error

        # feedback linearized force F_fl
        F_fl = P.k * z_hat
        # compute the state feedback controller
        F_tilde = -self.K @ x_hat + self.ki * self.integrator

        # compute the final force and saturate
        F_unsat = F_fl + F_tilde[0]
        F = saturate(F_unsat, P.F_max)

        self.F_d1 = F
        return F, x_hat
    
    def update_observer(self, y_m):
        # update the observer using RK4 integration
        F1 = self.observer_f(self.x_hat, y_m)
        F2 = self.observer_f(self.x_hat + P.Ts / 2 * F1, y_m)
        F3 = self.observer_f(self.x_hat + P.Ts / 2 * F2, y_m)
        F4 = self.observer_f(self.x_hat + P.Ts * F3, y_m)
        self.x_hat += P.Ts / 6 * (F1 + 2 * F2 + 2 * F3 + F4)
        return self.x_hat
    
    def observer_f(self, x_hat, y_m):
        # compute feedback linearizing force F_fl
        z_hat = x_hat[0][0]
        F_fl = P.k * z_hat
        # xhatdot = A*(xhat-xe) + B*(u-ue) + L(y-C*xhat)
        xhat_dot = self.A @ x_hat\
                   + self.B * (self.F_d1 - F_fl)\
                   + self.L * (y_m - self.C @ x_hat)
        return xhat_dot


def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u







