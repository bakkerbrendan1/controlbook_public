import numpy as np
from scipy import signal
import blockbeamParam as P
import control as cnt

class ctrlObserver:
    def __init__(self):
        # tuning parameters
        tr_th = 0.5              # Rise time for inner loop (z)
        tr_z = 1.2 # tr_th * 10.0       # 10 times rise time of theta
        zeta_th = 0.95           # inner loop Damping Coefficient
        zeta_z = 0.95            # outer loop Damping Coefficient
        integrator_pole = -4     # integrator pole
        tr_z_obs = tr_z/6.0       # rise time for position
        tr_theta_obs = tr_th/6.0  # rise time for angle


        self.A = P.A
        self.B = P.B
        self.C = P.C
        self.Cr = np.array([[1.0, 0.0, 0.0, 0.0]])

        # form the augmented system
        A1 = np.vstack((
                np.hstack((self.A, np.zeros((4,1)))),
                np.hstack((-self.Cr, np.zeros((1,1))))))
        B1 = np.vstack((self.B, np.zeros((1,1))))

        # gain calculation
        wn_th = 2.2 / tr_th
        wn_z = 2.2 / tr_z
        des_char_poly = np.convolve(
            np.convolve([1, 2*zeta_z*wn_z, wn_z**2],
                        [1, 2*zeta_th*wn_th, wn_th**2]),
            [1, -integrator_pole]
        )
        des_poles = np.roots(des_char_poly)
        # Compute the gains if the system is controllable
        if np.linalg.matrix_rank(cnt.ctrb(A1, B1)) != 5:
            print('The system is not controllable')
        else:
            K1 = cnt.acker(A1, B1, des_poles)
            self.K = K1[0][0:4]
            self.ki = K1[0][4]
        
        # compute observer gains
        wn_z_obs = 2.2 / tr_z_obs
        wn_th_obs = 2.2 / tr_theta_obs
        des_obs_char_poly = np.convolve(
                [1, 2 * zeta_z * wn_z_obs, wn_z_obs**2],
                [1, 2 * zeta_th * wn_th_obs, wn_th_obs**2])
        des_obs_poles = np.roots(des_obs_char_poly)
        # Compute the observer gains if the system is observable
        if np.linalg.matrix_rank(cnt.ctrb(self.A.T, self.C.T)) != 4:
            print('The system is not observable')
        else:
            self.L = cnt.place(self.A.T, self.C.T, des_obs_poles).T
            # self.L = signal.place_poles(self.A.T, self.C.T, 
            #                             des_obs_poles).gain_matrix.T

        # print the gains to the terminal
        print('K: ', self.K)
        print('ki: ', self.ki)
        print('L^T: ', self.L)

        #--------------------------------------------------
        # saturation limits
        theta_max = 30.0 * np.pi / 180.0  # Max theta, rads
        #--------------------------------------------------
        # variables to implement integrator
        self.integrator_z = 0.0  # integrator
        self.error_z_d1 = 0.0  # error signal delayed by 1 sample
        # estimated state variables
        self.x_hat = np.array([
            [0.0],  # initial estimate for z_hat
            [0.0],  # initial estimate for theta_hat
            [0.0],  # initial estimate for z_hat_dot
            [0.0],  # initial estimate for theta_hat_dot
        ])
        self.F_d1 = 0.0 # Computed force, delayed by one sample

    def update(self, z_r, y):
        x_hat = self.update_observer(y)
        # z_hat = x_hat[0][0]
        z_hat = self.Cr @ x_hat

        # integrate error
        error_z = z_r - z_hat

        # integrate error
        self.integrator_z = self.integrator_z \
                + (P.Ts / 2.0) * (error_z + self.error_z_d1)
        self.error_z_d1 = error_z

        # return saturated force
        F_fl = P.m2*P.g/2.0 + P.m1*P.g*P.ze / P.length # define equilibrium force
        
        xe = np.array([[P.ze], [0.0], [0.0], [0.0]])
        x_tilde = x_hat - xe
        
        F_tilde = -self.K @ x_tilde - self.ki * self.integrator_z

        

        F = F_tilde + F_fl
        F = saturate(F[0][0], P.Fmax)
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
        # xhatdot = A*xhat + B*u + L(y-C*xhat)
        xe = np.array([[P.ze], [0.0], [0.0], [0.0]])
        F_e = P.m2*P.g/2.0 + P.m1*P.g*P.ze / P.length

        xhat_dot = self.A @ (x_hat - xe) \
                   + self.B * (self.F_d1 - F_e) \
                   + self.L @ (y_m-self.C @ x_hat)
        return xhat_dot

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
