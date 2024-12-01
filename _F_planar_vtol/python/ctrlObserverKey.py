import numpy as np
import VTOLParam as P
from scipy import signal
import control as cnt

class ctrlObserver:
    def __init__(self):
        # tuning parameters
        # wn_h    = 1.0
        # zeta_h  = 0.95
        # wn_z    = 0.9905
        # zeta_z  = 0.95
        # wn_th   = 13.3803
        # zeta_th = 0.95
        # integrator_h = 1.0
        # integrator_z = 1.0

        # # observer gains
        # wn_h_obs    = 10.0*wn_h
        # wn_z_obs    = 10.0*wn_z
        # wn_th_obs   = 10.0*wn_th
        
        # tuning parameters
        tr_h = 2.0 # originally 4.0 # rise time for altitude - tuned for best rise time without saturation
        zeta_h = 0.707 # damping ratio for altitude
        tr_z = 1.0 # originally 4.0 # rise time for outer lateral loop (position) - tuned for best rise time w/o saturation
        zeta_z = 0.707
        tr_th = tr_z / 10.0
        zeta_th = 0.707
        integrator_h = 0.2 # integrator pole
        integrator_z = 0.38 # if tr_z= 2.0, int_pole_z=0.19
        # pick observer poles
        tr_h_obs = tr_h/6.0
        tr_z_obs = tr_z/6.0
        tr_th_obs = tr_th/6.0

        wn_h_obs = 2.2 / tr_h_obs
        wn_z_obs = 2.2 / tr_z_obs
        wn_th_obs = 2.2 / tr_th_obs

        # gain calculation
        wn_h = 2.2/tr_h
        wn_z = 2.2/tr_z
        wn_th = 2.2/tr_th

        

        # State Space Equations
        self.Fe = (P.mc + 2.0 * P.mr) * P.g  # equilibrium force 
        self.A_lon = np.array([[0.0, 1.0],
                        [0.0, 0.0]])
        self.B_lon = np.array([[0.0],
                        [1.0 / (P.mc + 2.0 * P.mr)]])
        self.C_lon = np.array([[1.0, 0.0]])
        self.A_lat = np.array([[0.0, 0.0, 1.0, 0.0],
                        [0.0, 0.0, 0.0, 1.0],
                        [0.0, -self.Fe / (P.mc + 2.0 * P.mr), -(P.mu / (P.mc + 2.0 * P.mr)), 0.0],
                        [0.0, 0.0, 0.0, 0.0]])
        self.B_lat = np.array([[0.0],
                        [0.0],
                        [0.0],
                        [1.0 / (P.Jc + 2 * P.mr * P.d ** 2)]])
        self.C_lat = np.array([[1.0, 0.0, 0.0, 0.0],
                          [0.0, 1.0, 0.0, 0.0]])
        
        # form augmented system
        A1_lon = np.vstack((
                np.hstack((self.A_lon, np.zeros((2,1)))),
                np.hstack((-self.C_lon, np.zeros((1,1))))))
        B1_lon = np.vstack((self.B_lon, np.zeros((1,1))))
        A1_lat = np.vstack((
                np.hstack((self.A_lat, np.zeros((4,1)))),
                np.hstack((-self.C_lat[0:1], np.zeros((1,1))))))
        B1_lat = np.vstack((self.B_lat, np.zeros((1,1))))

        # gain calculation
        des_char_poly_lon = np.convolve([1.0, 2.0 * zeta_h * wn_h, wn_h ** 2],
                                        [1, integrator_h])
        des_poles_lon = np.roots(des_char_poly_lon)
        des_char_poly_lat = np.convolve(
            np.convolve([1.0, 2.0 * zeta_z * wn_z, wn_z ** 2],
                        [1.0, 2.0 * zeta_th * wn_th, wn_th ** 2]),
            [1, integrator_z])
        des_poles_lat = np.roots(des_char_poly_lat)

        # Compute the gains if the system is controllable
        if np.linalg.matrix_rank(cnt.ctrb(A1_lon, B1_lon)) != 3:
            print("The longitudinal system is not controllable")
        else:
            K1_lon = cnt.place(A1_lon, B1_lon, des_poles_lon)
            self.K_lon = K1_lon[0][0:2]
            self.ki_lon = K1_lon[0][2]

        if np.linalg.matrix_rank(cnt.ctrb(A1_lat, B1_lat)) != 5:
            print("The lateral system is not controllable")
        else:
            K1_lat = cnt.place(A1_lat, B1_lat, des_poles_lat)
            self.K_lat = K1_lat[0][0:4]
            self.ki_lat = K1_lat[0][4]

        # compute observer poles
        des_obs_poles_lon = np.roots([1.0, 2.0*zeta_h*wn_h_obs, wn_h_obs**2])
        des_char_poly_lat = np.convolve([1.0, 2.0*zeta_z*wn_z_obs, wn_z_obs**2],
                                        [1.0, 2.0*zeta_th*wn_th_obs, wn_th_obs**2])
        des_obs_poles_lat = np.roots(des_char_poly_lat)

        # check if system is observable and calculate L
        if np.linalg.matrix_rank(cnt.ctrb(self.A_lon.T, self.C_lon.T)) != 2:
            print("The longitudinal system is not observable")
        else:
            self.L_lon = signal.place_poles(self.A_lon.T, self.C_lon.T, des_obs_poles_lon).gain_matrix.T
        if np.linalg.matrix_rank(cnt.ctrb(self.A_lat.T, self.C_lat.T)) != 4:
            print("The lateral system is not observable")
        else:
            self.L_lat = signal.place_poles(self.A_lat.T, self.C_lat.T, des_obs_poles_lat).gain_matrix.T

        print('K_lon: ', self.K_lon)
        print('ki_lon: ', self.ki_lon)
        print('L_lon^T: ', self.L_lon.T)
        print('K_lat: ', self.K_lat)
        print('ki_lat: ', self.ki_lat)
        print('L_lat^T: ', self.L_lat.T)

        self.xhat_lon = np.array([[0.0], [0.0]])
        self.xhat_lat = np.array([[0.0], [0.0], [0.0], [0.0]])

        self.integrator_z = 0.0  # integrator on position z
        self.error_z_d1 = 0.0  # error signal delayed by 1 sample
        self.integrator_h = 0.0  # integrator on altitude h
        self.error_h_d1 = 0.0  # error signal delayed by 1 sample
        self.F_d1 = 0.0  # Force signal delayed by 1 sample
        self.tau_d1 = 0.0  # torque signal delayed by 1 sample
        self.F_limit = 2*P.fmax
        self.tau_limit = P.tau_max
        self.Ts = P.Ts

    def update(self, r, y):
        z_r = r[0][0]
        h_r = r[1][0]
        y_lat = np.array([[y[0][0]],[y[2][0]]])
        y_lon = y[1][0]

        # update the observers
        xhat_lat = self.update_lat_observer(y_lat)
        xhat_lon = self.update_lon_observer(y_lon)
        z_hat = xhat_lat[0][0]
        h_hat = xhat_lon[0][0]
        theta_hat = xhat_lat[1][0]

        # integrate error
        error_z = z_r - z_hat
        self.integrator_z += (P.Ts/2.0)*(error_z + self.error_z_d1)
        self.error_z_d1 = error_z
        error_h = h_r - h_hat
        self.integrator_h += (P.Ts/2.0)*(error_h + self.error_h_d1)
        self.error_h_d1 = error_h

        # Compute the state feedback controllers
        F_tilde = -self.K_lon @ xhat_lon - self.ki_lon * self.integrator_h
        F = self.Fe / np.cos(theta_hat) + F_tilde[0]
        F_sat = saturate(F, self.F_limit)
        self.integratorAntiWindup(F_sat, F, self.ki_lon, self.integrator_h)

        tau = -self.K_lat @ xhat_lat - self.ki_lat*self.integrator_z
        tau_sat = saturate(tau[0], self.tau_limit)
        self.integratorAntiWindup(tau_sat, tau, self.ki_lat, self.integrator_z)

        u = np.array([[F_sat], [tau_sat]])
        self.F_d1 = F_sat
        self.tau_d1 = tau_sat

        return u, xhat_lat, xhat_lon
    
    def integratorAntiWindup(self, u_sat, u_unsat, ki, integrator):
        if ki != 0.0:
            integrator = integrator + P.Ts/ki*(u_sat-u_unsat)

    def update_lat_observer(self, y_m):
        # update the observer using RK4 integration
        F1 = self.observer_f_lat(self.xhat_lat, y_m)
        F2 = self.observer_f_lat(self.xhat_lat + self.Ts / 2 * F1, y_m)
        F3 = self.observer_f_lat(self.xhat_lat + self.Ts / 2 * F2, y_m)
        F4 = self.observer_f_lat(self.xhat_lat + self.Ts * F3, y_m)
        self.xhat_lat += self.Ts / 6 * (F1 + 2 * F2 + 2 * F3 + F4)

        return self.xhat_lat

    def observer_f_lat(self, x_hat, y_m):
        # xhatdot = A*xhat + B*u + L(y-C*xhat)
        xhat_dot = self.A_lat @ x_hat \
                   + self.B_lat * self.tau_d1 \
                   + self.L_lat @ (y_m - self.C_lat @ x_hat)
        
        return xhat_dot

    def update_lon_observer(self, y_m):
        # update the observer using RK4 integration
        F1 = self.observer_f_lon(self.xhat_lon, y_m)
        F2 = self.observer_f_lon(self.xhat_lon + self.Ts / 2 * F1, y_m)
        F3 = self.observer_f_lon(self.xhat_lon + self.Ts / 2 * F2, y_m)
        F4 = self.observer_f_lon(self.xhat_lon + self.Ts * F3, y_m)
        self.xhat_lon += self.Ts / 6 * (F1 + 2 * F2 + 2 * F3 + F4)

        return self.xhat_lon

    def observer_f_lon(self, x_hat, y_m):
        # xhatdot = A*xhat + B*u + L(y-C*xhat)
        xhat_dot = self.A_lon @ x_hat \
                   + self.B_lon * (self.F_d1 - self.Fe) \
                   + self.L_lon @ (y_m - self.C_lon @ x_hat)
        
        return xhat_dot


def saturate(u, limit):
    if abs(u) > limit:
        u = limit*np.sign(u)
    return u

