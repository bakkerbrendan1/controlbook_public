import numpy as np
from scipy import signal
import VTOLParam as P
import control as cnt


class ctrlObserver:
    def __init__(self):
        # saturation limits
        self.theta_max = 10.0 * np.pi / 180.0 # Max theta, rads
        ####################################################
        #       PD Control: Time Design Strategy
        ####################################################
        # tuning parameters
        tr_h = 2.0 # originally 4.0 # rise time for altitude - tuned for best rise time without saturation
        zeta_h = 0.707 # damping ratio for altitude
        tr_z = 1.0 # originally 4.0 # rise time for outer lateral loop (position) - tuned for best rise time w/o saturation
        zeta_z = 0.707
        tr_th = tr_z / 10.0
        zeta_th = 0.707
        int_pole_h = 0.2 # integrator pole
        int_pole_z = 0.38 # if tr_z= 2.0, int_pole_z=0.19
        # pick observer poles
        tr_h_obs = tr_h/6.0
        tr_z_obs = tr_z/6.0
        tr_th_obs = tr_th/6.0


        self.A_lat = P.A_lat
        self.B_lat = P.B_lat
        self.C_lat = P.C_lat
        Cr_lat = np.array([[1.0, 0.0, 0.0, 0.0]])
        # form the lateral augmented system
        A1_lat = np.vstack((np.hstack((self.A_lat, np.zeros((np.size(self.A_lat,1),1)))), 
                            np.hstack((-Cr_lat, np.array([[0.0]]))) ))
        B1_lat = np.vstack( (self.B_lat, 0.0) )

        self.A_lon = P.A_lon
        self.B_lon = P.B_lon
        self.C_lon = P.C_lon
        # form the longitudinal augmented system
        Cr_lon = np.array([[1.0, 0.0]])

        A1_lon = np.vstack((np.hstack((self.A_lon, np.zeros((np.size(self.A_lon,1),1)))), 
                        np.hstack((-Cr_lon, np.array([[0.0]]))) ))
        B1_lon = np.vstack( (self.B_lon, 0.0) )

        # gain calculation
        wn_h = 2.2/tr_h
        wn_z = 2.2/tr_z
        wn_th = 2.2/tr_th

        des_char_poly_lat = np.convolve(
                np.convolve([1, 2*zeta_th*wn_th, wn_th**2],
                            [1, 2*zeta_z*wn_z, wn_z**2]),
                [1, -int_pole_z])
        des_poles_lat = np.roots(des_char_poly_lat)

        des_char_poly_lon = np.convolve([1, 2*zeta_h*wn_h, wn_h**2],
                                        [1, -int_pole_h])
        des_poles_lon = np.roots(des_char_poly_lon)

        # Compute the gains if the system is controllable
        if np.linalg.matrix_rank(cnt.ctrb(A1_lat, B1_lat)) != 5:
            print('The latitude system is not controllable')
            print(A1_lat)
        else: 
            K1_lat = cnt.place(A1_lat, B1_lat, des_poles_lat)
            self.K_lat = K1_lat[0][0:4]
            self.ki_lat = K1_lat[0][4]

        if np.linalg.matrix_rank(cnt.ctrb(A1_lon, B1_lon)) != 3:
            print('The longitude system is not controllable')
        else: 
            K1_lon = cnt.place(A1_lon, B1_lon, des_poles_lon)
            self.K_lon = K1_lon[0][0:2]
            self.ki_lon = K1_lon[0][2]

        # compute observer gains
        wn_h_obs = 2.2 / tr_h_obs
        wn_z_obs = 2.2 / tr_z_obs
        wn_th_obs = 2.2 / tr_th_obs
        des_obs_char_poly_lat = np.convolve(
                [1, 2 * zeta_z * wn_z_obs, wn_z_obs**2],
                [1, 2 * zeta_th * wn_th_obs, wn_th_obs**2])
        des_obs_poles_lat = np.roots(des_obs_char_poly_lat)
        # Compute the observer gains if the system is observable
        if np.linalg.matrix_rank(cnt.ctrb(self.A_lat.T, self.C_lat.T)) != 4:
            print('The system is not observable')
        else:
            self.L_lat = signal.place_poles(self.A_lat.T, self.C_lat.T, 
                                        des_obs_poles_lat).gain_matrix.T
            
        des_obsv_char_poly_lon = [1, 2*zeta_h*wn_h_obs, wn_h_obs**2]
        des_obsv_poles_lon = np.roots(des_obsv_char_poly_lon)
        # Compute the gains if the system is controllable
        if np.linalg.matrix_rank(cnt.ctrb(self.A_lon.T, self.C_lon.T)) != 2:
            print('The system is not observable')
        else:
            self.L_lon = cnt.acker(self.A_lon.T, self.C_lon.T, des_obsv_poles_lon).T

        # print control gains to terminal        
        print('K_lat: ', self.K_lat)
        print('ki_lat: ', self.ki_lat)
        print('L^T_lat: ', self.L_lat)
        print('K_lon: ', self.K_lon)
        print('ki_lon: ', self.ki_lon)
        print('L^T_lon: ', self.L_lon)

        #--------------------------------------------------
        # variables to implement the integrator
        self.integrator_z = 0.0
        self.integrator_h = 0.0  # integrators
        self.error_z_d1 = 0.0
        self.error_h_d1 = 0.0  # errors delayed by 1 sample
        # estimated state variables
        self.xhat_lon = np.array([[0.0], [0.0]])
        self.xhat_lat = np.array([[0.0], [0.0], [0.0], [0.0]])
        self.F_d1 = 0.0
        self.tau_d1 = 0.0

    def update(self, r, y):
        z_r = r[0][0]
        h_r = r[1][0]
        y_lat = np.array([[y[0][0]], [y[2][0]]])
        y_lon = y[1][0]

        # update the observers
        x_hat_lat = self.update_lat_observer(y_lat)
        x_hat_lon = self.update_lon_observer(y_lon)
        z_hat = x_hat_lat[0][0]
        h_hat = x_hat_lon[0][0]
        theta_hat = x_hat_lat[1][0]

        error_z = z_r - z_hat
        error_h = h_r - h_hat

        self.integrator_z = self.integrator_z \
                            + (P.Ts / 2.0) * (error_z + self.error_z_d1)
        self.error_z_d1 = error_z
        self.integrator_h = self.integrator_h \
                            + (P.Ts / 2.0) * (error_h + self.error_h_d1)
        self.error_h_d1 = error_h

        F_tilde = -self.K_lon @ x_hat_lon - self.ki_lon * self.integrator_h
        F_unsat = F_tilde + P.Fe/np.cos(theta_hat)
        F = saturate(F_unsat, 2*P.fmax)
        self.integratorAntiWindup(F, F_unsat, self.ki_lon, self.integrator_h)

        tau_unsat = -self.K_lat @ x_hat_lat - self.ki_lat * self.integrator_z
        tau = saturate(tau_unsat[0], P.tau_max)
        self.integratorAntiWindup(tau, tau_unsat, self.ki_lat, self.integrator_z)


        motor_thrusts = P.mixing @ np.array([[F[0]], [tau]])

        return motor_thrusts, x_hat_lat, x_hat_lon
    
    def integratorAntiWindup(self, u_sat, u_unsat, ki, integrator):
        if ki != 0.0:
            integrator = integrator + P.Ts/ki*(u_sat-u_unsat)

    def update_lat_observer(self, y_m):
        # update the observer using RK4 integration
        F1 = self.observer_f_lat(self.xhat_lat, y_m)
        F2 = self.observer_f_lat(self.xhat_lat + P.Ts / 2 * F1, y_m)
        F3 = self.observer_f_lat(self.xhat_lat + P.Ts / 2 * F2, y_m)
        F4 = self.observer_f_lat(self.xhat_lat + P.Ts * F3, y_m)
        self.xhat_lat += P.Ts / 6 * (F1 + 2 * F2 + 2 * F3 + F4)

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
        F2 = self.observer_f_lon(self.xhat_lon + P.Ts / 2 * F1, y_m)
        F3 = self.observer_f_lon(self.xhat_lon + P.Ts / 2 * F2, y_m)
        F4 = self.observer_f_lon(self.xhat_lon + P.Ts * F3, y_m)
        self.xhat_lon += P.Ts / 6 * (F1 + 2 * F2 + 2 * F3 + F4)

        return self.xhat_lon

    def observer_f_lon(self, x_hat, y_m):
        # xhatdot = A*xhat + B*u + L(y-C*xhat)
        xhat_dot = self.A_lon @ x_hat \
                   + self.B_lon * (self.F_d1 - P.Fe) \
                   + self.L_lon @ (y_m - self.C_lon @ x_hat)
        
        return xhat_dot


def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u
