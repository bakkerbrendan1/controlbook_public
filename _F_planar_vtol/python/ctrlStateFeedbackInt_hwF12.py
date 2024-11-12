import numpy as np
import VTOLParam as P
import control as cnt


class ctrlStateFeedbackInt:
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
        int_pole_z = 0.38

        A_lat = P.A_lat
        B_lat = P.B_lat
        C_lat = P.C_lat
        Cr_lat = np.array([[1.0, 0.0, 0.0, 0.0]])
        A1_lat = np.vstack((np.hstack((A_lat, np.zeros((np.size(A_lat,1),1)))), 
                            np.hstack((-Cr_lat, np.array([[0.0]]))) ))
        B1_lat = np.vstack( (B_lat, 0.0) )

        A_lon = P.A_lon
        B_lon = P.B_lon
        C_lon = P.C_lon
        Cr_lon = np.array([[1.0, 0.0]])

        A1_lon = np.vstack((np.hstack((A_lon, np.zeros((np.size(A_lon,1),1)))), 
                        np.hstack((-Cr_lon, np.array([[0.0]]))) ))
        B1_lon = np.vstack( (B_lon, 0.0) )

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

        # print control gains to terminal        
        print('K_lat: ', self.K_lat)
        print('ki_lat: ', self.ki_lat)
        print('K_lon: ', self.K_lon)
        print('ki_lon: ', self.ki_lon)

        #--------------------------------------------------
        # variables to implement the integrator
        self.integrator_z = 0.0
        self.integrator_h = 0.0  # integrators
        self.error_z_d1 = 0.0
        self.error_h_d1 = 0.0  # errors delayed by 1 sample

    def update(self, reference, state):
        z_r = reference[0][0]
        h_r = reference[1][0]
        z = state[0][0]
        h = state[1][0]
        theta = state[2][0]

        state_lon = np.array([state[1][0], state[4][0]])

        state_lat = np.array([state[0][0],
                              state[2][0],
                              state[3][0],
                              state[5][0]])

        error_z = z_r - z
        error_h = h_r - h

        self.integrator_z = self.integrator_z \
                            + (P.Ts / 2.0) * (error_z + self.error_z_d1)
        self.error_z_d1 = error_z
        self.integrator_h = self.integrator_h \
                            + (P.Ts / 2.0) * (error_h + self.error_h_d1)
        self.error_h_d1 = error_h

        F_tilde = -self.K_lon @ state_lon + self.ki_lon * self.integrator_h
        F = saturate(F_tilde + P.Fe/np.cos(theta), 2*P.fmax)

        tau_unsat = -self.K_lat @ state_lat + self.ki_lat * self.integrator_z
        tau = saturate(tau_unsat, 2*P.fmax*P.d)

        # motor_thrusts = P.mixing @ np.array([[F], [tau]])
        motor_thrusts = P.mixing @ np.array([[F], [tau]])

        return motor_thrusts


def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u
