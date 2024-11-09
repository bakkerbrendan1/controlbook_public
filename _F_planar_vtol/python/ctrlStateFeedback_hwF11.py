import numpy as np
import VTOLParam as P
import control as cnt


class ctrlStateFeedback:
    def __init__(self):
        # saturation limits
        self.theta_max = 10.0 * np.pi / 180.0 # Max theta, rads
        # equilibrium force
        self.Fe = (P.mc + 2.0 * P.mr) * P.g
        ####################################################
        #       PD Control: Time Design Strategy
        ####################################################
        # tuning parameters
        tr_h = 2.0 # originally 4.0 # rise time for altitude - tuned for best rise time without saturation
        zeta_h = 0.707 # damping ratio for altitude
        tr_z = 2.0 # originally 4.0 # rise time for outer lateral loop (position) - tuned for best rise time w/o saturation
        zeta_z = 0.707
        tr_th = tr_z / 10.0
        zeta_th = 0.707

        A_lat = np.array([[0, 0, 1, 0],
                          [0, 0, 0, 1],
                          [0, -self.Fe/(P.mc + 2*P.mr), -P.mu/(P.mc + 2*P.mr), 0],
                          [0, 0, 0, 0]])
        B_lat = np.array([[0],
                          [0],
                          [0],
                          [1/(P.Jc + 2*P.mr*P.d**2)]])
        C_lat = np.array([[1, 0, 0, 0],
                          [0, 1, 0, 0]])
        
        A_lon = np.array([[0, 1],
                          [0, 0]])
        B_lon = np.array([[0],
                          [1/(P.mc + 2*P.mr)]])
        C_lon = np.array([1, 0])

        A = np.array([[0, 0, 0, 1, 0, 0],
                      [0, 0, 0, 0, 1, 0],
                      [0, 0, 0, 0, 0, 1],
                      [0, 0, -self.Fe / (P.mc + 2 * P.mr), -P.mu / (P.mc + 2 * P.mr), 0, 0],
                      [0, 0, 0, 0, 0, 0],
                      [0, 0, 0, 0, 0, 0]])
        
        B = np.array([[0, 0],
                      [0, 0],
                      [0, 0],
                      [0, 0],
                      [1 / (P.mc + 2 * P.mr), 0],
                      [0, 1 / (P.Jc + 2 * P.d**2 * P.mr)]])
        
        C = np.array([[1,0,0,0,0,0],
                      [0,1,0,0,0,0],
                      [0,0,1,0,0,0]])

        # gain calculation
        wn_h = 2.2/tr_h
        wn_z = 2.2/tr_z
        wn_th = 2.2/tr_th

        des_char_poly_lat = np.convolve([1, 2*zeta_th*wn_th, wn_th**2],
                                        [1, 2*zeta_z*wn_z, wn_z**2])
        des_poles_lat = np.roots(des_char_poly_lat)

        des_char_poly_lon = [1, 2*zeta_h*wn_h, wn_h**2]
        des_poles_lon = np.roots(des_char_poly_lon)

        des_char_poly = np.convolve([1, 2*zeta_th*wn_th, wn_th**2],
                                    [1, 2*zeta_z*wn_z, wn_z**2])
        des_char_poly = np.convolve(des_char_poly, 
                                    [1, 2*zeta_h*wn_h, wn_h**2])
        des_poles = np.roots(des_char_poly)

        # Compute the gains if the system is controllable
        if np.linalg.matrix_rank(cnt.ctrb(A_lat, B_lat)) != 4:
            print('The latitude system is not controllable')
        else: 
            self.K_lat = cnt.acker(A_lat, B_lat, des_poles_lat)
            Cr_lat = np.array([[1.0, 0.0, 0.0, 0.0]])
            self.kr_lat = -1.0 / (Cr_lat @ np.linalg.inv(A_lat - B_lat @ self.K_lat) @ B_lat)

        if np.linalg.matrix_rank(cnt.ctrb(A_lon, B_lon)) != 2:
            print('The longitude system is not controllable')
        else: 
            self.K_lon = cnt.acker(A_lon, B_lon, des_poles_lon)
            Cr_lon = np.array([[1.0, 0.0]])
            self.kr_lon = -1.0 / (Cr_lon @ np.linalg.inv(A_lon - B_lon @ self.K_lon) @ B_lon)

        # if np.linalg.matrix_rank(cnt.ctrb(A, B)) != 6:
        #     print('The system is not controllable')
        # else: 
        #     self.K = cnt.acker(A, B, des_poles)
        #     Cr = np.array([[1.0, 0.0, 0.0, 0.0, 0.0, 0.0]])
        #     self.kr = -1.0 / (Cr @ np.linalg.inv(A - B @ self.K) @ B)

        # print control gains to terminal        
        print('K_lat: ', self.K_lat)
        print('kr_lat: ', self.kr_lat)
        print('K_lon: ', self.K_lon)
        print('kr_lon: ', self.kr_lon)
        # print('K: ', self.K)
        # print('kr: ', self.kr)

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



        F_tilde = -self.K_lon @ state_lon + self.kr_lon * h_r
        F = saturate(F_tilde[0] + self.Fe/np.cos(theta), 2*P.fmax)

        tau_unsat = -self.K_lat @ state_lat + self.kr_lat * z_r
        tau = saturate(tau_unsat[0], 2*P.fmax*P.d)

        # motor_thrusts = P.mixing @ np.array([[F], [tau]])
        motor_thrusts = P.mixing @ np.array([F, tau])

        return motor_thrusts


def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u
