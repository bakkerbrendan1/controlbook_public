import numpy as np
import blockbeamParam as P
import control as cnt

class ctrlStateFeedbackInt:
    def __init__(self):
        # tuning parameters
        tr_th = 0.05           # Rise time for inner loop (z)
        tr_z = tr_th * 10.0   # 10 times rise time of theta
        zeta_th = 0.707       # inner loop Damping Coefficient
        zeta_z = 0.707        # outer loop Damping Coefficient
        integrator_pole = 0.7  # integrator pole

        A = P.A
        B = P.B
        C = P.C
        Cr = np.array([[1.0, 0.0, 0.0, 0.0]])

        A1 = np.vstack((
                np.hstack((A, np.zeros((4,1)))),
                np.hstack((-Cr, np.zeros((1,1))))))
        B1 = np.vstack((B, np.zeros((1,1))))

        # gain calculation
        wn_th = 2.2 / tr_th
        wn_z = 2.2 / tr_z
        des_char_poly = np.convolve(
            np.convolve([1, 2*zeta_z*wn_z, wn_z**2],
                        [1, 2*zeta_th*wn_th, wn_th**2]),
            np.poly([integrator_pole])
        )
        des_poles = np.roots(des_char_poly)
        # Compute the gains if the system is controllable
        if np.linalg.matrix_rank(cnt.ctrb(A1, B1)) != 5:
            print('The system is not controllable')
        else:
            K1 = cnt.acker(A1, B1, des_poles)
            self.K = K1[0][0:4]
            self.ki = K1[0][4]
        # print the gains to the terminal

        print('K: ', self.K)
        print('kr: ', self.ki)

        #--------------------------------------------------
        # saturation limits
        theta_max = 30.0 * np.pi / 180.0  # Max theta, rads
        #--------------------------------------------------
        # variables to implement integrator
        self.integrator_z = 0.0  # integrator
        self.error_z_d1 = 0.0  # error signal delayed by 1 sample


    def update(self, z_r, state):
        z = state[0][0]

        # integrate error
        error_z = z_r - z

        # integrate error
        self.integrator_z = self.integrator_z \
                + (P.Ts / 2.0) * (error_z + self.error_z_d1)
        self.error_z_d1 = error_z

        # Compute the state feedback controller
        F_fl = P.m2*P.g/2.0 + P.m1*P.g*z / P.length # define equilibrium force
        F_tilde = -self.K @ state + self.ki * self.integrator_z
        F = F_tilde + F_fl

        # return saturated force
        return saturate(F[0], P.Fmax)

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
