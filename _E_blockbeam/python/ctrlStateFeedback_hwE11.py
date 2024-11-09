import numpy as np
import blockbeamParam as P
import control as cnt

class ctrlStateFeedback:
    def __init__(self):
        # tuning parameters
        tr_th = 0.1           # Rise time for inner loop (z)
        tr_z = tr_th * 10.0   # 10 times rise time of theta
        zeta_th = 0.707       # inner loop Damping Coefficient
        zeta_z = 0.707        # outer loop Damping Coefficient

        A = np.array([[0.0, 0.0, 1.0, 0.0],
                      [0.0, 0.0, 0.0, 1.0],
                      [0.0, -P.g, 0.0, 0.0],
                      [0.0, 0.0, 0.0, 0.0]])
        B = np.array([[0.0],
                      [0.0],
                      [0.0],
                      [P.length/(P.length**2 * P.m2 / 3 + P.m1 * P.ze**2)]
                      ])
        C = np.array([[1.0, 0.0, 0.0, 0.0],
                      [0.0, 1, 0.0, 0.0]])

        # gain calculation
        wn_th = 2.2 / tr_th
        wn_z = 2.2 / tr_z
        des_char_poly = np.convolve(
            [1, 2*zeta_z*wn_z, wn_z**2],
            [1, 2*zeta_th*wn_th, wn_th**2]
        )
        des_poles = np.roots(des_char_poly)
        # Compute the gains if the system is controllable
        if np.linalg.matrix_rank(cnt.ctrb(A, B)) != 4:
            print('The system is not controllable')
        else:
            self.K = cnt.acker(A, B, des_poles)
            Cr = np.array([[1.0, 0.0, 0.0, 0.0]])
            self.kr = -1.0 / (Cr @ np.linalg.inv(A-B @ self.K) @ B)
        # print the gains to the terminal

        print('K: ', self.K)
        print('kr: ', self.kr)

    def update(self, z_r, state):
        z = state[0][0]

        F_fl = P.m2*P.g/2.0 + P.m1*P.g*z / P.length # define equilibrium force
        F_tilde = -self.K @ state + self.kr * z_r
        F = F_tilde + F_fl

        # return saturated force
        return saturate(F[0][0], P.Fmax)

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
