import numpy as np
import control as cnt
import pendulumParam as P

class ctrlStateFeedback:
    def __init__(self):
        #--------------------------------------------------
        # State Feedback Control Design
        #--------------------------------------------------
        # tuning parameters
        tr_z = 1.0        # rise time for position
        tr_theta = 0.5    # rise time for angle
        zeta_z = 0.707  # damping ratio position
        zeta_th = 0.707  # damping ratio angle
        # State Space Equations
        # xdot = A*x + B*u
        # y = C*x
        A = np.array([
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
            [0.0, -3 * P.m1 * P.g / 4 / (.25 * P.m1 + P.m2),
                -P.b / (.25 * P.m1 + P.m2), 0.0],
            [0.0, 
                3*(P.m1+P.m2) * P.g/2/(.25 * P.m1 + P.m2)/P.ell,
                3 * P.b / 2 / (.25 * P.m1 + P.m2) / P.ell, 0.0]
            ])
        B = np.array([[0.0],
                      [0.0],
                      [1 / (.25 * P.m1 + P.m2)],
                      [-3.0 / 2 / (.25 * P.m1 + P.m2) / P.ell]])
        C = np.array([[1.0, 0.0, 0.0, 0.0],
                      [0.0, 1.0, 0.0, 0.0]])
        # gain calculation
        Q = np.array([[10.0, 0, 0, 0],
                      [0, 1.0, 0, 0],
                      [0, 0, 0.1, 0],
                      [0, 0, 0, 0.1]])
        # Q = np.eye(4)
        R = np.array([[0.1]])

        # Compute the gains if the system is controllable
        if np.linalg.matrix_rank(cnt.ctrb(A, B)) != 4:
            print("The system is not controllable")
        else:
            # self.K = cnt.acker(A, B, des_poles)
            self.K, _S, _E = cnt.lqr(A, B, Q, R)
            Cr = np.array([[1.0, 0.0, 0.0, 0.0]])
            self.kr = -1.0 / (Cr @ np.linalg.inv(A-B @ self.K) @ B)
        # print gains to terminal
        print('K: ', self.K)
        print('kr: ', self.kr)

    def update(self, z_r, x):
        # Compute the state feedback controller
        F_unsat = -self.K @ x + self.kr * z_r
        F = saturate(F_unsat[0][0], P.F_max)
        return F


def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u

