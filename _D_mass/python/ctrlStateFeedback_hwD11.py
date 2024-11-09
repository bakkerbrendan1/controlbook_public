import numpy as np
import massParam as P
import control as cnt

class ctrlStateFeedback:
    def __init__(self):
        # tuning parameters
        tr = 1.0
        zeta = 0.707
        self.ki = 1.5 # integrator gain

        # State Space Equations
        # xdot = A*x + B*u
        # y = C*x
        A = np.array([[0.0, 1.0],
                      [-0.6, -0.1]])
        B = np.array([[0.0],
                      [0.2]])
        C = np.array([1, 0]).T

        # gain calculation
        wn = 2.2 / tr
        des_char_poly = [1, 2*zeta*wn, wn**2]
        des_poles = np.roots(des_char_poly)

        # compute the gains if the system is controllable
        if np.linalg.matrix_rank(cnt.ctrb(A, B)) != 2:
            print('The system is not controllable!')
        else:
            self.K = (cnt.acker(A, B, des_poles))
            self.kr = -1.0 / (C @ np.linalg.inv(A - B @ self.K) @ B)

        print('K: ', self.K)
        print('kr: ', self.kr)
        print(des_poles)

    def update(self, z_r, x):
        # update state variables, estimate zdot with dirty derivative
        z = x[0][0]

        # feedback linearized force
        F_fl = P.k * z
        # # compute the linearized force using PD control
        # F_tilde = self.kp * (z_r - z) - self.kd * self.zdot
        # PID control
        F_tilde = -self.K @ x + self.kr * z_r

        # compute the final force and saturate
        F_unsat = F_fl + F_tilde[0][0]
        F = saturate(F_unsat, P.F_max)
        return F


def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u







