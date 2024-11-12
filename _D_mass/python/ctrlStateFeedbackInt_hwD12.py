import numpy as np
import massParam as P
import control as cnt

class ctrlStateFeedbackInt:
    def __init__(self):
        # tuning parameters
        tr = 0.25
        zeta = 0.707
        integrator_pole = 1.2

        # State Space Equations
        # xdot = A*x + B*u
        # y = C*x
        A = P.A
        B = P.B
        C = P.C
        Cr = np.array([[1.0, 0.0]])

        # form augmented system
        A1 = np.vstack((np.hstack((A, np.zeros((np.size(A,1),1)))), 
                        np.hstack((-Cr, np.array([[0.0]]))) ))
        B1 = np.vstack( (B, 0.0) )

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

        print('K: ', self.K)
        print('kr: ', self.ki)
        print(des_poles)
        #--------------------------------------------------
        # variables to implement the integrator
        self.integrator = 0.0  # integrator
        self.error_d1 = 0.0  # error delayed by 1 sample

    def update(self, z_r, x):
        z = x[0][0]

        # integrate error
        error = z_r - z
        self.integrator = self.integrator \
                          + (P.Ts / 2.0) * (error + self.error_d1)
        self.error_d1 = error

        # feedback linearized force
        F_fl = P.k * z
        # PID control
        F_tilde = -self.K @ x + self.ki * self.integrator

        # compute the final force and saturate
        F_unsat = F_fl + F_tilde[0]
        F = saturate(F_unsat, P.F_max)
        return F


def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u







