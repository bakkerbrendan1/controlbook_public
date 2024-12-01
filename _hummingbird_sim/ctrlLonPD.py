import numpy as np
import hummingbirdParam as P

class ctrlLonPD:
    def __init__(self):
        # tuning parameters
        tr_pitch = 1.0
        zeta_pitch = 0.707
        # gain calculation
        b_theta = P.ellT/(P.m1 * P.ell1**2 + P.m2 * P.ell2**2 + P.J1y + P.J2y)
        #print('b_theta: ', b_theta)
        wn_pitch = 0.5 * (np.pi) / (tr_pitch * np.sqrt(1-zeta_pitch**2))
        self.kp_pitch = wn_pitch**2 / b_theta
        self.kd_pitch = 2 * wn_pitch * zeta_pitch / b_theta
        # print gains to terminal
        print('kp_pitch: ', self.kp_pitch)
        print('kd_pitch: ', self.kd_pitch) 
        # sample rate of the controller
        self.Ts = P.Ts
        pass 

    def update(self, ref, x):
        theta = x[1]
        thetadot = x[4]
        phi_r = ref[0]
        theta_r = ref[1]
        psi_r = ref[2]
        
        force_equilibrium = (P.m1*P.ell1 + P.m2*P.ell2)*P.g/P.ellT
        f_tilde = self.kp_pitch * (theta_r-theta) - self.kd_pitch * thetadot
        force = (force_equilibrium + f_tilde)[0]
        torque = 0.0


        # convert force and torque to pwm signals
        ul = 1/(2*P.km) * (force + torque/P.d)
        ur = 1/(2*P.km) * (force - torque/P.d)
        pwm = np.array([[ul], 
                        [ur]])
        pwm = saturate(pwm, 0, 1) 
        print(pwm)
        return pwm


def saturate(u, low_limit, up_limit):
    if isinstance(u, float) is True:
        u = np.max((np.min((u, up_limit)), low_limit))
    else:
        for i in range(0, u.shape[0]):
            print('u[i][0] = ', np.max((np.min((u[i][0], up_limit)), low_limit)))
            u[i][0] = np.max((np.min((u[i][0], up_limit)), low_limit))
    return u
