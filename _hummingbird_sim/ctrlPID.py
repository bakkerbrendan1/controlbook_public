import numpy as np
import hummingbirdParam as P

class ctrlPID:
    def __init__(self):
        # Roll, Yaw tuning parameters
        # Phi: inner loop
        # Psi: outer loop
        tr_psi = 2.
        tr_phi = tr_psi/10
        zeta_phi = 0.707
        zeta_psi = 0.707
        # self.ki_phi = 0.01
        self.ki_psi = 0.05

        # Pitch tuning parameters
        tr_pitch = 1.0
        zeta_pitch = 0.707
        self.ki_pitch = 0.15 # select integrator gain
        # DC gain = 1 for inner loop

        #---------------------------------------------------
        #              Inner Loop: Roll (phi ϕ)
        #---------------------------------------------------
        # gain calculation
        b_phi = 1/P.J1x
        wn_phi = 0.5 * (np.pi) / (tr_phi * np.sqrt(1-zeta_phi**2))
        self.kp_phi = wn_phi**2 / b_phi
        self.kd_phi = 2 * wn_phi * zeta_phi / b_phi
        # print gains to terminal
        print('kp_phi: ', self.kp_phi)
        print('kd_phi: ', self.kd_phi) 
        
        #---------------------------------------------------
        #              Outer Loop: Yaw (psi ψ)
        #---------------------------------------------------
        # gain calculation
        b_psi = P.ellT * P.km / (P.JT + P.J1z)
        wn_psi = 0.5 * (np.pi) / (tr_psi * np.sqrt(1-zeta_psi**2))
        self.kp_psi = wn_psi**2 / b_psi
        self.kd_psi = 2 * wn_psi * zeta_psi / b_psi
        # print gains to terminal
        print('kp_psi: ', self.kp_psi)
        print('kd_psi: ', self.kd_psi)

        #---------------------------------------------------
        #                 Pitch (theta θ)
        #---------------------------------------------------
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

        # sample rate of the controller
        self.Ts = P.Ts

        # initialize variables for integrator and differentiator
        self.integrator_pitch = 0.
        self.integrator_psi = 0.

        self.error_pitch_d1 = 0.
        self.error_psi_d1 = 0.


    def update(self, ref, x):
        phi = x[0]
        theta = x[1]
        psi = x[2]
        phidot = x[3]
        thetadot = x[4]
        psidot = x[5]

        phi_r = ref[0]
        theta_r = ref[1]
        psi_r = ref[2]
        
        #################################################
        #       Roll, Yaw PID control: psi, phi
        #################################################

        # ----------- Outer loop: psi -----------
        error_psi = psi_r - psi
        self.integrator_psi = self.integrator_psi \
                                + (P.Ts/2) * (error_psi + self.error_psi_d1)
        # phi_r = self.kp_psi * (psi_r - psi) - self.kd_psi * psidot
        phi_r = self.kp_psi * error_psi \
                + self.ki_psi * self.integrator_psi \
                - self.kd_psi * psidot
        

        # ----------- Inner loop: phi -----------
        torque = (self.kp_phi * (phi_r - phi) - self.kd_phi * phidot)[0]


        ######################################
        #      Pitch PID control: theta
        ######################################
        error_pitch = theta_r - theta
        self.integrator_pitch = self.integrator_pitch + (P.Ts/2) * (error_pitch + self.error_pitch_d1)
        # self.pitch_dot =  doesn't specify to use dirty derivative, use thetadot variable
        # Fe (eqn 4.2, page 24)
        force_equilibrium = (P.m1*P.ell1 + P.m2*P.ell2)*P.g*np.cos(theta)/P.ellT
        f_tilde = self.kp_pitch * error_pitch \
                    + self.ki_pitch * self.integrator_pitch \
                    - self.kd_pitch * thetadot
        force = (force_equilibrium + f_tilde)[0]

        
        # convert force and torque to pwm signals
        ul = (1/(2*P.km) * (force + torque/P.d))
        ur = (1/(2*P.km) * (force - torque/P.d))
        pwm = np.array([[ul], 
                        [ur]])
        pwm = saturate(pwm, 0, 1) 


        # update variables
        self.error_pitch_d1 = error_pitch
        self.error_psi_d1 = error_psi
        return pwm


def saturate(u, low_limit, up_limit):
    if isinstance(u, float) is True:
        u = np.max((np.min((u, up_limit)), low_limit))
    else:
        for i in range(0, u.shape[0]):
            u[i][0] = np.max((np.min((u[i][0], up_limit)), low_limit))
    return u
