import numpy as np
import VTOLParam as P


class ctrlPID:
    def __init__(self):
        self.sigma = 0.05
        self.beta = (2 * self.sigma - P.Ts) \
            / (2 * self.sigma + P.Ts) 
        self.ki_z = -0.18  # integral gain for outer loop
        self.ki_h = 0.2
        ####################################################
        #       PD Control: Time Design Strategy
        ####################################################
        # tuning parameters
        tr_h = 1.0 # originally 4.0 # rise time for altitude - tuned for best rise time without saturation
        zeta_h = 0.707 # damping ratio for altitude
        tr_z = 1.0 # originally 4.0 # rise time for outer lateral loop (position) - tuned for best rise time w/o saturation
        zeta_z = 0.707
        zeta_th = 0.707

        # saturation limits
        self.theta_max = 10.0 * np.pi / 180.0 # Max theta, rads

        # equilibrium force
        self.Fe = (P.mc + 2.0 * P.mr) * P.g

        #---------------------------------------------------
        #          Longitudinal control (altitude)
        #---------------------------------------------------
        wn_h = 2.2/tr_h  # natural frequency
        Delta_cl_d = [1, 2*zeta_h*wn_h, wn_h**2.0]  # desired closed loop
        self.kp_h = Delta_cl_d[2] * (P.mc+2.0*P.mr) # kp - altitude
        self.kd_h = Delta_cl_d[1] * (P.mc+2.0*P.mr) # kd = altitude

        #---------------------------------------------------
        #           Lateral control (left-right)
        #---------------------------------------------------
        # PD design for inner loop
        M = 10.0 # time separation between inner and outer lateral loops
        b0 = 1.0/(P.Jc + 2.0*P.mr*P.d**2)
        tr_th = tr_z/M
        wn_th = 2.2/tr_th
        self.kp_th = wn_th**2.0/b0
        self.kd_th = 2.0*zeta_th*wn_th/b0

        # PD design for outer loop
        b1 = -self.Fe/(P.mc+2.0*P.mr)
        a1 = P.mu/(P.mc+2.0*P.mr)
        wn_z = 2.2/tr_z
        self.kp_z = wn_z**2.0/b1
        self.kd_z = (2.0*zeta_z*wn_z - a1)/b1

        # print control gains to terminal        
        print('kp_z: ', self.kp_z)
        print('kd_z: ', self.kd_z)
        print('ki_z: ', self.ki_z)

        print('kp_h: ', self.kp_h)
        print('kd_h: ', self.kd_h)

        print('kp_th: ', self.kp_th)
        print('kd_th: ', self.kd_th)

        self.integrator_z = 0.
        self.integrator_h = 0.
        self.error_z_d1 = 0.
        self.error_h_d1 = 0.
        self.z_dot = 0.
        self.z_d1 = 0.
        self.theta_dot = 0.
        self.theta_d1 = 0.
        self.h_dot = 0.
        self.h_d1 = 0.

    def update(self, reference, state):
        z_r = reference[0][0]
        h_r = reference[1][0]
        z = state[0][0]
        h = state[1][0]
        theta = state[2][0]
        # self.z_dot = state[3][0]
        self.h_dot = state[4][0]
        # self.theta_dot = state[5][0]

        ######################################
        #      Altitude PID control: h
        ######################################
        # calculate control for the altitude loop
        error_h = h_r - h
        self.integrator_h = self.integrator_h + (P.Ts/2) * (error_h + self.error_h_d1)
        self.h_dot = self.beta * self.h_dot \
            + (1 - self.beta) * ((h - self.h_d1) / P.Ts)
        # F_tilde = self.kp_h * (error_h) - self.kd_h * self.h_dot
        F_tilde = self.kp_h * error_h \
            + self.ki_h * self.integrator_h \
            - self.kd_h * self.h_dot
        F = saturate( F_tilde * self.Fe, 2*P.fmax)

        ######################################
        #      Outer loop PID control: z
        ######################################
        error_z = z_r - z
        self.integrator_z = self.integrator_z \
            + (P.Ts / 2) * (error_z + self.error_z_d1)
        self.z_dot = self.beta * self.z_dot \
            + (1 - self.beta) * ((z - self.z_d1) / P.Ts)
        # PID control - unsaturated
        theta_r_unsat = self.kp_z * error_z \
                + self.ki_z * self.integrator_z \
                - self.kd_z * self.z_dot
        # saturate theta_r
        theta_r = saturate(theta_r_unsat, self.theta_max)
        # integrator anti-windup
        if self.ki_z != 0.0:
            self.integrator_z = self.integrator_z \
                + P.Ts / self.ki_z * (theta_r - theta_r_unsat)
        # calculate thera_r from the outer loop
        # theta_r = saturate( self.kp_z * (z_r-z) - self.kd_z * zdot, self.theta_max )

        ######################################
        #      Inner loop PID control: theta
        ######################################
        error_th = theta_r - theta
        # use dirty derivative
        self.theta_dot = self.beta * self.theta_dot \
            + (1 - self.beta) * ((theta - self.theta_d1) / P.Ts)
        # PD control on theta
        theta_r = self.kp_z * (error_z) - self.kd_z * self.z_dot
        # theta_r = self.filter.update(theta_r)
        # calculate torque from the inner loop
        tau = saturate( self.kp_th * (error_th) - self.kd_th*self.theta_dot, 2*P.fmax*P.d)

        motor_thrusts = P.mixing @ np.array([[F], [tau]])


        self.error_z_d1 = error_z
        self.error_h_d1 = error_h
        self.z_d1 = z
        self.theta_d1 = theta
        self.h_d1 = h
        return motor_thrusts


def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u
