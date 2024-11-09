import numpy as np
import blockbeamParam as P

class ctrlPID:
    def __init__(self):
        # dirty derivative parameters
        self.sigma = 0.05 # cutoff freq for dirty derivative
        self.beta = (2 * self.sigma - P.Ts) / (2 * self.sigma + P.Ts)
        self.theta_max = 5.0 * np.pi / 180.0  # Max theta

        ####################################################
        #       PD Control: Time Design Strategy
        ####################################################
        # tuning parameters
        tr_z = 1             # Rise time for inner loop (z)
        zeta_th = 0.707       # inner loop Damping Coefficient
        zeta_z = 0.707        # outer loop Damping Coefficient
        self.ki_z = -0.4 # select integrator gain

        #---------------------------------------------------
        #                    Inner Loop
        #---------------------------------------------------
        # parameters of the open loop transfer function
        ze = P.length/2.0 # equilibrium position - center of beam
        b0 = P.length/(P.m2*P.length**2/3.0+P.m1*ze**2)
        M = 10.0              # Time scale separation 
        tr_th = tr_z/M # rise time for inner loop
        # coefficients for desired inner loop
        wn_th = 2.2 / tr_th     # Natural frequency
        # compute gains
        self.kp_th = wn_th**2/b0 # 1.83
        self.kd_th = 2.0*zeta_th*wn_th/b0 # 1.17
        DC_gain = 1.0
        #---------------------------------------------------
        #                    Outer Loop
        #---------------------------------------------------
        # coefficients for desired outer loop
        # desired rise time, s, defined in "tuning parameters"
        wn_z = 2.2 / tr_z  # desired natural frequency
        # compute gains
        self.kp_z = -wn_z**2/P.g # -0.0317
        self.kd_z = -2.0*zeta_z*wn_z/P.g # -0.00494
        # print control gains to terminal
        print('DC_gain', DC_gain)
        print('kp_th: ', self.kp_th)
        print('kd_th: ', self.kd_th)
        print('kp_z: ', self.kp_z)
        print('kd_z: ', self.kd_z)
        print('ki_z: ', self.ki_z)
        #---------------------------------------------------
        #                    zero canceling filter
        #---------------------------------------------------
        self.filter = zeroCancelingFilter(DC_gain)

        #---------------------------------------------------
        # initialize variables for integrator and differentiators
        #---------------------------------------------------
        self.integrator_z = 0.
        self.error_z_d1 = 0.
        self.z_dot = 0.
        self.z_d1 = 0.
        self.theta_dot = 0.
        self.theta_d1 = 0.

    def update(self, z_r, state):
        z = state[0][0]
        theta = state[1][0]
        # self.z_dot = state[2][0]
        # self.theta_dot = state[3][0]
        # the reference angle for theta comes from the
        # outer loop PD control
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
        # theta_r = self.filter.update(theta_r)

        # Update inner loop (theta-control)
        
        error_th = theta_r - theta

        # use dirty derivative
        self.theta_dot = self.beta * self.theta_dot \
            + (1 - self.beta) * ((theta - self.theta_d1) / P.Ts)
        
        # PD control on theta
        theta_r = self.kp_z * (z_r - z) - self.kd_z * self.z_dot

        theta_r = self.filter.update(theta_r)

        F_fl = P.m2*P.g/2.0 + P.m1*P.g*z / P.length # define equilibrium force
        F_tilde = self.kp_th * error_th - self.kd_th * self.theta_dot
        F = F_tilde + F_fl

        self.error_z_d1 = error_z
        self.z_d1 = z
        self.theta_d1 = theta

        # return saturated force
        return saturate(F, P.Fmax)

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
