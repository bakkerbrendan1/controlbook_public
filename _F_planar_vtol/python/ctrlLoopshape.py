import numpy as np
from control import c2d, tf
import VTOLParam as P
import loopShaping_lat_in as L_lat_in
import loopShaping_lat_out as L_lat_out
import loopShaping_lon as L_lon

class ctrlLoopshape:
    def __init__(self, method="state_space"):
        if method == "state_space":
            self.prefilter_lat_out = transferFunction(L_lat_out.F_num, L_lat_out.F_den, P.Ts)
            self.prefilter_lon = transferFunction(L_lon.F_num, L_lon.F_den, P.Ts)

            self.control_lat_in = transferFunction(L_lat_in.C_num, L_lat_in.C_den, P.Ts)
            self.control_lat_out = transferFunction(L_lat_out.C_num, L_lat_out.C_den, P.Ts)
            self.control_lon = transferFunction(L_lon.C_num, L_lon.C_den, P.Ts)
        elif method == "digital_filter":
            self.prefilter_lat_out = digitalFilter(L_lat_out.F_num, L_lat_out.F_den, P.Ts)
            self.prefilter_lon = digitalFilter(L_lon.F_num, L_lon.F_den, P.Ts)
            
            self.control_lat_in = digitalFilter(L_lat_in.C_num, L_lat_in.C_den, P.Ts)
            self.control_lat_out = digitalFilter(L_lat_out.C_num, L_lat_out.C_den, P.Ts)
            self.control_lon = digitalFilter(L_lon.C_num, L_lon.C_den, P.Ts)
        self.method = method

    def update(self, reference, y):
        z_r = reference[0][0]
        h_r = reference[1][0]
        z = y[0][0]
        h = y[1][0]
        theta = y[2][0]
        # self.z_dot = state[3][0]
        self.h_dot = y[4][0]
        # self.theta_dot = state[5][0]

        ######################################
        #      Altitude PID control: h
        ######################################
        # calculate control for the altitude loop
        h_r_filtered = self.prefilter_lon.update(h_r)
        error_h = h_r_filtered - h
        F_tilde = self.control_lon.update(error_h)
        F = saturate( F_tilde + P.Fe, 2*P.fmax)

        ######################################
        #      Outer loop PID control: z
        ######################################
        # prefiltered for outer loop
        z_r_filtered = self.prefilter_lat_out.update(z_r)
        error_lat_out = z_r_filtered - z

        # outer loop control
        theta_r = self.control_lat_out.update(error_lat_out)

        ######################################
        #      Inner loop PID control: theta
        ######################################
        error_lat_in = theta_r - theta
        tau_unsat = self.control_lat_in.update(error_lat_in)
        tau = saturate( tau_unsat, 2*P.fmax*P.d)

        motor_thrusts = P.mixing @ np.array([[F], [tau]])
        return motor_thrusts

def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u


class transferFunction:
    def __init__(self, num, den, Ts):
        # expects num and den to be numpy arrays of
        # shape (1,m+1) and (1,n+1)
        m = num.shape[1]
        n = den.shape[1]
        # set initial conditions
        self.state = np.zeros((n-1, 1))
        self.Ts = Ts
        # make the leading coef of den == 1
        if den.item(0) != 1:
            tmp = den.item(0)
            num = num / tmp
            den = den / tmp
        self.num = num
        self.den = den
        # set up state space equations in control canonic form
        self.A = np.zeros((n-1, n-1))
        self.B = np.zeros((n-1, 1))
        self.C = np.zeros((1, n-1))
        for i in range(0, n-1):
            self.A[0][i] = - den.item(i + 1)
        for i in range(1, n-1):
            self.A[i][i - 1] = 1.0
        if n>1:
            self.B[0][0] = 1.0
        if m == n:
            self.D = num.item(0)
            for i in range(0, n-1):
                self.C[0][i] = num.item(i+1) \
                               - num.item(0)*den.item(i+1)
        else:
            self.D = 0.0
            for i in range(n-m-1, n-1):
                self.C[0][i] = num.item(i)

    def update(self, u):
        x = self.rk4(u)
        y = self.C @ x + self.D * u
        return y.item(0)

    def f(self, state, u):
        xdot = self.A @ state + self.B * u
        return xdot

    def rk4(self, u):
        # Integrate ODE using Runge-Kutta 4 algorithm
        F1 = self.f(self.state, u)
        F2 = self.f(self.state + self.Ts / 2 * F1, u)
        F3 = self.f(self.state + self.Ts / 2 * F2, u)
        F4 = self.f(self.state + self.Ts * F3, u)
        self.state += self.Ts / 6 * (F1 + 2 * F2 + 2 * F3 + F4)
        return self.state
    

class digitalFilter:
    def __init__(self, num, den, Ts):
        self.Ts = Ts
        sys = tf(num[0], den[0])
        sys_d = c2d(sys, Ts, method='tustin')
        self.den_d = sys_d.den[0][0]
        self.num_d = sys_d.num[0][0]
        self.prev_filt_output = np.zeros(len(self.num_d)-1)
        self.prev_filt_input = np.zeros(len(self.den_d))

    def update(self, u):
        # update vector with filter inputs (u)
        self.prev_filt_input = np.hstack(([u], self.prev_filt_input[0:-1]))
        # use filter coefficients to calculate new output (y)
        y = self.num_d @ self.prev_filt_input - self.den_d[1:] @ self.prev_filt_output
        # update vector with filter outputs
        self.prev_filt_output = np.hstack(([y], self.prev_filt_output[0:-1]))
        return y