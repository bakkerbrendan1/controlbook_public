# VTOL Parameter File
import numpy as np

# Physical parameters of the  VTOL known to the controller
mc = 1.0 # kg
mr = 0.25  # kg
Jc = 0.0042  # kg m^2
d = 0.3  # m
mu = 0.1  # kg/s
g = 9.81  # m/s^2
F_wind = 0.0 # record: 45 on hw12sim # wind disturbance force is zero in initial homeworks

# parameters for animation
length = 10.0

Fe = (mc + 2.0 * mr) * g

A_lat = np.array([[0.0, 0.0, 1.0, 0.0],
                  [0.0, 0.0, 0.0, 1.0],
                  [0.0, -Fe / (mc + 2.0 * mr), -(mu / (mc + 2.0 * mr)), 0.0],
                  [0.0, 0.0, 0.0, 0.0]])
B_lat = np.array([[0.0],
                  [0.0],
                  [0.0],
                  [1.0 / (Jc + 2.0*mr*d**2)]])
C_lat = np.array([[1.0, 0.0, 0.0, 0.0],
                  [0.0, 1.0, 0.0, 0.0]])

A_lon = np.array([[0.0, 1.0],
                  [0.0, 0.0]])
B_lon = np.array([[0.0],
                  [1.0 / (mc + 2*mr)]])
C_lon = np.array([[1.0, 0.0]])

A = np.array([[0, 0, 0, 1, 0, 0],
              [0, 0, 0, 0, 1, 0],
              [0, 0, 0, 0, 0, 1],
              [0, 0, -Fe / (mc + 2 * mr), -mu / (mc + 2 * mr), 0, 0],
              [0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0]])

B = np.array([[0, 0],
              [0, 0],
              [0, 0],
              [0, 0],
              [1 / (mc + 2 * mr), 0],
              [0, 1 / (Jc + 2 * d**2 * mr)]])

C = np.array([[1,0,0,0,0,0],
              [0,1,0,0,0,0],
              [0,0,1,0,0,0]])

# Initial Conditions
z0 = 0.0  # initial lateral position
h0 = 0.0  # initial altitude
theta0 = 0.0 # initial roll angle
zdot0 = 0.0  # initial lateral velocity
hdot0 = 0.0  # initial climb rate
thetadot0 = 0.0  # initial roll rate
target0 = 0.0

# Simulation Parameters
t_start = 0.0 # Start time of simulation
t_end = 30.0  # End time of simulation
Ts = 0.01  # sample time for simulation
t_plot = 0.1 # the plotting and animation is updated at this rate

# saturation limits
fmax = 500 # 20  # Max Force, N
tau_max = 2*fmax*d

# dirty derivative parameters
# sigma =   # cutoff freq for dirty derivative
# beta =  # dirty derivative gain

# equilibrium force
# Fe =

# mixing matrix
unmixing = np.array([[1.0, 1.0], [d, -d]]) # converts fl and fr (LR) to force and torque (FT)
mixing = np.linalg.inv(unmixing) # converts force and torque (FT) to fl and fr (LR) 

