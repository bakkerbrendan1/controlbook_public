# mass-spring-damper Parameter File
import numpy as np

# Physical parameters of the arm known to the controller
m = 5.0  # mass kg
k = 3.0  # spring constant Kg/s^2
b = 0.5  # damping coefficient Kg/s

# matrices for ctrl state feedback (chapter 11, 12)
A = np.array([[0.0, 1.0],
                      [-0.6, -0.1]])
B = np.array([[0.0],
                [0.2]])
C = np.array([[1.0, 0.0]])

# parameters for animation
length = 5.0
width = 1.0

# Initial Conditions
z0 = 0.0  # initial position of mass, m
zdot0 = 0.0  # initial velocity of mass m/s

# Simulation Parameters
t_start = 0.0 # Start time of simulation
t_end = 30  # End time of simulation
Ts = 0.01  # sample time for simulation
t_plot = 0.1  # the plotting and animation is updated at this rate

# dirty derivative parameters
sigma = 0.05 # cutoff freq for dirty derivative
# beta =   # dirty derivative gain

# saturation limits
F_max = 30  # Max force, N
