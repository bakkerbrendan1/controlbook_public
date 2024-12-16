import numpy as np

# one of these two forms for importing control library functions may be helpful.
#import control as cnt
#from control import tf, bode, etc.

# physical parameters
g = 9.8
theta = 45.0*np.pi/180.0   #this variable just defines the slope of the block, and may not be used explicitly
m = 0.5
k1 = 0.05
k2 = 0.02
F_max = 5.0
b = 0.1

# simulation parameters
t_start = 0.0
t_end = 100.0
Ts = 0.01
t_plot = 0.1
sigma = 0.05

