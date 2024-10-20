import matplotlib.pyplot as plt
import numpy as np
import blockbeamParam as P
from signalGenerator import signalGenerator
from blockbeamAnimation import blockbeamAnimation
from dataPlotter import dataPlotter
from blockbeamDynamics import blockbeamDynamics
from ctrlPDhw08 import ctrlPD
import matplotlib
matplotlib.use('tkagg')

# instantiate pendulum, controller, and reference classes
blockbeam = blockbeamDynamics(alpha=0.1)
controller = ctrlPD()

# HW asked for reference input frequency of 0.01
reference = signalGenerator(amplitude=0.5, frequency=0.05)
disturbance = signalGenerator(amplitude=0.0, frequency=0.0)

# instantiate the simulation plots and animation
dataPlot = dataPlotter()
animation = blockbeamAnimation()

t = P.t_start  # time starts at t_start
y = blockbeam.h()  # output of system at start of simulation

# for part e), we can uncomment below
# blockbeam.state[1,0] = 10.0*np.pi/180.0
#reference = signalGenerator(amplitude = 0.0, frequency=0.0)


while t < P.t_end:  # main simulation loop
    # Propagate dynamics in between plot samples
    t_next_plot = t + P.t_plot

    while t < t_next_plot:
        r = reference.square(t)  # reference input
        d = disturbance.step(t)  # input disturbance
        n = 0.0  #noise.random(t)  # simulate sensor noise
        x = blockbeam.state  # use state instead of output
        u = controller.update(r, x)  # update controller
        y = blockbeam.update(u + d)  # propagate system
        t = t + P.Ts  # advance time by Ts

    # update animation and data plots
    animation.update(blockbeam.state)
    dataPlot.update(t, r, blockbeam.state, u)
    plt.pause(0.0001)

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
