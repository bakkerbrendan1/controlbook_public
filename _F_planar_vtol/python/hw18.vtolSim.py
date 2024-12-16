import matplotlib.pyplot as plt
import numpy as np
import VTOLParam as P
from signalGenerator import signalGenerator
from VTOLAnimation import VTOLAnimation
from dataPlotter import dataPlotter
from VTOLDynamics import VTOLDynamics
from ctrlLoopshape import ctrlLoopshape
import matplotlib
matplotlib.use('tkagg')

# instantiate arm, controller, and reference classes
VTOL = VTOLDynamics(alpha=0.1)
controller = ctrlLoopshape()
z_reference = signalGenerator(amplitude=2.5, frequency=0.04, y_offset=3.0)
h_reference = signalGenerator(amplitude=3.0, frequency=0.03, y_offset=5.0)
disturbance = signalGenerator(amplitude=0.1)
noise = signalGenerator(amplitude=0.01)

# instantiate the simulation plots and animation
dataPlot = dataPlotter()
animation = VTOLAnimation()

t = P.t_start  # time starts at t_start
y = VTOL.h()  # output of system at start of simulation

while t < P.t_end:  # main simulation loop
    # Get referenced inputs from signal generators
    # Propagate dynamics in between plot samples
    t_next_plot = t + P.t_plot

    while t < t_next_plot:
        h_ref = h_reference.square(t)  # reference input
        z_ref = z_reference.square(t)
        r = np.array([[z_ref], [h_ref]])
        d_z = disturbance.square(t) # input disturbance
        d = np.array([[d_z], [d_z]])
        n = noise.random(t)  # simulate sensor noise
        x = VTOL.state
        u = controller.update(r, x + n)  # update controller
        y = VTOL.update(u + d)  # propagate system
        t = t + P.Ts  # advance time by Ts

    # update animation and data plots
    animation.update(VTOL.state)
    # print('u = ', u[0][0])
    dataPlot.update(t, VTOL.state, z_ref, h_ref, u[0], u[1])
    # plt.pause(0.0001)

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
