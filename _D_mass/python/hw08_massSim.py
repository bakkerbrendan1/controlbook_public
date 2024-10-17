import matplotlib.pyplot as plt
import massParam as P
import numpy as np
from signalGenerator import signalGenerator
from massAnimation import massAnimation
from dataPlotter import dataPlotter
from massDynamics import massDynamics
from ctrlPDhw08 import ctrlPD

import matplotlib
matplotlib.use('tkagg')

# instantiate arm, controller, and reference classes
mass = massDynamics(alpha=0)
controller = ctrlPD()
reference = signalGenerator(amplitude=1, frequency=0.05)
# force = signalGenerator(amplitude=2, frequency=0.05)
disturbance = signalGenerator(amplitude=0.01)
noise = signalGenerator(amplitude=0.01)

# instantiate the simulation plots and animation
dataPlot = dataPlotter()
animation = massAnimation()

t = P.t_start  # time starts at t_start
while t < P.t_end:  # main simulation loop
    # Propagate dynamics in between plot samples
    t_next_plot = t + P.t_plot
    # updates control and dynamics at faster simulation rate
    while t < t_next_plot:
        # Get referenced inputs from signal generators
        r = reference.square(t)
        d = disturbance.step(t)
        n = noise.random(t)
        x = mass.state
        u = controller.update(r, x)
        y = mass.update(u + d)
        t = t + P.Ts  # advance time by Ts

    # update animation and data plots
    animation.update(mass.state)
    #animation2.update(second_arm.state)
    dataPlot.update(t, r, mass.state, u)
    #dataPlot2.update(t, r, second_arm.state, u2)

    # the pause causes the figure to be displayed during the
    # simulation
    plt.pause(0.0001)  

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
