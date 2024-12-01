import matplotlib.pyplot as plt
import massParam as P
import numpy as np
from signalGenerator import signalGenerator
from massAnimation import massAnimation
from dataPlotter import dataPlotter
from massDynamics import massDynamics
from ctrlObserver_hwD13 import ctrlObserver
from dataPlotterObserver import dataPlotterObserver

import matplotlib
matplotlib.use('tkagg')

# instantiate arm, controller, and reference classes
mass = massDynamics(alpha=0.2)
controller = ctrlObserver()
reference = signalGenerator(amplitude=1, frequency=0.05)
disturbance = signalGenerator(amplitude=0.25)
noise = signalGenerator(amplitude=0.01)

# instantiate the simulation plots and animation
dataPlot = dataPlotter()
dataPlotObserver = dataPlotterObserver()
animation = massAnimation()

t = P.t_start  # time starts at t_start
y = mass.h()

while t < P.t_end:  # main simulation loop
    # Propagate dynamics in between plot samples
    t_next_plot = t + P.t_plot
    # updates control and dynamics at faster simulation rate
    while t < t_next_plot:
        # Get referenced inputs from signal generators
        r = reference.square(t)
        d = 0.0 # disturbance.step(t)
        n = noise.random(t)
        x = mass.state
        u, xhat = controller.update(r, y + n)
        y = mass.update(u + d)
        t = t + P.Ts  # advance time by Ts

    # update animation and data plots
    animation.update(mass.state)
    dataPlot.update(t, r, mass.state, u)
    dataPlotObserver.update(t, mass.state, xhat)

    # the pause causes the figure to be displayed during the
    # simulation
    plt.pause(0.0001)  

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
