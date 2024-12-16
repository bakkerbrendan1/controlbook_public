import matplotlib.pyplot as plt
import numpy as np
import massParam as P
from signalGenerator import signalGenerator
from massAnimation import massAnimation
from dataPlotter import dataPlotter
from massDynamics import massDynamics
from controllerPID import controllerPID

# instantiate system, controller, and reference classes
mass = massDynamics()
controller = controllerPID()
reference = signalGenerator(amplitude=0.5, frequency=0.05)
disturbance = signalGenerator(amplitude=0.1)
noise = signalGenerator(amplitude=0.01)

# instantiate the simulation plots and animation
dataPlot = dataPlotter()
animation = massAnimation()

t = P.t_start
y = mass.h()
while t < P.t_end:
    t_next_plot = t + P.t_plot
    while t < t_next_plot:
        r = reference.square(t)
        d = disturbance.step(t)
        n = 0.0  #noise.random(t)
        u = controller.update(r, y + n)
        y = mass.update(u + d)
        t = t + P.Ts
    # update animation and data plots
    animation.update(mass.state)
    dataPlot.update(t, r, mass.state, u)
    # plt.pause(0.0001)

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
