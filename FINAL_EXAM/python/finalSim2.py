import matplotlib.pyplot as plt
import massParam as P
from massAnimation import massAnimation
from dataPlotter import dataPlotter
from massDynamics import massDynamics

# instantiate arm, controller, and reference classes
mass = massDynamics()

# instantiate the simulation plots and animation
dataPlot = dataPlotter()
animation = massAnimation()

t = P.t_start
while t < P.t_end:
    t_next_plot = t + P.t_plot
    while t < t_next_plot:
        u = P.F_e
        y = mass.update(u)
        t = t + P.Ts
    # update animation and data plots
    animation.update(mass.state)
    dataPlot.update(t, 0, mass.state, u)
    # plt.pause(0.0001)

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
