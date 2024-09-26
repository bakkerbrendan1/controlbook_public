import matplotlib.pyplot as plt
import massParam as P
from signalGenerator import signalGenerator
from massAnimation import massAnimation
from dataPlotter import dataPlotter
from massDynamics import massDynamics

import matplotlib
matplotlib.use('tkagg')

# instantiate arm, controller, and reference classes
arm = massDynamics()
#second_arm = armDynamics()
reference = signalGenerator(amplitude=0.01, frequency=0.02)
force = signalGenerator(amplitude=0.2, frequency=0.05)
#torque2 = signalGenerator(amplitude=0.1, frequency=0.01)

# instantiate the simulation plots and animation
dataPlot = dataPlotter()
#dataPlot2 = dataPlotter()
animation = massAnimation()
#animation2 = armAnimation()

t = P.t_start  # time starts at t_start
while t < P.t_end:  # main simulation loop
    # Propagate dynamics in between plot samples
    t_next_plot = t + P.t_plot
    # updates control and dynamics at faster simulation rate
    while t < t_next_plot:
        # Get referenced inputs from signal generators
        r = reference.square(t)
        u = force.square(t)
        #u2 = torque2.sin(t)
        y = arm.update(u)  # Propagate the dynamics
        #y2 = second_arm.update(u2)
        t = t + P.Ts  # advance time by Ts
    # update animation and data plots
    animation.update(arm.state)
    #animation2.update(second_arm.state)
    dataPlot.update(t, r, arm.state, u)
    #dataPlot2.update(t, r, second_arm.state, u2)

    # the pause causes the figure to be displayed during the
    # simulation
    plt.pause(0.0001)  

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
