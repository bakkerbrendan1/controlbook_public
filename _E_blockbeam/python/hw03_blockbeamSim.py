import matplotlib.pyplot as plt
import blockbeamParam as P
from signalGenerator import signalGenerator
from blockbeamAnimation import blockbeamAnimation
from dataPlotter import dataPlotter
from blockbeamDynamics import blockbeamDynamics
import matplotlib
matplotlib.use('tkagg')

# instantiate pendulum, controller, and reference classes
blockbeam = blockbeamDynamics(alpha=0.0)
reference = signalGenerator(amplitude=0.0005, frequency=0.02)
force = signalGenerator(amplitude=0.00001, frequency=0.1)

# instantiate the simulation plots and animation
dataPlot = dataPlotter()
animation = blockbeamAnimation()

t = P.t_start  # time starts at t_start
while t < P.t_end:  # main simulation loop
    # Propagate dynamics at rate Ts
    t_next_plot = t + P.t_plot
    while t < t_next_plot:
        r = reference.square(t)
        u = force.sin(t) + 11.5
        y = blockbeam.update(u)  # Propagate the dynamics
        t = t + P.Ts  # advance time by Ts
    # update animation and data plots at rate t_plot
    animation.update(blockbeam.state)
    # def update(self, t, reference, states, ctrl):
    dataPlot.update(t, r, blockbeam.state, u)
    plt.pause(0.0001)  # allows time for animation to draw

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
