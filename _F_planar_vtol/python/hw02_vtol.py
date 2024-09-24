import matplotlib.pyplot as plt
import numpy as np
import VTOLParam as P
from signalGenerator import signalGenerator
from VTOLAnimation import VTOLAnimation
from dataPlotter import dataPlotter

# instantiate reference input classes
reference = signalGenerator(amplitude=0.5, frequency=0.1)
vRef = signalGenerator(amplitude=2.0, frequency=0.1)
zRef = signalGenerator(amplitude=0.25, frequency=0.5)
hRef = signalGenerator(amplitude=10, frequency=.5)
fRef = signalGenerator(amplitude=5, frequency=0.5)
tauRef = signalGenerator(amplitude=50, frequency=1)

# instantiate the simulation plots and animation
dataPlot = dataPlotter()
animation = VTOLAnimation()

t = P.t_start  # time starts at t_start
while t < P.t_end:  # main simulation loop
    # set variables
    r = reference.sin(t)
    v = vRef.sin(t) + 2
    z = zRef.sin(t) + 2
    h = hRef.sin(t) - 1
    f = fRef.sin(t) + 0.5
    tau = tauRef.sin(t) + 0.5

    # update animation
    state = np.array([[v], [z], [0.0], [0.0]])
    animation.update(state)
    dataPlot.update(t, state, z, h, f, tau)#, state, tau)
    # advance time by t_plot
    t = t + P.t_plot  
    plt.pause(0.001)

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
