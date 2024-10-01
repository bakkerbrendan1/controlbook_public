import matplotlib.pyplot as plt
import hummingbirdParam as P
from signalGenerator import SignalGenerator
from hummingbirdAnimation import HummingbirdAnimation
from dataPlotter import DataPlotter
from hummingbirdDynamics import HummingbirdDynamics
import matplotlib
import numpy as np
matplotlib.use('tkagg')

# instantiate satellite, controller, and reference classes
bird = HummingbirdDynamics()
reference = SignalGenerator(amplitude=0.5, frequency=0.1)
fl = SignalGenerator(amplitude=4, frequency=0.1)
fr = SignalGenerator(amplitude=4, frequency=0.1)

# instantiate the simulation plots and animation
dataPlot = DataPlotter()
animation = HummingbirdAnimation()

t = P.t_start  # time starts at t_start
while t < P.t_end:  # main simulation loop
    # Propagate dynamics in between plot samples
    t_next_plot = t + P.t_plot
    # updates control and dynamics at faster simulation rate
    while t < t_next_plot:  
        r = np.array([[reference.square(t), 0],
                      [reference.square(t), 0]])
        u = np.array([[fl.sin(t) + 5,0],
                      [fr.sin(t) + 5,0]])
        y = bird.update(u)  # Propagate the dynamics
        t = t + P.Ts  # advance time by Ts
    # update animation and data plots
    animation.update(t, bird.state)
    # def update(self, t, state, ref, force, torque):
    dataPlot.update(t, bird.state, r, u, [[0,0,0],[0,0,0]])
    plt.pause(0.0001)  

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
