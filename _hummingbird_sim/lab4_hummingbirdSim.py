import matplotlib.pyplot as plt
import hummingbirdParam as P
from signalGenerator import SignalGenerator
from hummingbirdAnimation import HummingbirdAnimation
from dataPlotter import DataPlotter
from hummingbirdDynamics import HummingbirdDynamics
from ctrlEquilibrium import ctrlEquilibrium
import matplotlib
import numpy as np
matplotlib.use('tkagg')

# instantiate satellite, controller, and reference classes
bird = HummingbirdDynamics()
reference = SignalGenerator(amplitude=0.5, frequency=0.1)
force = SignalGenerator(amplitude=4, frequency=0.1)
torque = SignalGenerator(amplitude=4, frequency=0.1)
ctrl = ctrlEquilibrium()

# instantiate the simulation plots and animation
dataPlot = DataPlotter()
animation = HummingbirdAnimation()

t = P.t_start  # time starts at t_start
while t < P.t_end:  # main simulation loop
    # Propagate dynamics in between plot samples
    t_next_plot = t + P.t_plot
    # updates control and dynamics at faster simulation rate
    while t < t_next_plot:  
        r = np.array([[reference.square(t)], # phi
                      [reference.square(t)], # theta
                      [reference.square(t)]]) # psi
        
        pwm = ctrl.update(bird.state)
        bird.update(pwm)
        t = t + P.Ts  # advance time by Ts
    # update animation and data plots
    animation.update(t, bird.state)
    # def update(self, t, state, ref, force, torque):
    force_tor = P.unmixing @ pwm
    dataPlot.update(t, bird.state, r, force_tor[0], force_tor[1])
    plt.pause(0.0001)  

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()

# km values:
# hardware: 0.354 (hummingbird 2)
# software: 0.23468322067605632