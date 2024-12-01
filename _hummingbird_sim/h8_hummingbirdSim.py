import matplotlib.pyplot as plt
import numpy as np
import hummingbirdParam as P
from signalGenerator import SignalGenerator
from hummingbirdAnimation import HummingbirdAnimation
from dataPlotter import DataPlotter
from hummingbirdDynamics import HummingbirdDynamics
from ctrlPD import ctrlPD

# instantiate pendulum, controller, and reference classes
hummingbird = HummingbirdDynamics(alpha=0.0)
controller = ctrlPD()
psi_ref = SignalGenerator(amplitude=30.*np.pi/180., frequency=0.05) # yaw
theta_ref = SignalGenerator(amplitude=15.*np.pi/180., frequency=0.05) # pitch
# doesn't make sense to have a phi ref because it's the roll of the bird
phi_ref = SignalGenerator(amplitude=0) # 15.*np.pi/180., frequency=0.05) 

# instantiate the simulation plots and animation
dataPlot = DataPlotter()
animation = HummingbirdAnimation()

t = P.t_start  # time starts at t_start
y = hummingbird.h()
while t < P.t_end:  # main simulation loop

    # Propagate dynamics at rate Ts
    t_next_plot = t + P.t_plot
    while t < t_next_plot:
        r = np.array([[phi_ref.square(t)],
                      [theta_ref.square(t)], 
                      [psi_ref.square(t)]])
        u = controller.update(r, hummingbird.state)
        y = hummingbird.update(u)  # Propagate the dynamics
        t = t + P.Ts  # advance time by Ts

    # update animation and data plots at rate t_plot
    animation.update(t, hummingbird.state)
    dataPlot.update(t, hummingbird.state, r, u[0], u[1])

    # the pause causes figure to be displayed during simulation
    plt.pause(0.0001)

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()


#######################################################
#        Values used to calibrate humming bird
#######################################################

# Humming bird 1

# Psi Ref: 20
# Psi kp: 0.34
# Psi kd: 0.14
# Psi ki: 0.20


# Theta: all zero

# Phi ref: 0
# Phi kp: 0.0275
# Phi kd: 0.0050

# Km: 0.355
