import matplotlib.pyplot as plt
import numpy as np
import VTOLParam as P
from signalGenerator import signalGenerator
from VTOLAnimation import VTOLAnimation
from dataPlotter import dataPlotter
from VTOLDynamics import VTOLDynamics
# from ctrlObserver_hwF13 import ctrlObserver
from ctrlObserverKey import ctrlObserver
from dataPlotterObserver import dataPlotterObserver
import matplotlib
matplotlib.use('tkagg')

# instantiate satellite, controller, and reference classes
VTOL = VTOLDynamics(alpha=0.0)
controller = ctrlObserver()
z_reference = signalGenerator(amplitude=2.5, frequency=0.04, y_offset=3.0)
h_reference = signalGenerator(amplitude=3.0, frequency=0.03, y_offset=5.0)
z_disturbance = signalGenerator(amplitude=0.1)
h_disturbance = signalGenerator(amplitude=0.1)

# instantiate the simulation plots and animation
dataPlot = dataPlotter()
dataPlotObserver = dataPlotterObserver()
animation = VTOLAnimation()

t = P.t_start  # time starts at t_start
y = VTOL.h()  # output of system at start of simulation

while t < P.t_end:  # main simulation loop
    # Propagate dynamics in between plot samples
    t_next_plot = t + P.t_plot

    # updates control and dynamics at faster simulation rate
    while t < t_next_plot:  
        h_ref = h_reference.square(t)  # reference input
        z_ref = z_reference.square(t)
        h_dis = 0.0 # h_disturbance.step(t)  # disturbance input
        z_dis = 0.0 # z_disturbance.step(t)
        r = np.array([[z_ref], [h_ref]])
        d = np.array([[z_dis], [h_dis]])
        n = 0.0 # noise.random(t)  # simulate sensor noise
        x = VTOL.state

        # controller here takes 2x1 nparray of zref and h
        u, xhat_lat, xhat_lon = controller.update(r, y)  # update controller
        y = VTOL.update(u) # + d)  # propagate system
        t = t + P.Ts  # advance time by Ts

    # update animation and data plots
    animation.update(VTOL.state, z_ref)
    dataPlot.update(t, VTOL.state, z_ref, h_ref, u[0], u[1])
    dataPlotObserver.update(t, VTOL.state, xhat_lat, xhat_lon)

    # the pause causes the figure to display for the simulation.
    # plt.pause(0.00001) # 0.0001)  

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
