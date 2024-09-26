import matplotlib.pyplot as plt
import VTOLParam as P
from signalGenerator import signalGenerator
from VTOLAnimation import VTOLAnimation
from dataPlotter import dataPlotter
from VTOLDynamics import VTOLDynamics
import matplotlib
matplotlib.use('tkagg')

# instantiate satellite, controller, and reference classes
VTOL = VTOLDynamics()
reference = signalGenerator(amplitude=0.5, frequency=0.1)
fl = signalGenerator(amplitude=4, frequency=0.1)
fr = signalGenerator(amplitude=4, frequency=0.1)

# instantiate the simulation plots and animation
dataPlot = dataPlotter()
animation = VTOLAnimation()

t = P.t_start  # time starts at t_start
while t < P.t_end:  # main simulation loop
    # Propagate dynamics in between plot samples
    t_next_plot = t + P.t_plot
    # updates control and dynamics at faster simulation rate
    while t < t_next_plot:  
        r = reference.square(t)
        u = [[fl.sin(t) + 5,0],
             [fr.sin(t) + 5,0]]
        y = VTOL.update(u)  # Propagate the dynamics
        t = t + P.Ts  # advance time by Ts
    # update animation and data plots
    animation.update(VTOL.state)
    # def update(self, t, states, z_ref, h_ref, force, torque):
    dataPlot.update(t, VTOL.state, 0, 0, 0, 0)
    plt.pause(0.0001)  

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
