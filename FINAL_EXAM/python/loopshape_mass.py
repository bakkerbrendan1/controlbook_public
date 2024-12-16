import massParam as P
import matplotlib.pyplot as plt
from control import tf, tf2ss
import control as cnt
import numpy as np

# in this file, you should will need to:
#### 1)define xfer functions for the plant, and the PID controller
#### 2)define specifications (for C_pid * P, and for the final controller C)
#### 3)make bode plots to help you in meeting specifications and responding to questions

# Compute plant transfer functions
Plant = tf()
C_pid = tf()

# compute final compensator C, and prefilter F
C =
F =

# convert them to state space models if you are going to use that method
# otherwise you can just use C and F directly in the controllerLoop.py file.
Css=cnt.tf2ss(C)
Fss=cnt.tf2ss(F)



