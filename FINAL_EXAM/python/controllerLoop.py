import massParam as P
import loopshape as L
import numpy as np

class controllerLoop:
    def __init__(self):
        # Either define your state space models for the compensator (C) and
        # prefilter (F) here, or define the variables you need for a digital filter.


        self.limit = P.F_max
        self.F_e = P.F_e
        self.Ts = P.Ts

    def update(self, z_r, y):
        z = y.item(0)

        #prefilter the reference command
        z_r_filtered =

        # error signal
        error = 

        # find the control
        F_tilde =

        # compute total force
        F =

        return F
