import massParam as P
from ctrlPIDhwD10 import ctrlPID
import hw15 as P15
from control import tf, bode
import matplotlib.pyplot as plt
import numpy as np

P10 = ctrlPID()

# flag to define if using dB or absolute scale for M(omega)
dB_flag = P15.dB_flag

# Assign plant from previous homework solution
Plant = P15.Plant

# Compute transfer function of controller from HW 10
C_pid = tf([(P10.kd+P10.kp*P.sigma), (P10.kp+P10.ki*P.sigma), P10.ki],
           [P.sigma, 1, 0])


if __name__=="__main__":

    # display bode plots of transfer functions
    fig = plt.figure()
    bode([Plant, Plant*C_pid], dB=dB_flag)
    plt.legend(['P(s)', 'C(s)P(s)'])
    fig.axes[0].set_title('HW 16 - Freq. Response Metrics for Mass Spring Damper')

    # the order here matters and is assumed to be from smallest to largest because
    # the function re-orders the omega values to be smallest to largest regardless
    # of the order in which they are passed in. The function will also plot things
    # dB, but returns all magnitudes in absolute.
    omegas = [0.1, 100.0]  # omega_d, omega_no
    mag_CP, phase, omegas = bode(Plant*C_pid, plot=False, omega = omegas)
    mag_P, phase, omegas = bode(Plant, plot=False, omega = omegas)

    # For part a), the transfer function C*P is shown in figure, this is a
    # shortcut to rendering latex from python and could be done in other ways.
    # We can use this result to find e_ss to any input.
    plt.figure()
    xfer_func = (Plant*C_pid)._repr_latex_()
    plt.text(0.1, 0.5,'$%s$'%xfer_func[2:-2], fontsize='xx-large')
    plt.tick_params(axis='both', which='both', bottom=False, top=False,
                    left=False, right=False, labelbottom=False, labelleft=False)

    # finding the necessary values for parts b) and c)
    print("\n\nvalues for metric calculation for parts b) and c):")

    if dB_flag == False:
        print("gamma_d = ", mag_P[0]/mag_CP[0])
        print("gamma_n = ", mag_CP[1])

    elif dB_flag == True:
        # this conversion from absolute to dB (given the output from the
        # bode function), is a little funny, but is provided to match
        # the equations from the book.
        mag_P_dB = 20.0*np.log10(mag_P)
        mag_CP_dB = 20.0*np.log10(mag_CP)
        print("gamma_d = ", 10.0**((mag_P_dB[0]-mag_CP_dB[0])/20.0))
        print("gamma_n = ", 10.0**(mag_CP_dB[1]/20.0))

    print('Close window to end program')
    plt.show()
