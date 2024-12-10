import massParam as P
from control import tf, bode
import matplotlib.pyplot as plt

# flag to define if using dB or absolute scale for M(omega)
dB_flag = False

# Compute plant transfer functions
Plant = tf([1.0/P.m],
           [1, P.b/P.m, P.k/P.m])


if __name__=="__main__":

    # Bode plot of the plant
    fig = plt.figure()
    bode(Plant, dB=dB_flag, omega_limits=[10**(-2), 10**(1)])
    fig.axes[0].set_title('HW 15 - P(s) for mass-spring-damper')

    print('Close window to end program')
    plt.show()
