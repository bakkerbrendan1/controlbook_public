# VTOL Parameter File
import VTOLParam as P
from control import tf, bode
import matplotlib.pyplot as plt

# flag to define if using dB or absolute scale for M(omega)
dB_flag = False 

# Compute open-loop transfer functions
P_lon = tf(1/(P.mc+2*P.mr), [1, 0, 0])
P_lat_in = tf(1/(P.Jc+2*P.mr*P.d**2), [1, 0, 0])
P_lat_out = tf(-P.Fe/(P.mc+2*P.mr), [1, P.mu/(P.mc+2*P.mr), 0])

if __name__=="__main__":

    # Plot the closed loop and open loop bode plots for the inner loop
    fig1 = plt.figure()
    bode(P_lon, dB=dB_flag)
    fig1.axes[0].set_title('$P(s)$ for longitudinal dynamics')

    fig2 = plt.figure()
    bode(P_lat_in, dB=dB_flag)
    fig2.axes[0].set_title('$P_{in}(s)$ for lateral dynamics')

    fig3 = plt.figure()
    bode(P_lat_out, dB=dB_flag)
    fig3.axes[0].set_title('$P_{out}(s)$ for lateral dynamics')

    print('Close window to end program')
    plt.show()
