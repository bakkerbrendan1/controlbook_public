# VTOL Parameter File
import VTOLParam as P
from ctrlPIDhwF10 import ctrlPID
import hw15 as P15
from control import tf, bode
import matplotlib.pyplot as plt
import numpy as np

P10 = ctrlPID()

# flag to define if using dB or absolute scale for M(omega)
dB_flag = P15.dB_flag

# Assign transfer functions from previous HW solution
P_lon = P15.P_lon
P_lat_in = P15.P_lat_in
P_lat_out = P15.P_lat_out

# Compute the controller transfer functions from HW10

# PID xfer function for longitudinal control
C_lon = tf([(P10.kd_h+P10.kp_h*P10.sigma), (P10.kp_h+P10.ki_h*P10.sigma), P10.ki_h],
           [P10.sigma, 1, 0])

# PD control xfer function for inner loop of lateral control
C_lat_in = tf([(P10.kd_th+P10.sigma*P10.kp_th), P10.kp_th], [P10.sigma, 1])

# based on performance from HW 10, we get best performance with
# no integrator term, therefore, this is just a PD controller.
# PD xfer function for outer loop of lateral control
C_lat_out = tf([(P10.kd_z+P10.sigma*P10.kp_z), P10.kp_z], [P10.sigma, 1])

if __name__=="__main__":

    # display bode plots of transfer functions
    fig1 = plt.figure()
    bode([P_lon, P_lon*C_lon], dB=dB_flag)
    plt.legend(['$P_{lon}(s)$', '$C_{lon}(s)P_{lon}(s)$'])
    fig1.axes[0].set_title('VTOL Longitudinal Loop')

    # display bode plots of transfer functions
    fig2 = plt.figure()
    bode([P_lat_in, P_lat_in*C_lat_in], dB=dB_flag)
    plt.legend(['$P_{lat,in}(s)$', '$C_{lat,in}(s)P_{lat,in}(s)$'])
    fig2.axes[0].set_title('VTOL Lateral Inner Loop')

    fig3 = plt.figure()
    bode([P_lat_out, P_lat_out*C_lat_out, tf([1.0], [1.0, 0.0])], 
         dB=dB_flag)
    plt.legend(['$P_{lat,out}(s)$', '$C_{lat,out}(s)P_{lat,out}(s)$', '$\\frac{1}{s}$'])
    fig3.axes[0].set_title('VTOL Lateral Outer Loop')

    ##########################
    # for parts a)-b)
    ##########################

    # For part a), the transfer function C*P is shown in figure, this is a
    # shortcut to rendering latex from python and could be done in other ways.
    # We can use this result to find e_ss to any input.
    plt.figure()
    xfer_func = (P_lon*C_lon)._repr_latex_()
    plt.text(0.1, 0.5,'$%s$'%xfer_func[2:-2], fontsize='xx-large')
    plt.tick_params(axis='both', which='both', bottom=False, top=False,
                    left=False, right=False, labelbottom=False, labelleft=False)
    plt.title("xfer function for $C_{lon}(s)P_{lon}(s)$")

    omegas = [30.0]  # omega_no
    mag_CP_lon, phase, omegas = bode(P_lon*C_lon, plot=False, omega = omegas)

    print("\n\n\nvalue for metric calculation (part b):")
    # part b)
    if dB_flag == False:
        print("gamma_n = ", mag_CP_lon[0])
    elif dB_flag == True:
        # this conversion from absolute to dB (given the output from the
        # bode function), is a little funny, but is provided to match
        # the equations from the book.
        mag_CP_lon_dB = 20.0*np.log10(mag_CP_lon)
        print("gamma_n = ", 10.0 ** (mag_CP_lon_dB[0] / 20.0))


    ##########################
    # for parts c)-d)
    ##########################
    omegas = [2.0]  # omega_d_in
    mag_CP_lat_in, phase, omegas = bode(P_lat_in*C_lat_in, plot=False, omega = omegas)
    mag_P_lat_in, phase, omegas = bode(P_lat_in, plot=False, omega = omegas)


    print("\n\n\nvalue for metric calculation (part c):")
    # part c)
    if dB_flag == False:
        print("gamma_d_in = ", mag_P_lat_in[0]/mag_CP_lat_in[0])
    elif dB_flag == True:
        # this conversion from absolute to dB (given the output from the
        # bode function), is a little funny, but is provided to match
        # the equations from the book.
        mag_P_lat_in_dB = 20.0*np.log10(mag_P_lat_in)
        mag_CP_lat_in_dB = 20.0*np.log10(mag_CP_lat_in)
        print("gamma_d_in = ", 10.0 ** ((mag_P_lat_in_dB[0] - mag_CP_lat_in_dB[0])/20.0))

    # part d)
    print('for part d), look at plot directly for C_{lat,in}P_{lat,in}')

    ##########################
    # for parts e)-f)
    ##########################
    omegas = [0.01, 0.1]  # omega_d_out, omega_r
    mag_CP_lat_out, phase, omegas = bode(P_lat_out*C_lat_out, plot=False, omega = omegas)

    print("\n\n\nvalue for metric calculation (parts e and f):")

    # part e)
    if dB_flag == False:
        # part e)
        print("gamma_r = ", 1.0 / mag_CP_lat_out[1] )
        # part f)
        print("gamma_d_out =", 1.0 / mag_CP_lat_out[0])
    elif dB_flag == True:
        # this conversion from absolute to dB (given the output from the
        # bode function), is a little funny, but is provided to match
        # the equations from the book.
        mag_CP_lat_out_dB = 20.0*np.log10(mag_CP_lat_out)

        # part e)
        print("gamma_r = ", 1.0 / 10 ** (mag_CP_lat_out_dB[1] / 20.0))
        # part f)
        print("gamma_d_out =", 1.0 / 10 ** (mag_CP_lat_out_dB[0] / 20.0))

    print('\n\nClose window(s) to end program')
    plt.show()

