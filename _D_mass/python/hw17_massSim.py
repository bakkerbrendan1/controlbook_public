# Single link arm Parameter File
from control import bode, margin, mag2db, tf
import matplotlib.pyplot as plt
from ctrlPIDhwD10 import ctrlPID

# flag to define if using dB or absolute scale for M(omega)
dB_flag = False

# assigning plant and controller from past HW (to make sure
# we don't introduce additional errors)

P10 = ctrlPID()

Plant = tf([0.2],
           [1, 0.1, 0.6])
C_pid = tf([(P10.kd+P10.kp*P10.sigma), 
            (P10.kp+P10.ki*P10.sigma), 
            P10.ki],
           [P10.sigma, 1, 0])
           

if __name__ == '__main__':
    # Calculate the phase and gain margins
    if dB_flag:
        gm, pm, Wcg, Wcp = margin(Plant*C_pid)
        gm = mag2db(gm)
        print("Inner Loop:", "gm: ",gm,
              " pm: ", pm," Wcg: ", Wcg, " Wcp: ", Wcp)

    else:
        gm, pm, Wcg, Wcp = margin(Plant*C_pid)
        print("Inner Loop:", "gm: ", gm,
              " pm: ", pm, " Wcg: ", Wcg, " Wcp: ", Wcp)

    # display bode plots of transfer functions
    fig1 = plt.figure()

    # this makes two bode plots for open and closed loop
    bode(Plant * C_pid, dB=dB_flag,
         label='C(s)P(s) - Open-loop')
    bode(Plant*C_pid/(1+Plant*C_pid), dB=dB_flag,
         label=r'$\frac{P(s)C(s)}{1+P(s)C(s)}$ - Closed-loop')

    # now we can add lines to show where we calculated the GM and PM
    gm_line = fig1.axes[0].plot([Wcg, Wcg],
                                plt.ylim(), 'k--', label='GM')
    gm_line[0].set_label('GM')
    fig1.axes[0].legend()
    fig1.axes[1].plot([Wcp, Wcp], plt.ylim(), 'b--', label='PM')
    plt.legend()

    # setting axis title
    fig1.axes[0].set_title('Single Link Robot Arm - '+
                           'GM:'+str(round(gm, 2))+
                           ', PM:'+str(round(pm, 2)))

    plt.show()
