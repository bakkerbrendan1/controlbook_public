import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np 
import massParam as P

# if you are having difficulty with the graphics, 
# try using one of the following backends 
# See https://matplotlib.org/stable/users/explain/backends.html
import matplotlib
# matplotlib.use('qtagg')  # requires pyqt or pyside
# matplotlib.use('ipympl')  # requires ipympl
# matplotlib.use('gtk3agg')  # requires pyGObject and pycairo
# matplotlib.use('gtk4agg')  # requires pyGObject and pycairo
# matplotlib.use('gtk3cairo')  # requires pyGObject and pycairo
# matplotlib.use('gtk4cairo')  # requires pyGObject and pycairo
matplotlib.use('tkagg')  # requires TkInter
# matplotlib.use('wxagg')  # requires wxPython


class massAnimation:
    def __init__(self):
        self.flag_init = True
        self.fig, self.ax = plt.subplots()
        self.handle = []
        self.L = 5.0
        self.w = 1.0
        self.theta = 45.0 * np.pi / 180.0

        # plot track
        R = np.array([[np.cos(self.theta), -np.sin(self.theta)],
                      [np.sin(self.theta), np.cos(self.theta)]])
        pts = np.array([[-3.0 * self.L - self.L / 5.0, 0.0], [0.0, 0.0]])
        pts = pts @ R
        self.handle.append(
            self.ax.plot(pts.T[0], pts.T[1], lw=2, c='black')
        )

        # plot wall
        pts = np.array([[-3.0 * self.L, 0],
                        [-3.0 * self.L, 2.0 * self.w]])
        pts = pts @ R
        self.handle.append(
            self.ax.plot(pts.T[0], pts.T[1], lw=1, c='black')
        )

        # set axis
        plt.axis([-3.0 * self.L - self.L / 5.0, self.L / 5.0,
                  - self.L / 5.0, 3.0 * self.L + 2.0 * self.L / 5.0])

    def update(self, state):
        z = state.item(0)
        self.draw_mass(z)
        self.draw_spring(z)
        self.ax.axis('equal')
        # Set initialization flag to False after first call
        if self.flag_init == True:
            self.flag_init = False

    def draw_mass(self, z):
        R = np.array([[np.cos(self.theta), -np.sin(self.theta)],
                      [np.sin(self.theta), np.cos(self.theta)]])
        pts = np.array([[z - self.w / 2.0 - 2.0 * self.L, 0],
                        [z + self.w / 2.0 - 2.0 * self.L, 0],
                        [z + self.w / 2.0 - 2.0 * self.L, self.w],
                        [z - self.w / 2.0 - 2.0 * self.L, self.w]])
        pts = pts @ R

        if self.flag_init == True:
            # Create the Rectangle patch and append its handle
            # to the handle list
            self.handle.append(
                mpatches.Polygon(pts, facecolor='blue', edgecolor='black')
            )
            self.ax.add_patch(self.handle[2]) # Add the patch to the axes
        else:
            self.handle[2].set_xy(pts)         # Update polygon

    def draw_spring(self, z):
        # specify x-y points of the spring
        R = np.array([[np.cos(self.theta), -np.sin(self.theta)],
                      [np.sin(self.theta), np.cos(self.theta)]])

        pts = np.array([[-3.0 * self.L, self.w / 2.0],
                        [z - self.w / 2.0 - 2.0 * self.L, self.w / 2.0]])
        pts = pts @ R
        # create spring on first call, update on subsequent calls
        if self.flag_init == True:
            # Create the line object and append its handle
            # to the handle list.
            line, =self.ax.plot(pts.T[0], pts.T[1], lw=1, c='blue')
            self.handle.append(line)
        else:
            self.handle[3].set_xdata(pts.T[0])
            self.handle[3].set_ydata(pts.T[1])


# Used see the animation from the command line
if __name__ == "__main__":

    simAnimation = massAnimation()
    z = 3
    state = np.array([[z],[0]])
    simAnimation.update(state)
    #plt.show()
    # Keeps the program from closing until the user presses a button.
    print('Press key to close')
    plt.waitforbuttonpress()
    plt.close()