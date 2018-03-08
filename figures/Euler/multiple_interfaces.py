import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import random

import grid_plot as gp

mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['mathtext.rm'] = 'serif'

# font sizes
mpl.rcParams['font.size'] = 12
mpl.rcParams['legend.fontsize'] = 'large'
mpl.rcParams['figure.titlesize'] = 'medium'


class RiemannWaves(object):
    """draw a diagram with a the three waves """

    def __init__(self, x0, dx):
        """ create a random collection of Riemann waves """

        states = ["shock", "rarefaction"]

        self.x0 = x0
        self.dx = dx

        # first pick if the angle of the contact, from a normalize distribution centered on 0
        _a = np.degrees(0.4*np.random.randn())
        while _a > 45:
            _a = np.degrees(0.4*np.random.randn())
        self.contact_pos = _a
        print(_a)

        # now the left state -- first set how far it is from the contact
        _a = np.degrees(0.4*np.random.randn())
        while _a > 45:
            _a = np.degrees(0.4*np.random.randn())
        print(_a)
        self.left_pos = self.contact_pos - np.abs(_a)
        self.left_type = random.choice(states)

        # and the right
        _a = np.degrees(0.4*np.random.randn())
        while _a > 45:
            _a = np.degrees(0.4*np.random.randn())
        print(_a)
        self.right_pos = self.contact_pos + np.abs(_a)
        self.right_type = random.choice(states)

    def draw(self, output="plot.png", label=None):
        """ draw the wave structure """

        print(self.left_pos, self.contact_pos, self.right_pos)

        # each wave has a point on the origin.  The maximum height we
        # want to draw to is self.dx -- we'll do the center and right
        x1 = self.x0 + self.dx * np.tan(np.radians(self.left_pos))
        xc = self.x0 + self.dx * np.tan(np.radians(self.contact_pos))
        x3 = self.x0 + self.dx * np.tan(np.radians(self.right_pos))

        # draw the waves
        if self.left_type == "shock":
            plt.plot([self.x0, x1], [0, self.dx], color="C0", lw=3)
        else:
            for n in range(5):
                x1 = self.x0 + self.dx * np.tan(np.radians(self.left_pos - 2*n))
                plt.plot([self.x0, x1], [0, self.dx], color="C0", lw=1)

        plt.plot([self.x0, xc], [0, self.dx], color="C1")

        if self.right_type == "shock":
            plt.plot([self.x0, x3], [0, self.dx], color="C0", lw=3)
        else:
            for n in range(5):
                x3 = self.x0 + self.dx * np.tan(np.radians(self.right_pos + 2*n))
                plt.plot([self.x0, x3], [0, self.dx], color="C0", lw=1)


if __name__ == "__main__":


    xmin = 0.0
    xmax = 8.0

    nzones = 4
    ng = 0

    gr = gp.FVGrid(nzones, xmin=xmin, xmax=xmax)

    gr.draw_grid()

    gr.label_edge(0, r"$i-\mythreehalf$", fontsize="medium")
    gr.label_edge(1, r"$i-\myhalf$", fontsize="medium")
    gr.label_edge(2, r"$i+\myhalf$", fontsize="medium")
    gr.label_edge(3, r"$i+\mythreehalf$", fontsize="medium")
    gr.label_edge(3, r"$i+\myfivehalf$", fontsize="medium", right_edge=True)

    gr.label_center(0, r"$i-1$", fontsize="medium")
    gr.label_center(1, r"$i$", fontsize="medium")
    gr.label_center(2, r"$i+1$", fontsize="medium")
    gr.label_center(3, r"$i+2$", fontsize="medium")

    riemann = []
    for n in range(nzones):
        riemann.append(RiemannWaves(gr.xl[n], 0.5*gr.dx))
    riemann.append(RiemannWaves(gr.xr[nzones-1], 0.5*gr.dx))

    # draw
    for r in riemann:
        r.draw()

    gr.clean_axes()

    f = plt.gcf()
    f.set_size_inches(8.0, 2.5)

    plt.tight_layout()

    plt.savefig("multiple_interfaces.pdf", dpi=150)
