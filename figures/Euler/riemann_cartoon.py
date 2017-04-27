import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['mathtext.rm'] = 'serif'

# font sizes
mpl.rcParams['font.size'] = 12
mpl.rcParams['legend.fontsize'] = 'large'
mpl.rcParams['figure.titlesize'] = 'medium'


class RiemannWaves(object):
    """draw a diagram with the 3 waves for the Euler Riemann problem"""

    def __init__(self, left=-1, contact=-1, right=1):
        """for each of the waves, -1 is left of the interface, +1 is right of
	the interface"""

        self.left = left
        self.contact = contact
        self.right = right

        # these control the spacing of the waves in the figure
        self.theta = np.radians(25)
        self.h = 1.0

    def draw(self, output="plot.png", label=None):
        """ draw the wave structure """

        # by default, we'll make the waves theta degrees apart
        if self.right < 0:
            angle3 = -0.75*self.theta
        elif self.right > 0 and self.contact < 0:
            angle3 = 0.534*self.theta
        elif self.right > 0 and self.contact > 0 and self.left < 0:
            angle3 = 0.966*self.theta
        else:
            angle3 = 2.25*self.theta

        anglec = angle3 - 0.75*self.theta
        angle1 = anglec - 0.75*self.theta

        print(np.degrees(angle1), np.degrees(anglec), np.degrees(angle3))

        # each wave has a point on the origin.  The maximum height we want to draw to is self.h
        x1 = self.h*np.tan(angle1)
        xc = self.h*np.tan(anglec)
        x3 = self.h*np.tan(angle3)

        print(x1, xc, x3)

        plt.clf()

        # draw the waves
        plt.plot([0, x1], [0, self.h], color="C0")
        plt.plot([0, xc], [0, self.h], color="C0")
        plt.plot([0, x3], [0, self.h], color="C0")

        # label regions
        x1h = 0.75*self.h*np.tan(angle1)
        xch = 0.75*self.h*np.tan(anglec)
        x3h = 0.75*self.h*np.tan(angle3)

        dx = max(xch - x1h, x3h - xch)
        xL = x1h - 0.6*dx
        xLstar = 0.5*(x1h + xch)
        xRstar = 0.5*(xch + x3h)
        xR = x3h + 0.6*dx

        plt.text(xL, 0.65*self.h, r"$L$", horizontalalignment="right", color="C0")
        plt.text(xLstar, 0.75*self.h, r"$L_\star$", horizontalalignment="center", color="C0")
        plt.text(xRstar, 0.75*self.h, r"$R_\star$", horizontalalignment="center", color="C0")
        plt.text(xR, 0.65*self.h, r"$R$", horizontalalignment="left", color="C0")

        # draw the grid
        L = 1.1*max(abs(x1), abs(x3))

        plt.plot([-L, L], [0, 0], color="k")
        plt.plot([0, 0], [0, 1.2*self.h], color="k")

        f = plt.gcf()

        # label?
        if label is not None:
            plt.text(0.1, 0.9, label, transform=f.transFigure)

        plt.axis("off")
        plt.subplots_adjust(left=0.02, right=0.95, bottom=0.05, top=0.95)


        f.set_size_inches(5.0, 3.5)

        plt.tight_layout()
        plt.savefig(output, bbox_inches="tight", pad_inches=0)


if __name__ == "__main__":

    rw = RiemannWaves(left=-1, contact=-1, right=-1)
    rw.draw("riemann_waves_ifc_R.pdf", label="(a)")

    rw = RiemannWaves(left=-1, contact=-1, right=1)
    rw.draw("riemann_waves_ifc_Rstar.pdf", label="(b)")

    rw = RiemannWaves(left=-1, contact=1, right=1)
    rw.draw("riemann_waves_ifc_Lstar.pdf", label="(c)")

    rw = RiemannWaves(left=1, contact=1, right=1)
    rw.draw("riemann_waves_ifc_L.pdf", label="(d)")
