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
    """draw a diagram with a left rarefaction and 2 other waves """

    def __init__(self, rpos=-1):
        """rpos is position of rarefaction, -1 is left of the interface, +1 is
	right of the interface

        """

        self.rpos = rpos

        # these control the spacing of the waves in the figure
        self.h = 1.0

        # this controls the spacing of the structure in the rarefaction
        self.rare_theta = np.radians(2.5)

        if rpos < 0:
            self.head = np.radians(-45)
            self.center = np.radians(5)
            self.right = np.radians(40)
        elif rpos == 0:
            self.head = 0 - 2.5*self.rare_theta
            self.center = np.radians(30)
            self.right = np.radians(45)
        else:
            self.head = np.radians(5)
            self.center = np.radians(35)
            self.right = np.radians(50)

        self.tail = self.head + 5*self.rare_theta

    def draw(self, output="plot.png", label=None):
        """ draw the wave structure """

        # each wave has a point on the origin.  The maximum height we
        # want to draw to is self.h -- we'll do the center and right
        xc = self.h*np.tan(self.center)
        x3 = self.h*np.tan(self.right)

        plt.clf()

        # draw the waves
        plt.plot([0, xc], [0, self.h], color="C0")
        plt.plot([0, x3], [0, self.h], color="C0")

        # rarefaction
        for t in np.arange(self.head, self.tail+self.rare_theta, 
                           self.rare_theta):
            xr = self.h*np.tan(t)
            plt.plot([0, xr], [0, self.h], color="C0")

        xhead = self.h*np.tan(self.head)
        xtail = self.h*np.tan(self.tail) #+self.rare_theta)

        # draw the grid
        L = 1.1*max(abs(xhead), abs(x3))

        plt.plot([-L, L], [0, 0], color="k", zorder=-100)
        plt.plot([0, 0], [0, 1.2*self.h], color="k", zorder=-100)

        plt.text(xhead, self.h, r"$\lambda_\mathrm{head}$", 
                 horizontalalignment="right", color="C0")

        plt.text(xtail, self.h, r"$\lambda_\mathrm{tail}$", 
                 horizontalalignment="left", color="C0")


        dx = abs(0.75*self.h*np.tan(0.5*(self.head + self.tail)))

        if self.rpos == -1:
            fac = 0.4
        elif self.rpos == 0:
            fac = 1.0
        else:
            fac = 1.0

        xh = 0.75*self.h*np.tan(self.head)
        xt = 0.75*self.h*np.tan(self.tail)
        xc = 0.75*self.h*np.tan(self.center)
        xr = 0.75*self.h*np.tan(self.right)

        if dx == 0:
            dx = abs(xh)

        xL = xh - fac*dx
        xLstar = 0.5*(xt + xc)
        xRstar = 0.5*(xc + xr)
        xR = xr + fac*dx

        plt.text(xL, 0.65*self.h, r"$L$", horizontalalignment="right", color="C0")
        plt.text(xLstar, 0.75*self.h, r"$L_\star$", horizontalalignment="center", color="C0")
        plt.text(xRstar, 0.75*self.h, r"$R_\star$", horizontalalignment="center", color="C0")
        plt.text(xR, 0.65*self.h, r"$R$", horizontalalignment="left", color="C0")


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

    rw = RiemannWaves(rpos=-1)
    rw.draw("rarefaction_left.pdf", label="(a)")

    rw = RiemannWaves(rpos=0)
    rw.draw("rarefaction_center.pdf", label="(b)")

    rw = RiemannWaves(rpos=1)
    rw.draw("rarefaction_right.pdf", label="(c)")
