# plot the Hugoniot loci for a compressible Riemann problem

import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as optimize

class State(object):
    def __init__(self, p=1.0, u=0.0, rho=1.0):
         self.p = p
         self.u = u
         self.rho = rho

    def __str__(self):
        return "rho: {}; u: {}; p: {}".format(self.rho, self.u, self.p)

class RiemannProblem(object):
    def __init__(self, left_state, right_state, gamma=1.4):
        self.left = left_state
        self.right = right_state
        self.gamma = gamma

    def u_hugoniot(self, p, side):

        if side == "left":
            state = self.left
            s = 1.0
        elif side == "right":
            state = self.right
            s = -1.0

        c = np.sqrt(self.gamma*state.p/state.rho)

        if p < state.p:
            # rarefaction
            u = state.u + s*(2.0*c/(self.gamma-1.0))* \
                (1.0 - (p/state.p)**((self.gamma-1.0)/(2.0*self.gamma)))
        else:
            # shock
            beta = (self.gamma+1.0)/(self.gamma-1.0)
            u = state.u + s*(2.0*c/np.sqrt(2.0*self.gamma*(self.gamma-1.0)))* \
                (1.0 - p/state.p)/np.sqrt(1.0 + beta*p/state.p)

        return u

    def find_star_state(self, p_min=0.001, p_max=1000.0):
        # we need to root-find on
        pstar = optimize.brentq(lambda p: self.u_hugoniot(p, "left") - self.u_hugoniot(p, "right"),
                               p_min, p_max)
        ustar = self.u_hugoniot(pstar, "left")

        return pstar, ustar

    def plot_hugoniot(self, p_min = 0.0, p_max=1.5, N=200):

        p = np.linspace(p_min, p_max, num=N)
        u_left = np.zeros_like(p)
        u_right = np.zeros_like(p)

        for n in range(N):
            u_left[n] = self.u_hugoniot(p[n], "left")
        ish = np.where(p > self.left.p)
        ir = np.where(p < self.left.p)

        plt.plot(p[ish], u_left[ish], c="b", ls=":", lw=2)
        plt.plot(p[ir], u_left[ir], c="b", ls="-", lw=2)
        plt.scatter([self.left.p], [self.left.u], marker="x", c="b", s=40)

        for n in range(N):
            u_right[n] = self.u_hugoniot(p[n], "right")
        ish = np.where(p > self.right.p)
        ir = np.where(p < self.right.p)

        plt.plot(p[ish], u_right[ish], c="r", ls=":", lw=2)
        plt.plot(p[ir], u_right[ir], c="r", ls="-", lw=2)
        plt.scatter([self.right.p], [self.right.u], marker="x", c="r", s=40)

        du = 0.025*(max(np.max(u_left), np.max(u_right)) -
                    min(np.min(u_left), np.min(u_right)))

        plt.text(self.left.p, self.left.u+du, "left",
                 horizontalalignment="center", color="b")

        plt.text(self.right.p, self.right.u+du, "right",
                 horizontalalignment="center", color="r")

        plt.xlim(p_min, p_max)

        plt.xlabel(r"$p$", fontsize="large")
        plt.ylabel(r"$u$", fontsize="large")

        legs = []
        legnames = []

        legs.append(plt.Line2D((0,1),(0,0), color="k", ls=":", marker=None))
        legnames.append("shock")

        legs.append(plt.Line2D((0,1),(0,0), color="k", ls="-", marker=None))
        legnames.append("rarefaction")

        plt.legend(legs, legnames, frameon=False, loc="best")

        plt.tight_layout()

        plt.savefig("riemann-phase.pdf")


if __name__ == "__main__":

    # setup the problem -- Sod
    left = State(p = 1.0, u = 0.0, rho = 1.0)
    right = State(p = 0.1, u = 0.0, rho = 0.125)

    rp = RiemannProblem(left, right)

    rp.plot_hugoniot()

    pstar, ustar = rp.find_star_state()

    print pstar, ustar
