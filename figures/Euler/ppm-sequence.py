import math
import numpy as np
import matplotlib.pyplot as plt
import grid_plot as gp
import riemann

class RiemannProblem(object):
    def __init__(self, q_l, q_r, gamma=5./3., N=3):

        self.q_l = q_l
        self.q_r = q_r
        self.nx = 2*N

        self.gr = gp. FVGrid(self.nx)

        # allocate space -- we just care about N zones to the left and
        # right of the interface
        self.vars = {}

        rho = self.gr.scratch_array()
        rho[0:N] = q_l.rho
        rho[N:2**N] = q_r.rho

        u = self.gr.scratch_array()
        u[0:N] = q_l.u
        u[N:2*N] = q_r.u

        p = self.gr.scratch_array()
        p[0:N] = q_l.p
        p[N:2*N] = q_r.p

        rscale = max(q_l.rho, q_r.rho)
        uscale = max(q_l.u, q_r.u)
        pscale = max(q_l.p, q_r.p)

        # do the ppm reconstruction on all of these variables
        self.vars["rho"] = gp.PiecewiseParabolic(self.gr, rho, scale=rscale)
        self.vars["u"] = gp.PiecewiseParabolic(self.gr, u, scale=uscale)
        self.vars["p"] = gp.PiecewiseParabolic(self.gr, p, scale=pscale)

    def draw_cubic_points(self, var, color="r"):
        # these are defined on the i+1/2 interface
        # note that we require 2 zones on either side of the interface, so
        # we cannot draw points for the first 2 or last 2 zones
        aint = self.vars[var].aint
        plt.scatter(self.gr.xr[1:-2], aint[1:-2], marker="x", s=40, color=color, zorder=10)

    def draw_parabola(self, var):
        ap, am = self.vars[var].ap, self.vars[var].am
        for n in range(2,self.nx-2):
            self.vars[var].draw_parabola(n)

    def draw_grid(self):
        self.gr.draw_grid()

        self.gr.label_center(self.nx/2,   r"$i$")
        self.gr.label_center(self.nx/2-1, r"$i-1$")
        self.gr.label_center(self.nx/2+1, r"$i+1$")
        self.gr.label_center(self.nx/2-2, r"$i-2$")
        self.gr.label_center(self.nx/2+2, r"$i+2$")

    def draw_var_avg(self, var, scale=1.0):
        for n in range(self.gr.nx):
            self.vars[var].draw_cell_avg(n, color="r")


#-----------------------------------------------------------------------------



# multipage PDF
# http://matplotlib.org/examples/pylab_examples/multipage_pdf.html


# left state
q_l = riemann.State(rho=1.0, u=0.0, p=1.0)

# right state
q_r = riemann.State(rho=0.125, u=0.0, p=0.1)

r = RiemannProblem(q_l, q_r)


subidx = [311, 312, 313]
states = ["rho", "u", "p"]
state_labels = [r"$\rho$", r"$u$", r"$p$"]


#-----------------------------------------------------------------------------
#  plot 1: initial state
plt.clf()

for n, s in enumerate(states):

    plt.subplot(subidx[n])
    r.draw_grid()
    r.draw_var_avg(s)
    r.gr.clean_axes()

    plt.title(state_labels[n])

f = plt.gcf()
f.set_size_inches(8.0,7.0)
plt.tight_layout()

plt.savefig("ppm-seq-1.png")


#-----------------------------------------------------------------------------
#  plot 2: cubic points (three vertical)
plt.clf()

for n, s in enumerate(states):

    plt.subplot(subidx[n])
    r.draw_grid()
    r.draw_var_avg(s)
    r.draw_cubic_points(s)

    r.gr.clean_axes()

    plt.title(state_labels[n])

f = plt.gcf()
f.set_size_inches(8.0,7.0)
plt.tight_layout()

plt.savefig("ppm-seq-2.png")


#-----------------------------------------------------------------------------
#  plot 3: parabola (three vertical)
plt.clf()

for n, s in enumerate(states):

    plt.subplot(subidx[n])
    r.draw_grid()
    r.draw_var_avg(s)
    r.draw_cubic_points(s)
    r.draw_parabola(s)
    r.gr.clean_axes()

    plt.title(state_labels[n])

f = plt.gcf()
f.set_size_inches(8.0,7.0)
plt.tight_layout()

plt.savefig("ppm-seq-3.png")



#-----------------------------------------------------------------------------
#  plot 4: cell-center Riemann waves to show what hits the interface


#-----------------------------------------------------------------------------
#  plot 5: tracing (zoom in, 3 horizontal)



#-----------------------------------------------------------------------------
#  plot 6: final interface state (zoom in, 3 horizontal)



#-----------------------------------------------------------------------------
#  plot 7: Riemann phase



#-----------------------------------------------------------------------------
#  plot 8: Riemann solution



