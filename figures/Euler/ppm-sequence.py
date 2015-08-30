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

    def get_cubic_points(self, var):
        return self.vars[var].aint

    def get_parabola_coefficients(self, var):
        return self.vars[var].ap, self.vars[var].am, self.vars[var].a6

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
q_l = riemann.State(rho=1.0, u=0.0, p=0.1)

# right state
q_r = riemann.State(rho=0.1, u=0.0, p=0.125)


r = RiemannProblem(q_l, q_r)


#-----------------------------------------------------------------------------
#  plot 1: initial state
plt.clf()


# rho
plt.subplot(311)

r.draw_grid()
r.draw_var_avg("rho")
r.gr.clean_axes()

plt.title(r"$\rho$")

# u
plt.subplot(312)

r.draw_grid()
r.draw_var_avg("u")
r.gr.clean_axes()

plt.title("$u$")

# p
plt.subplot(313)

r.draw_grid()
r.draw_var_avg("p")
r.gr.clean_axes()

plt.title("$p$")

f = plt.gcf()
f.set_size_inches(8.0,7.0)

plt.tight_layout()

plt.savefig("ppm-seq-1.png")



# #------------- PLM -------------
# plt.clf()

# gpu.drawGrid(gr)

# gpu.labelCenter(gr, nzones/2,   r"$i$")
# gpu.labelCenter(gr, nzones/2-1, r"$i-1$")
# gpu.labelCenter(gr, nzones/2+1, r"$i+1$")
# gpu.labelCenter(gr, nzones/2-2, r"$i-2$")
# gpu.labelCenter(gr, nzones/2+2, r"$i+2$")


# n = 0
# while (n < nzones):
#     gpu.drawCellAvg(gr, n, a[n], color="0.5")
#     n += 1

# # compute the slopes
# da = gpu.lslopes(a, nolimit=1)
# lda = gpu.lslopes(a)

# n = 2
# while (n < nzones-2):
#     gpu.drawSlope(gr, n, da[n], a[n], color="r", ls=":")
#     gpu.drawSlope(gr, n, lda[n], a[n], color="r")
#     n += 1


# plt.axis([gr.xmin-0.5*gr.dx,gr.xmax+0.5*gr.dx, -0.25, 1.2])
# plt.axis("off")

# plt.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95)

# f = plt.gcf()
# f.set_size_inches(8.0,2.0)

# plt.savefig("piecewise-linear.eps")
# plt.savefig("piecewise-linear.png")


# #------------- PPM -------------
# plt.clf()

# gpu.drawGrid(gr)

# gpu.labelCenter(gr, nzones/2,   r"$i$")
# gpu.labelCenter(gr, nzones/2-1, r"$i-1$")
# gpu.labelCenter(gr, nzones/2+1, r"$i+1$")
# gpu.labelCenter(gr, nzones/2-2, r"$i-2$")
# gpu.labelCenter(gr, nzones/2+2, r"$i+2$")


# n = 0
# while (n < nzones):
#     gpu.drawCellAvg(gr, n, a[n], color="0.5")
#     n += 1


# # compute the parabolic coefficients
# ap, am, a6 = gpu.ppm(a, nolimit=1)
# lap, lam, la6 = gpu.ppm(a)

# n = 2
# while (n < nzones-2):
#     gpu.drawParabola(gr, n, ap[n], am[n], a6[n], color="r", ls=":")
#     gpu.drawParabola(gr, n, lap[n], lam[n], la6[n], color="r")
#     #plt.scatter([gr.xl[n], gr.xr[n]], [am[n], ap[n]])
#     n += 1




# plt.axis([gr.xmin-0.5*gr.dx,gr.xmax+0.5*gr.dx, -0.25, 1.2])
# plt.axis("off")

# plt.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95)

# f = plt.gcf()
# f.set_size_inches(8.0,2.0)


# plt.savefig("piecewise-parabolic.eps")
# plt.savefig("piecewise-parabolic.png")
