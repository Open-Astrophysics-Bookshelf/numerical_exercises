import math
import numpy as np
import matplotlib.pyplot as plt
import grid_plot_util as gpu


class State(object):
    def __init__(self, rho, u, p):
        self.rho = rho
        self.u = u
        self.p = p

class RiemannProblem(object):
    def __init__(self, q_l, q_r, gamma=5./3., N=3):

        self.q_l = q_l
        self.q_r = q_r
        self.nx = 2*N

        self.gr = gpu.grid(self.nx)

        # allocate space -- we just care about N zones to the left and
        # right of the interface
        self.vars = {}

        self.vars["rho"] = np.zeros(2*N, dtype=np.float64)
        self.vars["u"] = np.zeros(2*N, dtype=np.float64)
        self.vars["p"] = np.zeros(2*N, dtype=np.float64)

        self.vars["rho"][0:N] = q_l.rho
        self.vars["rho"][N:2*N] = q_r.rho

        self.vars["u"][0:N] = q_l.u
        self.vars["u"][N:2*N] = q_r.u

        self.vars["p"][0:N] = q_l.p
        self.vars["p"][N:2*N] = q_r.p

    def get_cubic_points(self, var):
        q = self.vars[var]
        return gpu.ppm_cubic(q)

    def get_parabola_coefficients(self, var):
        q = self.vars[var]
        return gpu.ppm(q)
        
    def draw_grid(self):
        gpu.drawGrid(self.gr)

        gpu.labelCenter(self.gr, self.nx/2,   r"$i$")
        gpu.labelCenter(self.gr, self.nx/2-1, r"$i-1$")
        gpu.labelCenter(self.gr, self.nx/2+1, r"$i+1$")
        gpu.labelCenter(self.gr, self.nx/2-2, r"$i-2$")
        gpu.labelCenter(self.gr, self.nx/2+2, r"$i+2$")
        
    def draw_var_avg(self, var):
        q = self.vars[var]        
        
        for n in range(r.nx):
            gpu.drawCellAvg(self.gr, n, q[n], color="r")



#-----------------------------------------------------------------------------



# multipage PDF
# http://matplotlib.org/examples/pylab_examples/multipage_pdf.html


# left state
q_l = State(1.0, 0.0, 0.1)

# right state
q_r = State(0.1, 0.0, 0.125)


r = RiemannProblem(q_l, q_r)

plt.clf()



r.draw_grid()
r.draw_var_avg("rho")


plt.axis([r.gr.xmin-0.5*r.gr.dx,r.gr.xmax+0.5*r.gr.dx, -0.25, 1.2])
plt.axis("off")

plt.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95)

f = plt.gcf()
f.set_size_inches(8.0,2.0)

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

               


