# This is a refactorization of the grid_plot_util routines that allow one to
# plot grids

import math
import numpy as np
import matplotlib.pyplot as plt
import sys

class FVGrid(object):

    def __init__(self, nx, ng = 0, xmin=0.0, xmax=1.0, voff=0.0):

        # finite-volume
        self.nx = nx
        self.ng = ng
        self.xmin = xmin
        self.xmax = xmax

        self.ilo = ng
        self.ihi = ng+nx-1

        self.dx = (xmax - xmin)/float(nx)

        self.xl = (np.arange(2*ng+nx)-ng)*self.dx + xmin
        self.xr = (np.arange(2*ng+nx)+1-ng)*self.dx + xmin
        self.xc = 0.5*(self.xl + self.xr)

        # vertical offset (if we want to stack grids)
        self.voff = voff


    def draw_grid(self, center_only=0, draw_ghost=0,
                  emphasize_end=0, edge_ticks=1, color="k"):

        if center_only and emphasize_end:
            sys.exit("center_only and emphasize_end are incompatible")

        grid_top = 1.0

        if not draw_ghost:
            if center_only == 1:
                nstart = self.ng + self.nx/2-1
                nstop = self.ng + self.nx/2
            else:
                nstart = self.ilo
                nstop = self.ihi
        else:
            nstart = self.ilo -self.ng
            nstop = self.ihi + self.ng


        if emphasize_end:
            # horizontal line
            plt.plot([self.xl[nstart], self.xr[nstop]],
                       [self.voff,self.voff], color=color, lw=2)

        else:
            # horizontal line
            plt.plot([self.xl[nstart]-0.5*self.dx, self.xr[nstop]+0.5*self.dx],
                       [self.voff,self.voff], color=color, lw=2)

        # draw first left edge
        plt.plot([self.xl[nstart], self.xl[nstart]],
                   [self.voff, grid_top+self.voff], color=color, lw=2)


        for n in range(nstart, nstop+1):

            # emphasize?
            if emphasize_end and n == self.ilo:
                plt.plot([self.xl[n], self.xl[n]],
                           [self.voff, grid_top+self.voff], color=color, lw=4)

            # draw right edge
            if emphasize_end and n == self.ihi:
                plt.plot([self.xr[n], self.xr[n]], [self.voff, grid_top+self.voff],
                           color=color, lw=4)
            else:
                plt.plot([self.xr[n], self.xr[n]], [self.voff, grid_top+self.voff],
                           color=color, lw=2)

            # draw center marker
            plt.plot([self.xc[n], self.xc[n]], [-0.05+self.voff, self.voff], color=color)

            # draw edge marker
            if n == nstart and edge_ticks:
                plt.plot([self.xl[nstart], self.xl[nstart]], [-0.05+self.voff, self.voff],
                           color=color)

            if edge_ticks:
                plt.plot([self.xr[n], self.xr[n]], [-0.05+self.voff, self.voff], color=color)


    def label_center(self, idx, string, fontsize="small"):

        plt.text(self.xc[idx], self.voff-0.1, string,
                 horizontalalignment='center', verticalalignment='top',
                 fontsize=fontsize)

    def label_edge(self, idx, string, fontsize="small"):

        plt.text(self.xl[idx], self.voff-0.075, string,
                 horizontalalignment='center', verticalalignment='top',
                 fontsize=fontsize)

    def label_cell_center(self, idx, string):

        plt.text(self.xc[idx], self.voff+0.5, string,
                 horizontalalignment='center', verticalalignment='center',
                 fontsize="large")

    def mark_cell_left_state(self, idx, string, color="k", value=0.5,
                             vertical="center"):

        plt.scatter(self.xl[idx]+0.05*self.dx, self.voff+value, marker="x", color=color)

        plt.text(self.xl[idx]+0.075*self.dx, self.voff+value, string,
                 horizontalalignment='left', verticalalignment=vertical, color=color)

    def mark_cell_right_state(self, idx, string, color="k", value=0.5,
                              vertical="center"):

        plt.scatter(self.xr[idx]-0.05*self.dx, self.voff+value, marker="x", color=color)

        plt.text(self.xr[idx]-0.075*self.dx, self.voff+value, string,
                 horizontalalignment='right', verticalalignment=vertical, color=color)

    def clean_axes(self):
        plt.xlim(self.xmin-0.5*self.dx, self.xmax+0.5*self.dx)
        plt.axis("off")


class PiecewiseConstant(object):
    def __init__(self, gr, a):
        if not len(a) == len(gr.xc):
            sys.exit("ERROR: grid length != data length")

        self.gr = gr
        self.a = a

    def label_cell_avg(self, idx, string, color="k"):
        plt.text(self.gr.xc[idx], self.gr.voff+self.a[idx]+0.1, string,
                 horizontalalignment='center', verticalalignment='bottom',
                 fontsize="large", color=color)

    def draw_cell_avg(self, idx, color="0.5", ls="-"):
        plt.plot([self.gr.xl[idx], self.gr.xr[idx]], 
                 [self.gr.voff+self.a[idx], self.gr.voff+self.a[idx]], color=color, ls=ls)


class PLM(PiecewiseConstant):
    def __init__(self, gr, a, nolimit=0):

        PiecewiseConstant.__init__(self, gr, a)
        
        self.slope = np.zeroslike(self.a)

        # calculate the slopes
        for n in range(1, len(self.a)-1):
            test = (self.a[n+1] - self.a[n])*(self.a[n] - self.a[n-1])
            da = 0.5*(self.a[n+1] - self.a[n-1])

            if not nolimit:
                if test > 0.0:
                    self.slope[n] = min(math.fabs(da), 
                                        min(2.0*math.fabs(self.a[n+1] - self.a[n]),
                                            2.0*math.fabs(self.a[n] - self.a[n-1]))) * \
                                            np.sign(self.a[n+1] - self.a[n-1])
                else:
                    self.slope[n] = 0.0
            else:
                self.slope[n] = da


    def draw_slope(self, idx, color="r", ls="-"):

        yl = self.slope[idx]*(self.gr.xl[idx] - self.gr.xc[idx])/self.gr.dx + \
             self.gr.voff+self.a[idx]
        yr = self.slope[idx]*(self.gr.xr[idx] - self.gr.xc[idx])/self.gr.dx + \
             self.gr.voff+value

        plt.plot([self.gr.xl[idx], self.gr.xr[idx]], [yl, yr],
                 color=color, ls=ls, lw=1, zorder=10)


    def slope_trace_left(self, idx, sigma, color="0.5"):

        # sigma is the fraction of the domain -- the CFL number
        x = np.linspace(self.gr.xr[idx]-sigma*self.gr.dx, self.gr.xr[idx], 50)

        a = self.gr.voff+self.a[idx] + (self.slope[idx]/self.gr.dx) * (x - self.gr.xc[idx])

        # vertices of a polygon
        xx = np.zeros(len(x) + 3, dtype=np.float64)
        yy = np.zeros(len(x) + 3, dtype=np.float64)

        xx[0:len(x)] = x
        xx[len(x):] = [self.gr.xr[idx], 
                       self.gr.xr[idx]-sigma*self.gr.dx, 
                       self.gr.xr[idx]-sigma*self.gr.dx]

        yy[0:len(x)] = a
        yy[len(x):] = [self.gr.voff, self.gr.voff, a[0]]

        plt.fill(xx, yy, color=color, lw=1, zorder=-1)


    def slope_trace_right(self, idx, slope, color="0.5"):

        # sigma is the fraction of the domain -- the CFL number
        x = np.linspace(self.gr.xl[idx], self.gr.xl[idx]+sigma*self.dx, 50)

        a = self.gr.voff+self.a[idx] + (self.slope[idx]/self.gr.dx) * (x - self.gr.xc[idx])

        # vertices of a polygon
        xx = np.zeros(len(x) + 3, dtype=np.float64)
        yy = np.zeros(len(x) + 3, dtype=np.float64)

        xx[0:len(x)] = x
        xx[len(x):] = [self.gr.xl[idx]+sigma*self.gr.dx, 
                       self.gr.xl[idx], 
                       self.gr.xl[idx]]

        yy[0:len(x)] = a
        yy[len(x):] = [self.gr.voff, self.gr.voff, a[0]]

        plt.fill(xx, yy, color=color, lw=1, zorder=-1)


    def evolve_to_right(self, idx, sigma, color="0.5", ls="-"):

        # sigma is the fraction of the domain -- the CFL number
        # show the reconstructed profile as we evolve to the right

        xm = np.linspace(self.gr.xr[idx-1]-sigma*self.gr.dx, self.gr.xr[idx-1], 50)
        am = self.a[idx-1] + (self.slopes[idx-1]/self.gr.dx) * (xm - self.gr.xc[idx-1])
        xm = xm + sigma*self.gr.dx

        xp = np.linspace(self.gr.xl[idx], self.gr.xl[idx]+(1.0-sigma)*self.gr.dx, 50)
        ap = self.a[idx] + (self.slopes[idx]/gr.dx) * (xp - gr.xc[idx])
        xp = xp + sigma*self.gr.dx

        plt.plot(xm, am, color=color, lw=1, ls=ls, zorder=10)
        plt.plot(xp, ap, color=color, lw=1, ls=ls, zorder=10)



#-----------------------------------------------------------------------------
def ppm_cubic(a, nolimit=0):
    # this computes the (limited) interface states just from
    # cubic interpolation through the 4 zones centered on an
    # interface
    aint = np.zeros(len(a), dtype=np.float64)

    n = 1
    while n < len(a)-2:

        da0 = 0.5*(a[n+1] - a[n-1])
        dap = 0.5*(a[n+2] - a[n])

        if not nolimit:

            if (a[n+1] - a[n])*(a[n] - a[n-1]) > 0.0:
                da0 = np.sign(da0)*min(math.fabs(da0),
                                          min(2.0*math.fabs(a[n] - a[n-1]),
                                              2.0*math.fabs(a[n+1] - a[n])) )
            else:
                da0 = 0.0


            if ( (a[n+2] - a[n+1])*(a[n+1] - a[n]) > 0.0):
                dap = np.sign(dap)*min(math.fabs(dap),
                                          min(2.0*math.fabs(a[n+1] - a[n]),
                                              2.0*math.fabs(a[n+2] - a[n+1])) )
            else:
                dap = 0.0


        # cubic
        aint[n] = 0.5*(a[n] + a[n+1]) - (1.0/6.0)*(dap - da0)

        n += 1

    return aint


def ppm(a, nolimit=0):

    ap = np.zeros(len(a), dtype=np.float64)
    am = np.zeros(len(a), dtype=np.float64)
    a6 = np.zeros(len(a), dtype=np.float64)

    # parabola of form:
    #    a(xi) = aminus + xi*(aplus - aminus + a6 * (1-xi) )
    # with  xi = (x - xl)/dx

    ap[:] = ppm_cubic(a, nolimit=nolimit)
    am[1:] = ap[:-1]

    if not nolimit:

        n = 2
        while n < len(a)-2:

            if ( (ap[n] - a[n])*(a[n] - am[n]) <= 0.0):
                am[n] = a[n]
                ap[n] = a[n]

            elif ( (ap[n] - am[n])*(a[n] - 0.5*(am[n] + ap[n])) >
                   (ap[n] - am[n])**2/6.0 ):
                am[n] = 3.0*a[n] - 2.0*ap[n]

            elif ( -(ap[n] - am[n])**2/6.0 >
                    (ap[n] - am[n])*(a[n] - 0.5*(am[n] + ap[n])) ):
                ap[n] = 3.0*a[n] - 2.0*am[n]

            n += 1

    n = 2
    while n < len(a)-2:
        a6[n] = 6.0*a[n] - 3.0*(am[n] + ap[n])
        n += 1

    return ap, am, a6


def drawParabola(gr, idx, ap, am, a6, color="r", ls="-"):

    x = np.linspace(gr.xl[idx], gr.xr[idx], 50)
    xi = (x - gr.xl[idx])/gr.dx
    print xi
    a = am + xi*(ap - am + a6 * (1.0-xi) )

    plt.plot(x, a,
               color=color, ls=ls, lw=1, zorder=10)


def ppmTraceLeft(gr, idx, ap, am, a6, sigma, color="0.5"):

    # sigma is the fraction of the domain from the right interface
    # (since we are tracing to create the left state there).

    x = np.linspace(gr.xr[idx]-sigma*gr.dx, gr.xr[idx], 50)
    xi = (x - gr.xl[idx])/gr.dx

    a = am + xi*(ap - am + a6 * (1.0-xi) )

    xx = np.zeros(len(x) + 3, dtype=np.float64)
    yy = np.zeros(len(x) + 3, dtype=np.float64)

    xx[0:len(x)] = x
    xx[len(x):] = [gr.xr[idx], gr.xr[idx]-sigma*gr.dx, gr.xr[idx]-sigma*gr.dx]

    yy[0:len(x)] = a
    yy[len(x):] = [0.0, 0.0, a[0]]

    plt.fill(xx, yy, color=color, lw=1, zorder=-1)
