# grid_plot is a collection of classes for defining and drawing
# finite-volume and finite-difference grids in 1- and 2-d.

import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# Use LaTeX for rendering
mpl.rcParams["text.usetex"] = True
# load the xfrac package
mpl.rcParams["text.latex.preamble"].append(r'\usepackage{xfrac}')
mpl.rcParams["text.latex.preamble"].append(r'\newcommand{\myhalf}{\sfrac{1}{2}}')
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['mathtext.rm'] = 'serif'

# font sizes
mpl.rcParams['font.size'] = 12
mpl.rcParams['legend.fontsize'] = 'large'
mpl.rcParams['figure.titlesize'] = 'medium'

import sys

class FDGrid(object):
    """ a 1-d finite-difference grid """

    def __init__(self, nx, ng=0, xmin=0.0, xmax=1.0, voff=0.0):

        # finite-difference -- the xc will now be the node (where
        # the data lives)
        self.nx = nx
        self.ng = ng
        self.xmin = xmin
        self.xmax = xmax
        
        self.ilo = ng
        self.ihi = ng+nx-1

        self.dx = (xmax-xmin)/float(nx-1)

        self.xl = None
        self.xr = None
        self.xc = xmin + (np.arange(2*ng+nx)-ng)*self.dx

        # vertical offset (if we want to stack grids)
        self.voff = voff

    def draw_grid(self, draw_ghost=0, emphasize_end=0, edge_ticks=1, color="k"):

        grid_top = 0.25

        if not draw_ghost:
            nstart = self.ilo
            nstop = self.ihi
        else:
            nstart = self.ilo-self.ng
            nstop = self.ihi+self.ng

        if not draw_ghost:
            if emphasize_end:
                plt.plot([self.xmin, self.xmax], 
                         [self.voff, self.voff], color=color, lw=2)
            else:
                plt.plot([self.xmin-0.5*self.dx, self.xmax+0.5*self.dx], 
                         [self.voff,self.voff], color=color, lw=2)
        else:
            plt.plot([self.xmin-self.ng*self.dx, self.xmin],
                     [self.voff,self.voff], color=color, lw=2, ls=":")
            plt.plot([self.xmax, self.xmax+self.ng*self.dx], 
                     [self.voff,self.voff], color=color, lw=2, ls=":")
            plt.plot([self.xmin, self.xmax], 
                     [self.voff,self.voff], color=color, lw=2)

        for n in range(nstart, nstop+1):

            # draw center (node) indicator line
            if n < self.ilo or n > self.ihi:
                plt.plot([self.xc[n], self.xc[n]], 
                         [-0.05+self.voff, grid_top+self.voff], 
                         color=color, ls=":", lw=2)
            else:
                plt.plot([self.xc[n], self.xc[n]], 
                         [-0.05+self.voff, grid_top+self.voff], color=color, lw=2)
      

        if emphasize_end:
            plt.plot([self.xc[self.ilo], self.xc[self.ilo]], 
                     [-0.05+self.voff, grid_top+self.voff], color=color, lw=4)

            plt.plot([self.xc[self.ihi], self.xc[self.ihi]], 
                     [-0.05+self.voff, grid_top+self.voff], color=color, lw=4)

    def label_node(self, idx, string, fontsize="small"):

        plt.text(self.xc[idx], self.voff-0.1, string,
                 horizontalalignment='center', verticalalignment='top',
                 fontsize=fontsize)

    def label_node_data(self, idx, string):

        plt.text(self.xc[idx], self.voff+0.5, string,
                 horizontalalignment='center', verticalalignment='center',
                 fontsize="large")

    def label_value(self, idx, value, string, color="k", fontsize="large"):

        plt.text(self.xc[idx], self.voff+value+0.1, string,
                 horizontalalignment='center', verticalalignment='bottom',
                 fontsize=fontsize, color=color)

    def draw_data(self, idx, value, color="0.5", marker="o"):
        plt.scatter([self.xc[idx]], [self.voff+value], 
                    color=color, marker=marker, zorder=100)


    def label_dx(self, idx):
        # idx is the right edge of the dx interval drawn
        plt.plot([self.xc[idx-1], self.xc[idx-1]], [-0.35,-0.25], color="k")
        plt.plot([self.xc[idx], self.xc[idx]], [-0.35,-0.25], color="k")
        plt.plot([self.xc[idx-1], self.xc[idx]], [-0.3,-0.3], color="k")
        plt.text(0.5*(self.xc[idx-1] + self.xc[idx]), -0.45, r"$\Delta x$",
                 horizontalalignment="center")

            
    def clean_axes(self, show_ghost=False, padding=True, ylim=None, pad_fac=1.0):
        xmin = self.xmin
        xmax = self.xmax
        if show_ghost:
            xmin -= self.ng*self.dx
            xmax += self.ng*self.dx
        if padding:
            xmin -= pad_fac*self.dx
            xmax += pad_fac*self.dx

        plt.xlim(xmin, xmax)

        if ylim is not None:
            plt.ylim(ylim)
        plt.axis("off")
        plt.subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.95)


class FVGrid(object):
    """ a 1-d finite-volume grid """

    def __init__(self, nx, ng=0, xmin=0.0, xmax=1.0, voff=0.0):

        # finite-volume or cell-centered finite-difference
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

    def scratch_array(self):
        return np.zeros(2*self.ng+self.nx, dtype=np.float64)

    def draw_grid(self, center_only=0, draw_ghost=0,
                  emphasize_end=0, draw_end=True,
                  edge_ticks=True, color="k"):

        if center_only and emphasize_end:
            sys.exit("center_only and emphasize_end are incompatible")

        grid_top = 1.0

        if not draw_ghost:
            if center_only == 1:
                nstart = self.ng + self.nx//2-1
                nstop = self.ng + self.nx//2
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
        if draw_end:
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
            elif n < nstop or (n == nstop and draw_end):
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
                             vertical="center", fontsize="medium"):

        plt.scatter(self.xl[idx]+0.05*self.dx, self.voff+value, marker="x", color=color)

        plt.text(self.xl[idx]+0.075*self.dx, self.voff+value, string,
                 horizontalalignment='left', verticalalignment=vertical, color=color,
                 fontsize=fontsize)

    def mark_cell_right_state(self, idx, string, color="k", value=0.5,
                              vertical="center", fontsize="medium"):

        plt.scatter(self.xr[idx]-0.05*self.dx, self.voff+value, marker="x", color=color)

        plt.text(self.xr[idx]-0.075*self.dx, self.voff+value, string,
                 horizontalalignment='right', verticalalignment=vertical, color=color,
                 fontsize=fontsize)

    def label_dx(self, idx):
        # idx is the right edge of the dx interval drawn
        plt.plot([self.xr[idx-1], self.xr[idx-1]], [-0.35,-0.25], color="k")
        plt.plot([self.xr[idx], self.xr[idx]], [-0.35,-0.25], color="k")
        plt.plot([self.xr[idx-1], self.xr[idx]], [-0.3,-0.3], color="k")
        plt.text(self.xc[idx], -0.45, r"$\Delta x$", horizontalalignment="center")

    def label_center_dx(self, idx):
        # idx is the right edge of the dx interval drawn
        plt.plot([self.xc[idx-1], self.xc[idx-1]], [-0.35,-0.25], color="k")
        plt.plot([self.xc[idx], self.xc[idx]], [-0.35,-0.25], color="k")
        plt.plot([self.xc[idx-1], self.xc[idx]], [-0.3,-0.3], color="k")
        plt.text(0.5*(self.xc[idx-1] + self.xc[idx]), -0.45, r"$\Delta x$",
                 horizontalalignment="center")

    def clean_axes(self, show_ghost=False, padding=True, ylim=None, pad_fac=0.5):
        xmin = self.xmin
        xmax = self.xmax
        if show_ghost:
            xmin -= self.ng*self.dx
            xmax += self.ng*self.dx
        if padding:
            xmin -= pad_fac*self.dx
            xmax += pad_fac*self.dx

        plt.xlim(xmin, xmax)

        if ylim is not None:
            plt.ylim(ylim)

        plt.axis("off")
        plt.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95)


class CellCentered(object):
    """ cell-centered data defined on a 1-d finite-volume grid """

    def __init__(self, gr, a):
        if not len(a) == len(gr.xc):
            sys.exit("ERROR: grid length != data length")

        self.gr = gr
        self.a = a

    def label_data_point(self, idx, string, color="k", fontsize="large"):
        plt.text(self.gr.xc[idx], self.gr.voff+self.a[idx]+0.1, string,
                 horizontalalignment='center', verticalalignment='bottom',
                 fontsize=fontsize, color=color)

    def draw_data_point(self, idx, color="0.5", marker="o"):
        plt.scatter([self.gr.xc[idx]], [self.gr.voff+self.a[idx]],
                    color=color, marker=marker, zorder=100)


class PiecewiseConstant(object):
    """ piecewise constant data defined on a 1-d finite-volume grid """

    def __init__(self, gr, a, scale=1.0):
        if not len(a) == len(gr.xc):
            sys.exit("ERROR: grid length != data length")

        self.gr = gr
        self.a = a

        # scale is used for plotting only -- it is the normalization
        # factor for a
        if scale <= 0.0: scale = 1.0

        self.scale = scale


    def fill_zero_gradient(self):
        self.a[0:self.gr.ilo] = self.a[self.gr.ilo]
        self.a[self.gr.ihi:2*self.gr.ng+self.gr.nx] = self.a[self.gr.ihi]


    def label_cell_avg(self, idx, string, color="k"):
        plt.text(self.gr.xc[idx], self.gr.voff+self.a[idx]/self.scale+0.1, string,
                 horizontalalignment='center', verticalalignment='bottom',
                 fontsize="large", color=color)

    def draw_cell_avg(self, idx, color="0.5", ls="-"):
        plt.plot([self.gr.xl[idx], self.gr.xr[idx]],
                 [self.gr.voff+self.a[idx]/self.scale, 
                  self.gr.voff+self.a[idx]/self.scale], color=color, ls=ls)


class PiecewiseLinear(PiecewiseConstant):
    """ piecewise linear data defined on a 1-d finite-volume grid """

    def __init__(self, gr, a, nolimit=0, scale=1.0):

        PiecewiseConstant.__init__(self, gr, a, scale=scale)

        self.slope = np.zeros_like(self.a)
        self.nolimit = nolimit

        self.calculate_slopes()
 
    def calculate_slopes(self):
        # calculate the slopes
        for n in range(1, len(self.a)-1):
            test = (self.a[n+1] - self.a[n])*(self.a[n] - self.a[n-1])
            da = 0.5*(self.a[n+1] - self.a[n-1])

            if not self.nolimit:
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
             self.gr.voff+self.a[idx]

        plt.plot([self.gr.xl[idx], self.gr.xr[idx]], [yl/self.scale, yr/self.scale],
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

        plt.fill(xx, yy/self.scale, color=color, lw=1, zorder=-1)


    def slope_trace_right(self, idx, sigma, color="0.5"):

        # sigma is the fraction of the domain -- the CFL number
        x = np.linspace(self.gr.xl[idx], self.gr.xl[idx]+sigma*self.gr.dx, 50)

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

        plt.fill(xx, yy/self.scale, color=color, lw=1, zorder=-1)


    def evolve_to_right(self, idx, sigma, color="0.5", ls="-"):

        # sigma is the fraction of the domain -- the CFL number
        # show the reconstructed profile as we evolve to the right

        xm = np.linspace(self.gr.xr[idx-1]-sigma*self.gr.dx, self.gr.xr[idx-1], 50)
        am = self.a[idx-1] + (self.slope[idx-1]/self.gr.dx) * (xm - self.gr.xc[idx-1])
        xm = xm + sigma*self.gr.dx

        xp = np.linspace(self.gr.xl[idx], self.gr.xl[idx]+(1.0-sigma)*self.gr.dx, 50)
        ap = self.a[idx] + (self.slope[idx]/self.gr.dx) * (xp - self.gr.xc[idx])
        xp = xp + sigma*self.gr.dx

        plt.plot(xm, am/self.scale, color=color, lw=1, ls=ls, zorder=10)
        plt.plot(xp, ap/self.scale, color=color, lw=1, ls=ls, zorder=10)


class PiecewiseParabolic(PiecewiseConstant):
    """ piecewise parabolic data defined on a 1-d finite-volume grid """

    def __init__(self, gr, a, nolimit=0, scale=1.0):

        PiecewiseConstant.__init__(self, gr, a, scale=scale)

        # this computes the (limited) interface states just from
        # cubic interpolation through the 4 zones centered on an
        # interface
        self.aint = np.zeros_like(a)

        for n in range(1, len(a)-2):

            da0 = 0.5*(self.a[n+1] - self.a[n-1])
            dap = 0.5*(self.a[n+2] - self.a[n])

            if not nolimit:

                if (self.a[n+1] - self.a[n])*(self.a[n] - self.a[n-1]) > 0.0:
                    da0 = np.sign(da0)*min(math.fabs(da0),
                                           min(2.0*math.fabs(self.a[n] - self.a[n-1]),
                                               2.0*math.fabs(self.a[n+1] - self.a[n])) )
                else:
                    da0 = 0.0


                if (self.a[n+2] - self.a[n+1])*(self.a[n+1] - self.a[n]) > 0.0:
                    dap = np.sign(dap)*min(math.fabs(dap),
                                           min(2.0*math.fabs(self.a[n+1] - self.a[n]),
                                               2.0*math.fabs(self.a[n+2] - self.a[n+1])) )
                else:
                    dap = 0.0


            # cubic
            self.aint[n] = 0.5*(self.a[n] + self.a[n+1]) - (1.0/6.0)*(dap - da0)


        self.ap = np.zeros_like(a)
        self.am = np.zeros_like(a)
        self.a6 = np.zeros_like(a)

        # parabola of form:
        #    a(xi) = aminus + xi*(aplus - aminus + a6 * (1-xi) )
        # with  xi = (x - xl)/dx

        self.ap[:] = self.aint[:]
        self.am[1:] = self.ap[:-1]

        if not nolimit:

            for n in range(2, len(self.a)-2):

                if (self.ap[n] - self.a[n])*(self.a[n] - self.am[n]) <= 0.0:
                    self.am[n] = self.a[n]
                    self.ap[n] = self.a[n]

                elif ( (self.ap[n] - self.am[n])*(self.a[n] - 0.5*(self.am[n] + self.ap[n])) >
                       (self.ap[n] - self.am[n])**2/6.0 ):
                    self.am[n] = 3.0*self.a[n] - 2.0*self.ap[n]

                elif ( -(self.ap[n] - self.am[n])**2/6.0 >
                       (self.ap[n] - self.am[n])*(self.a[n] - 0.5*(self.am[n] + self.ap[n])) ):
                    self.ap[n] = 3.0*self.a[n] - 2.0*self.am[n]

        for n in range(2, len(self.a)-2):
            self.a6[n] = 6.0*self.a[n] - 3.0*(self.am[n] + self.ap[n])



    def draw_parabola(self, idx, color="r", ls="-"):
        x = np.linspace(self.gr.xl[idx], self.gr.xr[idx], 50)
        xi = (x - self.gr.xl[idx])/self.gr.dx

        a = self.am[idx] + xi*(self.ap[idx] - self.am[idx] + self.a6[idx] * (1.0-xi) )

        plt.plot(x, a/self.scale, color=color, ls=ls, lw=1, zorder=10)


    def ppm_trace_left(self, idx, sigma, color="0.5"):

        # sigma is the fraction of the domain from the right interface
        # (since we are tracing to create the left state there).

        x = np.linspace(self.gr.xr[idx]-sigma*self.gr.dx, self.gr.xr[idx], 50)
        xi = (x - self.gr.xl[idx])/self.gr.dx

        a = self.am[idx] + xi*(self.ap[idx] - self.am[idx] + self.a6[idx] * (1.0-xi) )

        xx = np.zeros(len(x) + 3, dtype=np.float64)
        yy = np.zeros(len(x) + 3, dtype=np.float64)

        xx[0:len(x)] = x
        xx[len(x):] = [self.gr.xr[idx],
                       self.gr.xr[idx]-sigma*self.gr.dx,
                       self.gr.xr[idx]-sigma*self.gr.dx]

        yy[0:len(x)] = a
        yy[len(x):] = [0.0, 0.0, a[0]]

        plt.fill(xx, yy/self.scale, color=color, lw=1, zorder=-1)


class Grid2d(object):
    """ the base 2-d grid """

    def __init__(self, nx, ny, ng = 0, 
                 xmin=0.0, xmax=1.0, ymin=0.0, ymax=1.0):

        # finite-volume or cell-centered finite-difference
        self.nx = nx
        self.ny = ny
        self.ng = ng
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax

        self.ilo = ng
        self.ihi = ng+nx-1
        self.jlo = ng
        self.jhi = ng+ny-1


        self.dx = (xmax - xmin)/float(nx)

        self.xl = (np.arange(2*ng+nx)-ng)*self.dx + xmin
        self.xr = (np.arange(2*ng+nx)+1-ng)*self.dx + xmin
        self.xc = 0.5*(self.xl + self.xr)

        self.dy = (ymax - ymin)/float(ny)

        self.yl = (np.arange(2*ng+ny)-ng)*self.dy + ymin
        self.yr = (np.arange(2*ng+ny)+1-ng)*self.dy + ymin
        self.yc = 0.5*(self.yl + self.yr)

    def draw_grid(self, color="k"):
        # x lines
        for n in range(self.ny):
            plt.plot([self.xmin-0.25*self.dx, self.xmax+0.25*self.dx], 
                     [self.yl[self.ng+n], self.yl[self.ng+n]], 
                     color=color, lw=2)
        
        plt.plot([self.xmin-0.25*self.dx, self.xmax+0.25*self.dx], 
                 [self.yr[self.ng+self.ny-1], self.yr[self.ng+self.ny-1]], 
                 color=color, lw=2)

        # y lines
        for n in range(self.nx):
            plt.plot([self.xl[self.ng+n], self.xl[self.ng+n]], 
                     [self.ymin-0.25*self.dy, self.ymax+0.25*self.dy], 
                     color=color, lw=2)
        
        plt.plot([self.xr[self.ng+self.nx-1], self.xr[self.ng+self.nx-1]], 
                 [self.ymin-0.25*self.dy, self.ymax+0.25*self.dy], 
                 color=color, lw=2)


    def label_center_x(self, idx, string):
        plt.text(self.xc[idx], self.yl[0]-0.35*self.dy, string,
                   horizontalalignment='center', fontsize="medium")

    def label_center_y(self, jdx, string):
        plt.text(self.xl[0]-0.35*self.dx, self.yc[jdx], string,
                   verticalalignment='center', fontsize="medium")

    def clean_axes(self):
        plt.xlim(self.xmin-0.5*self.dx, self.xmax+0.5*self.dx)
        plt.ylim(self.ymin-0.5*self.dy, self.ymax+0.5*self.dy)
        plt.axis("off")


class FVGrid2d(Grid2d):
    """ a 2-d finite-volume grid """

    def label_cell_center(self, idx, jdx, string, fontsize="medium", color="k"):
        plt.text(self.xc[idx], self.yc[jdx], 
                 string, fontsize=fontsize, color=color,
                 horizontalalignment='center', verticalalignment='center')

    def shade_cell(self, idx, jdx):
        xl = self.xl[idx]
        xr = self.xr[idx]
        yl = self.yl[jdx]
        yr = self.yr[jdx]
        plt.fill([xl, xl, xr, xr, xl], [yl, yr, yr, yl, yl], "0.75")

    def mark_cell_left_state_x(self, idx, jdx, string, color="k", 
                               fontsize="medium"):
        plt.scatter(self.xr[idx]-0.05*self.dx, self.yc[jdx], 
                      marker="x", s=50, color=color)
        plt.text(self.xr[idx]-0.075*self.dx, self.yc[jdx], string, 
                   fontsize=fontsize, rotation="270", color=color,
                   horizontalalignment='right', verticalalignment='center')

    def mark_cell_right_state_x(self, idx, jdx, string, color="k", 
                                fontsize="medium"):
        plt.scatter(self.xl[idx]+0.05*self.dx, self.yc[jdx], 
                      marker="x", s=50, color=color)
        plt.text(self.xl[idx]+0.075*self.dx, self.yc[jdx], string, 
                   fontsize=fontsize, rotation="270", color=color,
                   horizontalalignment='left', verticalalignment='center')

    def mark_cell_state_y(self, idx, jdx, string, color="k", 
                          fontsize="medium", off_sign=1.0):
        plt.scatter(self.xc[idx], self.yr[jdx], 
                      marker="x", s=50, color=color)
        if off_sign > 0:
            align = "bottom"
        else:
            align = "top"

        plt.text(self.xc[idx], self.yr[jdx]+off_sign*0.05*self.dy, string, 
                   fontsize=fontsize, rotation="0", color=color,
                   horizontalalignment='center', verticalalignment=align)

    def mark_cell_left_state_y(self, idx, jdx, string, color="k", 
                               fontsize="medium"):
        plt.scatter(self.xc[idx], self.yr[jdx]-0.05*self.dy, 
                      marker="x", s=50, color=color)
        plt.text(self.xc[idx], self.yr[jdx]-0.075*self.dy, string, 
                   fontsize=fontsize, rotation="0", color=color,
                   horizontalalignment='center', verticalalignment='top')

    def mark_cell_right_state_y(self, idx, jdx, string, color="k", 
                                fontsize="medium"):
        plt.scatter(self.xc[idx], self.yl[jdx]+0.05*self.dy, 
                      marker="x", s=50, color=color)
        plt.text(self.xc[idx], self.yl[jdx]+0.075*self.dy, string, 
                   fontsize=fontsize, rotation="0", color=color,
                   horizontalalignment='center', verticalalignment='bottom')



class FDGrid2d(Grid2d):
    """ a 2-d finite-difference grid """

    def label_cell_center(self, idx, jdx, string, fontsize="medium", color="k"):
        plt.scatter([self.xc[idx]], [self.yc[jdx]], marker="x", color=color)
        plt.text(self.xc[idx]+0.075*self.dx, self.yc[jdx]+0.075*self.dy, 
                 string, fontsize=fontsize, color=color,
                 horizontalalignment='left', verticalalignment='center')

