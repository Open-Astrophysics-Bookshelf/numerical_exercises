# 2nd-order accurate finite-volume implementation of linear advection with 
# piecewise linear slope reconstruction
# 
# We are solving a_t + u a_x = 0
#
# M. Zingale (2013-03-24)

import numpy
import pylab
import math

class ccFVgrid:

    def __init__(self, nx, ng, xmin=0.0, xmax=1.0):

        self.xmin = xmin
        self.xmax = xmax
        self.ng = ng
        self.nx = nx

        # python is zero-based.  Make easy intergers to know where the
        # real data lives
        self.ilo = ng
        self.ihi = ng+nx-1

        # physical coords -- cell-centered, left and right edges
        self.dx = (xmax - xmin)/(nx)
        self.x = xmin + (numpy.arange(nx+2*ng)-ng+0.5)*self.dx
        self.xl = xmin + (numpy.arange(nx+2*ng)-ng)*self.dx
        self.xr = xmin + (numpy.arange(nx+2*ng)-ng+1.0)*self.dx

        # storage for the solution
        self.a = numpy.zeros((nx+2*ng), dtype=numpy.float64)

    def period(self, u):
        """ return the period for advection with velocity u """
        return (self.xmax - self.xmin)/u

    def scratchArray(self):
        """ return a scratch array dimensioned for our grid """
        return numpy.zeros((self.nx+2*self.ng), dtype=numpy.float64)

    def fillBCs(self):
        """ fill all single ghostcell with periodic boundary conditions """

        # left boundary
        n = 0
        while (n < self.ng):
            self.a[self.ilo-1-n] = self.a[self.ihi-n]
            n += 1

        # right boundary
        n = 0
        while (n < self.ng):
            self.a[self.ihi+1+n] = self.a[self.ilo+n]
            n += 1

    def initCond(self, type="tophat"):

        if type == "tophat":
            self.a[numpy.logical_and(self.x >= 0.333, self.x <= 0.666)] = 1.0
        elif type == "sine":
            self.a[:] = numpy.sin(2.0*math.pi*self.x/(self.xmax-self.xmin))
        elif type == "gaussian":
            self.a[:] = 1.0 + numpy.exp(-60.0*(self.x - 0.5)**2)

        self.ainit = self.a.copy()

    def norm(self, e):
        """ return the norm of quantity e which lives on the grid """
        if not len(e) == (2*self.ng + self.nx):
            return None

        return numpy.sqrt(self.dx*numpy.sum(e[self.ilo:self.ihi+1]**2))



#-----------------------------------------------------------------------------
# advection-specific routines

def timestep(g, C, u):
    return C*g.dx/u


def states(g, dt, u, slopeType):
    """ compute the left and right interface states """

    # compute the piecewise linear slopes
    slope = g.scratchArray()


    if slopeType == "godunov":

        # piecewise constant = 0 slopes
        slope[:] = 0.0

    elif slopeType == "centered":

        # unlimited centered difference slopes

        i = g.ilo-1
        while (i <= g.ihi+1):
            slope[i] = 0.5*(g.a[i+1] - g.a[i-1])/g.dx
            i += 1

    elif slopeType == "minmod":

        # minmod limited slope

        i = g.ilo-1
        while (i <= g.ihi+1):
            slope[i] = minmod( (g.a[i] - g.a[i-1])/g.dx, 
                               (g.a[i+1] - g.a[i])/g.dx )
            i += 1
        
    elif slopeType == "MC":

        # MC limiter

        i = g.ilo-1
        while (i <= g.ihi+1):
            slope[i] = minmod(minmod( 2.0*(g.a[i] - g.a[i-1])/g.dx, 
                                      2.0*(g.a[i+1] - g.a[i])/g.dx ),
                              0.5*(g.a[i+1] - g.a[i-1])/g.dx)
            i += 1

    elif slopeType == "superbee":

        # superbee limiter

        i = g.ilo-1
        while (i <= g.ihi+1):
            A = minmod( (g.a[i+1] - g.a[i])/g.dx,
                        2.0*(g.a[i] - g.a[i-1])/g.dx )

            B = minmod( (g.a[i] - g.a[i-1])/g.dx,
                        2.0*(g.a[i+1] - g.a[i])/g.dx )
            
            slope[i] = maxmod(A, B)
            i += 1



    # loop over all the interfaces.  Here, i refers to the left
    # interface of the zone.  Note that thre are 1 more interfaces
    # than zones
    al = g.scratchArray()
    ar = g.scratchArray()

    i = g.ilo
    while (i <= g.ihi+1):

        # left state on the current interface comes from zone i-1
        al[i] = g.a[i-1] + 0.5*g.dx*(1.0 - u*dt/g.dx)*slope[i-1]

        # right state on the current interface comes from zone i
        ar[i] = g.a[i] - 0.5*g.dx*(1.0 + u*dt/g.dx)*slope[i]

        i += 1

    return al, ar


def riemann(u, al, ar):
    """ Riemann problem for advection -- this is simply upwinding,
        but we return the flux """

    if u > 0.0:
        return u*al
    else:
        return u*ar


def update(g, dt, flux):
    """ conservative update """

    anew = g.scratchArray()

    anew[g.ilo:g.ihi+1] = g.a[g.ilo:g.ihi+1] + \
        dt/g.dx * (flux[g.ilo:g.ihi+1] - flux[g.ilo+1:g.ihi+2])

    return anew


def evolve(nx, C, u, numPeriods, ICname, slopeType="centered"):

    ng = 2

    # create the grid
    g = ccFVgrid(nx, ng)

    t = 0.0
    tmax = numPeriods*g.period(u)

    # initialize the data
    g.initCond(ICname)


    # main evolution loop
    while (t < tmax):

        # fill the boundary conditions
        g.fillBCs()

        # get the timestep
        dt = timestep(g, C, u)

        if (t + dt > tmax):
            dt = tmax - t

        # get the interface states
        al, ar = states(g, dt, u, slopeType)

        # solve the Riemann problem at all interfaces
        flux = riemann(u, al, ar)
        
        # do the conservative update
        anew = update(g, dt, flux)

        g.a[:] = anew[:]

        t += dt

    return g


def minmod(a, b):
    if (abs(a) < abs(b) and a*b > 0.0):
        return a
    elif (abs(b) < abs(a) and a*b > 0.0):
        return b
    else:
        return 0.0

def maxmod(a, b):
    if (abs(a) > abs(b) and a*b > 0.0):
        return a
    elif (abs(b) > abs(a) and a*b > 0.0):
        return b
    else:
        return 0.0



#-----------------------------------------------------------------------------
# tophat

u = 1.0
nx = 64
C = 0.7

gc = evolve(nx, C, u, 5, "tophat", slopeType="centered")
gm = evolve(nx, C, u, 5, "tophat", slopeType="minmod")

pylab.plot(gc.x[gc.ilo:gc.ihi+1], gc.a[gc.ilo:gc.ihi+1], color="r", label="unlimited")
pylab.plot(gm.x[gm.ilo:gm.ihi+1], gm.a[gm.ilo:gm.ihi+1], color="b", label="minmod limiter")
pylab.plot(gc.x[gc.ilo:gc.ihi+1], gc.ainit[gc.ilo:gc.ihi+1], ls=":", color="0.5", label="exact")
pylab.legend(frameon=False, loc="best")

pylab.xlabel(r"$x$")
pylab.ylabel(r"$a$")

pylab.savefig("fv-advect.eps")


#-----------------------------------------------------------------------------
# convergence test
problem = "gaussian"
N = [32, 64, 128, 256, 512]
u = 1.0
C = 0.8

err = []

for nx in N:

    g = evolve(nx, C, u, 5, problem)

    # compute the error
    err.append(g.norm(g.a - g.ainit))
    print g.dx, nx, err[-1]


pylab.clf()

N = numpy.array(N, dtype=numpy.float64)
err = numpy.array(err)

pylab.scatter(N, err, color="r")
pylab.plot(N, err[len(N)-1]*(N[len(N)-1]/N)**2, color="k", label=r"$\mathcal{O}(\Delta x^2)$")

ax = pylab.gca()
ax.set_xscale('log')
ax.set_yscale('log')

pylab.xlim(20,600)

pylab.xlabel("number of zones", fontsize=16)
pylab.ylabel(r"$\|\| a^\mathrm{final} - a^\mathrm{init} \|\|_2$", fontsize=16)

pylab.legend(frameon=False, loc="best")

pylab.savefig("plm-converge.eps")


