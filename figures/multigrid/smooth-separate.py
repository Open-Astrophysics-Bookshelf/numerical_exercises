#!/usr/bin/env python

"""

an example of solving Poisson's equation via smoothing only.  Here, we
solve

u_xx = sin(x)
u = 0 on the boundary [0,1]

The analytic solution is u(x) = -sin(x) + x sin(1)

M. Zingale (2013-03-31)

"""
#from io import *
import numpy
import pylab
import sys

# the analytic solution
def true(x):
    return -numpy.sin(x) + x*numpy.sin(1.0)


# the L2 error norm
def error(ilo, ihi, dx, r):

    # L2 norm of elements in r, multiplied by dx to
    # normalize
    return numpy.sqrt(dx*numpy.sum((r[ilo:ihi+1]**2)))


# the righthand side
def f(x):
    return numpy.sin(x)



def computeResidual(ilo, ihi, dx, phi, frhs):

    r = numpy.zeros(len(phi))

    r[ilo:ihi+1] = frhs[ilo:ihi+1] - \
        (phi[ilo+1:ihi+2] - 2.0*phi[ilo:ihi+1] + phi[ilo-1:ihi])/dx**2


    return r


def smoothRun(nx):

    xmin = 0.0
    xmax = 1.0

    ng = 1

    print nx

    # initialize the solution to zero.  Put one ghost cell on either end
    phi = numpy.zeros(nx + 2*ng, dtype=numpy.float64)
    phinew = numpy.zeros(nx + 2*ng, dtype=numpy.float64)

    ilo = ng
    ihi = ng + nx - 1

    # coordinates of centers
    dx = (xmax - xmin)/nx
    x = (numpy.arange(nx+2*ng) - ng + 0.5)*dx + xmin

    # initialize the RHS using the function f
    frhs = numpy.zeros(nx + 2*ng, dtype=numpy.float64)
    frhs[ilo:ihi+1] = f(x[ilo:ihi+1])

    # smooth 
    n = numpy.arange(20000) + 1
    e = []
    r = []

    # fill the ghost cells
    phi[ilo-1] = -phi[ilo]
    phi[ihi+1] = -phi[ihi]

    print "source norm: ", error(ilo, ihi, dx, frhs)
    print numpy.sum(frhs[ilo:ihi+1])

    for i in n:

        # do Jacobi
        #phinew[ilo:ihi+1] = \
        #    (-dx*dx*frhs[ilo:ihi+1] + phi[ilo+1:ihi+2] + phi[ilo-1:ihi])/2.0
        # phi[:] = phinew[:]

        # this also does Jacobi, since phi on the RHS is evaluated before
        # the LHS is updated
        #phi[ilo:ihi+1] = \
        #    0.5*(-dx*dx*frhs[ilo:ihi+1] + phi[ilo+1:ihi+2] + phi[ilo-1:ihi])


        # red-black Gauss-Seidel -- first do the odd, then even points
        phi[ilo:ihi+1:2] = \
            0.5*(-dx*dx*frhs[ilo:ihi+1:2] + \
                      phi[ilo+1:ihi+2:2] + phi[ilo-1:ihi:2])

        # fill the ghost cells
        phi[ilo-1] = -phi[ilo]
        phi[ihi+1] = -phi[ihi]

        phi[ilo+1:ihi+1:2] = \
            0.5*(-dx*dx*frhs[ilo+1:ihi+1:2] + \
                      phi[ilo+2:ihi+2:2] + phi[ilo:ihi:2])


        # fill the ghost cells
        phi[ilo-1] = -phi[ilo]
        phi[ihi+1] = -phi[ihi]

        
        # compute the true error (wrt the analytic solution)
        e.append(error(ilo, ihi, dx, phi - true(x)))
        
        # compute the residual
        resid = computeResidual(ilo, ihi, dx, phi, frhs)

        r.append(error(ilo, ihi, dx, resid))


    r = numpy.array(r)
    e = numpy.array(e)

    return n, r, e


# test the multigrid solver
N = [16, 32, 64]

c = ["r", "g", "b"]

for nx in N:

    n, r, e = smoothRun(nx)
    color = c.pop()
    pylab.plot(n, e, color=color, label = `nx`)
    pylab.plot(n, r, color=color, ls=":")

ax = pylab.gca()
ax.set_xscale('log')
ax.set_yscale('log')

pylab.xlabel("# of iterations")
pylab.ylabel("L2 norm of true error (solid) and residual (dotted)")
pylab.legend(frameon=False, fontsize="small")

pylab.savefig("smooth-error.png")


