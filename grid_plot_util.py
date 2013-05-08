import math
import numpy
import pylab
import sys

class grid:

    def __init__(self, nx, ng = 0, xmin=0.0, xmax=1.0, fd=0, voff=0.0):

        # grid stuff

        if not fd:
            # finite-volume
            self.nx = nx
            self.ng = ng
            self.xmin = xmin
            self.xmax = xmax
        
            self.ilo = ng
            self.ihi = ng+nx-1

            self.dx = (xmax - xmin)/float(nx)
        
            self.xl = (numpy.arange(2*ng+nx)-ng)*self.dx + xmin
            self.xr = (numpy.arange(2*ng+nx)+1-ng)*self.dx + xmin
            self.xc = 0.5*(self.xl + self.xr)
        
        else:
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
            self.xc = xmin + (numpy.arange(2*ng+nx)-ng)*self.dx


        # if this a finite-difference (node-centered) grid?
        self.fd = fd

        # vertical offset (if we want to stack grids)
        self.voff = voff


def drawGrid(gr, centerOnly=0, drawGhost=0, emphasizeEnd=0, edgeTicks=1, color="k"):

    if centerOnly and emphasizeEnd:
        sys.exit("centerOnly and emphasizeEnd are incompatible")

    if gr.fd:
        gridTop = 0.25

        if not drawGhost:
            nstart = gr.ilo
            nstop = gr.ihi
        else:
            print "grad ghost"
            nstart = gr.ilo-gr.ng
            nstop = gr.ihi+gr.ng

        print nstart, nstop

        if not drawGhost:
            if (emphasizeEnd):
                pylab.plot([gr.xmin, gr.xmax], 
                           [0,0], color=color, lw=2)
            else:
                pylab.plot([gr.xmin-0.5*gr.dx, gr.xmax+0.5*gr.dx], 
                           [gr.voff,gr.voff], color=color, lw=2)
        else:
            pylab.plot([gr.xmin-gr.ng*gr.dx, gr.xmin],
                       [gr.voff,gr.voff], color=color, lw=2, ls=":")
            pylab.plot([gr.xmax, gr.xmax+gr.ng*gr.dx], 
                       [gr.voff,gr.voff], color=color, lw=2, ls=":")
            pylab.plot([gr.xmin, gr.xmax], 
                       [gr.voff,gr.voff], color=color, lw=2)

        n = nstart
        while (n <= nstop):

            # draw center (node) indicator line
            if (n < gr.ilo or n > gr.ihi):
                pylab.plot([gr.xc[n], gr.xc[n]], 
                           [-0.05+gr.voff, gridTop+gr.voff], color=color, ls=":", lw=2)
            else:
                pylab.plot([gr.xc[n], gr.xc[n]], 
                           [-0.05+gr.voff, gridTop+gr.voff], color=color, lw=2)
      
            n += 1

        if (emphasizeEnd):
            pylab.plot([gr.xc[gr.ilo], gr.xc[gr.ilo]], 
                       [-0.05+gr.voff, gridTop+gr.voff], color=color, lw=4)

            pylab.plot([gr.xc[gr.ihi], gr.xc[gr.ihi]], 
                       [-0.05+gr.voff, gridTop+gr.voff], color=color, lw=4)



    else:
        # finite volume grid
        gridTop = 1.0

        if not drawGhost:
            if (centerOnly == 1):
                nstart = gr.ng + gr.nx/2-1
                nstop = gr.ng + gr.nx/2
            else:
                nstart = gr.ilo
                nstop = gr.ihi
        else:
            print "grad ghost"
            nstart = gr.ilo -gr.ng
            nstop = gr.ihi + gr.ng

        print nstart, nstop

        if (emphasizeEnd):
            # horizontal line
            print "drawing line", nstart, nstop

            pylab.plot([gr.xl[nstart], gr.xr[nstop]], 
                       [gr.voff,gr.voff], color=color, lw=2)

        else:
            # horizontal line
            pylab.plot([gr.xl[nstart]-0.5*gr.dx, gr.xr[nstop]+0.5*gr.dx], 
                       [gr.voff,gr.voff], color=color, lw=2)

        # draw first left edge
        pylab.plot([gr.xl[nstart], gr.xl[nstart]], 
                   [gr.voff, gridTop+gr.voff], color=color, lw=2)


        n = nstart
        while (n <= nstop):

            # emphasize?
            if (emphasizeEnd and n == gr.ilo):
                pylab.plot([gr.xl[n], gr.xl[n]], 
                           [gr.voff, gridTop+gr.voff], color=color, lw=4)                


            # draw right edge
            if (emphasizeEnd and n == gr.ihi):
                pylab.plot([gr.xr[n], gr.xr[n]], [gr.voff, gridTop+gr.voff], 
                           color=color, lw=4)        
            else:
                pylab.plot([gr.xr[n], gr.xr[n]], [gr.voff, gridTop+gr.voff], 
                           color=color, lw=2)        

            # draw center marker
            pylab.plot([gr.xc[n], gr.xc[n]], [-0.05+gr.voff, gr.voff], color=color)      

            # draw edge marker
            if (n == nstart and edgeTicks):
                pylab.plot([gr.xl[nstart], gr.xl[nstart]], [-0.05+gr.voff, gr.voff],
                           color=color)

            if edgeTicks:
                pylab.plot([gr.xr[n], gr.xr[n]], [-0.05+gr.voff, gr.voff], color=color)
                
            n += 1


def labelCenter(gr, idx, string):

    pylab.text(gr.xc[idx], gr.voff-0.1, string, 
               horizontalalignment='center', verticalalignment='top', 
               fontsize="small")


def labelEdge(gr, idx, string):

    pylab.text(gr.xl[idx], gr.voff-0.075, string,
               horizontalalignment='center', verticalalignment='top', 
               fontsize="small")


def labelCellAvg(gr, idx, value, string, color="k"):

    pylab.text(gr.xc[idx], gr.voff+value+0.1, string,
               horizontalalignment='center', verticalalignment='bottom',
               fontsize="large", color=color)

def labelCellCenter(gr, idx, string):

    pylab.text(gr.xc[idx], gr.voff+0.5, string,
               horizontalalignment='center', verticalalignment='center',
               fontsize="large")


def labelFD(gr, idx, value, string, color="k"):

    pylab.text(gr.xc[idx], gr.voff+value+0.1, string,
               horizontalalignment='center', verticalalignment='bottom',
               fontsize="large", color=color)


#-----------------------------------------------------------------------------
def markCellLeftState(gr, idx, string, color="k"):

    pylab.scatter(gr.xl[idx]+0.05*gr.dx, gr.voff+0.5, marker="x", color=color)

    pylab.text(gr.xl[idx]+0.075*gr.dx, gr.voff+0.5, string,
               horizontalalignment='left', verticalalignment='center', color=color)


def markCellRightState(gr, idx, string, color="k"):

    pylab.scatter(gr.xr[idx]-0.05*gr.dx, gr.voff+0.5, marker="x", color=color)

    pylab.text(gr.xr[idx]-0.075*gr.dx, gr.voff+0.5, string,
               horizontalalignment='right', verticalalignment='center', color=color)



#-----------------------------------------------------------------------------
def drawFDData(gr, idx, value, color="0.5", marker="o"):
    pylab.scatter([gr.xc[idx]], [gr.voff+value], color=color, marker=marker, zorder=100)



#-----------------------------------------------------------------------------
def drawCellAvg(gr, idx, value, color="0.5", ls="-"):
    print "drawing average: ", idx
    pylab.plot([gr.xl[idx], gr.xr[idx]], [gr.voff+value, gr.voff+value], color=color, ls=ls)



#-----------------------------------------------------------------------------
def lslopes(a, nolimit=0):

    lda = numpy.zeros(len(a), dtype=numpy.float64)

    n = 1
    while (n < len(a)-1):
        test = (a[n+1] - a[n])*(a[n] - a[n-1])
        da = 0.5*(a[n+1] - a[n-1])

        if (not nolimit):
            if (test > 0.0):
                lda[n] = min(math.fabs(da), min(2.0*math.fabs(a[n+1] - a[n]),
                                                2.0*math.fabs(a[n] - a[n-1]))) * \
                                                numpy.sign(a[n+1] - a[n-1])

            else:
                lda[n] = 0.0

        else:
            lda[n] = da

        n += 1

    return lda


def drawSlope(gr, idx, slope, value, color="r", ls="-"):

    yl = slope*(gr.xl[idx] - gr.xc[idx])/gr.dx + gr.voff+value
    yr = slope*(gr.xr[idx] - gr.xc[idx])/gr.dx + gr.voff+value

    pylab.plot([gr.xl[idx], gr.xr[idx]], [yl, yr], 
               color=color, ls=ls, lw=1, zorder=10)
    

def slopeTraceLeft(gr, idx, slope, value, sigma, color="0.5"):

    # sigma is the fraction of the domain -- the CFL number
    x = numpy.linspace(gr.xr[idx]-sigma*gr.dx, gr.xr[idx], 50)

    a = gr.voff+value + (slope/gr.dx) * (x - gr.xc[idx])

    xx = numpy.zeros(len(x) + 3, dtype=numpy.float64)
    yy = numpy.zeros(len(x) + 3, dtype=numpy.float64)

    xx[0:len(x)] = x
    xx[len(x):] = [gr.xr[idx], gr.xr[idx]-sigma*gr.dx, gr.xr[idx]-sigma*gr.dx]

    yy[0:len(x)] = a
    yy[len(x):] = [gr.voff, gr.voff, a[0]]

    pylab.fill(xx, yy, color=color, lw=1, zorder=-1)


def slopeTraceRight(gr, idx, slope, value, sigma, color="0.5"):

    # sigma is the fraction of the domain -- the CFL number
    x = numpy.linspace(gr.xl[idx], gr.xl[idx]+sigma*gr.dx, 50)

    a = gr.voff+value + (slope/gr.dx) * (x - gr.xc[idx])

    xx = numpy.zeros(len(x) + 3, dtype=numpy.float64)
    yy = numpy.zeros(len(x) + 3, dtype=numpy.float64)

    xx[0:len(x)] = x
    xx[len(x):] = [gr.xl[idx]+sigma*gr.dx, gr.xl[idx], gr.xl[idx]]

    yy[0:len(x)] = a
    yy[len(x):] = [gr.voff, gr.voff, a[0]]

    pylab.fill(xx, yy, color=color, lw=1, zorder=-1)


def evolveToRight(gr, idx, slopes, values, sigma, color="0.5", ls="-"):

    # sigma is the fraction of the domain -- the CFL number
    # show the reconstructed profile as we evolve to the right

    xm = numpy.linspace(gr.xr[idx-1]-sigma*gr.dx, gr.xr[idx-1], 50)
    am = values[idx-1] + (slopes[idx-1]/gr.dx) * (xm - gr.xc[idx-1])
    xm = xm + sigma*gr.dx

    xp = numpy.linspace(gr.xl[idx], gr.xl[idx]+(1.0-sigma)*gr.dx, 50)
    ap = values[idx] + (slopes[idx]/gr.dx) * (xp - gr.xc[idx])
    xp = xp + sigma*gr.dx

    pylab.plot(xm, am, color=color, lw=1, ls=ls, zorder=10)
    pylab.plot(xp, ap, color=color, lw=1, ls=ls, zorder=10)



#-----------------------------------------------------------------------------
def ppm(a, nolimit=0):

    ap = numpy.zeros(len(a), dtype=numpy.float64)
    am = numpy.zeros(len(a), dtype=numpy.float64)
    a6 = numpy.zeros(len(a), dtype=numpy.float64)

    # parabola of form: 
    #    a(xi) = aminus + xi*(aplus - aminus + a6 * (1-xi) )  
    # with  xi = (x - xl)/dx

    n = 1
    while (n < len(a)-2):

        da0 = 0.5*(a[n+1] - a[n-1])
        dap = 0.5*(a[n+2] - a[n])

        if (not nolimit):

            if ( (a[n+1] - a[n])*(a[n] - a[n-1]) > 0.0):
                da0 = numpy.sign(da0)*min(math.fabs(da0),
                                          min(2.0*math.fabs(a[n] - a[n-1]),
                                              2.0*math.fabs(a[n+1] - a[n])) )
            else:
                da0 = 0.0


            if ( (a[n+2] - a[n+1])*(a[n+1] - a[n]) > 0.0):
                dap = numpy.sign(dap)*min(math.fabs(dap),
                                          min(2.0*math.fabs(a[n+1] - a[n]),
                                              2.0*math.fabs(a[n+2] - a[n+1])) )
            else:
                dap = 0.0


        # cubic
        ap[n] = 0.5*(a[n] + a[n+1]) - (1.0/6.0)*(dap - da0)
        am[n+1] = ap[n]

        n += 1

    if (not nolimit):

        n = 2
        while (n < len(a)-2):

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
    while (n < len(a)-2):
        a6[n] = 6.0*a[n] - 3.0*(am[n] + ap[n])
        n += 1

    return ap, am, a6


def drawParabola(gr, idx, ap, am, a6, color="r", ls="-"):

    x = numpy.linspace(gr.xl[idx], gr.xr[idx], 50)
    xi = (x - gr.xl[idx])/gr.dx
    print xi
    a = am + xi*(ap - am + a6 * (1.0-xi) )  

    pylab.plot(x, a,
               color=color, ls=ls, lw=1, zorder=10)
    

def ppmTraceLeft(gr, idx, ap, am, a6, sigma, color="0.5"):

    # sigma is the fraction of the domain from the right interface
    # (since we are tracing to create the left state there).

    x = numpy.linspace(gr.xr[idx]-sigma*gr.dx, gr.xr[idx], 50)
    xi = (x - gr.xl[idx])/gr.dx

    a = am + xi*(ap - am + a6 * (1.0-xi) )  

    xx = numpy.zeros(len(x) + 3, dtype=numpy.float64)
    yy = numpy.zeros(len(x) + 3, dtype=numpy.float64)

    xx[0:len(x)] = x
    xx[len(x):] = [gr.xr[idx], gr.xr[idx]-sigma*gr.dx, gr.xr[idx]-sigma*gr.dx]

    yy[0:len(x)] = a
    yy[len(x):] = [0.0, 0.0, a[0]]

    pylab.fill(xx, yy, color=color, lw=1, zorder=-1)

