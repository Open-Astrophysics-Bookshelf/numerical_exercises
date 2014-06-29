import math
import numpy
import pylab

def simplegrid():

    # grid info
    xmin = 0.0
    xmax = 1.0

    nzones = 7
    ng = 1
    
    dx = (xmax - xmin)/float(nzones)

    xl = (numpy.arange(2*ng+nzones) - ng)*dx
    xr = (numpy.arange(2*ng+nzones)+1 - ng)*dx

    xc = 0.5*(xl + xr)


    #------------------------------------------------------------------------
    # plot a domain without ghostcells
    pylab.plot([xmin,xmax], [0,0], color="k", lw=2)

    # domain left edge
    pylab.plot([xl[ng], xl[ng]], [0, 0.5], color="k", lw=4)

    n = ng
    while (n < ng+nzones):

        # draw right edge
        pylab.plot([xr[n], xr[n]], [0, 0.5], color="k", lw=2)        

        # draw center marker
        pylab.plot([xc[n], xc[n]], [-0.05, 0], color="k")                
        n += 1

    # domain right edge
    pylab.plot([xl[ng+nzones], xl[ng+nzones]], [0, 0.5], color="k", lw=4)
       

    # label a few
    pylab.text(xc[ng+nzones/2], -0.1, r"$i$", 
               horizontalalignment='center', verticalalignment='top')

    pylab.text(xc[ng+nzones/2-1], -0.1, r"$i-1$", 
               horizontalalignment='center', verticalalignment='top')

    pylab.text(xc[ng+nzones/2+1], -0.1, r"$i+1$", 
               horizontalalignment='center', verticalalignment='top')

    pylab.xlim(xl[0]-0.5*dx,xr[2*ng+nzones-1]+0.5*dx)
    pylab.ylim(-0.25, 0.6)
    pylab.axis("off")

    pylab.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95)

    f = pylab.gcf()
    f.set_size_inches(8.0,2.0)


    pylab.savefig("simplegrid.png")
    pylab.savefig("simplegrid.eps")
               


    #------------------------------------------------------------------------
    # now draw one with ghostcells
    pylab.plot([xmin-ng*dx,xmax+ng*dx], [0,0], color="k", lw=2)


    # domain (w/ ghostcells) left edge
    pylab.plot([xl[0], xl[0]], [0, 0.5], color="k", lw=2)

    n = 0
    while (n < 2*ng+nzones):

        # draw right edge
        pylab.plot([xr[n], xr[n]], [0, 0.5], color="k", lw=2)        

        # draw center marker
        pylab.plot([xc[n], xc[n]], [-0.05, 0], color="k")                
        n += 1


    # domain left edge
    pylab.plot([xl[ng], xl[ng]], [0, 0.5], color="k", lw=4)

    # domain right edge
    pylab.plot([xl[ng+nzones], xl[ng+nzones]], [0, 0.5], color="k", lw=4)
       

    # label a few
    pylab.text(xc[ng+nzones/2], -0.1, r"$i$", 
               horizontalalignment='center', verticalalignment='top')

    pylab.text(xc[ng+nzones/2-1], -0.1, r"$i-1$", 
               horizontalalignment='center', verticalalignment='top')

    pylab.text(xc[ng+nzones/2+1], -0.1, r"$i+1$", 
               horizontalalignment='center', verticalalignment='top')

    pylab.text(xc[ng-1], -0.1, r"$-1$", 
               horizontalalignment='center', verticalalignment='top')

    pylab.text(xc[ng], -0.1, r"$0$", 
               horizontalalignment='center', verticalalignment='top')

    pylab.text(xc[ng+nzones-1], -0.1, r"$N-1$", 
               horizontalalignment='center', verticalalignment='top')

    pylab.text(xc[ng+nzones], -0.1, r"$N$", 
               horizontalalignment='center', verticalalignment='top')

    pylab.xlim(xl[0]-0.5*dx,xr[2*ng+nzones-1]+0.5*dx)
    pylab.ylim(-0.25, 0.6)
    pylab.axis("off")

    pylab.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95)

    f = pylab.gcf()
    f.set_size_inches(8.0,2.0)


    pylab.savefig("simplegrid_gc.png")
    pylab.savefig("simplegrid_gc.eps")
               



if __name__== "__main__":
    simplegrid()
