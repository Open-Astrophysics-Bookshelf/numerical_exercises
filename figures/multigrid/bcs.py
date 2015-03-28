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

    pylab.text(xc[ng-1], -0.1, r"$\mathrm{lo}-1$", 
               horizontalalignment='center', verticalalignment='top')

    pylab.text(xc[ng], -0.1, r"$\mathrm{lo}$", 
               horizontalalignment='center', verticalalignment='top')

    pylab.text(xc[ng+nzones-1], -0.1, r"$\mathrm{hi}$", 
               horizontalalignment='center', verticalalignment='top')

    pylab.text(xc[ng+nzones], -0.1, r"$\mathrm{hi+1}$", 
               horizontalalignment='center', verticalalignment='top')


    # label boundaries
    pylab.text(xl[ng], -0.2, r"$\mathrm{left}$",
               horizontalalignment='center', verticalalignment='top')

    pylab.text(xl[ng], -0.3, r"$\mathrm{BC}$",
               horizontalalignment='center', verticalalignment='top')

    pylab.text(xl[ng+nzones], -0.2, r"$\mathrm{right}$",
               horizontalalignment='center', verticalalignment='top')

    pylab.text(xl[ng+nzones], -0.3, r"$\mathrm{BC}$",
               horizontalalignment='center', verticalalignment='top')
               

    pylab.plot([xl[ng], xl[ng]], [-0.19, 0], 
               linestyle=":", color="0.5",
               zorder=-1)

    pylab.plot([xl[ng+nzones], xl[ng+nzones]], [-0.19, 0], 
               linestyle=":", color="0.5",
               zorder=-1)

    pylab.xlim(xl[0]-0.5*dx,xr[2*ng+nzones-1]+0.5*dx)
    pylab.ylim(-0.4, 0.6)
    pylab.axis("off")


    pylab.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95)

    f = pylab.gcf()
    f.set_size_inches(8.0,2.0)


    pylab.savefig("mg-bcs.png")
    pylab.savefig("mg-bcs.eps")
               



if __name__== "__main__":
    simplegrid()
