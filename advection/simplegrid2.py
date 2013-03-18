import math
import numpy
import pylab

def simplegrid():

    xmin = 0.0
    xmax = 1.0

    nzones = 7

    dx = (xmax - xmin)/float(nzones)

    xl = numpy.arange(nzones)*dx
    xr = (numpy.arange(nzones)+1)*dx

    xc = 0.5*(xl + xr)

    pylab.plot([xmin-0.5*dx,xmax+0.5*dx], [0,0], color="k", lw=2)

    pylab.plot([xl[0], xl[0]], [0, 0.5], color="k", lw=2)
    n = 0
    while (n < nzones):

        # draw right edge
        pylab.plot([xr[n], xr[n]], [0, 0.5], color="k", lw=2)        

        # draw center marker
        pylab.plot([xc[n], xc[n]], [-0.05, 0], color="k")                

        # draw edge marker
        if (n == 0):
            pylab.plot([xl[0], xl[0]], [-0.05, 0], color="k")                

        pylab.plot([xr[n], xr[n]], [-0.05, 0], color="k")                

        n += 1


    # label a few cell-centers
    pylab.text(xc[nzones/2], -0.1, r"$i$", 
               horizontalalignment='center', verticalalignment='top', 
               fontsize="small")

    pylab.text(xc[nzones/2-1], -0.1, r"$i-1$", 
               horizontalalignment='center', verticalalignment='top', 
               fontsize="small")

    pylab.text(xc[nzones/2+1], -0.1, r"$i+1$", 
               horizontalalignment='center', verticalalignment='top',
               fontsize="small")


    # label a few edges
    pylab.text(xl[nzones/2], -0.075, r"$i-1/2$", 
               horizontalalignment='center', verticalalignment='top', 
               fontsize="small")

    pylab.text(xr[nzones/2], -0.075, r"$i+1/2$", 
               horizontalalignment='center', verticalalignment='top', 
               fontsize="small")



    pylab.axis([xmin-0.5*dx,xmax+0.5*dx, -0.25, 0.6])
    pylab.axis("off")

    pylab.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95)

    f = pylab.gcf()
    f.set_size_inches(8.0,2.0)


    pylab.savefig("simplegrid2.png")
    pylab.savefig("simplegrid2.eps")
               


if __name__== "__main__":
    simplegrid()
