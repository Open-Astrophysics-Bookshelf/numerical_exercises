import math
import numpy
import pylab

def simplegrid():

    # grid info
    xmin = 0.0
    xmax = 1.0

    nzones = 2
    ng = 0
    
    dx = (xmax - xmin)/float(nzones)

    xl = (numpy.arange(2*ng+nzones) - ng)*dx
    xr = (numpy.arange(2*ng+nzones)+1 - ng)*dx

    xc = 0.5*(xl + xr)


    #------------------------------------------------------------------------
    # plot a domain without ghostcells
    pylab.plot([xmin,xmax], [0,0], lw=2, color="0.5")

    # draw line separating zones
    pylab.plot([xr[0], xr[0]], [0, 0.5], lw=2, color="0.5")        
       

    # label a few
    pylab.text(xc[ng+nzones/2], -0.1, r"$i$", 
               horizontalalignment='center', verticalalignment='top')

    pylab.text(xc[ng+nzones/2-1], -0.1, r"$i-1$", 
               horizontalalignment='center', verticalalignment='top')


    # draw waves
    # u - c
    pylab.plot([xr[0], xr[0]-0.75*dx], [0,0.6], color="k", ls="--", lw=2)
    pylab.text(xr[0]-0.75*dx, 0.6+0.05, "$\lambda = u - c$", 
               horizontalalignment="center")

    # u
    pylab.plot([xr[0], xr[0]-0.2*dx], [0,0.6], color="k", ls="--", lw=2)
    pylab.text(xr[0]-0.2*dx, 0.6+0.05, "$\lambda = u$", 
               horizontalalignment="center")

    # u + c
    pylab.plot([xr[0], xr[0]+0.4*dx], [0,0.6], color="k", ls="--", lw=2)    
    pylab.text(xr[0]+0.4*dx, 0.6+0.05, "$\lambda = u + c$", 
               horizontalalignment="center")

    # label regions
    pylab.text(xr[0]-0.5*dx, 0.2, r"$L$")
    pylab.text(xr[0]-0.33*dx, 0.4, r"$L^*$")
    pylab.text(xr[0]+0.05*dx, 0.4, r"$R^*$")
    pylab.text(xr[0]+0.3*dx, 0.2, r"$R$")
    


    pylab.xlim(xl[0]-0.5*dx,xr[2*ng+nzones-1]+0.5*dx)
    pylab.ylim(-0.15, 0.75)
    pylab.axis("off")

    pylab.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95)

    f = pylab.gcf()
    f.set_size_inches(5,2.5)


    pylab.savefig("riemann-waves.png")
    pylab.savefig("riemann-waves.eps", bbox_inches="tight")
               





if __name__== "__main__":
    simplegrid()
