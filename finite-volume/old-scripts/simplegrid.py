import math
import numpy
import pylab
import grid_plot_util as gpu

def simplegrid():

    # grid info
    nzones = 7
    ng = 1

    gr = gpu.grid(nzones, ng, xmin=0.0, xmax=1.0)


    #------------------------------------------------------------------------
    gpu.drawGrid(gr, emphasizeEnd=1, edgeTicks=0)


    # label a few
    gpu.labelCenter(gr, ng+nzones/2, r"$i$")
    gpu.labelCenter(gr, ng+nzones/2-1, r"$i-1$")
    gpu.labelCenter(gr, ng+nzones/2+1, r"$i+1$")


    pylab.xlim(gr.xl[0]-0.5*gr.dx,gr.xr[2*ng+nzones-1]+0.5*gr.dx)

    # label dx
    pylab.plot([gr.xl[ng+nzones/2-2], gr.xl[ng+nzones/2-2]], 
               [-0.35,-0.25], color="k")

    pylab.plot([gr.xr[ng+nzones/2-2], gr.xr[ng+nzones/2-2]], 
               [-0.35,-0.25], color="k")

    pylab.plot([gr.xl[ng+nzones/2-2], gr.xr[ng+nzones/2-2]], 
               [-0.3,-0.3], color="k")

    pylab.text(gr.xc[ng+nzones/2-2], -0.5, r"$\Delta x$", 
               horizontalalignment="center")

    pylab.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95)

    pylab.ylim(-0.5, 1.5)
    pylab.axis("off")


    f = pylab.gcf()
    f.set_size_inches(10.0,2.5)


    pylab.savefig("simplegrid.png")
    pylab.savefig("simplegrid.eps")
               
if __name__== "__main__":
    simplegrid()
