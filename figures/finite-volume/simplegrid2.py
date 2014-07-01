import math
import numpy
import pylab
import grid_plot_util as gpu

def simplegrid():

    nzones = 7

    gr = gpu.grid(nzones, xmin=0, xmax=1)

    gpu.drawGrid(gr, edgeTicks=0)

    # label a few cell-centers
    gpu.labelCenter(gr, nzones/2, r"$i$")
    gpu.labelCenter(gr, nzones/2-1, r"$i-1$")
    gpu.labelCenter(gr, nzones/2+1, r"$i+1$")

    # label a few edges
    gpu.labelEdge(gr, nzones/2, r"$i-1/2$")
    gpu.labelEdge(gr, nzones/2+1, r"$i+1/2$")


    # draw an average quantity
    gpu.drawCellAvg(gr, nzones/2, 0.4, color="r")
    gpu.labelCellAvg(gr, nzones/2, 0.4, r"$\,\langle a \rangle_i$", color="r")

    pylab.axis([gr.xmin-1.5*gr.dx,gr.xmax+1.5*gr.dx, -0.25, 1.5])
    pylab.axis("off")

    pylab.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95)

    f = pylab.gcf()
    f.set_size_inches(10.0,2.5)


    pylab.savefig("simplegrid2.png")
    pylab.savefig("simplegrid2.eps")
               


if __name__== "__main__":
    simplegrid()
