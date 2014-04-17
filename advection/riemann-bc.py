import math
import numpy
import pylab
import grid_plot_util as gpu

def riemann():

    # grid info
    xmin = 0.0
    xmax = 1.0

    nzones = 4
    ng = 2

    gr = gpu.grid(nzones, ng=ng)

    # interior
    atemp = numpy.array([0.8, 0.7, 0.4, 0.5])

    a = numpy.zeros(2*gr.ng + gr.nx, dtype=numpy.float64)

    # fill interior and ghost cells
    a[gr.ilo:gr.ihi+1] = atemp[:]
    a[0:gr.ilo] = a[gr.ihi-1:gr.ihi+1]
    a[gr.ihi:2*gr.ng+gr.nx] = a[gr.ihi]



    #------------------------------------------------------------------------
    # plot a domain without ghostcells
    gpu.drawGrid(gr, emphasizeEnd=1, drawGhost=1)

    gpu.labelCenter(gr, gr.ng-2, r"$\mathrm{lo-2}$")
    gpu.labelCenter(gr, gr.ng-1, r"$\mathrm{lo-1}$")
    gpu.labelCenter(gr, gr.ng, r"$\mathrm{lo}$")
    gpu.labelCenter(gr, gr.ng+1, r"$\mathrm{lo+1}$")

    gpu.labelEdge(gr, gr.ng, r"$\mathrm{lo}-1/2$")
    
    # draw cell averages
    n = 0
    while n < gr.ng+gr.nx:
        gpu.drawCellAvg(gr, n, a[n], color="0.5", ls=":")
        n += 1

    # get slopes
    lda = gpu.lslopes(a, nolimit=1)

    n = gr.ilo-1
    while (n <= gr.ihi):
        gpu.drawSlope(gr, n, lda[n], a[n], color="r")
        n += 1

    # compute the states to the left and right of lo-1/2
    C = 0.7 # CFL
    al = a[gr.ilo-1] + 0.5*gr.dx*(1.0 - C)*lda[gr.ilo-1]
    ar = a[gr.ilo] - 0.5*gr.dx*(1.0 + C)*lda[gr.ilo]

    # L
    gpu.markCellRightState(gr, ng-1, r"$a_{\mathrm{lo}+1/2,L}^{n+1/2}$", 
                           value=al, vertical="top", color="b")

    # R
    gpu.markCellLeftState(gr, ng, r"$a_{\mathrm{lo}+1/2,R}^{n+1/2}$", 
                          value=ar, vertical="top", color="b")



    pylab.xlim(gr.xl[0]-0.15*gr.dx,gr.xr[ng+1]+0.15*gr.dx)
    pylab.ylim(-0.25, 1.1)
    pylab.axis("off")

    pylab.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95)

    f = pylab.gcf()
    f.set_size_inches(8.0,2.0)

    pylab.tight_layout()

    pylab.savefig("riemann-bc.png")
    pylab.savefig("riemann-bc.eps")


if __name__== "__main__":
    riemann()
