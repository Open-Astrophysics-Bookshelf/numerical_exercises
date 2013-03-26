# reconstruct - evolve - average: demonstrate what happens when we don't
# limit

import math
import numpy
import pylab
import grid_plot_util as gpu


#-----------------------------------------------------------------------------
def fillgc(gr, a):

    # simple zero-gradient
    a[0:gr.ilo] = a[gr.ilo]
    a[gr.ihi:2*gr.ng+gr.nx] = a[gr.ihi]


#-----------------------------------------------------------------------------
def evolve(gr, a, C, num, nolimit=1):
    
    #-------------------------------------------------------------------------
    # first frame -- the original cell-averages

    pylab.clf()

    gpu.drawGrid(gr)

    gpu.labelCenter(gr, gr.ng + nzones/2,   r"$i$")
    gpu.labelCenter(gr, gr.ng + nzones/2-1, r"$i-1$")
    gpu.labelCenter(gr, gr.ng + nzones/2+1, r"$i+1$")
    gpu.labelCenter(gr, gr.ng + nzones/2-2, r"$i-2$")
    gpu.labelCenter(gr, gr.ng + nzones/2+2, r"$i+2$")


    # draw cell averages
    n = gr.ilo
    while (n <= gr.ihi):
        gpu.drawCellAvg(gr, n, a[n], color="r")
        n += 1

    pylab.axis([gr.xmin-0.5*gr.dx,gr.xmax+0.5*gr.dx, -0.25, 1.2])
    pylab.axis("off")

    print gr.xmin-0.5*gr.dx, gr.xmax+0.5*gr.dx

    pylab.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95)

    f = pylab.gcf()
    f.set_size_inches(8.0,2.0)

    if (nolimit):
        pylab.savefig("rea-nolimit-start_%3.3d.png" % (num))
        pylab.savefig("rea-nolimit-start_%3.3d.eps" % (num))
    else:
        pylab.savefig("rea-start_%3.3d.png" % (num))
        pylab.savefig("rea-start_%3.3d.eps" % (num))

    #-------------------------------------------------------------------------
    # second frame -- reconstruction

    # compute the slopes
    lda = gpu.lslopes(a, nolimit=nolimit)

    # draw
    pylab.clf()

    gpu.drawGrid(gr)

    gpu.labelCenter(gr, gr.ng + nzones/2,   r"$i$")
    gpu.labelCenter(gr, gr.ng + nzones/2-1, r"$i-1$")
    gpu.labelCenter(gr, gr.ng + nzones/2+1, r"$i+1$")
    gpu.labelCenter(gr, gr.ng + nzones/2-2, r"$i-2$")
    gpu.labelCenter(gr, gr.ng + nzones/2+2, r"$i+2$")


    # draw cell averages
    n = gr.ilo
    while (n <= gr.ihi):
        gpu.drawCellAvg(gr, n, a[n], color="0.5", ls=":")
        n += 1

    n = gr.ilo
    while (n <= gr.ihi):
        gpu.drawSlope(gr, n, lda[n], a[n], color="r")
        n += 1

    pylab.axis([gr.xmin-0.5*gr.dx,gr.xmax+0.5*gr.dx, -0.25, 1.2])
    pylab.axis("off")


    pylab.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95)

    f = pylab.gcf()
    f.set_size_inches(8.0,2.0)

    if (nolimit):
        pylab.savefig("rea-nolimit-reconstruction_%3.3d.png" % (num))
        pylab.savefig("rea-nolimit-reconstruction_%3.3d.eps" % (num))
    else:
        pylab.savefig("rea-reconstruction_%3.3d.png" % (num))
        pylab.savefig("rea-reconstruction_%3.3d.eps" % (num))


    #-------------------------------------------------------------------------
    # third frame -- evolve

    # draw

    pylab.clf()

    gpu.drawGrid(gr)

    gpu.labelCenter(gr, gr.ng + nzones/2,   r"$i$")
    gpu.labelCenter(gr, gr.ng + nzones/2-1, r"$i-1$")
    gpu.labelCenter(gr, gr.ng + nzones/2+1, r"$i+1$")
    gpu.labelCenter(gr, gr.ng + nzones/2-2, r"$i-2$")
    gpu.labelCenter(gr, gr.ng + nzones/2+2, r"$i+2$")


    # draw cell slopes
    n = gr.ilo
    while (n <= gr.ihi):
        gpu.drawSlope(gr, n, lda[n], a[n], color="0.75", ls=":")
        n += 1

    # evolve
    n = gr.ilo
    while (n <= gr.ihi):
        gpu.evolveToRight(gr, n, lda, a, C, color="r")
        n += 1


    pylab.axis([gr.xmin-0.5*gr.dx,gr.xmax+0.5*gr.dx, -0.25, 1.2])
    pylab.axis("off")

    print gr.xmin-0.5*gr.dx, gr.xmax+0.5*gr.dx
    
    pylab.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95)

    f = pylab.gcf()
    f.set_size_inches(8.0,2.0)

    if (nolimit):
        pylab.savefig("rea-nolimit-evolve_%3.3d.png" % (num))
        pylab.savefig("rea-nolimit-evolve_%3.3d.eps" % (num))
    else:
        pylab.savefig("rea-evolve_%3.3d.png" % (num))
        pylab.savefig("rea-evolve_%3.3d.eps" % (num))


    #-------------------------------------------------------------------------
    # fourth frame -- re-average

    # left states (we don't need the right state when u > 0)
    al = numpy.zeros(2*gr.ng + gr.nx, dtype=numpy.float64)

    n = gr.ilo
    while (n <= gr.ihi+1):
        al[n] = a[n-1] + 0.5*(1 - C)*lda[n-1]
        n += 1


    # the Riemann problem just picks the left state.  Do a conservative
    # update
    anew = numpy.zeros(2*gr.ng + gr.nx, dtype=numpy.float64)

    anew[gr.ilo:gr.ihi+1] = a[gr.ilo:gr.ihi+1] + \
        C*(al[gr.ilo:gr.ihi+1] - al[gr.ilo+1:gr.ihi+2])


    pylab.clf()

    gpu.drawGrid(gr)

    gpu.labelCenter(gr, gr.ng + nzones/2,   r"$i$")
    gpu.labelCenter(gr, gr.ng + nzones/2-1, r"$i-1$")
    gpu.labelCenter(gr, gr.ng + nzones/2+1, r"$i+1$")
    gpu.labelCenter(gr, gr.ng + nzones/2-2, r"$i-2$")
    gpu.labelCenter(gr, gr.ng + nzones/2+2, r"$i+2$")


    # show the evolved profiles from the old time
    n = gr.ilo
    while (n <= gr.ihi):
        gpu.evolveToRight(gr, n, lda, a, C, color="0.5", ls=":")
        n += 1

    # draw new averages
    n = gr.ilo
    while (n <= gr.ihi):
        gpu.drawCellAvg(gr, n, anew[n], color="red")
        n += 1

    pylab.axis([gr.xmin-0.5*gr.dx,gr.xmax+0.5*gr.dx, -0.25, 1.2])
    pylab.axis("off")

    print gr.xmin-0.5*gr.dx, gr.xmax+0.5*gr.dx

    pylab.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95)

    f = pylab.gcf()
    f.set_size_inches(8.0,2.0)

    if (nolimit):
        pylab.savefig("rea-nolimit-final_%3.3d.png" % (num))
        pylab.savefig("rea-nolimit-final_%3.3d.eps" % (num))
    else:
        pylab.savefig("rea-final_%3.3d.png" % (num))
        pylab.savefig("rea-final_%3.3d.eps" % (num))

    return anew



#-----------------------------------------------------------------------------
ainit = numpy.array([1.0, 1.0, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1])
nzones = len(ainit)

# CFL number
C = 0.9

gr = gpu.grid(nzones, ng=4)

a = numpy.zeros(2*gr.ng + gr.nx, dtype=numpy.float64)
a[gr.ilo:gr.ihi+1] = ainit[:]


# loop
for i in range(1,9):

    fillgc(gr, a)

    anew = evolve(gr, a, C, i, nolimit=1)

    a[:] = anew[:]

