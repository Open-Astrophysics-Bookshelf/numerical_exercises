import math
import numpy
import pylab
import grid_plot_util as gpu


#-----------------------------------------------------------------------------

nzones = 8

# data that lives on the grid
#a = numpy.array([0.3, 1.0, 0.9, 0.8, 0.25, 0.15, 0.5, 0.55])
a = numpy.array([0.3, 1.0, 0.9, 0.8, 0.25, 0.1, 0.5, 0.55])

gr = gpu.grid(nzones)


pylab.clf()

gpu.drawGrid(gr)

gpu.labelCenter(gr, nzones/2,   r"$i$")
gpu.labelCenter(gr, nzones/2-1, r"$i-1$")
gpu.labelCenter(gr, nzones/2+1, r"$i+1$")
gpu.labelCenter(gr, nzones/2-2, r"$i-2$")
gpu.labelCenter(gr, nzones/2+2, r"$i+2$")


n = 0
while (n < nzones):
    gpu.drawCellAvg(gr, n, a[n], color="r")
    n += 1


pylab.axis([gr.xmin-0.5*gr.dx,gr.xmax+0.5*gr.dx, -0.25, 1.2])
pylab.axis("off")

pylab.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95)

f = pylab.gcf()
f.set_size_inches(8.0,2.0)

pylab.savefig("piecewise-constant.eps")
pylab.savefig("piecewise-constant.png")


#------------- PLM -------------
pylab.clf()

gpu.drawGrid(gr)

gpu.labelCenter(gr, nzones/2,   r"$i$")
gpu.labelCenter(gr, nzones/2-1, r"$i-1$")
gpu.labelCenter(gr, nzones/2+1, r"$i+1$")
gpu.labelCenter(gr, nzones/2-2, r"$i-2$")
gpu.labelCenter(gr, nzones/2+2, r"$i+2$")


n = 0
while (n < nzones):
    gpu.drawCellAvg(gr, n, a[n], color="0.5")
    n += 1

# compute the slopes
da = gpu.lslopes(a, nolimit=1)
lda = gpu.lslopes(a)

n = 2
while (n < nzones-2):
    gpu.drawSlope(gr, n, da[n], a[n], color="r", ls=":")
    gpu.drawSlope(gr, n, lda[n], a[n], color="r")
    n += 1


pylab.axis([gr.xmin-0.5*gr.dx,gr.xmax+0.5*gr.dx, -0.25, 1.2])
pylab.axis("off")

pylab.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95)

f = pylab.gcf()
f.set_size_inches(8.0,2.0)

pylab.savefig("piecewise-linear.eps")
pylab.savefig("piecewise-linear.png")


#------------- PPM -------------
pylab.clf()

gpu.drawGrid(gr)

gpu.labelCenter(gr, nzones/2,   r"$i$")
gpu.labelCenter(gr, nzones/2-1, r"$i-1$")
gpu.labelCenter(gr, nzones/2+1, r"$i+1$")
gpu.labelCenter(gr, nzones/2-2, r"$i-2$")
gpu.labelCenter(gr, nzones/2+2, r"$i+2$")


n = 0
while (n < nzones):
    gpu.drawCellAvg(gr, n, a[n], color="0.5")
    n += 1


# compute the parabolic coefficients
ap, am, a6 = gpu.ppm(a, nolimit=1)
lap, lam, la6 = gpu.ppm(a)

n = 2
while (n < nzones-2):
    gpu.drawParabola(gr, n, ap[n], am[n], a6[n], color="r", ls=":")
    gpu.drawParabola(gr, n, lap[n], lam[n], la6[n], color="r")
    #pylab.scatter([gr.xl[n], gr.xr[n]], [am[n], ap[n]])
    n += 1




pylab.axis([gr.xmin-0.5*gr.dx,gr.xmax+0.5*gr.dx, -0.25, 1.2])
pylab.axis("off")

pylab.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95)

f = pylab.gcf()
f.set_size_inches(8.0,2.0)


pylab.savefig("piecewise-parabolic.eps")
pylab.savefig("piecewise-parabolic.png")

               


