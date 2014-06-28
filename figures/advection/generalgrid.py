import math
import numpy
import pylab
import grid_plot_util as gpu


nzones = 7

# data that lives on the grid
a = numpy.array([0.3, 1.0, 0.9, 0.8, 0.25, 0.15, 0.5])

gr = gpu.grid(nzones)


gpu.drawGrid(gr)

gpu.labelCenter(gr, nzones/2,   r"$i$")
gpu.labelCenter(gr, nzones/2-1, r"$i-1$")
gpu.labelCenter(gr, nzones/2+1, r"$i+1$")
gpu.labelCenter(gr, nzones/2-2, r"$i-2$")
gpu.labelCenter(gr, nzones/2+2, r"$i+2$")

#labelEdge(gr, nzones/2,   r"$i-1/2$")
#labelEdge(gr, nzones/2+1, r"$i+1/2$")

n = 0
while (n < nzones):
    gpu.drawCellAvg(gr, n, a[n])
    n += 1


# compute the slopes as simple centered-differences
da = gpu.lslopes(a, nolimit=1)

n = 1
while (n < nzones-1):
    gpu.drawSlope(gr, n, da[n], a[n], color="blue")
    n += 1


# compute the limited slopes
lda = gpu.lslopes(a)

n = 1
while (n < nzones-1):
    gpu.drawSlope(gr, n, lda[n], a[n], color="r")
    n += 1




pylab.axis([gr.xmin-0.5*gr.dx,gr.xmax+0.5*gr.dx, -0.25, 1.2])
pylab.axis("off")

pylab.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95)

f = pylab.gcf()
f.set_size_inches(8.0,2.0)


pylab.savefig("generalgrid.png")
pylab.savefig("generalgrid.eps")
               


