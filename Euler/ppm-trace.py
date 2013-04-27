import math
import numpy
import pylab
import grid_plot_util as gpu


#-----------------------------------------------------------------------------

nzones = 8

# data that lives on the grid
#a = numpy.array([0.3, 1.0, 0.9, 0.8, 0.25, 0.15, 0.5, 0.55])
a = numpy.array([0.3, 1.0, 1.0, 0.8, 0.2, 0.15, 0.5, 0.55])

gr = gpu.grid(nzones)


pylab.clf()

gpu.drawGrid(gr, centerOnly=1)

gpu.labelCenter(gr, nzones/2-1,   r"$i$")
gpu.labelCenter(gr, nzones/2, r"$i+1$")


# compute the parabolic coefficients
ap, am, a6 = gpu.ppm(a, nolimit=1)
lap, lam, la6 = gpu.ppm(a)

n = gr.nx/2-1
while (n <= gr.nx/2):
    gpu.drawParabola(gr, n, lap[n], lam[n], la6[n], color="r")
    n += 1

nn = gr.nx/2-1
sigma = 0.6
gpu.ppmTraceLeft(gr, nn, lap[nn], lam[nn], la6[nn], sigma, color="0.75")

pylab.plot([gr.xr[gr.nx/2-1]-sigma*gr.dx, gr.xr[gr.nx/2-1]], [1.2, 1.2], color="k")
pylab.plot([gr.xr[gr.nx/2-1]-sigma*gr.dx, gr.xr[gr.nx/2-1]-sigma*gr.dx], [1.15, 1.25], color="k")
pylab.plot([gr.xr[gr.nx/2-1], gr.xr[gr.nx/2-1]], [1.15, 1.25], color="k")
pylab.text(gr.xr[gr.nx/2-1]-0.5*sigma*gr.dx, 1.3, r"$\sigma_{i}^{(\nu)} \Delta x$", horizontalalignment="center")

pylab.text(gr.xr[gr.nx/2-1]-0.5*sigma*gr.dx, 0.4, r"$\mathcal{I}_+^{(\nu)}$",
           color="r", horizontalalignment="center")

pylab.axis([gr.xl[gr.nx/2-1]-0.5*gr.dx,gr.xr[gr.nx/2]+0.5*gr.dx, -0.25, 1.4])
pylab.axis("off")

pylab.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95)

f = pylab.gcf()
f.set_size_inches(3.0,2.0)


pylab.savefig("ppm-trace.eps")
pylab.savefig("ppm-trace.png")
               


