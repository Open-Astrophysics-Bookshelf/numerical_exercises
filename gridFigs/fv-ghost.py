import math
import numpy
import pylab
import grid_plot_util as gpu

# plot a simple finite-difference grid

#-----------------------------------------------------------------------------

nzones = 8
ng = 2

# data that lives on the grid
#a = numpy.array([0.3, 1.0, 0.9, 0.8, 0.25, 0.15, 0.5, 0.55])
a = numpy.array([0.3, 1.0, 0.9, 0.8, 0.25, 0.1, 0.5, 0.55])

gr = gpu.grid(nzones, ng)


pylab.clf()

gpu.drawGrid(gr, drawGhost=1, emphasizeEnd=1)

gpu.labelCenter(gr, ng+nzones/2,   r"$i$")
gpu.labelCenter(gr, ng+nzones/2-1, r"$i-1$")
gpu.labelCenter(gr, ng+nzones/2+1, r"$i+1$")

gpu.labelCenter(gr, gr.ilo, r"$\mathrm{lo}$")
gpu.labelCenter(gr, gr.ilo-1, r"$\mathrm{lo-1}$")
gpu.labelCenter(gr, gr.ilo-2, r"$\mathrm{lo-2}$")

gpu.labelCenter(gr, gr.ihi, r"$\mathrm{hi}$")
gpu.labelCenter(gr, gr.ihi+1, r"$\mathrm{hi+1}$")
gpu.labelCenter(gr, gr.ihi+2, r"$\mathrm{hi+2}$")

gpu.labelEdge(gr, ng+nzones/2,   r"$i-1/2$")
gpu.labelEdge(gr, ng+nzones/2+1,   r"$i+1/2$")

# draw the data
i = 0
while i < nzones:
    gpu.drawCellAvg(gr, ng+i, a[i], color="r")    
    i += 1
    
gpu.labelCellAvg(gr, ng+nzones/2, a[nzones/2], r"$\langle a\rangle_i$", color="r")

# label dx
pylab.plot([gr.xr[gr.ng+nzones/2-1], gr.xr[gr.ng+nzones/2-1]], [-0.35,-0.25], color="k")
pylab.plot([gr.xr[gr.ng+nzones/2], gr.xr[gr.ng+nzones/2]], [-0.35,-0.25], color="k")
pylab.plot([gr.xr[gr.ng+nzones/2-1], gr.xr[gr.ng+nzones/2]], [-0.3,-0.3], color="k")
pylab.text(gr.xc[gr.ng+nzones/2], -0.45, r"$\Delta x$", 
           horizontalalignment="center")



pylab.axis([gr.xmin-2.1*gr.dx,gr.xmax+2.1*gr.dx, -0.5, 1.6])
pylab.axis("off")

pylab.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95)

f = pylab.gcf()
f.set_size_inches(10.0,2.5)

pylab.savefig("fv_ghost.png")
pylab.savefig("fv_ghost.eps")

