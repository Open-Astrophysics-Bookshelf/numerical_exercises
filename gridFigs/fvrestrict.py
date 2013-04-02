import math
import numpy
import pylab
import grid_plot_util as gpu

# plot two stacked fv grids of different (2x) resolution to show prolongation

#-----------------------------------------------------------------------------

nf = 4
nc = nf/2

grf = gpu.grid(nf, voff=2.0)
grc = gpu.grid(nc)


pylab.clf()

gpu.drawGrid(grf)
gpu.drawGrid(grc)

gpu.labelCenter(grf, nf/2-2,   r"$i-2$")
gpu.labelCenter(grf, nf/2-1,   r"$i-1$")
gpu.labelCenter(grf, nf/2,   r"$i$")
gpu.labelCenter(grf, nf/2+1,   r"$i+1$")


gpu.labelCenter(grc, nc/2-1,   r"$j-1$")
gpu.labelCenter(grc, nc/2,   r"$j$")

gpu.labelCellCenter(grf, nf/2-2, r"$\phi_{i-2}^h$")
gpu.labelCellCenter(grf, nf/2-1, r"$\phi_{i-1}^h$")
gpu.labelCellCenter(grf, nf/2, r"$\phi_i^h$")
gpu.labelCellCenter(grf, nf/2+1, r"$\phi_{i+1}^h$")

gpu.labelCellCenter(grc, nc/2-1, r"$\phi_{j-1}^{2h}$")
gpu.labelCellCenter(grc, nc/2,   r"$\phi_{j}^{2h}$")
    

# connect the dots...

pylab.plot([grf.xl[nf/2-2], grf.xl[nf/2-2]], [-0.25, 3.25], ":", color="0.5")
pylab.plot([grf.xl[nf/2], grf.xl[nf/2]], [-0.25, 3.25], ":", color="0.5")
pylab.plot([grf.xr[nf/2+1], grf.xr[nf/2+1]], [-0.25, 3.25], ":", color="0.5")


pylab.axis([grf.xmin-0.5*grf.dx,grf.xmax+0.5*grf.dx, -0.5, 3.5])
pylab.axis("off")

pylab.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95)

f = pylab.gcf()
f.set_size_inches(6.0,5.0)

pylab.savefig("fvrestrict.png")
pylab.savefig("fvrestrict.eps")

