import math
import numpy
import pylab
import grid_plot_util as gpu

# plot two stacked fv grids of different (2x) resolution to show prolongation

#-----------------------------------------------------------------------------

gr = []

nf = 2
while (nf <= 16):

    gr.append(gpu.grid(nf, ng=1, voff=2.0*len(gr)))
    nf = nf*2

pylab.clf()

for g in gr:
    gpu.drawGrid(g, emphasizeEnd=1, drawGhost=1, edgeTicks=0)

f = pylab.gcf()
f.set_size_inches(7.0,5.0)

grf = gr[0]
pylab.xlim(grf.xmin-1.1*grf.dx,grf.xmax+1.1*grf.dx)

pylab.axis("off")

pylab.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95)


pylab.savefig("mgtower.png")
pylab.savefig("mgtower.eps")

