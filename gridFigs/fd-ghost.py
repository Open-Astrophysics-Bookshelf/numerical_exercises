import math
import numpy
import pylab
import grid_plot_util as gpu

# plot a simple finite-difference grid

#-----------------------------------------------------------------------------

nzones = 9

# data that lives on the grid
#a = numpy.array([0.3, 1.0, 0.9, 0.8, 0.25, 0.15, 0.5, 0.55])
a = numpy.array([0.55, 0.3, 1.0, 0.9, 0.8, 0.25, 0.1, 0.5, 0.55])

gr = gpu.grid(nzones, ng=1, fd=1)


pylab.clf()

gpu.drawGrid(gr, drawGhost=1)

labels = ["-1", "0", "1", "", "i-1", "i", "i+1", "", "N-2", "N-1", "N"]

i = gr.ilo-gr.ng
while (i < gr.ng+gr.nx+1):

    if not labels[i] == "":
        gpu.labelCenter(gr, i,   r"$%s$" % (labels[i]), fontsize="medium")
    i += 1

    
# draw the data
i = gr.ilo
while i < gr.ihi+1:
    gpu.drawFDData(gr, i, a[i-gr.ng], color="r")    
    i += 1
    

gpu.labelFD(gr, gr.ilo+4, a[gr.ilo+4-gr.ng], r"$a_i$", color="r")

# label dx
pylab.plot([gr.xc[gr.ng+nzones/2-1], gr.xc[gr.ng+nzones/2-1]], [-0.35,-0.25], color="k")
pylab.plot([gr.xc[gr.ng+nzones/2], gr.xc[gr.ng+nzones/2]], [-0.35,-0.25], color="k")
pylab.plot([gr.xc[gr.ng+nzones/2-1], gr.xc[gr.ng+nzones/2]], [-0.3,-0.3], color="k")
pylab.text(0.5*(gr.xc[gr.ng+nzones/2-1] + gr.xc[gr.ng+nzones/2]), -0.45, 
           r"$\Delta x$", 
           horizontalalignment="center", fontsize=16)


pylab.axis([gr.xmin-1.1*gr.dx,gr.xmax+1.1*gr.dx, -0.5, 1.3])
pylab.axis("off")

pylab.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95)

f = pylab.gcf()
f.set_size_inches(10.0,3.0)


pylab.savefig("fd_ghost.png")
pylab.savefig("fd_ghost.eps")

