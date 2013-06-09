import math
import numpy
import pylab
import grid_plot_util as gpu

# compare a finite-difference, cc-finite-difference, and finite-volume
# grid showing where the data falls at the boundaries

# optionally show a ghost cell

#-----------------------------------------------------------------------------

ncells = 8
nnodes = ncells+1
ng = 1

# data that lives on the grid nodes (this will be the f-d data)
#anode = numpy.random.rand(nnodes)
anode = numpy.array([ 0.22510379,  0.25017587,  0.09887572,  0.05568329,  0.55446147,  
                      0.56273893,  0.61099057,  0.63142573,  0.45309529])

acell = 0.5*(anode[0:nnodes-1] + anode[1:nnodes])


# periodic ghost cells
if (ng == 1):
    anode = numpy.insert(anode, 0, anode[nnodes-1])
    anode = numpy.append(anode, anode[ng])

    acell = numpy.insert(acell, 0, acell[ncells-1])
    acell = numpy.append(acell, acell[ng])


print "anode: ", len(anode), nnodes, ng


#-----------------------------------------------------------------------------
# finite difference
gr = gpu.grid(nnodes, ng=ng, fd=1)


pylab.clf()

gpu.drawGrid(gr, emphasizeEnd=1, drawGhost=1)

gpu.labelCenter(gr, ng+nnodes/2,   r"$i$")
gpu.labelCenter(gr, ng+nnodes/2-1, r"$i-1$")
gpu.labelCenter(gr, ng+nnodes/2+1, r"$i+1$")
#gpu.labelCenter(gr, ng+nnodes/2-2, r"$i-2$")
#gpu.labelCenter(gr, ng+nnodes/2+2, r"$i+2$")


# draw the data
i = 0
while i < nnodes+2*ng:
    gpu.drawFDData(gr, i, anode[i], color="r")    
    i += 1
    

gpu.labelFD(gr, ng+nnodes/2, anode[ng+nnodes/2], r"$f_i$", color="r")

# label dx
pylab.plot([gr.xc[gr.ng+nnodes/2-1], gr.xc[gr.ng+nnodes/2-1]], 
           [-0.35,-0.25], color="k")
pylab.plot([gr.xc[gr.ng+nnodes/2], gr.xc[gr.ng+nnodes/2]], 
           [-0.35,-0.25], color="k")
pylab.plot([gr.xc[gr.ng+nnodes/2-1], gr.xc[gr.ng+nnodes/2]], 
           [-0.3,-0.3], color="k")
pylab.text(0.5*(gr.xc[gr.ng+nnodes/2-1] + gr.xc[gr.ng+nnodes/2]), 
           -0.45, r"$\Delta x$", 
           horizontalalignment="center")


pylab.axis([gr.xmin-ng*gr.dx-0.25*gr.dx,gr.xmax+ng*gr.dx+0.25*gr.dx, 
            -0.5, 1.2])
pylab.axis("off")

pylab.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95)

f = pylab.gcf()
f.set_size_inches(10.0,2.5)

pylab.savefig("fd_grid_bnd.png")
pylab.savefig("fd_grid_bnd.eps")



#-----------------------------------------------------------------------------
# cell-centered finite-differences

pylab.clf()

gr = gpu.grid(ncells, ng=ng)

gpu.drawGrid(gr, emphasizeEnd=1, drawGhost=1)

gpu.labelCenter(gr, ng+ncells/2,   r"$i$")
gpu.labelCenter(gr, ng+ncells/2-1, r"$i-1$")
gpu.labelCenter(gr, ng+ncells/2+1, r"$i+1$")
#gpu.labelCenter(gr, ng+ncells/2-2, r"$i-2$")
#gpu.labelCenter(gr, ng+ncells/2+2, r"$i+2$")


# draw the data
i = 0
while i < ncells+2*ng:
    gpu.drawFDData(gr, i, acell[i], color="r")    
    i += 1
    

gpu.labelFD(gr, ng+ncells/2, acell[ng+ncells/2], r"$f_i$", color="r")

# label dx
pylab.plot([gr.xc[gr.ng+ncells/2-1], gr.xc[gr.ng+ncells/2-1]], 
           [-0.35,-0.25], color="k")
pylab.plot([gr.xc[gr.ng+ncells/2], gr.xc[gr.ng+ncells/2]], 
           [-0.35,-0.25], color="k")
pylab.plot([gr.xc[gr.ng+ncells/2-1], gr.xc[gr.ng+ncells/2]], 
           [-0.3,-0.3], color="k")
pylab.text(0.5*(gr.xc[gr.ng+ncells/2-1] + gr.xc[gr.ng+ncells/2]), 
           -0.45, r"$\Delta x$", 
           horizontalalignment="center")



pylab.axis([gr.xmin-ng*gr.dx-0.25*gr.dx,gr.xmax+ng*gr.dx+0.25*gr.dx, 
            -0.5, 1.2])
pylab.axis("off")

pylab.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95)

f = pylab.gcf()
f.set_size_inches(10.0,2.5)

pylab.savefig("ccfd_grid_bnd.png")
pylab.savefig("ccfd_grid_bnd.eps")

#-----------------------------------------------------------------------------
# finite-volume

gr = gpu.grid(ncells, ng=ng)


pylab.clf()

gpu.drawGrid(gr, emphasizeEnd=1, drawGhost=1)

gpu.labelCenter(gr, ng+ncells/2,   r"$i$")
gpu.labelCenter(gr, ng+ncells/2-1, r"$i-1$")
gpu.labelCenter(gr, ng+ncells/2+1, r"$i+1$")
#gpu.labelCenter(gr, ng+ncells/2-2, r"$i-2$")
#gpu.labelCenter(gr, ng+ncells/2+2, r"$i+2$")

gpu.labelCenter(gr, ng+0, r"$\mathrm{lo}$")
gpu.labelCenter(gr, ng-1, r"$\mathrm{lo-1}$")

gpu.labelCenter(gr, ng+gr.nx-1, r"$\mathrm{hi}$")
gpu.labelCenter(gr, ng+gr.nx, r"$\mathrm{hi+1}$")


gpu.labelEdge(gr, ng+ncells/2,   r"$i-1/2$")
gpu.labelEdge(gr, ng+ncells/2+1,   r"$i+1/2$")

gpu.labelEdge(gr, ng+gr.nx,   r"$b$")
gpu.labelEdge(gr, ng+0,   r"$a$")


# draw the data
i = 0
while i < ncells+2*ng:
    gpu.drawCellAvg(gr, i, acell[i], color="r")    
    i += 1
    
gpu.labelCellAvg(gr, ng+ncells/2, acell[ng+ncells/2], 
                 r"$\langle f\rangle_i$", color="r")

# label dx
pylab.plot([gr.xr[gr.ng+ncells/2-1], gr.xr[gr.ng+ncells/2-1]], 
           [-0.35,-0.25], color="k")
pylab.plot([gr.xr[gr.ng+ncells/2], gr.xr[gr.ng+ncells/2]], 
           [-0.35,-0.25], color="k")
pylab.plot([gr.xr[gr.ng+ncells/2-1], gr.xr[gr.ng+ncells/2]], 
           [-0.3,-0.3], color="k")
pylab.text(gr.xc[gr.ng+ncells/2], -0.45, r"$\Delta x$", 
           horizontalalignment="center")



pylab.axis([gr.xmin-ng*gr.dx-0.25*gr.dx,gr.xmax+ng*gr.dx+0.25*gr.dx, 
            -0.5, 1.2])
pylab.axis("off")

pylab.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95)

f = pylab.gcf()
f.set_size_inches(10.0,2.5)

pylab.savefig("fv_grid_bnd.png")
pylab.savefig("fv_grid_bnd.eps")

