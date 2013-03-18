import math
import numpy
import pylab
import grid_plot_util as gpu

def riemann():

    
    # grid info
    xmin = 0.0
    xmax = 1.0

    nzones = 1
    ng = 0
    
    gr = gpu.grid(nzones, xmin=xmin, xmax=xmax)


    #------------------------------------------------------------------------
    # plot a domain without ghostcells
    gpu.drawGrid(gr)

    gpu.labelCenter(gr, 0, r"$i$")

    gpu.labelCellCenter(gr, 0, r"$q_i$")
    
    gpu.markCellLeftState(gr, 0, r"$q_{i-1/2,R}^{n+1/2}$", color="r")
    gpu.markCellRightState(gr, 0, r"$q_{i+1/2,L}^{n+1/2}$", color="r")
    

    pylab.arrow(gr.xc[0]-0.05*gr.dx, 0.5, -0.13*gr.dx, 0, 
                shape='full', head_width=0.075, head_length=0.05, 
                lw=1, width=0.01,
                edgecolor="none", facecolor="r",
                length_includes_head=True, zorder=100)

    pylab.arrow(gr.xc[0]+0.05*gr.dx, 0.5, 0.13*gr.dx, 0, 
                shape='full', head_width=0.075, head_length=0.05, 
                lw=1, width=0.01,
                edgecolor="none", facecolor="r",
                length_includes_head=True, zorder=100)
    

    # pylab.xlim(xl[0]-0.5*dx,xr[2*ng+nzones-1]+0.5*dx)
    # pylab.ylim(-0.25, 0.75)
    pylab.axis("off")

    pylab.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95)

    f = pylab.gcf()
    f.set_size_inches(8.0,2.5)


    pylab.savefig("states.png")
    pylab.savefig("states.eps")


if __name__== "__main__":
    riemann()
