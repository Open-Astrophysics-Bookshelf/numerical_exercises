import numpy as np
import matplotlib.pyplot as plt
import grid_plot as gp

def riemann():

    # grid info
    xmin = 0.0
    xmax = 1.0

    nzones = 4
    ng = 2

    gr = gp.FVGrid(nzones, ng=ng, xmin=xmin, xmax=xmax)
    
    # interior and ghost cell initialization
    a = gr.scratch_array()

    a[gr.ilo:gr.ihi+1] = np.array([0.8, 0.7, 0.4, 0.5])
    a[0:gr.ilo] = a[gr.ihi-1:gr.ihi+1]
    a[gr.ihi:2*gr.ng+gr.nx] = a[gr.ihi]

    pc = gp.PiecewiseConstant(gr, a)
    pl = gp.PiecewiseLinear(gr, a, nolimit=1)



    #------------------------------------------------------------------------
    # plot a domain without ghostcells
    gr.draw_grid(draw_ghost=0, emphasize_end=0)

    zi = gr.ilo+1

    gr.label_center(zi-1, r"${i-1}$", fontsize="medium")
    gr.label_center(zi, r"${i}$", fontsize="medium")
    gr.label_center(zi+1, r"${i+1}$", fontsize="medium")

    gr.label_edge(zi+1, r"${i}+\myhalf$", fontsize="medium")

    # draw cell averages
    for n in range(gr.ilo, gr.ihi+1):
        pc.draw_cell_avg(n, color="0.5", ls=":")


    # draw slopes
    for n in range(gr.ilo, gr.ihi+1):
        pl.draw_slope(n, color="r")

    # compute the states to the left and right of lo-1/2
    al = a[zi] + 0.5*pl.slope[zi]
    ar = a[zi+1] - 0.5*pl.slope[zi+1]

    # L
    gr.mark_cell_right_state(zi, r"$a_{i+\myhalf,L}$",
                             value=al, vertical="top", color="b", zorder=100)

    # R
    gr.mark_cell_left_state(zi+1, r"$a_{i+\myhalf,R}$",
                            value=ar, vertical="top", color="b", zorder=100)


    plt.xlim(gr.xl[gr.ilo]-0.025*gr.dx, gr.xr[gr.ihi]+0.15*gr.dx)
    plt.ylim(-0.25, 1.1)
    plt.axis("off")

    plt.subplots_adjust(left=0.025,right=0.95,bottom=0.05,top=0.95)

    f = plt.gcf()
    f.set_size_inches(8.0,2.0)

    plt.tight_layout()

    plt.savefig("riemann-mol.pdf")

if __name__== "__main__":
    riemann()
