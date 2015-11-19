import math
import numpy as np
import matplotlib.pylab as plt
import grid_plot as gp

def riemann():

    # grid info
    xmin = 0.0
    xmax = 1.0

    nzones = 2
    ng = 0

    gr = gp.FVGrid(nzones, xmin=xmin, xmax=xmax)


    #------------------------------------------------------------------------
    # plot a domain without ghostcells
    gr.draw_grid()

    gr.label_center(0, r"$i$", fontsize="medium")
    gr.label_center(1, r"$i+1$", fontsize="medium")

    gr.label_edge(1, r"$i+1/2$", fontsize="medium")


    plt.arrow(gr.xc[0]+0.05*gr.dx, 0.5, 0.12*gr.dx, 0,
                shape='full', head_width=0.05, head_length=0.025,
                lw=1, width=0.02,
                edgecolor="none", facecolor="r",
                length_includes_head=True, zorder=100)

    plt.arrow(gr.xc[1]-0.1*gr.dx, 0.5, -0.12*gr.dx, 0,
                shape='full', head_width=0.05, head_length=0.025,
                lw=1, width=0.02,
                edgecolor="none", facecolor="r",
                length_includes_head=True, zorder=100)


    gr.mark_cell_left_state(1, r"$a_{i+1/2,L}^{n+1/2}$", fontsize="large",
                            color="b")
    gr.mark_cell_right_state(0, r"$a_{i+1/2,R}^{n+1/2}$", fontsize="large",
                             color="b")

    gr.label_cell_center(0, r"$a_i$")
    gr.label_cell_center(1, r"$a_{i+1}$")


    plt.xlim(gr.xl[0]-0.125*gr.dx,gr.xr[2*ng+nzones-1]+0.125*gr.dx)

    plt.ylim(-0.25, 1.0)
    plt.axis("off")

    plt.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95)

    f = plt.gcf()
    f.set_size_inches(7.0,2.0)

    plt.tight_layout()

    plt.savefig("riemann.pdf")

if __name__== "__main__":
    riemann()
