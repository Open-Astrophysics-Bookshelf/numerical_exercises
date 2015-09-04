import math
import numpy as np
import matplotlib.pyplot as plt
import grid_plot as gp

def riemann():

    # grid info
    xmin = 0.0
    xmax = 1.0

    nzones = 1
    ng = 0

    gr = gp.FVGrid(nzones, xmin=xmin, xmax=xmax)


    #------------------------------------------------------------------------
    # plot a domain without ghostcells
    gr.draw_grid()

    gr.label_center(0, r"$i$")

    gr.label_cell_center(0, r"$q_i$")

    gr.mark_cell_left_state(0, r"$q_{i-1/2,R}^{n+1/2}$", color="r")
    gr.mark_cell_right_state(0, r"$q_{i+1/2,L}^{n+1/2}$", color="r")


    plt.arrow(gr.xc[0]-0.05*gr.dx, 0.5, -0.13*gr.dx, 0,
                shape='full', head_width=0.075, head_length=0.05,
                lw=1, width=0.01,
                edgecolor="none", facecolor="r",
                length_includes_head=True, zorder=100)

    plt.arrow(gr.xc[0]+0.05*gr.dx, 0.5, 0.13*gr.dx, 0,
                shape='full', head_width=0.075, head_length=0.05,
                lw=1, width=0.01,
                edgecolor="none", facecolor="r",
                length_includes_head=True, zorder=100)


    plt.xlim(gr.xl[0]-0.25*gr.dx,gr.xr[2*ng+nzones-1]+0.25*gr.dx)
    # plt.ylim(-0.25, 0.75)
    plt.axis("off")

    plt.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95)

    f = plt.gcf()
    f.set_size_inches(4.0,2.5)


    plt.savefig("states.png")
    plt.savefig("states.pdf")


if __name__== "__main__":
    riemann()
