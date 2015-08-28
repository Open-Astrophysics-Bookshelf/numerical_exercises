import math
import matplotlib.pyplot as plt
import numpy as np
import grid_plot as gp

def simplegrid():

    nzones = 7

    gr = gp.FVGrid(nzones, xmin=0, xmax=1)

    gr.draw_grid(edge_ticks=0)

    # label a few cell-centers
    gr.label_center(nzones/2, r"$i$")
    gr.label_center(nzones/2-1, r"$i-1$")
    gr.label_center(nzones/2+1, r"$i+1$")

    # label a few edges
    gr.label_edge(nzones/2, r"$i-1/2$")
    gr.label_edge(nzones/2+1, r"$i+1/2$")

    # sample data
    A = 0.4
    a = A*np.ones(nzones, dtype=np.float64)

    pc = gp.PiecewiseConstant(gr, a)

    # draw an average quantity
    pc.draw_cell_avg(nzones/2, color="r")
    pc.label_cell_avg(nzones/2, r"$\,\langle a \rangle_i$", color="r")

    gr.clean_axes()
    plt.ylim(-0.25, 1.5)

    plt.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95)

    f = plt.gcf()
    f.set_size_inches(10.0,2.5)

    plt.savefig("simplegrid2.png")
    plt.savefig("simplegrid2.pdf")

if __name__== "__main__":
    simplegrid()
