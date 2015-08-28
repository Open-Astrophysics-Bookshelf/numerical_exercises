import math
import numpy as np
import matplotlib.pyplot as plt
import grid_plot as gp

def simplegrid():

    # grid info
    nzones = 7
    ng = 1

    gr = gp.FVGrid(nzones, ng, xmin=0.0, xmax=1.0)

    gr.draw_grid(emphasize_end=1, edge_ticks=0, draw_ghost=1)

    # label a few
    gr.label_center(ng+nzones/2, r"$i$", fontsize="medium")
    gr.label_center(ng+nzones/2-1, r"$i-1$", fontsize="medium")
    gr.label_center(ng+nzones/2+1, r"$i+1$", fontsize="medium")

    gr.label_center(ng-1, r"$\mathrm{lo}-1$", fontsize="medium")
    gr.label_center(ng, r"$\mathrm{lo}$", fontsize="medium")
    gr.label_center(ng+nzones-1, r"$\mathrm{hi}$", fontsize="medium")
    gr.label_center(ng+nzones, r"$\mathrm{hi+1}$", fontsize="medium")


    # label dx
    plt.plot([gr.xl[ng+nzones/2-2], gr.xl[ng+nzones/2-2]], 
               [-0.35,-0.25], color="k")

    plt.plot([gr.xr[ng+nzones/2-2], gr.xr[ng+nzones/2-2]], 
               [-0.35,-0.25], color="k")

    plt.plot([gr.xl[ng+nzones/2-2], gr.xr[ng+nzones/2-2]], 
               [-0.3,-0.3], color="k")

    plt.text(gr.xc[ng+nzones/2-2], -0.5, r"$\Delta x$", 
               horizontalalignment="center")


    plt.xlim(gr.xl[0]-0.05*gr.dx,gr.xr[2*ng+nzones-1]+0.25*gr.dx)
    plt.ylim(-0.5, 1.5)
    plt.axis("off")

    plt.subplots_adjust(left=0.025,right=0.95,bottom=0.05,top=0.95)

    f = plt.gcf()
    f.set_size_inches(10.0,2.5)

    plt.tight_layout()

    plt.savefig("simplegrid_gc.png")
    plt.savefig("simplegrid_gc.pdf", bbox_inches="tight")
               

if __name__== "__main__":
    simplegrid()
