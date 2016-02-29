import math
import numpy as np
import matplotlib.pyplot as plt
import grid_plot as gp

def simplegrid():

    nzones = 5
    gr = gp.FVGrid2d(nzones, nzones, ng=0)

    # plot a domain without ghostcells
    gr.draw_grid()

    #------------------------------------------------------------------------
    # label
    gr.label_cell_center(nzones/2, nzones/2, r"$a_{i,j}$", 
                         fontsize="x-large", color="r")
    gr.shade_cell(nzones/2, nzones/2)

    # grid labels
    gr.label_center_x(nzones/2-2, r"$i-2$")
    gr.label_center_x(nzones/2-1, r"$i-1$")
    gr.label_center_x(nzones/2, r"$i$")
    gr.label_center_x(nzones/2+1, r"$i+1$")
    gr.label_center_x(nzones/2+2, r"$i+2$")

    gr.label_center_y(nzones/2-2, r"$j-2$")
    gr.label_center_y(nzones/2-1, r"$j-1$")
    gr.label_center_y(nzones/2, r"$j$")
    gr.label_center_y(nzones/2+1, r"$j+1$")
    gr.label_center_y(nzones/2+2, r"$j+2$")

    # axes
    gr.clean_axes()
    plt.subplots_adjust(left=0.025,right=0.98,bottom=0.1,top=0.98)

    f = plt.gcf()
    f.set_size_inches(7.0,7.0)

    plt.savefig("2dFV.png")

if __name__== "__main__":
    simplegrid()
