import math
import numpy
import matplotlib.pyplot as plt
import grid_plot as gp

def simplegrid():

    nzones = 3
    gr = gp.FVGrid2d(nzones, nzones, ng=0)


    # plot a domain without ghostcells
    gr.draw_grid()

    #------------------------------------------------------------------------
    # label
    gr.label_cell_center(nzones//2, nzones//2, r"$a_{i,j}$", fontsize="large")
    gr.label_cell_center(nzones//2+1, nzones//2, r"$a_{i+1,j}$", fontsize="large")
    gr.label_cell_center(nzones//2, nzones//2+1, r"$a_{i,j+1}$", fontsize="large")
    gr.label_cell_center(nzones//2, nzones//2-1, r"$a_{i,j-1}$", fontsize="large")

    # i+1/2,j interface
    gr.mark_cell_left_state_x(nzones//2, nzones//2, r"$\hat{a}^{n+\myhalf}_{i+\myhalf,j,L}$", color="b")

    # i,j+1/2 interface
    gr.mark_cell_state_y(nzones//2, nzones//2, r"${a}^T_{i,j+\myhalf}$")

    # i,j-1/2 interface
    gr.mark_cell_state_y(nzones//2, nzones//2-1, r"${a}^T_{i,j-\myhalf}$", off_sign=-1)

    # helpful line showing the transverse bits
    plt.plot([gr.xc[nzones//2], gr.xc[nzones//2]], 
               [gr.yl[nzones//2]+0.025*gr.dy, gr.yr[nzones//2]-0.025*gr.dy], linestyle=":",
               color="0.5", zorder=1000)

    plt.plot([gr.xc[nzones//2], gr.xr[nzones//2]-0.26*gr.dx],
               [gr.yc[nzones//2], gr.yc[nzones//2]], linestyle=":",
               color="0.5", zorder=1000)

    plt.plot([gr.xr[nzones//2]-0.34*gr.dx, gr.xr[nzones//2]-0.26*gr.dx],
               [gr.yc[nzones//2]+0.04*gr.dy, gr.yc[nzones//2]], linestyle=":",
               color="0.5", zorder=1000)

    plt.plot([gr.xr[nzones//2]-0.34*gr.dx, gr.xr[nzones//2]-0.26*gr.dx],
               [gr.yc[nzones//2]-0.04*gr.dy, gr.yc[nzones//2]], linestyle=":",
               color="0.5", zorder=1000)


    # grid labels
    gr.label_center_x(nzones//2-1, r"$i-1$")
    gr.label_center_x(nzones//2, r"$i$")
    gr.label_center_x(nzones//2+1, r"$i+1$")

    gr.label_center_y(nzones//2-1, r"$j-1$")
    gr.label_center_y(nzones//2, r"$j$")
    gr.label_center_y(nzones//2+1, r"$j+1$")


    # axes
    gr.clean_axes()

    plt.subplots_adjust(left=0.02,right=0.98,bottom=0.02,top=0.98)

    f = plt.gcf()
    f.set_size_inches(7.0,7.0)

    plt.savefig("2dgrid-transverse.pdf")
               

if __name__== "__main__":
    simplegrid()
