import numpy as np
import matplotlib.pyplot as plt
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
    gr.draw_grid(emphasize_end=True)

    # label a few
    gr.label_center(0, r"$i$", fontsize="medium")
    gr.label_center(1, r"$i+1$", fontsize="medium")
    gr.label_edge(1, r"$i+\myhalf$", fontsize="medium")

    gr.mark_cell_left_state(1, r"$U_{i+\myhalf,R}^{n+\myhalf}$", fontsize="large",
                            color="b")
    gr.mark_cell_right_state(0, r"$U_{i+\myhalf,L}^{n+\myhalf}$", fontsize="large",
                             color="b")

    gr.label_cell_center(0, r"$U_i$")
    gr.label_cell_center(1, r"$U_{i+1}$")



    # flux
    plt.arrow(gr.xl[ng+nzones//2]-0.25*gr.dx, 1.05, 0.5*gr.dx, 0, 
                shape='full', head_width=0.075, head_length=0.05, 
                lw=1, width=0.03,
                edgecolor="none", facecolor="red",
                length_includes_head=True, zorder=100)
    
    plt.text(gr.xl[ng+nzones//2], 1.15, r"$F(U_{i+\myhalf}^{n+\myhalf})$", color="red",
               horizontalalignment="center")


    gr.clean_axes(padding=False)
    plt.ylim(-0.25, 1.25)

    plt.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95)

    f = plt.gcf()
    f.set_size_inches(6.0,2.25)


    plt.tight_layout()

    plt.savefig("riemann_comp.pdf")


if __name__== "__main__":
    riemann()
