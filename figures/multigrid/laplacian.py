import numpy as np
import matplotlib.pyplot as plt
import grid_plot as gp

def laplace():

    # grid info
    xmin = 0.0
    xmax = 1.0

    nzones = 3
    ng = 0

    gr = gp.FVGrid(nzones, xmin=xmin, xmax=xmax)


    #------------------------------------------------------------------------
    # plot a domain without ghostcells
    gr.draw_grid(emphasize_end=True)

    # label a few
    gr.label_center(0, r"$i-1$", fontsize="medium")
    gr.label_center(1, r"$i$", fontsize="medium")
    gr.label_center(2, r"$i+1$", fontsize="medium")
    gr.label_edge(1, r"$i-\myhalf$", fontsize="medium")
    gr.label_edge(2, r"$i+\myhalf$", fontsize="medium")

    gr.label_cell_center(0, r"$\phi_{i-1}$")
    gr.label_cell_center(1, r"$\phi_{i}$", value=0.7)
    gr.label_cell_center(2, r"$\phi_{i+1}$")

    gr.mark_cell_edge(1, r"$\left .\frac{d\phi}{dx} \right |_{i-\myhalf}$", color="C0")
    gr.mark_cell_edge(2, r"$\left .\frac{d\phi}{dx} \right |_{i+\myhalf}$", color="C0")
    
    gr.label_cell_center(1, r"$\left .\frac{d^2\phi}{dx^2} \right |_{i}$", 
                         value=0.3, color="C1")

    gr.clean_axes(padding=False)
    plt.ylim(-0.25, 1.25)

    plt.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95)

    f = plt.gcf()
    f.set_size_inches(6.0,2.25)


    plt.tight_layout()

    plt.savefig("laplacian.pdf")


if __name__== "__main__":
    laplace()
