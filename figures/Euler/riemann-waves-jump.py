import matplotlib.pyplot as plt
import grid_plot as gp

def simplegrid():

    # grid info
    xmin = 0.0
    xmax = 1.0

    nzones = 2
    ng = 0

    gr = gp.FVGrid(nzones, xmin=xmin, xmax=xmax)


    #------------------------------------------------------------------------
    # plot a domain without ghostcells
    gr.draw_grid(draw_end=False, edge_ticks=False)

    gr.label_center(0, r"$i$")
    gr.label_center(1, r"$i+1$")

    gr.label_edge(1, r"$q_{i+\myhalf}$")


    # draw waves
    # u - c
    plt.plot([gr.xr[0], gr.xr[0]-0.75*gr.dx], [0,1.0], color="C0", ls="-")
    plt.text(gr.xr[0]-0.75*gr.dx, 1.0+0.05, "$\lambda^{(-)} =\, u - c$", 
               horizontalalignment="center")

    # u
    plt.plot([gr.xr[0], gr.xr[0]-0.2*gr.dx], [0,1.0], color="C0", ls="-")
    plt.text(gr.xr[0]-0.2*gr.dx, 1.0+0.05, "$\lambda^{(\circ)} =\, u$", 
               horizontalalignment="center")

    # u + c
    plt.plot([gr.xr[0], gr.xr[0]+0.4*gr.dx], [0,1.0], color="C0", ls="-")
    plt.text(gr.xr[0]+0.4*gr.dx, 1.0+0.05, "$\lambda^{(+)} =\, u + c$", 
               horizontalalignment="center")


    plt.plot([gr.xl[0], gr.xr[0]], [0.3, 0.3], color="C1", linewidth=2)
    plt.text(gr.xc[0], 0.33, r"$\langle q \rangle_i$", color="C1")
    
    plt.plot([gr.xl[1], gr.xr[1]], [0.6, 0.6], color="C1", linewidth=2)
    plt.text(gr.xc[1], 0.63, r"$\langle q \rangle_{i+1}$", color="C1")

    gr.clean_axes(padding=False)
    plt.ylim(-0.2,1.2)

    plt.tight_layout()

    f = plt.gcf()
    f.set_size_inches(6, 3.5)


    plt.savefig("riemann-waves-jump.png", dpi=150)

if __name__== "__main__":
    simplegrid()
