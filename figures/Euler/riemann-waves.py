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

    # label regions
    plt.text(gr.xr[0]-0.5*gr.dx, 0.3, r"$L$", horizontalalignment="center", color="C1")
    plt.text(gr.xr[0]-0.5*gr.dx, 0.22, r"$(\rho_L, u_L, p_L)$",
             horizontalalignment="center", color="C1", fontsize="small")

    plt.text(gr.xr[0]-0.33*gr.dx, 0.75, r"$L_*$",
             horizontalalignment="center", color="C1")
    plt.text(gr.xr[0]-0.33*gr.dx, 0.67, r"$(\rho_{\star L}, u_\star, p_\star)$",
             horizontalalignment="center", color="C1", fontsize="small")

    plt.text(gr.xr[0]+0.1*gr.dx, 0.75, r"$R_*$", horizontalalignment="center", color="C1")
    plt.text(gr.xr[0]+0.1*gr.dx, 0.67, r"$(\rho_{\star R}, u_\star, p_\star)$",
             horizontalalignment="center", color="C1", fontsize="small", zorder=1000)

    plt.text(gr.xr[0]+0.3*gr.dx, 0.3, r"$R$", horizontalalignment="center", color="C1")
    plt.text(gr.xr[0]+0.3*gr.dx, 0.22, r"$(\rho_R, u_R, p_R)$",
             horizontalalignment="center", color="C1", fontsize="small")



    gr.clean_axes(padding=False)
    plt.ylim(-0.2,1.2)

    plt.tight_layout()

    f = plt.gcf()
    f.set_size_inches(6,3.5)


    plt.savefig("riemann-waves.pdf")

if __name__== "__main__":
    simplegrid()
