import numpy as np
import matplotlib.pylab as plt
import matplotlib as mpl

mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['mathtext.rm'] = 'serif'

mpl.rcParams['font.size'] = 12
mpl.rcParams['legend.fontsize'] = 'large'
mpl.rcParams['figure.titlesize'] = 'medium'

def rh():

    # grid info
    xl = 0.0
    xr = 1.0

    dx = xr - xl

    xc = 0.5*(xl + xr)

    t0 = 0.0
    t1 = 1.0

    dt = t1 - t0

    # plot a square representing [x, x+dx] x [t, t+dt] 

    # x-axis
    plt.arrow(0, 0, 1.3*dx, 0,
              shape="full", head_width=0.04, head_length=0.06,
              lw=1, width=0.005,
              facecolor="k",
              length_includes_head=True, zorder=100)

    plt.text(1.35*dx, 0, r"$x$", fontsize=20, verticalalignment="center")

    plt.text(dx, -0.1, r"$x_r$", fontsize=20)
    plt.text(0, -0.1, r"$x_l$", fontsize=20)

    # time axis
    plt.arrow(0, 0, 0, 1.3*dt,
              shape="full", head_width=0.04, head_length=0.06,
              lw=1, width=0.005,
              facecolor="k",
              length_includes_head=True, zorder=100)

    plt.text(0, 1.35*dt, r"$t$", fontsize=20, horizontalalignment="center")

    plt.text(-0.17, dt, r"$t^{n+1}$", fontsize=20)
    plt.text(-0.17, 0, r"$t^n$", fontsize=20)

    
    # space-time volume
    plt.plot([dx,dx], [0,dt], ls="--", color="0.5", lw=2)
    plt.plot([0,dx], [dt,dt], ls="--", color="0.5", lw=2)


    # diagonal representing S
    plt.plot([0,dx], [0,dt], color="k", lw=5, solid_capstyle="butt")
    
    plt.annotate(r"shock: $S = \Delta x / \Delta t$", xy=(0.74*dx,0.76*dt),
                 xytext=(0.2*dx, 1.2*dt), textcoords="data",
                 fontsize=18,
                 arrowprops=dict(arrowstyle="->", 
                                 connectionstyle="arc3,rad=.2"))

    # states
    plt.text(0.66*dx, 0.33*dt, r"$u_r$", color="k", fontsize=20)
    plt.text(0.33*dx, 0.66*dt, r"$u_l$", color="k", fontsize=20)

    
    # fluxes
    plt.arrow(-0.1*dx, 0.5*dt, 0.2*dx, 0,
              shape="full", head_width=0.04, head_length=0.06,
              lw=1, width=0.005,
              edgecolor="r", facecolor="r",
              length_includes_head=True, zorder=100)

    plt.text(-0.12*dx, 0.5*dt, r"$f = f(u_l)$",
             horizontalalignment="right",
             verticalalignment="center", color="r", fontsize=18)

    plt.arrow(0.9*dx, 0.5*dt, 0.2*dx, 0,
              shape="full", head_width=0.04, head_length=0.06,
              lw=1, width=0.005, 
              edgecolor="r", facecolor="r",
              length_includes_head=True, zorder=100)

    plt.text(1.12*dx, 0.5*dt, r"$f = f(u_r)$",
             horizontalalignment="left",
             verticalalignment="center", color="r", fontsize=18)


    plt.axis("off")

    plt.subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.95)

    plt.xlim(-0.4,1.4*dx)
    plt.ylim(-0.1,1.4*dx)

    f = plt.gcf()
    f.set_size_inches(7.0,6.0)

    plt.savefig("rh.pdf")


if __name__== "__main__":
    rh()
