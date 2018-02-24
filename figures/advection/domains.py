import numpy as np
import matplotlib.pylab as plt
import matplotlib as mpl

mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['mathtext.rm'] = 'serif'

mpl.rcParams['font.size'] = 12
mpl.rcParams['legend.fontsize'] = 'large'
mpl.rcParams['figure.titlesize'] = 'medium'

def domains(method="FTCS"):

    # grid info
    xl = 0.0
    xc = 1.0
    xr = 2.0

    dx = xr - xc

    t0 = 0.0
    t1 = 1.0

    dt = t1 - t0

    C = 0.6
    xp = xc - C*dx

    plt.clf()

    # plot a square representing [x, x+dx] x [t, t+dt]

    # x-axis
    plt.arrow(0, 0, 2.5*dx, 0,
              shape="full", head_width=0.04, head_length=0.06,
              lw=1, width=0.005,
              facecolor="k",
              length_includes_head=True, zorder=100)

    plt.text(2.55*dx, 0, r"$x$", fontsize=16, verticalalignment="center")

    # x points labeled
    for p, l in [((xc, 0), r"$x_i$"), ((xl, 0), r"$x_{i-1}$"), ((xr, 0), r"$x_{i+1}$")]:
        plt.text(p[0], p[1]-0.07, l, fontsize=16,
             horizontalalignment="center", verticalalignment="top")
        plt.plot([p[0], p[0]], [p[1], p[1]-0.04], color="k")

    plt.text(xp, -0.07, r"$x_i - C \Delta x$", fontsize=16,
             horizontalalignment="center", verticalalignment="top")
    plt.plot([xp, xp], [0, -0.04], color="k")

    # points
    plt.scatter([xl, xc, xr, xc], [0, 0, 0, dt], color="C3",
                marker="o", s=40, zorder=1000, alpha=1.0)

    # data points
    eps = 0.02*dx
    for p, l in [((xl, 0), r"$a_{i-1}^n$"), ((xc, 0), r"$a_{i}^n$"),
                 ((xr, 0), r"$a_{i+1}^n$"), ((xc, dt), r"$a_{i}^{n+1}$")]:
        plt.text(p[0]+eps, p[1]+eps, l,
                 horizontalalignment="left", verticalalignment="bottom",
                 fontsize=16, color="C3", zorder=1000)

    # time axis
    plt.arrow(0, 0, 0, 1.3*dt,
              shape="full", head_width=0.04, head_length=0.06,
              lw=1, width=0.005,
              facecolor="k",
              length_includes_head=True, zorder=100)

    plt.text(0, 1.35*dt, r"$t$", fontsize=16, horizontalalignment="center")

    plt.text(-0.25, dt, r"$t^{n+1}$", fontsize=16)
    plt.plot([-0.04, 0], [dt, dt], color="k")

    plt.text(-0.25, 0, r"$t^n$", fontsize=16)
    plt.plot([-0.04, 0], [0, 0], color="k")

    plt.plot([0, 2.0*dx], [dt, dt], ls=":", color="0.5")
    plt.plot([dx, dx], [0, dt], ls=":", color="0.5")
    plt.plot([2.0*dx, 2.0*dx], [0, dt], ls=":", color="0.5")

    # domain of dependence
    if method == "upwind":
        plt.fill([0, xc, xc, 0], [0, dt, 0, 0], color="C0", lw=2, alpha=0.5)
    elif method == "FTCS":
        plt.fill([0, xc, xr, 0], [0, dt, 0, 0], color="C0", lw=2, alpha=0.5)
    elif method == "downwind":
        plt.fill([xc, xc, xr, 0], [0, dt, 0, 0], color="C0", lw=2, alpha=0.5)

    # true domain of dependence
    plt.fill([xp, xc, xc, 0], [0, dt, 0, 0], color="C1", lw=2, alpha=0.5)


    # label
    plt.text(xr+0.1*dx, 0.5*dt, method, fontsize=16,
             horizontalalignment="left")

    plt.axis("off")

    plt.subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.95)

    plt.xlim(-0.4,2.8*dx)
    plt.ylim(-0.1,1.4*dx)

    f = plt.gcf()
    f.set_size_inches(9.0,5.0)

    plt.tight_layout()
    plt.savefig("domains_{}.pdf".format(method), bbox_inches="tight")


if __name__== "__main__":
    domains(method="upwind")
    domains(method="FTCS")
    domains(method="downwind")
