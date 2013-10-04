import math
import numpy
import pylab

def rh():

    # grid info
    xl = 0.0
    xr = 1.0

    dx = xr - xl

    xc = 0.5*(xl + xr)

    t0 = 0.0
    t1 = 1.0

    dt = t1 - t0


    #------------------------------------------------------------------------
    # plot a square representing [x, x+dx] x [t, t+dt] 

    pylab.arrow(0, 0, 1.3*dx, 0,
                shape="full", head_width=0.04, head_length=0.06,
                lw=1, width=0.005,
                edgecolor="none", facecolor="k",
                length_includes_head=True, zorder=100)

    pylab.text(1.35*dx, 0, r"$x$", fontsize=20, verticalalignment="center")

    pylab.arrow(0, 0, 0, 1.3*dt,
                shape="full", head_width=0.04, head_length=0.06,
                lw=1, width=0.005,
                edgecolor="none", facecolor="k",
                length_includes_head=True, zorder=100)

    pylab.text(0, 1.35*dt, r"$t$", fontsize=20, horizontalalignment="center")

    pylab.plot([dx,dx], [0,dt], ls=":", color="k", lw=2)
    pylab.text(dx, -0.1, r"$x_r$", fontsize=20)
    pylab.text(0, -0.1, r"$x_l$", fontsize=20)
    
    pylab.plot([0,dx], [dt,dt], ls=":", color="k", lw=2)
    pylab.text(-0.15, dt, r"$t^{n+1}$", fontsize=20)
    pylab.text(-0.15, 0, r"$t^n$", fontsize=20)

    pylab.plot([0,dx], [0,dt], color="k", lw=5, solid_capstyle="butt")
    
    pylab.annotate(r"shock: $S = \Delta x / \Delta t$", xy=(0.74*dx,0.76*dt),
                   xytext=(0.2*dx, 1.2*dt), textcoords="data",
                   arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2"))

    pylab.text(0.66*dx, 0.33*dt, r"$u_r$", color="k", fontsize=20)
    pylab.text(0.33*dx, 0.66*dt, r"$u_l$", color="k", fontsize=20)

    pylab.axis("off")

    pylab.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95)

    pylab.xlim(-0.1,1.7*dx)
    pylab.ylim(-0.1,1.7*dx)

    f = pylab.gcf()
    f.set_size_inches(8.0,8.0)


    pylab.savefig("rh.png")
    pylab.savefig("rh.eps")
               



if __name__== "__main__":
    rh()
