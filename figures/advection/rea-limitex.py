# reconstruct - evolve - average: demonstrate what happens when we don't
# limit

from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
import grid_plot as gp

#-----------------------------------------------------------------------------
def evolve(pl, C, num, nolimit=1):
    
    # first frame -- the original cell-averages

    plt.clf()

    pl.gr.draw_grid()

    pl.gr.label_center(pl.gr.ng + pl.gr.nx//2,   r"$i$", fontsize="medium")
    pl.gr.label_center(pl.gr.ng + pl.gr.nx//2-1, r"$i-1$", fontsize="medium")
    pl.gr.label_center(pl.gr.ng + pl.gr.nx//2+1, r"$i+1$", fontsize="medium")
    pl.gr.label_center(pl.gr.ng + pl.gr.nx//2-2, r"$i-2$", fontsize="medium")
    pl.gr.label_center(pl.gr.ng + pl.gr.nx//2+2, r"$i+2$", fontsize="medium")

    # draw cell averages
    for n in range(pl.gr.ilo, gr.ihi+1):
        pl.draw_cell_avg(n, color="r")

    pl.gr.clean_axes(ylim=(-0.25, 1.2))

    f = plt.gcf()
    f.set_size_inches(8.0,2.0)

    if nolimit:
        plt.savefig("rea-nolimit-start_%3.3d.pdf" % (num))
    else:
        plt.savefig("rea-start_%3.3d.pdf" % (num))


    # second frame -- reconstruction

    # draw
    plt.clf()

    pl.gr.draw_grid()

    pl.gr.label_center(pl.gr.ng + pl.gr.nx//2,   r"$i$", fontsize="medium")
    pl.gr.label_center(pl.gr.ng + pl.gr.nx//2-1, r"$i-1$", fontsize="medium")
    pl.gr.label_center(pl.gr.ng + pl.gr.nx//2+1, r"$i+1$", fontsize="medium")
    pl.gr.label_center(pl.gr.ng + pl.gr.nx//2-2, r"$i-2$", fontsize="medium")
    pl.gr.label_center(pl.gr.ng + pl.gr.nx//2+2, r"$i+2$", fontsize="medium")

    # draw cell averages and slopes
    for n in range(pl.gr.ilo, gr.ihi+1):
        pl.draw_cell_avg(n, color="0.5", ls=":")
        pl.draw_slope(n, color="r")

    plt.axis([pl.gr.xmin-0.5*pl.gr.dx,pl.gr.xmax+0.5*pl.gr.dx, -0.25, 1.2])
    plt.axis("off")

    plt.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95)

    f = plt.gcf()
    f.set_size_inches(8.0,2.0)

    if nolimit:
        plt.savefig("rea-nolimit-reconstruction_%3.3d.pdf" % (num))
    else:
        plt.savefig("rea-reconstruction_%3.3d.pdf" % (num))


    # third frame -- evolve

    # draw
    plt.clf()

    pl.gr.draw_grid()

    pl.gr.label_center(pl.gr.ng + pl.gr.nx//2,   r"$i$", fontsize="medium")
    pl.gr.label_center(pl.gr.ng + pl.gr.nx//2-1, r"$i-1$", fontsize="medium")
    pl.gr.label_center(pl.gr.ng + pl.gr.nx//2+1, r"$i+1$", fontsize="medium")
    pl.gr.label_center(pl.gr.ng + pl.gr.nx//2-2, r"$i-2$", fontsize="medium")
    pl.gr.label_center(pl.gr.ng + pl.gr.nx//2+2, r"$i+2$", fontsize="medium")

    # draw cell averages and slopes
    for n in range(pl.gr.ilo, gr.ihi+1):
        #pl.draw_cell_avg(n, color="r")
        pl.draw_slope(n, color="0.5", ls=":")


    # evolve
    for n in range(pl.gr.ilo, pl.gr.ihi+1):
        pl.evolve_to_right(n, C, color="r")

    pl.gr.clean_axes(ylim=(-0.25, 1.2))

    f = plt.gcf()
    f.set_size_inches(8.0,2.0)

    if nolimit:
        plt.savefig("rea-nolimit-evolve_%3.3d.pdf" % (num))
    else:
        plt.savefig("rea-evolve_%3.3d.pdf" % (num))


    #-------------------------------------------------------------------------
    # fourth frame -- re-average

    # left states (we don't need the right state when u > 0)
    al = pl.gr.scratch_array()

    for n in range(pl.gr.ilo, pl.gr.ihi+2):
        al[n] = pl.a[n-1] + 0.5*(1.0 - C)*pl.slope[n-1]

    # the Riemann problem just picks the left state.  Do a conservative
    # update
    anew = pl.gr.scratch_array()

    anew[pl.gr.ilo:pl.gr.ihi+1] = pl.a[pl.gr.ilo:pl.gr.ihi+1] + \
            C*(al[pl.gr.ilo:pl.gr.ihi+1] - al[pl.gr.ilo+1:pl.gr.ihi+2])


    plt.clf()

    pl.gr.draw_grid()

    pl.gr.label_center(pl.gr.ng + pl.gr.nx//2,   r"$i$", fontsize="medium")
    pl.gr.label_center(pl.gr.ng + pl.gr.nx//2-1, r"$i-1$", fontsize="medium")
    pl.gr.label_center(pl.gr.ng + pl.gr.nx//2+1, r"$i+1$", fontsize="medium")
    pl.gr.label_center(pl.gr.ng + pl.gr.nx//2-2, r"$i-2$", fontsize="medium")
    pl.gr.label_center(pl.gr.ng + pl.gr.nx//2+2, r"$i+2$", fontsize="medium")


    # show old evolved profiles and new averages
    for n in range(pl.gr.ilo, pl.gr.ihi+1):
        pl.evolve_to_right(n, C, color="0.5", ls=":")

    # draw new averages
    pl.a[:] = anew[:]
    pl.fill_zero_gradient()
    pl.calculate_slopes()
    for n in range(pl.gr.ilo, pl.gr.ihi+1):
        pl.draw_cell_avg(n, color="r")

    pl.gr.clean_axes(ylim=(-0.25, 1.2))

    f = plt.gcf()
    f.set_size_inches(8.0,2.0)

    if nolimit:
        plt.savefig("rea-nolimit-final_%3.3d.pdf" % (num))
    else:
        plt.savefig("rea-final_%3.3d.pdf" % (num))

    return anew



#-----------------------------------------------------------------------------
ainit = np.array([1.0, 1.0, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1])
nzones = len(ainit)

# CFL number
C = 0.7
for nolimit in range(2):

    gr = gp.FVGrid(nzones, ng=4)

    a = gr.scratch_array()
    a[gr.ilo:gr.ihi+1] = ainit[:]

    pl = gp.PiecewiseLinear(gr, a, nolimit=nolimit)

    # loop
    for i in range(1,9):

        pl.fill_zero_gradient()
        print(i, pl.a[:])
        evolve(pl, C, i, nolimit=nolimit)

