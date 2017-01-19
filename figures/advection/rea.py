# reconstruct - evolve - average

import numpy as np
import matplotlib.pyplot as plt
import grid_plot as gp


#-----------------------------------------------------------------------------

atemp = np.array([1.0, 1.0, 0.9, 0.8, 0.25, 0.1, 0.1, 0.1])
nzones = len(atemp)

# CFL number
C = 0.6


gr = gp.FVGrid(nzones, ng=4)

a = gr.scratch_array()
a[gr.ilo:gr.ihi+1] = atemp[:]

pl = gp.PiecewiseLinear(gr, a)
pl.fill_zero_gradient()


#-----------------------------------------------------------------------------
# first frame -- the original cell-averages

plt.clf()

gr.draw_grid()

labels = ["$i-2$", "$i-1$", "$i$", "$i+1$", "$i+2$"]
indices = [gr.ng+nzones//2-2, gr.ng+nzones//2-1, 
           gr.ng+nzones//2, gr.ng+nzones//2+1, gr.ng+nzones//2+2]

for i, l in zip(indices, labels):
    gr.label_center(i, l)

# draw cell averages
for n in range(gr.ilo, gr.ihi+1):
    pl.draw_cell_avg(n, color="r")

gr.clean_axes(ylim=(-0.25, 1.2))

f = plt.gcf()
f.set_size_inches(8.0,2.0)

plt.savefig("rea-start.pdf")


#-----------------------------------------------------------------------------
# second frame -- reconstruction

plt.clf()

gr.draw_grid()

for i, l in zip(indices, labels):
    gr.label_center(i, l)

# draw cell averages and slopes
for n in range(gr.ilo, gr.ihi+1):
    pl.draw_cell_avg(n, color="0.5", ls=":")
    pl.draw_slope(n, color="r")


gr.clean_axes(ylim=(-0.25, 1.2))

f = plt.gcf()
f.set_size_inches(8.0,2.0)

plt.savefig("rea-reconstruction.pdf")


#-----------------------------------------------------------------------------
# third frame -- shade the portion that will come into the cell

plt.clf()

gr.draw_grid()

for i, l in zip(indices, labels):
    gr.label_center(i, l)

# draw cell averages and slopes
for n in range(gr.ilo, gr.ihi+1):
    pl.draw_slope(n, color="r")


# shade regions
ii = gr.ng + nzones//2-1
pl.slope_trace_left(ii, C, color="0.75")
pl.slope_trace_right(ii+1, 1.0-C, color="0.75")

gr.clean_axes(ylim=(-0.25, 1.2))

f = plt.gcf()
f.set_size_inches(8.0,2.0)

plt.savefig("rea-trace.pdf")


#-----------------------------------------------------------------------------
# fourth frame -- evolve

plt.clf()

gr.draw_grid()

for i, l in zip(indices, labels):
    gr.label_center(i, l)

# draw cell averages and slopes
for n in range(gr.ilo, gr.ihi+1):
    pl.draw_slope(n, color="0.75", ls=":")
    pl.evolve_to_right(n, C, color="r")

gr.clean_axes(ylim=(-0.25, 1.2))

f = plt.gcf()
f.set_size_inches(8.0,2.0)

plt.savefig("rea-evolve.pdf")


#-----------------------------------------------------------------------------
# fifth frame -- re-average

# left states (we don't need the right state when u > 0)
al = np.zeros(2*gr.ng + gr.nx, dtype=np.float64)

for n in range(gr.ilo, gr.ihi+2):
    al[n] = a[n-1] + 0.5*(1 - C)*pl.slope[n-1]

# the Riemann problem just picks the right state.  Do a conservative
# update
anew = gr.scratch_array()

anew[gr.ilo:gr.ihi+1] = a[gr.ilo:gr.ihi+1] + \
    C*(al[gr.ilo:gr.ihi+1] - al[gr.ilo+1:gr.ihi+2])


plt.clf()

gr.draw_grid()

for i, l in zip(indices, labels):
    gr.label_center(i, l)


# show the evolved profiles from the old time
for n in range(gr.ilo, gr.ihi+1):
    pl.evolve_to_right(n, C, color="0.5", ls=":")


pl.a[:] = anew[:]
pl.fill_zero_gradient()
pl.calculate_slopes()

for n in range(pl.gr.ilo, pl.gr.ihi+1):
    pl.draw_cell_avg(n, color="r")


gr.clean_axes(ylim=(-0.25, 1.2))

f = plt.gcf()
f.set_size_inches(8.0,2.0)

plt.savefig("rea-final.pdf")



