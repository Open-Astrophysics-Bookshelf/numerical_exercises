import math
import numpy as np
import matplotlib.pyplot as plt
import grid_plot as gp


nzones = 7

# data that lives on the grid
ainit = np.array([0.3, 1.0, 0.9, 0.8, 0.25, 0.15, 0.5])

gr = gp.FVGrid(nzones)
a = gr.scratch_array()
a[gr.ilo:gr.ihi+1] = ainit[:]

pl_nolim = gp.PiecewiseLinear(gr, a, nolimit=1)
pl = gp.PiecewiseLinear(gr, a)

gr.draw_grid()

gr.label_center(nzones//2,   r"$i$")
gr.label_center(nzones//2-1, r"$i-1$")
gr.label_center(nzones//2+1, r"$i+1$")
gr.label_center(nzones//2-2, r"$i-2$")
gr.label_center(nzones//2+2, r"$i+2$")

#labelEdge(gr, nzones//2,   r"$i-1/2$")
#labelEdge(gr, nzones//2+1, r"$i+1/2$")

for n in range(nzones):
    pl.draw_cell_avg(n, ls=":", color="0.5")
    pl_nolim.draw_slope(n, color="0.5")
    pl.draw_slope(n, color="r")

gr.clean_axes(ylim=(-0.25, 1.2))

f = plt.gcf()
f.set_size_inches(8.0,2.0)

plt.savefig("generalgrid.pdf")
               


