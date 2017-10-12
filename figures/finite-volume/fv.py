import math
import numpy as np
import matplotlib.pyplot as plt
import grid_plot as gp



# plot a simple finite-difference grid

#-----------------------------------------------------------------------------

nzones = 8

# data that lives on the grid
#a = np.array([0.3, 1.0, 0.9, 0.8, 0.25, 0.15, 0.5, 0.55])
a = np.array([0.3, 1.0, 0.9, 0.8, 0.25, 0.1, 0.5, 0.55])

gr = gp.FVGrid(nzones)

pc = gp.PiecewiseConstant(gr, a)

plt.clf()

gr.draw_grid()

gr.label_center(nzones//2,   r"$i$", fontsize="medium")
gr.label_center(nzones//2-1, r"$i-1$", fontsize="medium")
gr.label_center(nzones//2+1, r"$i+1$", fontsize="medium")
gr.label_center(nzones//2-2, r"$i-2$", fontsize="medium")
gr.label_center(nzones//2+2, r"$i+2$", fontsize="medium")

gr.label_edge(nzones//2,   r"$i-\sfrac{1}{2}$", fontsize="small")
gr.label_edge(nzones//2+1, r"$i+\sfrac{1}{2}$", fontsize="small")

# draw the data
for i in range(nzones):
    pc.draw_cell_avg(i, color="r", filled=True)

pc.label_cell_avg(nzones//2, r"$\langle f\rangle_i$", color="r")

# label dx
gr.label_dx(gr.ng+nzones//2)

gr.clean_axes(ylim=(-0.5, 1.2))

f = plt.gcf()
f.set_size_inches(10.0,3.0)

plt.savefig("fv_grid.pdf")
