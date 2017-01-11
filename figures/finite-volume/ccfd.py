import math
import numpy as np

import matplotlib
# Use LaTeX for rendering
matplotlib.rcParams["text.usetex"] = True
# load the xfrac package
matplotlib.rcParams["text.latex.preamble"].append(r'\usepackage{xfrac}')

import matplotlib.pyplot as plt
import grid_plot as gp

# plot a simple finite-difference grid

#-----------------------------------------------------------------------------

nzones = 8

# data that lives on the grid
#a = np.array([0.3, 1.0, 0.9, 0.8, 0.25, 0.15, 0.5, 0.55])
a = np.array([0.3, 1.0, 0.9, 0.8, 0.25, 0.1, 0.5, 0.55])

gr = gp.FVGrid(nzones)


plt.clf()

gr.draw_grid()

gr.label_center(nzones/2,   r"$i$", fontsize="medium")
gr.label_center(nzones/2-1, r"$i-1$", fontsize="medium")
gr.label_center(nzones/2+1, r"$i+1$", fontsize="medium")
gr.label_center(nzones/2-2, r"$i-2$", fontsize="medium")
gr.label_center(nzones/2+2, r"$i+2$", fontsize="medium")


cc = gp.CellCentered(gr, a)

# draw the data
for i in range(nzones):
    cc.draw_data_point(i, color="r")


cc.label_data_point(nzones/2, r"$f_i$", color="r")

# label dx
plt.plot([gr.xc[gr.ng+nzones/2-1], gr.xc[gr.ng+nzones/2-1]], [-0.35,-0.25], color="k")
plt.plot([gr.xc[gr.ng+nzones/2], gr.xc[gr.ng+nzones/2]], [-0.35,-0.25], color="k")
plt.plot([gr.xc[gr.ng+nzones/2-1], gr.xc[gr.ng+nzones/2]], [-0.3,-0.3], color="k")
plt.text(0.5*(gr.xc[gr.ng+nzones/2-1] + gr.xc[gr.ng+nzones/2]), -0.45, r"$\Delta x$",
           horizontalalignment="center")


gr.clean_axes()
plt.ylim(-0.5, 1.2)


plt.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95)

f = plt.gcf()
f.set_size_inches(10.0,3.0)

plt.savefig("ccfd_grid.pdf")
