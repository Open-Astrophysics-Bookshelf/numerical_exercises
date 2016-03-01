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

nzones = 7
ng = 2

# data that lives on the grid
#a = np.array([0.3, 1.0, 0.9, 0.8, 0.25, 0.15, 0.5, 0.55])
#a = np.array([0.3, 1.0, 0.9, 0.8, 0.25, 0.1, 0.5, 0.55])
a = np.array([1.0, 0.9, 0.8, 0.25, 0.1, 0.5, 0.55])

gr = gp.FVGrid(nzones, ng)

aa = gr.scratch_array()
aa[gr.ilo:gr.ihi+1] = a

cc = gp.PiecewiseConstant(gr, aa)

plt.clf()

gr.draw_grid(draw_ghost=1, emphasize_end=1)

gr.label_center(ng+nzones/2,   r"$i$")
gr.label_center(ng+nzones/2-1, r"$i-1$")
gr.label_center(ng+nzones/2+1, r"$i+1$")

gr.label_center(gr.ilo, r"$\mathrm{lo}$")
gr.label_center(gr.ilo-1, r"$\mathrm{lo-1}$")
gr.label_center(gr.ilo-2, r"$\mathrm{lo-2}$")

gr.label_center(gr.ihi, r"$\mathrm{hi}$")
gr.label_center(gr.ihi+1, r"$\mathrm{hi+1}$")
gr.label_center(gr.ihi+2, r"$\mathrm{hi+2}$")

gr.label_edge(ng+nzones/2,   r"$i-\sfrac{1}{2}$")
gr.label_edge(ng+nzones/2+1,   r"$i+\sfrac{1}{2}$")

# draw the data
for i in range(nzones):
    cc.draw_cell_avg(ng+i, color="r")    
    
cc.label_cell_avg(ng+nzones/2, r"$\langle a\rangle_i$", color="r")

# label dx
plt.plot([gr.xr[gr.ng+nzones/2-1], gr.xr[gr.ng+nzones/2-1]], [-0.35,-0.25], color="k")
plt.plot([gr.xr[gr.ng+nzones/2], gr.xr[gr.ng+nzones/2]], [-0.35,-0.25], color="k")
plt.plot([gr.xr[gr.ng+nzones/2-1], gr.xr[gr.ng+nzones/2]], [-0.3,-0.3], color="k")
plt.text(gr.xc[gr.ng+nzones/2], -0.55, r"$\Delta x$", 
           horizontalalignment="center", fontsize=16)



plt.axis([gr.xmin-2.02*gr.dx,gr.xmax+2.02*gr.dx, -0.5, 1.6])
plt.axis("off")

plt.subplots_adjust(left=0.025,right=0.975,bottom=0.05,top=0.95)

f = plt.gcf()
f.set_size_inches(10.0,2.5)

plt.savefig("fv_ghost.pdf")
