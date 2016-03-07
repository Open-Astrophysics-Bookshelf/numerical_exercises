import matplotlib.pyplot as plt
import numpy as np
import grid_plot as gp

# plot a simple finite-difference grid showing the boundary data
# and a finite-volume grid with the same data, showing it's boundary

#-----------------------------------------------------------------------------


# finite-difference

# data that lives on the grid
#a = np.array([0.3, 1.0, 0.9, 0.8, 0.25, 0.15, 0.5, 0.55])
a = np.array([0.3, 1.0, 0.9, 0.8, 0.25, 0.1, 0.2, 0.5, 0.55])
nzones = len(a)

voff = 2.0
gr = gp.FDGrid(nzones, voff=voff)

gr.draw_grid(emphasize_end=True)

gr.label_node(nzones/2,   r"$i$",   fontsize="large")
gr.label_node(nzones/2-1, r"$i-1$", fontsize="large")
gr.label_node(nzones/2+1, r"$i+1$", fontsize="large")
gr.label_node(nzones/2-2, r"$i-2$", fontsize="large")
gr.label_node(nzones/2+2, r"$i+2$", fontsize="large")

# draw the data
for i in range(nzones):
    gr.draw_data(i, a[i], color="r")


gr.label_value(nzones/2, a[nzones/2], r"$f_i$", color="r")

# label dx
plt.plot([gr.xc[gr.ng+nzones/2-1], gr.xc[gr.ng+nzones/2-1]], [-0.35+voff,-0.25+voff], color="k")
plt.plot([gr.xc[gr.ng+nzones/2], gr.xc[gr.ng+nzones/2]], [-0.35+voff,-0.25+voff], color="k")
plt.plot([gr.xc[gr.ng+nzones/2-1], gr.xc[gr.ng+nzones/2]], [-0.3+voff,-0.3+voff], color="k")
plt.text(0.5*(gr.xc[gr.ng+nzones/2-1] + gr.xc[gr.ng+nzones/2]), -0.45+voff, r"$\Delta x$",
           horizontalalignment="center")



# finite-volume
av = 0.5*(a[0:nzones-1] + a[1:])
nzones = len(av)
ng = 1

gr = gp.FVGrid(nzones, ng=1)

gr.draw_grid(emphasize_end=True, draw_ghost=True)

gr.label_center(ng+nzones/2,   r"$i$",   fontsize="large")
gr.label_center(ng+nzones/2-1, r"$i-1$", fontsize="large")
gr.label_center(ng+nzones/2+1, r"$i+1$", fontsize="large")

gr.label_center(ng+nzones-1, r"$\mathrm{hi}$", fontsize="large")
gr.label_center(ng+nzones, r"$\mathrm{hi+1}$", fontsize="large")

gr.label_center(gr.ilo, r"$\mathrm{lo}$", fontsize="large")
gr.label_center(gr.ilo-1, r"$\mathrm{lo-1}$", fontsize="large")

a = gr.scratch_array()
a[gr.ilo:gr.ihi+1] = av

cc = gp.CellCentered(gr, a)

# draw the data
for i in range(gr.ilo, gr.ihi+1):
    cc.draw_data_point(i, color="r")

cc.label_data_point(ng+nzones/2, r"$f_i$", color="r")

# label dx
plt.plot([gr.xl[gr.ng+nzones/2], gr.xl[gr.ng+nzones/2]], [-0.35,-0.25], color="k")
plt.plot([gr.xr[gr.ng+nzones/2], gr.xr[gr.ng+nzones/2]], [-0.35,-0.25], color="k")
plt.plot([gr.xl[gr.ng+nzones/2], gr.xr[gr.ng+nzones/2]], [-0.3,-0.3], color="k")
plt.text(gr.xc[gr.ng+nzones/2], -0.45, r"$\Delta x$",
           horizontalalignment="center")



# illustrate the boundaries
plt.plot([gr.xl[gr.ilo], gr.xl[gr.ilo]], [-0.5, 2.0], ls=":", color="0.5")
plt.plot([gr.xr[gr.ihi], gr.xr[gr.ihi]], [-0.5, 2.0], ls=":", color="0.5")

plt.text(gr.xl[gr.ilo], -0.5, "left BC", horizontalalignment="center",
         verticalalignment="top")

plt.text(gr.xr[gr.ihi], -0.5, "right BC", horizontalalignment="center",
         verticalalignment="top")

plt.axis([gr.xmin-gr.dx,gr.xmax+gr.dx, -0.5, 3.2])
plt.axis("off")

plt.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95)

f = plt.gcf()
f.set_size_inches(10.0,6.0)

plt.savefig("fv-fd_grid_bc.pdf")
