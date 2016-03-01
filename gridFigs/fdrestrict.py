import numpy as np
import matplotlib.pyplot as plt
import grid_plot as gp

# plot two stacked fv grids of different (2x) resolution to show prolongation

#-----------------------------------------------------------------------------

nf = 5
nc = 3

grf = gp.FDGrid(nf, voff=2.0)
grc = gp.FDGrid(nc)

plt.clf()

grf.draw_grid()
grc.draw_grid()

grf.label_node(nf/2-2, r"$i-2$")
grf.label_node(nf/2-1, r"$i-1$")
grf.label_node(nf/2,   r"$i$")
grf.label_node(nf/2+1, r"$i+1$")
grf.label_node(nf/2+2, r"$i+2$")


grc.label_node(nc/2-1, r"$j-1$")
grc.label_node(nc/2,   r"$j$")
grc.label_node(nc/2+1, r"$j+1$")

grf.label_node_data(nf/2-2, r"$\phi_{i-2}^h$")
grf.label_node_data(nf/2-1, r"$\phi_{i-1}^h$")
grf.label_node_data(nf/2,   r"$\phi_i^h$")
grf.label_node_data(nf/2+1, r"$\phi_{i+1}^h$")
grf.label_node_data(nf/2+2, r"$\phi_{i+2}^h$")

grc.label_node_data(nc/2-1, r"$\phi_{j-1}^{2h}$")
grc.label_node_data(nc/2,   r"$\phi_{j}^{2h}$")
grc.label_node_data(nc/2+1, r"$\phi_{j+1}^{2h}$")
    

# connect the dots...

plt.plot([grf.xc[nf/2-2], grf.xc[nf/2-2]], [-0.25, 3.25], ":", color="0.5")
plt.plot([grf.xc[nf/2],   grf.xc[nf/2]],   [-0.25, 3.25], ":", color="0.5")
plt.plot([grf.xc[nf/2+2], grf.xc[nf/2+2]], [-0.25, 3.25], ":", color="0.5")


plt.axis([grf.xmin-0.5*grf.dx,grf.xmax+0.5*grf.dx, -0.5, 3.5])
plt.axis("off")

plt.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95)

f = plt.gcf()
f.set_size_inches(6.0,5.0)

plt.savefig("fdrestrict.pdf")


