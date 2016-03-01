import numpy as np
import matplotlib.pyplot as plt
import grid_plot as gp

# plot a simple finite-difference grid

#-----------------------------------------------------------------------------

nzones = 9

# data that lives on the grid
#a = np.array([0.3, 1.0, 0.9, 0.8, 0.25, 0.15, 0.5, 0.55])
a = np.array([0.55, 0.3, 1.0, 0.9, 0.8, 0.25, 0.1, 0.5, 0.55])

gr = gp.FDGrid(nzones, ng=1)

plt.clf()

gr.draw_grid(draw_ghost=1)

labels = ["-1", "0", "1", "", "i-1", "i", "i+1", "", "N-2", "N-1", "N"]

for i in range(gr.ilo-gr.ng, gr.ng+gr.nx+1):
    if not labels[i] == "":
        gr.label_node(i, r"$%s$" % (labels[i]), fontsize="medium")
    
# draw the data
for i in range(gr.ilo, gr.ihi+1):
    gr.draw_data(i, a[i-gr.ng], color="r")    
    

gr.label_value(gr.ilo+4, a[gr.ilo+4-gr.ng], r"$a_i$", color="r")

# label dx
plt.plot([gr.xc[gr.ng+nzones/2-1], gr.xc[gr.ng+nzones/2-1]], [-0.35,-0.25], color="k")
plt.plot([gr.xc[gr.ng+nzones/2], gr.xc[gr.ng+nzones/2]], [-0.35,-0.25], color="k")
plt.plot([gr.xc[gr.ng+nzones/2-1], gr.xc[gr.ng+nzones/2]], [-0.3,-0.3], color="k")
plt.text(0.5*(gr.xc[gr.ng+nzones/2-1] + gr.xc[gr.ng+nzones/2]), -0.45, 
           r"$\Delta x$", 
           horizontalalignment="center", fontsize=16)


plt.axis([gr.xmin-1.1*gr.dx,gr.xmax+1.1*gr.dx, -0.5, 1.3])
plt.axis("off")

plt.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95)

f = plt.gcf()
f.set_size_inches(10.0,3.0)


plt.savefig("fd_ghost.pdf")

