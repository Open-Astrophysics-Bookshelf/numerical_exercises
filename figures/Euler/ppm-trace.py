import numpy as np
import matplotlib.pyplot as plt
import grid_plot as gp


#-----------------------------------------------------------------------------

nzones = 8

# data that lives on the grid
#a = np.array([0.3, 1.0, 0.9, 0.8, 0.25, 0.15, 0.5, 0.55])
a = np.array([0.3, 1.0, 1.0, 0.8, 0.2, 0.15, 0.5, 0.55])

gr = gp.FVGrid(nzones)


plt.clf()

gr.draw_grid(center_only=1)

gr.label_center(nzones//2-1,   r"$i$")
gr.label_center(nzones//2, r"$i+1$")


ppm = gp.PiecewiseParabolic(gr, a)

for n in range(gr.nx//2-1, gr.nx//2+1):
    ppm.draw_parabola(n, color="r")

nn = gr.nx//2-1
sigma = 0.6
ppm.ppm_trace_left(nn, sigma, color="0.75")

plt.plot([gr.xr[gr.nx//2-1]-sigma*gr.dx, 
          gr.xr[gr.nx//2-1]], [1.2, 1.2], color="k")
plt.plot([gr.xr[gr.nx//2-1]-sigma*gr.dx, 
          gr.xr[gr.nx//2-1]-sigma*gr.dx], [1.15, 1.25], color="k")
plt.plot([gr.xr[gr.nx//2-1], 
          gr.xr[gr.nx//2-1]], [1.15, 1.25], color="k")

plt.text(gr.xr[gr.nx//2-1]-0.5*sigma*gr.dx, 1.3, 
         r"$\sigma_{i}^{(\nu)} \Delta x$", 
         horizontalalignment="center")

plt.text(gr.xr[gr.nx//2-1]-0.5*sigma*gr.dx, 0.4, 
         r"$\mathcal{I}_+^{(\nu)}$",
           color="r", horizontalalignment="center")

plt.axis([gr.xl[gr.nx//2-1]-0.5*gr.dx,gr.xr[gr.nx//2]+0.5*gr.dx, -0.25, 1.4])
plt.axis("off")

plt.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95)

f = plt.gcf()
f.set_size_inches(3.0,2.0)

plt.savefig("ppm-trace.pdf")
plt.savefig("ppm-trace.png")
               


