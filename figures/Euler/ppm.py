import math
import numpy as np
import matplotlib.pyplot as plt
import grid_plot as gp


#-----------------------------------------------------------------------------

nzones = 8

# data that lives on the grid
a = np.array([0.3, 1.0, 0.9, 0.8, 0.25, 0.1, 0.5, 0.55])

gr = gp.FVGrid(nzones)

plt.clf()

gr.draw_grid()

gr.label_center(nzones//2,   r"$i$")
gr.label_center(nzones//2-1, r"$i-1$")
gr.label_center(nzones//2+1, r"$i+1$")
gr.label_center(nzones//2-2, r"$i-2$")
gr.label_center(nzones//2+2, r"$i+2$")


pc = gp.PiecewiseConstant(gr, a)

for n in range(nzones):
    pc.draw_cell_avg(n, color="r")


gr.clean_axes()
plt.ylim(-0.25, 1.2)

plt.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95)

f = plt.gcf()
f.set_size_inches(8.0,2.0)

plt.savefig("piecewise-constant.pdf")
plt.savefig("piecewise-constant.png")


#------------- PLM -------------
plt.clf()

gr.draw_grid()

gr.label_center(nzones//2,   r"$i$")
gr.label_center(nzones//2-1, r"$i-1$")
gr.label_center(nzones//2+1, r"$i+1$")
gr.label_center(nzones//2-2, r"$i-2$")
gr.label_center(nzones//2+2, r"$i+2$")


# not limited and limited
pl_n = gp.PiecewiseLinear(gr, a, nolimit=1)
pl_y = gp.PiecewiseLinear(gr, a)

for n in range(nzones):
    pc.draw_cell_avg(n, color="0.5")

for n in range(2, nzones-2):
    pl_n.draw_slope(n, color="r", ls=":")
    pl_y.draw_slope(n, color="r")


gr.clean_axes()
plt.ylim(-0.25, 1.2)

plt.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95)

f = plt.gcf()
f.set_size_inches(8.0,2.0)

plt.savefig("piecewise-linear.pdf")
plt.savefig("piecewise-linear.png")


#------------- PPM -------------
plt.clf()

gr.draw_grid()

gr.label_center(nzones//2,   r"$i$")
gr.label_center(nzones//2-1, r"$i-1$")
gr.label_center(nzones//2+1, r"$i+1$")
gr.label_center(nzones//2-2, r"$i-2$")
gr.label_center(nzones//2+2, r"$i+2$")

for n in range(nzones):
    pc.draw_cell_avg(n, color="0.5")


# compute the parabolic coefficients -- limited and not limited
pm_n = gp.PiecewiseParabolic(gr, a, nolimit=1)
pm_y = gp.PiecewiseParabolic(gr, a)

for n in range(2, nzones-2):
    pm_n.draw_parabola(n, color="r", ls=":")
    pm_y.draw_parabola(n, color="r")


gr.clean_axes()
plt.ylim(-0.25, 1.2)

plt.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95)

f = plt.gcf()
f.set_size_inches(8.0,2.0)


plt.savefig("piecewise-parabolic.pdf")
plt.savefig("piecewise-parabolic.png")

               


