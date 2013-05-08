import math
import numpy
import pylab
import grid_plot_util as gpu

# plot two stacked fv grids of different (2x) resolution to show prolongation

#-----------------------------------------------------------------------------

gr = []

nf = 8
gr.append(gpu.grid(nf, ng=1, voff=0.0))
gr.append(gpu.grid(nf, ng=1, voff=0.0, xmin=0.25, xmax=0.75))
    

pylab.clf()

gpu.drawGrid(gr[0], emphasizeEnd=1, drawGhost=0, edgeTicks=0, color="0.75")
gpu.drawGrid(gr[1], emphasizeEnd=1, drawGhost=0, edgeTicks=0)

f = pylab.gcf()
f.set_size_inches(7.0,1.0)


grf = gr[0]
pylab.xlim(grf.xmin-0.75*grf.dx,grf.xmax+0.25*grf.dx)

pylab.axis("off")

pylab.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95)


pylab.savefig("nested1.png")
pylab.savefig("nested1.eps")

#-----------------------------------------------------------------------------

pylab.clf()

gr = []
gr.append(gpu.grid(nf, ng=1, voff=0.0))
gr.append(gpu.grid(nf, ng=1, voff=4.0))


gpu.drawGrid(gr[0], emphasizeEnd=1, drawGhost=0, edgeTicks=0)
gpu.drawGrid(gr[1], emphasizeEnd=1, drawGhost=0, edgeTicks=0)


# labels
pylab.text(gr[0].xmin-0.7*gr[0].dx, gr[0].voff+0.5, r"$\phi^{c,n}$", fontsize="medium")
pylab.text(gr[1].xmin-0.7*gr[1].dx, gr[1].voff+0.5, r"$\phi^{c,n+1}$", fontsize="medium")

pylab.arrow(gr[0].xmin-0.5*gr[0].dx, 0.5+3.0*gr[0].dx, 0.001, gr[1].voff-4.0*gr[0].dx,
            shape='full', head_width=0.025, head_length=0.25, 
            lw=1, width=0.01,
            edgecolor="none", facecolor="r",
            length_includes_head=True, zorder=100)


f = pylab.gcf()
f.set_size_inches(7.0,5.0)


grf = gr[0]
pylab.xlim(grf.xmin-0.75*grf.dx,grf.xmax+0.25*grf.dx)

pylab.axis("off")

pylab.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95)


pylab.savefig("nested2.png")
pylab.savefig("nested2.eps")



#-----------------------------------------------------------------------------

pylab.clf()

gr1 = []
gr1.append(gpu.grid(nf, ng=1, voff=0.0))
gr1.append(gpu.grid(nf, ng=1, voff=0.0, xmin=0.25, xmax=0.75))

gr2 = []
gr2.append(gpu.grid(nf, ng=1, voff=2.0))
gr2.append(gpu.grid(nf, ng=1, voff=2.0, xmin=0.25, xmax=0.75))

gr3 = []
gr3.append(gpu.grid(nf, ng=1, voff=4.0))
gr3.append(gpu.grid(nf, ng=1, voff=4.0, xmin=0.25, xmax=0.75))


gpu.drawGrid(gr1[0], emphasizeEnd=1, drawGhost=0, edgeTicks=0, color="0.75")
gpu.drawGrid(gr1[1], emphasizeEnd=1, drawGhost=0, edgeTicks=0)

gpu.drawGrid(gr2[0], emphasizeEnd=1, drawGhost=0, edgeTicks=0, color="0.75")
gpu.drawGrid(gr2[1], emphasizeEnd=1, drawGhost=0, edgeTicks=0)

gpu.drawGrid(gr3[0], emphasizeEnd=1, drawGhost=0, edgeTicks=0, color="0.75")
gpu.drawGrid(gr3[1], emphasizeEnd=1, drawGhost=0, edgeTicks=0)

pylab.text(gr1[1].xmin-0.95*gr1[0].dx, gr1[1].voff+0.5, r"$\phi^{f,n}$", fontsize="11", horizontalalignment="left")
pylab.text(gr2[1].xmin-0.95*gr2[0].dx, gr2[1].voff+0.5, r"$\phi^{f,n+1/2}$", fontsize="11", horizontalalignment="left")
pylab.text(gr3[1].xmin-0.95*gr3[0].dx, gr3[1].voff+0.5, r"$\phi^{f,n+1}$", fontsize="11", horizontalalignment="left")


pylab.annotate(r"Dirichlet BCs $\phi^{c,n}$", xy=(gr1[1].xmin, gr1[1].voff), xycoords='data',
               xytext=(30,-25), textcoords="offset points",
               arrowprops=dict(arrowstyle="->", connectionstyle="angle,angleA=180,angleB=-90,rad=10", ec="0.5"),
               fontsize="small",color="0.5")

pylab.annotate(r"Dirichlet BCs $(\phi^{c,n} + \phi^{c,n+1})/2$", xy=(gr2[1].xmin, gr2[1].voff), xycoords='data',
               xytext=(30,-25), textcoords="offset points",
               arrowprops=dict(arrowstyle="->", connectionstyle="angle,angleA=180,angleB=-90,rad=10", ec="0.5"),
               fontsize="small",color="0.5")

pylab.annotate(r"Dirichlet BCs $\phi^{c,n+1}$", xy=(gr3[1].xmin, gr3[1].voff), xycoords='data',
               xytext=(30,-25), textcoords="offset points",
               arrowprops=dict(arrowstyle="->", connectionstyle="angle,angleA=180,angleB=-90,rad=10", ec="0.5"),
               fontsize="small",color="0.5")


pylab.arrow(gr1[1].xmin-0.75*gr1[0].dx, 0.5+3.0*gr1[0].dx, 0.001, gr2[0].voff-4.0*gr2[0].dx,
            shape='full', head_width=0.025, head_length=0.25, 
            lw=1, width=0.01,
            edgecolor="none", facecolor="r",
            length_includes_head=True, zorder=100)

pylab.arrow(gr2[1].xmin-0.75*gr2[0].dx, gr2[1].voff+0.5+3.0*gr2[0].dx, 0.001, gr3[1].voff-gr2[1].voff-4.0*gr3[0].dx,
            shape='full', head_width=0.025, head_length=0.25, 
            lw=1, width=0.01,
            edgecolor="none", facecolor="r",
            length_includes_head=True, zorder=100)




f = pylab.gcf()
f.set_size_inches(7.0,5.0)


grf = gr[0]
pylab.xlim(grf.xmin-0.75*grf.dx,grf.xmax+0.25*grf.dx)

pylab.axis("off")

pylab.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95)


pylab.savefig("nested3.png")
pylab.savefig("nested3.eps")



