import matplotlib.pyplot as plt
import grid_plot as gp

# plot two stacked fv grids of different (2x) resolution to show prolongation

#-----------------------------------------------------------------------------

gr = []

nf = 8
gr.append(gp.FVGrid(nf, ng=1, voff=0.0))
gr.append(gp.FVGrid(nf, ng=1, voff=0.0, xmin=0.25, xmax=0.75))
    

plt.clf()

gr[0].draw_grid(emphasize_end=1, draw_ghost=0, edge_ticks=0, color="0.75")
gr[1].draw_grid(emphasize_end=1, draw_ghost=0, edge_ticks=0)

f = plt.gcf()
f.set_size_inches(7.0,1.0)


grf = gr[0]
plt.xlim(grf.xmin-0.75*grf.dx,grf.xmax+0.25*grf.dx)

plt.axis("off")
plt.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95)

plt.savefig("nested1.pdf")

#-----------------------------------------------------------------------------

plt.clf()

gr = []
gr.append(gp.FVGrid(nf, ng=1, voff=0.0))
gr.append(gp.FVGrid(nf, ng=1, voff=4.0))

gr[0].draw_grid(emphasize_end=1, draw_ghost=0, edge_ticks=0)
gr[1].draw_grid(emphasize_end=1, draw_ghost=0, edge_ticks=0)


# labels
plt.text(gr[0].xmin-0.7*gr[0].dx, gr[0].voff+0.5, r"$\phi^{c,n}$", fontsize="medium")
plt.text(gr[1].xmin-0.7*gr[1].dx, gr[1].voff+0.5, r"$\phi^{c,n+1}$", fontsize="medium")

plt.arrow(gr[0].xmin-0.5*gr[0].dx, 0.5+3.0*gr[0].dx, 0.001, gr[1].voff-4.0*gr[0].dx,
            shape='full', head_width=0.025, head_length=0.25, 
            lw=1, width=0.01,
            edgecolor="none", facecolor="r",
            length_includes_head=True, zorder=100)


f = plt.gcf()
f.set_size_inches(7.0,5.0)


grf = gr[0]
plt.xlim(grf.xmin-0.75*grf.dx,grf.xmax+0.25*grf.dx)

plt.axis("off")

plt.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95)

plt.savefig("nested2.pdf")


#-----------------------------------------------------------------------------

plt.clf()

gr1 = []
gr1.append(gp.FVGrid(nf, ng=1, voff=0.0))
gr1.append(gp.FVGrid(nf, ng=1, voff=0.0, xmin=0.25, xmax=0.75))

gr2 = []
gr2.append(gp.FVGrid(nf, ng=1, voff=2.0))
gr2.append(gp.FVGrid(nf, ng=1, voff=2.0, xmin=0.25, xmax=0.75))

gr3 = []
gr3.append(gp.FVGrid(nf, ng=1, voff=4.0))
gr3.append(gp.FVGrid(nf, ng=1, voff=4.0, xmin=0.25, xmax=0.75))


gr1[0].draw_grid(emphasize_end=1, draw_ghost=0, edge_ticks=0, color="0.75")
gr1[1].draw_grid(emphasize_end=1, draw_ghost=0, edge_ticks=0)

gr2[0].draw_grid(emphasize_end=1, draw_ghost=0, edge_ticks=0, color="0.75")
gr2[1].draw_grid(emphasize_end=1, draw_ghost=0, edge_ticks=0)

gr3[0].draw_grid(emphasize_end=1, draw_ghost=0, edge_ticks=0, color="0.75")
gr3[1].draw_grid(emphasize_end=1, draw_ghost=0, edge_ticks=0)

plt.text(gr1[1].xmin-0.95*gr1[0].dx, gr1[1].voff+0.5, r"$\phi^{f,n}$", 
         fontsize="11", horizontalalignment="left")
plt.text(gr2[1].xmin-0.95*gr2[0].dx, gr2[1].voff+0.5, r"$\phi^{f,n+1/2}$", 
         fontsize="11", horizontalalignment="left")
plt.text(gr3[1].xmin-0.95*gr3[0].dx, gr3[1].voff+0.5, r"$\phi^{f,n+1}$", 
         fontsize="11", horizontalalignment="left")


plt.annotate(r"Dirichlet BCs $\phi^{c,n}$", xy=(gr1[1].xmin, gr1[1].voff), xycoords='data',
               xytext=(30,-25), textcoords="offset points",
               arrowprops=dict(arrowstyle="->", 
                               connectionstyle="angle,angleA=180,angleB=-90,rad=10", 
                               ec="0.5"),
               fontsize="small",color="0.5")

plt.annotate(r"Dirichlet BCs $(\phi^{c,n} + \phi^{c,n+1})/2$", xy=(gr2[1].xmin, gr2[1].voff), xycoords='data',
               xytext=(30,-25), textcoords="offset points",
               arrowprops=dict(arrowstyle="->", 
                               connectionstyle="angle,angleA=180,angleB=-90,rad=10", 
                               ec="0.5"),
               fontsize="small",color="0.5")

plt.annotate(r"Dirichlet BCs $\phi^{c,n+1}$", xy=(gr3[1].xmin, gr3[1].voff), xycoords='data',
               xytext=(30,-25), textcoords="offset points",
               arrowprops=dict(arrowstyle="->", 
                               connectionstyle="angle,angleA=180,angleB=-90,rad=10", 
                               ec="0.5"),
               fontsize="small",color="0.5")


plt.arrow(gr1[1].xmin-0.75*gr1[0].dx, 0.5+3.0*gr1[0].dx, 0.001, gr2[0].voff-4.0*gr2[0].dx,
            shape='full', head_width=0.025, head_length=0.25, 
            lw=1, width=0.01,
            edgecolor="none", facecolor="r",
            length_includes_head=True, zorder=100)

plt.arrow(gr2[1].xmin-0.75*gr2[0].dx, gr2[1].voff+0.5+3.0*gr2[0].dx, 0.001, gr3[1].voff-gr2[1].voff-4.0*gr3[0].dx,
            shape='full', head_width=0.025, head_length=0.25, 
            lw=1, width=0.01,
            edgecolor="none", facecolor="r",
            length_includes_head=True, zorder=100)

f = plt.gcf()
f.set_size_inches(7.0,5.0)


grf = gr[0]
plt.xlim(grf.xmin-0.75*grf.dx,grf.xmax+0.25*grf.dx)

plt.axis("off")

plt.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95)

plt.savefig("nested3.pdf")



