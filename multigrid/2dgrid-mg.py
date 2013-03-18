import math
import numpy
import pylab

def simplegrid():

    #-------------------------------------------------------------------------
    # prolongation
    #-------------------------------------------------------------------------

    # grid info
    xmin = 0.0
    xmax = 1.0

    ymin = 0.0
    ymax = 1.0

    nzones = 3
    ng = 0
    
    dx = (xmax - xmin)/float(nzones)
    dy = (ymax - ymin)/float(nzones)

    xl = (numpy.arange(2*ng+nzones) - ng)*dx
    xr = (numpy.arange(2*ng+nzones)+1 - ng)*dx

    xc = 0.5*(xl + xr)

    yl = (numpy.arange(2*ng+nzones) - ng)*dy
    yr = (numpy.arange(2*ng+nzones)+1 - ng)*dy

    yc = 0.5*(yl + yr)


    #------------------------------------------------------------------------
    # plot a domain without ghostcells
    
    # x lines
    n = 0
    while (n < nzones):
        pylab.plot([xmin-0.25*dx,xmax+0.25*dx], [yl[n],yl[n]], color="k", lw=2)
        n += 1
        
    pylab.plot([xmin-0.25*dx,xmax+0.25*dx], [yr[nzones-1],yr[nzones-1]], color="k", lw=2)


    # y lines
    n = 0
    while (n < nzones):
        pylab.plot([xl[n],xl[n]], [ymin-0.25*dy,ymax+0.25*dy], color="k", lw=2)
        n += 1
        
    pylab.plot([xr[nzones-1],xr[nzones-1]], [ymin-0.25*dy,ymax+0.25*dy], color="k", lw=2)


    #------------------------------------------------------------------------
    # label
    pylab.text(xc[nzones/2], yc[nzones/2], r"$\phi_{i,j}^c$", fontsize="18",
               horizontalalignment='center', verticalalignment='center', 
               zorder=100, color="b")

    pylab.text(xc[nzones/2+1], yc[nzones/2], r"$\phi_{i+1,j}^c$", fontsize="18",
               horizontalalignment='center', verticalalignment='center',
               color="b")

    pylab.text(xc[nzones/2], yc[nzones/2+1], r"$\phi_{i,j+1}^c$", fontsize="18",
               horizontalalignment='center', verticalalignment='center',
               color="b")

    pylab.text(xc[nzones/2-1], yc[nzones/2], r"$\phi_{i-1,j}^c$", fontsize="18",
               horizontalalignment='center', verticalalignment='center',
               color="b")

    pylab.text(xc[nzones/2], yc[nzones/2-1], r"$\phi_{i,j-1}^c$", fontsize="18",
               horizontalalignment='center', verticalalignment='center',
               color="b")


    # shading
    ii = nzones/2; jj = nzones/2
    pylab.fill([xl[ii], xl[ii], xr[ii], xr[ii], xl[ii]],
               [yl[jj], yr[jj], yr[jj], yl[jj], yl[jj]], color="0.85", 
               zorder=-1)

    ii = nzones/2+1; jj = nzones/2
    pylab.fill([xl[ii], xl[ii], xr[ii], xr[ii], xl[ii]],
               [yl[jj], yr[jj], yr[jj], yl[jj], yl[jj]], color="0.85", 
               zorder=-1)

    ii = nzones/2-1; jj = nzones/2
    pylab.fill([xl[ii], xl[ii], xr[ii], xr[ii], xl[ii]],
               [yl[jj], yr[jj], yr[jj], yl[jj], yl[jj]], color="0.85", 
               zorder=-1)

    ii = nzones/2; jj = nzones/2+1
    pylab.fill([xl[ii], xl[ii], xr[ii], xr[ii], xl[ii]],
               [yl[jj], yr[jj], yr[jj], yl[jj], yl[jj]], color="0.85", 
               zorder=-1)

    ii = nzones/2; jj = nzones/2-1
    pylab.fill([xl[ii], xl[ii], xr[ii], xr[ii], xl[ii]],
               [yl[jj], yr[jj], yr[jj], yl[jj], yl[jj]], color="0.85", 
               zorder=-1)


    # fine cells
    ii = nzones/2; jj = nzones/2
    pylab.plot([xc[ii], xc[ii]], [yl[ii], yr[ii]], linestyle="--", color="0.3")
    pylab.plot([xl[ii], xr[ii]], [yc[ii], yc[ii]], linestyle="--", color="0.3")


    pylab.text(xc[ii]-dx/4, yc[jj]-dy/4, r"$\phi_{--}^{f}$", fontsize="18",
               horizontalalignment='center', verticalalignment='center', 
               zorder=100, color="r")

    pylab.text(xc[ii]-dx/4, yc[jj]+dy/4, r"$\phi_{-+}^{f}$", fontsize="18",
               horizontalalignment='center', verticalalignment='center', 
               zorder=100, color="r")

    pylab.text(xc[ii]+dx/4, yc[jj]-dy/4, r"$\phi_{+-}^{f}$", fontsize="18",
               horizontalalignment='center', verticalalignment='center', 
               zorder=100, color="r")

    pylab.text(xc[ii]+dx/4, yc[jj]+dy/4, r"$\phi_{++}^{f}$", fontsize="18",
               horizontalalignment='center', verticalalignment='center', 
               zorder=100, color="r")


    # grid labels
    pylab.text(xc[nzones/2-1], yl[0]-0.35*dy, r"$i-1$",
               horizontalalignment='center', fontsize="16")

    pylab.text(xc[nzones/2], yl[0]-0.35*dy, r"$i$",
               horizontalalignment='center', fontsize="16")

    pylab.text(xc[nzones/2+1], yl[0]-0.35*dy, r"$i+1$",
               horizontalalignment='center', fontsize="16")


    pylab.text(xl[0]-0.35*dx, yc[nzones/2-1], r"$j-1$",
               verticalalignment='center', fontsize="16")

    pylab.text(xl[0]-0.35*dx, yc[nzones/2], r"$j$",
               verticalalignment='center', fontsize="16")

    pylab.text(xl[0]-0.35*dx, yc[nzones/2+1], r"$j+1$",
               verticalalignment='center', fontsize="16")



    # axes
    pylab.xlim(xl[0]-0.5*dx,xr[2*ng+nzones-1]+0.25*dx)
    pylab.ylim(yl[0]-0.5*dy,yr[2*ng+nzones-1]+0.25*dy)
    pylab.axis("off")

    pylab.subplots_adjust(left=0.02,right=0.98,bottom=0.02,top=0.98)

    f = pylab.gcf()
    f.set_size_inches(8.0,8.0)


    pylab.savefig("2dgrid-prolong.png")
    pylab.savefig("2dgrid-prolong.eps")


               

if __name__== "__main__":
    simplegrid()
