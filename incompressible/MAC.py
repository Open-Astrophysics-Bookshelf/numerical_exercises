import math
import numpy
import pylab

def simplegrid():

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
    pylab.text(xc[nzones/2], yc[nzones/2]+0.05*dy, r"$u_{i,j}$", fontsize="15",
               horizontalalignment='center', verticalalignment='center')

    pylab.text(xc[nzones/2], yc[nzones/2]-0.05*dy, r"$v_{i,j}$", fontsize="15",
               horizontalalignment='center', verticalalignment='center')


    pylab.text(xc[nzones/2+1], yc[nzones/2]+0.05*dy, r"$u_{i+1,j}$", 
               fontsize="15",
               horizontalalignment='center', verticalalignment='center')

    pylab.text(xc[nzones/2+1], yc[nzones/2]-0.05*dy, r"$v_{i+1,j}$", 
               fontsize="15",
               horizontalalignment='center', verticalalignment='center')


    pylab.text(xc[nzones/2-1], yc[nzones/2]+0.05*dy, r"$u_{i-1,j}$", 
               fontsize="15",
               horizontalalignment='center', verticalalignment='center')

    pylab.text(xc[nzones/2-1], yc[nzones/2]-0.05*dy, r"$v_{i-1,j}$", 
               fontsize="15",
               horizontalalignment='center', verticalalignment='center')


    pylab.text(xc[nzones/2], yc[nzones/2+1]+0.05*dy, r"$u_{i,j+1}$", 
               fontsize="15",
               horizontalalignment='center', verticalalignment='center')

    pylab.text(xc[nzones/2], yc[nzones/2+1]-0.05*dy, r"$v_{i,j+1}$", 
               fontsize="15",
               horizontalalignment='center', verticalalignment='center')


    pylab.text(xc[nzones/2], yc[nzones/2-1]+0.05*dy, r"$u_{i,j-1}$", 
               fontsize="15",
               horizontalalignment='center', verticalalignment='center')

    pylab.text(xc[nzones/2], yc[nzones/2-1]-0.05*dy, r"$v_{i,j-1}$", 
               fontsize="15",
               horizontalalignment='center', verticalalignment='center')


    # i+1/2,j interface
    pylab.scatter(xr[nzones/2], yc[nzones/2], marker="x", s=50)
    pylab.text(xr[nzones/2]+0.05*dx, yc[nzones/2], 
               r"${u}^\mathrm{adv}_{i+1/2,j}$", 
               fontsize="15", rotation=270,
               horizontalalignment='left', verticalalignment='center')

    # i-1/2,j interface
    pylab.scatter(xr[nzones/2-1], yc[nzones/2], marker="x", s=50)
    pylab.text(xr[nzones/2-1]-0.05*dx, yc[nzones/2], 
               r"${u}^\mathrm{adv}_{i-1/2,j}$", 
               fontsize="15", rotation=270, 
               horizontalalignment='right', verticalalignment='center')



    # i,j+1/2 interface
    pylab.scatter(xc[nzones/2], yr[nzones/2], marker="x", s=50)
    pylab.text(xc[nzones/2], yr[nzones/2]+0.05*dy, 
               r"${v}^\mathrm{adv}_{i,j+1/2}$", 
               fontsize="15", 
               horizontalalignment='center', verticalalignment='bottom')


    # i,j-1/2 interface
    pylab.scatter(xc[nzones/2], yr[nzones/2-1], marker="x", s=50)
    pylab.text(xc[nzones/2], yr[nzones/2-1]-0.05*dy, 
               r"${v}^\mathrm{adv}_{i,j-1/2}$", 
               fontsize="15", 
               horizontalalignment='center', verticalalignment='top')




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


    pylab.savefig("MAC.png")
    pylab.savefig("MAC.eps")
               

if __name__== "__main__":
    simplegrid()
