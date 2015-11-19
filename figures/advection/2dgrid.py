import math
import numpy as np
import matplotlib.pyplot as plt

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

    xl = (np.arange(2*ng+nzones) - ng)*dx
    xr = (np.arange(2*ng+nzones)+1 - ng)*dx

    xc = 0.5*(xl + xr)

    yl = (np.arange(2*ng+nzones) - ng)*dy
    yr = (np.arange(2*ng+nzones)+1 - ng)*dy

    yc = 0.5*(yl + yr)


    #------------------------------------------------------------------------
    # plot a domain without ghostcells
    
    # x lines
    n = 0
    while (n < nzones):
        plt.plot([xmin-0.25*dx,xmax+0.25*dx], [yl[n],yl[n]], color="k", lw=2)
        n += 1
        
    plt.plot([xmin-0.25*dx,xmax+0.25*dx], [yr[nzones-1],yr[nzones-1]], color="k", lw=2)


    # y lines
    n = 0
    while (n < nzones):
        plt.plot([xl[n],xl[n]], [ymin-0.25*dy,ymax+0.25*dy], color="k", lw=2)
        n += 1
        
    plt.plot([xr[nzones-1],xr[nzones-1]], [ymin-0.25*dy,ymax+0.25*dy], color="k", lw=2)


    #------------------------------------------------------------------------
    # label
    plt.text(xc[nzones/2], yc[nzones/2], r"$a_{i,j}$", fontsize="large",
               horizontalalignment='center', verticalalignment='center')

    plt.text(xc[nzones/2+1], yc[nzones/2], r"$a_{i+1,j}$", fontsize="large",
               horizontalalignment='center', verticalalignment='center')

    plt.text(xc[nzones/2], yc[nzones/2+1], r"$a_{i,j+1}$", fontsize="large",
               horizontalalignment='center', verticalalignment='center')


    # i+1/2,j interface
    plt.scatter(xr[nzones/2]-0.05*dx, yc[nzones/2], marker="x", s=50, color="b")
    plt.text(xr[nzones/2]-0.075*dx, yc[nzones/2], r"$a^{n+1/2}_{i+1/2,j,L}$", 
               fontsize="medium", rotation="270",
               horizontalalignment='right', verticalalignment='center', color="b")

    plt.scatter(xl[nzones/2+1]+0.05*dx, yc[nzones/2], marker="x", s=50, color="b")
    plt.text(xl[nzones/2+1]+0.075*dx, yc[nzones/2], r"$a^{n+1/2}_{i+1/2,j,R}$", 
               fontsize="medium", rotation="270",
               horizontalalignment='left', verticalalignment='center', color="b")


    # i,j+1/2 interface
    plt.scatter(xc[nzones/2], yr[nzones/2]-0.05*dy, marker="x", s=50, color="b")
    plt.text(xc[nzones/2], yr[nzones/2]-0.075*dy, r"$a^{n+1/2}_{i,j+1/2,L}$", 
               fontsize="medium", rotation="0",
               horizontalalignment='center', verticalalignment='top', color="b")

    plt.scatter(xc[nzones/2], yl[nzones/2+1]+0.05*dx, marker="x", s=50, color="b")
    plt.text(xc[nzones/2], yl[nzones/2+1]+0.075*dy, r"$a^{n+1/2}_{i,j+1/2,R}$", 
               fontsize="medium", rotation="0",
               horizontalalignment='center', verticalalignment='bottom', color="b")


    # grid labels
    plt.text(xc[nzones/2-1], yl[0]-0.35*dy, r"$i-1$",
               horizontalalignment='center', fontsize="large")

    plt.text(xc[nzones/2], yl[0]-0.35*dy, r"$i$",
               horizontalalignment='center', fontsize="large")

    plt.text(xc[nzones/2+1], yl[0]-0.35*dy, r"$i+1$",
               horizontalalignment='center', fontsize="large")


    plt.text(xl[0]-0.35*dx, yc[nzones/2-1], r"$j-1$",
               verticalalignment='center', fontsize="large")

    plt.text(xl[0]-0.35*dx, yc[nzones/2], r"$j$",
               verticalalignment='center', fontsize="large")

    plt.text(xl[0]-0.35*dx, yc[nzones/2+1], r"$j+1$",
               verticalalignment='center', fontsize="large")



    # axes
    plt.xlim(xl[0]-0.5*dx,xr[2*ng+nzones-1]+0.25*dx)
    plt.ylim(yl[0]-0.5*dy,yr[2*ng+nzones-1]+0.25*dy)
    plt.axis("off")

    plt.subplots_adjust(left=0.02,right=0.98,bottom=0.02,top=0.98)

    f = plt.gcf()
    f.set_size_inches(6.0,6.0)

    plt.savefig("2dgrid.pdf")

if __name__== "__main__":
    simplegrid()
