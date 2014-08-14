import math
import numpy
import pylab



class Grid2d:

    def __init__(self, nx, ny,
                 xmin=0.0, xmax=1.0, ymin=0.0, ymax=1.0):
                 
        self.xmin = xmin
        self.xmax = xmax

        self.ymin = ymin
        self.ymax = ymax

        self.nx = nx
        self.ny = ny
    
        self.dx = (self.xmax - self.xmin)/float(self.nx)
        self.dy = (self.ymax - self.ymin)/float(self.ny)

        self.xl = (numpy.arange(self.nx)  )*self.dx
        self.xr = (numpy.arange(self.nx)+1)*self.dx

        self.xc = 0.5*(self.xl + self.xr)

        self.yl = (numpy.arange(self.ny)  )*self.dy
        self.yr = (numpy.arange(self.ny)+1)*self.dy

        self.yc = 0.5*(self.yl + self.yr)


    def draw_grid(self):
        """ plot a domain without ghostcells """
    
        # lines parallel to the x-axis
        n = 0
        while n < self.ny:
            pylab.plot([self.xmin-0.25*self.dx,self.xmax+0.25*self.dx], 
                       [self.yl[n],self.yl[n]], color="k", lw=2)
            n += 1
        
        pylab.plot([self.xmin-0.25*self.dx,self.xmax+0.25*self.dx], 
                   [self.yr[self.ny-1],self.yr[self.ny-1]], color="k", lw=2)


        # lines parallel to the y-axis
        n = 0
        while n < self.nx:
            pylab.plot([self.xl[n],self.xl[n]], 
                       [self.ymin-0.25*self.dy,self.ymax+0.25*self.dy], color="k", lw=2)
            n += 1
        
        pylab.plot([self.xr[self.nx-1],self.xr[self.nx-1]], 
                   [self.ymin-0.25*self.dy,self.ymax+0.25*self.dy], color="k", lw=2)


        # grid labels
        pylab.text(self.xc[self.nx/2], self.yl[0]-0.35*self.dy, r"$i$",
                   horizontalalignment='center', fontsize="16")

        if self.nx/2-1 >= 0:
            pylab.text(self.xc[self.nx/2-1], self.yl[0]-0.35*self.dy, r"$i-1$",
                       horizontalalignment='center', fontsize="16")
        
        if self.nx/2+1 < self.nx:
            pylab.text(self.xc[self.nx/2+1], self.yl[0]-0.35*self.dy, r"$i+1$",
                       horizontalalignment='center', fontsize="16")


        pylab.text(self.xl[0]-0.35*self.dx, self.yc[self.ny/2], r"$j$",
                   verticalalignment='center', fontsize="16")

        if self.ny/2-1 >= 0:
            pylab.text(self.xl[0]-0.35*self.dx, self.yc[self.ny/2-1], r"$j-1$",
                       verticalalignment='center', fontsize="16")

        if self.ny/2+1 < self.ny:
            pylab.text(self.xl[0]-0.35*self.dx, self.yc[self.ny/2+1], r"$j+1$",
                       verticalalignment='center', fontsize="16")



def simplegrid():

    g = Grid2d(3,3)

    pylab.subplot(121)

    g.draw_grid()

    #------------------------------------------------------------------------
    # label
    pylab.text(g.xc[g.nx/2], g.yc[g.ny/2], 
               r"$\phi_{i,j}$", fontsize="15",
               horizontalalignment='center', verticalalignment='center')


    # phi
    pylab.text(g.xc[g.nx/2+1], g.yc[g.ny/2]+0.05*g.dy, 
               r"$\phi_{i+1,j}$", 
               fontsize="15",
               horizontalalignment='center', verticalalignment='center')

    pylab.text(g.xc[g.nx/2-1], g.yc[g.ny/2]+0.05*g.dy, 
               r"$\phi_{i-1,j}$", 
               fontsize="15",
               horizontalalignment='center', verticalalignment='center')

    pylab.text(g.xc[g.nx/2], g.yc[g.ny/2+1]+0.05*g.dy, 
               r"$\phi_{i,j+1}$", 
               fontsize="15",
               horizontalalignment='center', verticalalignment='center')

    pylab.text(g.xc[g.nx/2], g.yc[g.ny/2-1]+0.05*g.dy, 
               r"$\phi_{i,j-1}$", 
               fontsize="15",
               horizontalalignment='center', verticalalignment='center')


    # MAC velocities

    # i+1/2,j interface
    pylab.scatter(g.xr[g.nx/2], g.yc[g.ny/2], marker="x", s=50)
    pylab.text(g.xr[g.nx/2]+0.05*g.dx, g.yc[g.ny/2], 
               r"${u}^\mathrm{adv}_{i+1/2,j}$", 
               fontsize="15", rotation=270,
               horizontalalignment='left', verticalalignment='center')

    # i-1/2,j interface
    pylab.scatter(g.xr[g.nx/2-1], g.yc[g.ny/2], marker="x", s=50)
    pylab.text(g.xr[g.nx/2-1]-0.05*g.dx, g.yc[g.ny/2], 
               r"${u}^\mathrm{adv}_{i-1/2,j}$", 
               fontsize="15", rotation=270, 
               horizontalalignment='right', verticalalignment='center')

    # i,j+1/2 interface
    pylab.scatter(g.xc[g.nx/2], g.yr[g.ny/2], marker="x", s=50)
    pylab.text(g.xc[g.nx/2], g.yr[g.ny/2]+0.05*g.dy, 
               r"${v}^\mathrm{adv}_{i,j+1/2}$", 
               fontsize="15", 
               horizontalalignment='center', verticalalignment='bottom')

    # i,j-1/2 interface
    pylab.scatter(g.xc[g.nx/2], g.yr[g.ny/2-1], marker="x", s=50)
    pylab.text(g.xc[g.nx/2], g.yr[g.ny/2-1]-0.05*g.dy, 
               r"${v}^\mathrm{adv}_{i,j-1/2}$", 
               fontsize="15", 
               horizontalalignment='center', verticalalignment='top')



    # axes
    pylab.xlim(g.xl[0]-0.5*g.dx,g.xr[g.nx-1]+0.25*g.dx)
    pylab.ylim(g.yl[0]-0.5*g.dy,g.yr[g.ny-1]+0.25*g.dy)
    pylab.axis("off")

    pylab.subplots_adjust(left=0.02,right=0.98,bottom=0.02,top=0.98)

    ax = pylab.gca()
    ax.set_aspect("equal", "datalim")



    pylab.subplot(122)

    g.draw_grid()

    #------------------------------------------------------------------------
    # label
    pylab.text(g.xc[g.nx/2], g.yc[g.ny/2]+0.1*g.dy, 
               r"$[L\phi]_{i,j}$", fontsize="15",
               horizontalalignment='center', verticalalignment='center')

    pylab.text(g.xc[g.nx/2], g.yc[g.ny/2]-0.1*g.dy, 
               r"$[DU^\mathrm{\,adv}]_{i,j}$", fontsize="15",
               horizontalalignment='center', verticalalignment='center')



    # Gphi velocities

    # i+1/2,j interface
    pylab.scatter(g.xr[g.nx/2], g.yc[g.ny/2], marker="x", s=50)
    pylab.text(g.xr[g.nx/2]+0.05*g.dx, g.yc[g.ny/2], 
               r"${(G\phi)}_{i+1/2,j}$", 
               fontsize="15", rotation=270,
               horizontalalignment='left', verticalalignment='center')

    # i-1/2,j interface
    pylab.scatter(g.xr[g.nx/2-1], g.yc[g.ny/2], marker="x", s=50)
    pylab.text(g.xr[g.nx/2-1]-0.05*g.dx, g.yc[g.ny/2], 
               r"${(G\phi)}_{i-1/2,j}$", 
               fontsize="15", rotation=270, 
               horizontalalignment='right', verticalalignment='center')

    # i,j+1/2 interface
    pylab.scatter(g.xc[g.nx/2], g.yr[g.ny/2], marker="x", s=50)
    pylab.text(g.xc[g.nx/2], g.yr[g.ny/2]+0.05*g.dy, 
               r"${(G\phi)}_{i,j+1/2}$", 
               fontsize="15", 
               horizontalalignment='center', verticalalignment='bottom')

    # i,j-1/2 interface
    pylab.scatter(g.xc[g.nx/2], g.yr[g.ny/2-1], marker="x", s=50)
    pylab.text(g.xc[g.nx/2], g.yr[g.ny/2-1]-0.05*g.dy, 
               r"${(G\phi)}_{i,j-1/2}$", 
               fontsize="15", 
               horizontalalignment='center', verticalalignment='top')



    # axes
    pylab.xlim(g.xl[0]-0.5*g.dx,g.xr[g.nx-1]+0.25*g.dx)
    pylab.ylim(g.yl[0]-0.5*g.dy,g.yr[g.ny-1]+0.25*g.dy)
    pylab.axis("off")

    pylab.subplots_adjust(left=0.02,right=0.98,bottom=0.02,top=0.98)

    ax = pylab.gca()
    ax.set_aspect("equal", "datalim")


    f = pylab.gcf()
    f.set_size_inches(12.0,6.0)


    pylab.savefig("MAC_solve.png")
    pylab.savefig("MAC_solve.eps")
               

if __name__== "__main__":
    simplegrid()
