import math
import numpy
import pylab
import random

# red black gauss seidel pattern

class marker:

    def __init__(self, xc, yc, color):

        # a marker is indicated by its center (xc,yc).
        self.xc = xc
        self.yc = yc
        
        # keep track of the color
        self.color = color

        
def RB():

    # define the number of markers in x and y
    nx = 10
    ny = 10

    # define the length of a marker side
    L = 0.8

    # create a list of marker objects, one at each grid location
    markers = []
    

    color = 0
    while (color <= 1):

        j = 0
        while (j < ny):

            if (color == 0):
                ioff = j % 2
            else:
                ioff = 1 - (j % 2)
                
            i = ioff
            while (i < nx):
                markers.append(marker(i, j, color))
                i += 2

            j += 1

        color += 1


    pylab.clf()

    # the margins are funny -- we pick them to ensure that the
    # plot size is an integer multiple of the number of markers in
    # each dimension
    pylab.subplots_adjust(left=0.1,  right=0.9,
                          bottom=0.1,top=0.9)

    # draw the current state
    n = 0
    while (n < len(markers)):

        if (markers[n].color == 1):
            c = "r"
        else:
            c = "k"

        pylab.fill([markers[n].xc-L/2, markers[n].xc-L/2,
                    markers[n].xc+L/2, markers[n].xc+L/2,
                    markers[n].xc-L/2],
                   [markers[n].yc-L/2, markers[n].yc+L/2,
                    markers[n].yc+L/2, markers[n].yc-L/2,
                    markers[n].yc-L/2], 
                   c)

        n += 1

            
    ax = pylab.axis([-0.5,nx+0.5,-0.5,ny+0.5])
    pylab.axis("off")

    f = pylab.gcf()
    f.set_size_inches(6.0,6.0)

    pylab.tight_layout()

    pylab.savefig("rb.png")
    pylab.savefig("rb.eps", bbox_inches="tight", pad_inches=0)


if __name__== "__main__":
    RB()


    
        
