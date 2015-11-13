import pylab
import numpy


def drawBox(ll, uu, nx, ny, gridColor="0.5", ng=0):

    # draw the frame
    pylab.plot([ll[0], ll[0], uu[0], uu[0], ll[0]],
               [ll[1], uu[1], uu[1], ll[1], ll[1]], color="k", lw=2)

    # draw the x grid lines
    dx = (uu[0] - ll[0])/nx
    n = 1
    while (n < nx):
        pylab.plot([ll[0]+n*dx, ll[0]+n*dx],
                   [ll[1], uu[1]], color=gridColor, ls=":")

        n += 1

    # draw the y grid lines
    dy = (uu[1] - ll[1])/ny
    n = 1
    while (n < ny):
        pylab.plot([ll[0], uu[0]],
                   [ll[1]+n*dy, ll[1]+n*dy], color=gridColor, ls=":")

        n += 1
    
    # ghostcells?  
    if (ng > 0):
        print "here"
        xmin = ll[0]-ng*dx
        xmax = uu[0]+ng*dx
        ymin = ll[1]-ng*dy
        ymax = uu[1]+ng*dy
        pylab.plot([xmin, xmin, xmax, xmax, xmin], 
                   [ymin, ymax, ymax, ymin, ymin], 
                   ls="--", color="r")


pylab.clf()

drawBox([0., 0.], [1., 1.], 5, 5, ng=1)
drawBox([1., 0.], [2., 1.], 5, 5)
drawBox([2., 0.], [3., 1.], 5, 5)

drawBox([0., 1.], [1., 2.], 5, 5)
drawBox([1., 1.], [2., 2.], 5, 5)
drawBox([2., 1.], [3., 2.], 5, 5)

pylab.xlim(-0.25,3.25)
pylab.ylim(-0.25,2.25)

a = pylab.gca()
a.set_aspect("equal", "datalim")
pylab.axis("off")

f = pylab.gcf()
f.set_size_inches(7.0,4.25)


pylab.savefig("domain.pdf", bbox_inches="tight")


