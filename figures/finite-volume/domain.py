import matplotlib.pyplot as plt

def draw_box(ll, uu, nx, ny, gridColor="0.5", ng=0):

    # draw the frame
    plt.plot([ll[0], ll[0], uu[0], uu[0], ll[0]],
             [ll[1], uu[1], uu[1], ll[1], ll[1]], color="k", lw=2)

    # draw the x grid lines
    dx = (uu[0] - ll[0])/nx
    for n in range(1, nx):
        plt.plot([ll[0]+n*dx, ll[0]+n*dx],
                   [ll[1], uu[1]], color=gridColor, ls=":", lw=1)

    # draw the y grid lines
    dy = (uu[1] - ll[1])/ny
    for n in range(1, ny):
        plt.plot([ll[0], uu[0]],
                   [ll[1]+n*dy, ll[1]+n*dy], color=gridColor, ls=":", lw=1)
    
    # ghostcells?  
    if ng > 0:
        xmin = ll[0]-ng*dx
        xmax = uu[0]+ng*dx
        ymin = ll[1]-ng*dy
        ymax = uu[1]+ng*dy
        plt.plot([xmin, xmin, xmax, xmax, xmin], 
                 [ymin, ymax, ymax, ymin, ymin], 
                 ls="--", color="r")


plt.clf()

draw_box([0., 0.], [1., 1.], 5, 5, ng=1)
draw_box([1., 0.], [2., 1.], 5, 5)
draw_box([2., 0.], [3., 1.], 5, 5)

draw_box([0., 1.], [1., 2.], 5, 5)
draw_box([1., 1.], [2., 2.], 5, 5)
draw_box([2., 1.], [3., 2.], 5, 5)

plt.xlim(-0.25,3.25)
plt.ylim(-0.25,2.25)

a = plt.gca()
a.set_aspect("equal", "datalim")
plt.axis("off")

f = plt.gcf()
f.set_size_inches(7.0,4.25)


plt.savefig("domain.pdf", bbox_inches="tight")


