# make plots using the root finding methods in roots.py
#
# M. Zingale (2013-02-14)

import math
import numpy
import matplotlib.pyplot as plt
from roots import *
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset


r = root(f, 1.e-5, fprime=fprime)

xmin = 0.0                                                                                                                      
xmax = 5.0 

xfine = numpy.linspace(xmin, xmax, 200)                                                                                        

#---------------------------------------------------------------------------
# plot Newton
xmin = 0.0
xmax = 5.0

root, xeval = r.newton(xmax)

for n, x in enumerate(xeval):
    plt.clf()
    plt.plot(xfine, r.f(xfine))

    plt.scatter(numpy.array(xeval[0:n+1]), 
                  r.f(numpy.array(xeval[0:n+1])), 
                  marker="x", s=25, color="r")
 
    plt.plot(xfine, r.fprime(x)*(xfine-x) + r.f(x), color="0.5")

    xintercept = -r.f(x)/r.fprime(x) + x
    print "xintercept = ", xintercept
    plt.plot([xintercept, xintercept], [0, r.f(xintercept)], color="0.5", ls=":")

    if n%2 == 0:
        plt.text(x, r.f(x)-0.3, `n`, color="r", fontsize="10",
                   verticalalignment="top", horizontalalignment="center") 
    else:
        plt.text(x, r.f(x)+0.3, `n`, color="r", fontsize="10",
                   verticalalignment="bottom", horizontalalignment="center") 

    F = plt.gcf()
    plt.text(0.5, 0.02, "root approx = %s" % (`x`),
               transform = F.transFigure, color="k")

    # axes through origin
    ax = plt.gca()
    ax.spines['left'].set_position('zero')
    ax.spines['right'].set_color('none')
    ax.spines['bottom'].set_position('zero')
    ax.spines['top'].set_color('none')
    ax.spines['left'].set_smart_bounds(True)
    ax.spines['bottom'].set_smart_bounds(True)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    plt.ylim(-10.0,15)

    # inset near the root (x=3)
    axins = zoomed_inset_axes(ax, 6, loc=9) # zoom = 6
    x1, x2, y1, y2 = 2.8, 3.2, -0.5, 0.5
    axins.set_xlim(x1, x2)
    axins.set_ylim(y1, y2)

    plt.xticks(visible=False)
    plt.yticks(visible=False)

    axins.plot(xfine, r.f(xfine))

    axins.scatter(numpy.array(xeval[0:n+1]), 
                  r.f(numpy.array(xeval[0:n+1])), 
                  marker="x", s=25, color="r")

    axins.plot(xfine, r.fprime(x)*(xfine-x) + r.f(x), color="0.5")

    xintercept = -r.f(x)/r.fprime(x) + x
    axins.plot([xintercept, xintercept], [0, r.f(xintercept)], color="0.5", ls=":")

    if n%2 == 0:
        axins.text(x, r.f(x)-0.05, `n`, color="r", fontsize="10",
                   verticalalignment="top", horizontalalignment="center",
                   clip_on=True) 
    else:
        axins.text(x, r.f(x)+0.05, `n`, color="r", fontsize="10",
                   verticalalignment="bottom", horizontalalignment="center",
                   clip_on=True) 

    axins.plot(xfine, xfine*0, ls=":", color="0.5")

    mark_inset(ax, axins, loc1=3, loc2=4, fc="none", ec="0.5")

    
    f = plt.gcf()
    f.set_size_inches(7.20,7.20)


    plt.tight_layout()

    plt.savefig("newton_%2.2d.eps" % (n))



