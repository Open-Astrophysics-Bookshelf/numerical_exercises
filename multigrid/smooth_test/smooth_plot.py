import math
import numpy
import pylab

def plot_convergence():

    data1    = numpy.loadtxt("relax.mixed.00001")
    data10   = numpy.loadtxt("relax.mixed.00010")
    data100  = numpy.loadtxt("relax.mixed.00100")
    data1000 = numpy.loadtxt("relax.mixed.01000")

    nx1     = data1[:,0]
    phi1    = data1[:,1]

    nx10    = data10[:,0]
    phi10   = data10[:,1]

    nx100   = data100[:,0]
    phi100  = data100[:,1]

    nx1000  = data1000[:,0]
    phi1000 = data1000[:,1]


    ax = pylab.subplot(111)

    #pylab.scatter(nx1, phi1, marker="x", color="r")
    pylab.plot(nx1,    phi1,    color="r", label="N = 1")
    pylab.plot(nx10,   phi10,   color="g", label="N = 10")
    pylab.plot(nx100,  phi100,  color="b", label="N = 100")
    pylab.plot(nx1000, phi1000, color="c", label="N = 1000")

    pylab.xlabel("x")
    pylab.ylabel("solution error")

    pylab.xlim(0,1)

    leg = pylab.legend()
    leg.draw_frame(0)

    ltext = leg.get_texts()
    pylab.setp(ltext, fontsize='small')

    ax.spines['left'].set_position('zero')
    ax.spines['right'].set_color('none')
    ax.spines['bottom'].set_position('zero')
    ax.spines['top'].set_color('none')
    ax.spines['left'].set_smart_bounds(True)
    ax.spines['bottom'].set_smart_bounds(True)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    f = pylab.gcf()
    f.set_size_inches(5.0,5.0)

    pylab.savefig("smooth_error.png", bbox_inches="tight")
    pylab.savefig("smooth_error.eps", bbox_inches="tight")

    

if __name__== "__main__":
    plot_convergence()

