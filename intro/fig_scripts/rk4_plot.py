# make a plot of what the 4th-order Runge-Kutta is doing using the
# example dy/dt = -y (this comes from Garcia).

import numpy
import pylab

def rhs(y):
    """ return dy/dt """
    return -y

def exact(t):
    """ analytic solution """
    return numpy.exp(-t)


# frame 1 -- show the initial condition (y^n)

y0 = 1.0
t0 = 0.0

dt = 0.75

tt = numpy.linspace(t0, t0+2.0*dt, 100)


def start():
    """ default starting point """

    pylab.plot(tt, exact(tt), label="analytic solution", color="k")

    # draw the timeline
    pylab.plot(tt, 0*tt, color="k")

    # label the current point
    pylab.scatter([t0], [y0], color="r")
    pylab.text(t0, y0+0.03, r"$y^n$", 
               horizontalalignment="left", color="r", fontsize=18) 

    pylab.plot([t0,t0], [0, y0], ls=":", color="0.5")
    pylab.text(t0, -0.05, r"$t^n$", 
               verticalalignment="top", horizontalalignment="center",
               fontsize=18)


    # illustrate where t^n+1 is
    pylab.plot([t0+dt,t0+dt], [0, y0], ls=":", color="0.5")
    pylab.text(t0+dt, -0.05, r"$t^{n+1}$", 
               verticalalignment="top", horizontalalignment="center",
               fontsize=18)


start()

pylab.axis("off")

leg = pylab.legend()
ltext = leg.get_texts()
pylab.setp(ltext, fontsize='medium')
leg.draw_frame(0)

pylab.xlim(t0-0.05*dt, t0+2.1*dt)

pylab.tight_layout()
pylab.savefig("rk4_initial.eps", bbox_inches='tight')


# now draw the solution Euler would get
slope = rhs(y0)
tEuler = numpy.linspace(t0, t0+dt, 2)
yEuler = y0 + slope*(tEuler-t0)

pylab.plot(tEuler, yEuler, label="Euler step")

pylab.scatter([tEuler[1]], [yEuler[1]], color="r")
pylab.text(tEuler[1]+0.015, yEuler[1], r"$y^{n+1}$", 
           horizontalalignment="left", verticalalignment="center",
           color="r", fontsize=18) 


leg = pylab.legend()
ltext = leg.get_texts()
pylab.setp(ltext, fontsize='medium')
leg.draw_frame(0)

pylab.xlim(t0-0.05*dt, t0+2.1*dt)

pylab.tight_layout()
pylab.savefig("rk4_Euler.eps")



# k1 slope uses y^n position
pylab.clf()
start()

pylab.axis("off")

k1 = rhs(y0)
thalf = numpy.linspace(t0, t0+0.5*dt, 2)
yEuler = y0 + k1*(thalf-t0)

pylab.plot(thalf, yEuler, label="half-dt k1 step", color="b")

pylab.scatter([thalf[1]], [yEuler[1]], color="b")

# illustrate the slope at that position
k2 = rhs(yEuler[1])
tSlope = numpy.linspace(thalf[1]-0.2*dt, thalf[1]+0.2*dt, 2)
ySlope = yEuler[1] + k2*(tSlope - thalf[1])

pylab.plot(tSlope, ySlope, label=r"slope", color="g", ls="--", lw=2)

leg = pylab.legend()
ltext = leg.get_texts()
pylab.setp(ltext, fontsize='medium')
leg.draw_frame(0)

pylab.xlim(t0-0.05*dt, t0+2.1*dt)

pylab.tight_layout()
pylab.savefig("rk4_k1.eps")


# k2 slope goes from y^n to 1/2 time to get k3 slope
ytmp = y0 + k2*0.5*dt
k3 = rhs(ytmp)

# draw k2 half step
pylab.plot([t0, t0+0.5*dt], [y0, ytmp], color="c", label="half-dt k2 step")
pylab.scatter(t0+0.5*dt, ytmp, color="c")

# draw slope there
ySlope = ytmp + k3*(tSlope - (t0 + 0.5*dt))
pylab.plot(tSlope, ySlope, color="g", ls="--", lw=2)


leg = pylab.legend()
ltext = leg.get_texts()
pylab.setp(ltext, fontsize='medium')
leg.draw_frame(0)

pylab.xlim(t0-0.05*dt, t0+2.1*dt)

pylab.tight_layout()
pylab.savefig("rk4_k2.eps")


# k3 slope goes from y^n to full time to get k4 slope
ytmp = y0 + k3*dt
k4 = rhs(ytmp)

# draw k3 full step
pylab.plot([t0, t0+dt], [y0, ytmp], color="0.5", label="full-dt k3 step")
pylab.scatter(t0+dt, ytmp, color="0.5")

# draw slope there
tSlope2 = numpy.linspace(t0+dt-0.2*dt, t0+dt+0.2*dt, 2)
ySlope = ytmp + k4*(tSlope2 - (t0 + dt))
pylab.plot(tSlope2, ySlope, color="g", ls="--", lw=2)


leg = pylab.legend()
ltext = leg.get_texts()
pylab.setp(ltext, fontsize='medium')
leg.draw_frame(0)

pylab.xlim(t0-0.05*dt, t0+2.1*dt)

pylab.tight_layout()
pylab.savefig("rk4_k3.eps")



# final RK-4 step
ynew = y0 + (dt/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4)


# draw full RK-4 step
pylab.plot([t0, t0+dt], [y0, ynew], color="r", label="full-dt R-K 4 step")
pylab.scatter(t0+dt, ynew, color="r")
pylab.text(t0+1.05*dt, ynew+0.015, r"$y^{n+1}$", color = "r", 
           horizontalalignment="left", fontsize=18)


leg = pylab.legend()
ltext = leg.get_texts()
pylab.setp(ltext, fontsize='medium')
leg.draw_frame(0)

pylab.xlim(t0-0.05*dt, t0+2.1*dt)

pylab.tight_layout()
pylab.savefig("rk4_final.eps")

