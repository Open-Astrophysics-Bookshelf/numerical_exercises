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
               horizontalalignment="left", color="r", fontsize=20) 

    pylab.plot([t0,t0], [0, y0], ls=":", color="0.5")
    pylab.text(t0, -0.05, r"$t^n$", 
               verticalalignment="top", horizontalalignment="center",
               fontsize=20)


    # illustrate where t^n+1 is
    pylab.plot([t0+dt,t0+dt], [0, y0], ls=":", color="0.5")
    pylab.text(t0+dt, -0.05, r"$t^{n+1}$", 
               verticalalignment="top", horizontalalignment="center",
               fontsize=20)


start()

pylab.axis("off")

leg = pylab.legend()
ltext = leg.get_texts()
pylab.setp(ltext, fontsize='large')
leg.draw_frame(0)

pylab.xlim(t0-0.05*dt, t0+2.1*dt)

pylab.tight_layout()
pylab.savefig("rk4_initial.png", bbox_inches='tight')


# now draw the solution Euler would get
slope = rhs(y0)
tEuler = numpy.linspace(t0, t0+dt, 2)
yEuler = y0 + slope*(tEuler-t0)

pylab.plot(tEuler, yEuler, label="Euler step")

pylab.scatter([tEuler[1]], [yEuler[1]], color="r")
pylab.text(tEuler[1]+0.015, yEuler[1], r"$y^{n+1}$", 
           horizontalalignment="left", verticalalignment="center",
           color="r", fontsize=20) 


leg = pylab.legend()
ltext = leg.get_texts()
pylab.setp(ltext, fontsize='large')
leg.draw_frame(0)

pylab.xlim(t0-0.25*dt, t0+2.1*dt)

pylab.tight_layout()
pylab.savefig("rk4_Euler.png", bbox_inches='tight')


#----------------------------------------------------------------------------
# show k1 slope 
pylab.clf()
start()

pylab.axis("off")

k1 = rhs(y0)
tSlope = numpy.linspace(t0-0.2*dt, t0+0.2*dt, 2)
ySlope = y0 + k1*(tSlope - t0)

pylab.plot(tSlope, ySlope, label=r"slope", color="g", ls="--", lw=2, zorder=10)
pylab.text(t0-0.1*dt, y0, r"$k_1$", color="g")

leg = pylab.legend()
ltext = leg.get_texts()
pylab.setp(ltext, fontsize='large')
leg.draw_frame(0)

pylab.xlim(t0-0.25*dt, t0+2.1*dt)

pylab.tight_layout()
pylab.savefig("rk4_k1.png", bbox_inches='tight')


#----------------------------------------------------------------------------
# follow k1 to define k2

ytmp = y0 + k1*0.5*dt
k2 = rhs(ytmp)

# draw the k1 half step
pylab.plot([t0, t0+0.5*dt], [y0, ytmp], color="b", label="half-dt k1 step")
pylab.scatter(t0+0.5*dt, ytmp, color="b")

# draw slope there

tSlope = numpy.linspace(t0+0.5*dt-0.2*dt, t0+0.5*dt+0.2*dt, 2)
ySlope = ytmp + k2*(tSlope - (t0 + 0.5*dt))

pylab.plot(tSlope, ySlope, color="g", ls="--", lw=2)
pylab.text(t0+0.5*dt-0.1*dt, ytmp, r"$k_2$", color="g", verticalalignment="top")

leg = pylab.legend()
ltext = leg.get_texts()
pylab.setp(ltext, fontsize='large')
leg.draw_frame(0)

pylab.xlim(t0-0.25*dt, t0+2.1*dt)

pylab.tight_layout()
pylab.savefig("rk4_k2.png", bbox_inches='tight')


#----------------------------------------------------------------------------
# follow k2 to define k3
ytmp = y0 + k2*0.5*dt
k3 = rhs(ytmp)

# draw k2 half step
pylab.plot([t0, t0+0.5*dt], [y0, ytmp], color="c", label="half-dt k2 step")
pylab.scatter(t0+0.5*dt, ytmp, color="c")

# draw slope there
ySlope = ytmp + k3*(tSlope - (t0 + 0.5*dt))
pylab.plot(tSlope, ySlope, color="g", ls="--", lw=2)
pylab.text(t0+0.5*dt+0.05*dt, ytmp, r"$k_3$", color="g", verticalalignment="bottom")

leg = pylab.legend()
ltext = leg.get_texts()
pylab.setp(ltext, fontsize='large')
leg.draw_frame(0)

pylab.xlim(t0-0.25*dt, t0+2.1*dt)

pylab.tight_layout()
pylab.savefig("rk4_k3.png", bbox_inches='tight')


#----------------------------------------------------------------------------
# follow k3 to define k4
ytmp = y0 + k3*dt
k4 = rhs(ytmp)

# draw k3 full step
pylab.plot([t0, t0+dt], [y0, ytmp], color="0.5", label="full-dt k3 step")
pylab.scatter(t0+dt, ytmp, color="0.5")

# draw slope there
tSlope2 = numpy.linspace(t0+dt-0.2*dt, t0+dt+0.2*dt, 2)
ySlope = ytmp + k4*(tSlope2 - (t0 + dt))
pylab.plot(tSlope2, ySlope, color="g", ls="--", lw=2)
pylab.text(t0+dt-0.1*dt, ytmp, r"$k_4$", color="g", verticalalignment="top")

leg = pylab.legend()
ltext = leg.get_texts()
pylab.setp(ltext, fontsize='large')
leg.draw_frame(0)

pylab.xlim(t0-0.25*dt, t0+2.1*dt)

pylab.tight_layout()
pylab.savefig("rk4_k4.png", bbox_inches='tight')



#----------------------------------------------------------------------------
# final RK-4 step
ynew = y0 + (dt/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4)


# draw full RK-4 step
pylab.plot([t0, t0+dt], [y0, ynew], color="r", label="full 4th-order RK step")
pylab.scatter(t0+dt, ynew, color="r")
pylab.text(t0+1.05*dt, ynew+0.015, r"$y^{n+1}$", color = "r", 
           horizontalalignment="left", fontsize=20)


leg = pylab.legend()
ltext = leg.get_texts()
pylab.setp(ltext, fontsize='large')
leg.draw_frame(0)

pylab.xlim(t0-0.25*dt, t0+2.1*dt)

pylab.tight_layout()
pylab.savefig("rk4_final.png", bbox_inches='tight')

