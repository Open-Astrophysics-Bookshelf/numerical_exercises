import numpy
import pylab

x_smooth = numpy.linspace(0.0, numpy.pi, 500)
f_smooth = numpy.sin(x_smooth)


x = numpy.linspace(0.0, numpy.pi, 10)
f = numpy.sin(x)

i = 5
x_0 = x[i]
f_0 = f[i]

pylab.plot(x_smooth, f_smooth)
pylab.scatter(x, f)

pylab.scatter(x_0, numpy.sin(x_0), color="r", zorder=10)

dx = x[1] - x[0]
x_d = numpy.linspace(x[i]-1.5*dx, x[i]+1.5*dx, 2)

d_exact = numpy.cos(x_0)

d_l = (numpy.sin(x[i]) - numpy.sin(x[i-1]))/dx
d_r = (numpy.sin(x[i+1]) - numpy.sin(x[i]))/dx
d_c = 0.5*(numpy.sin(x[i+1]) - numpy.sin(x[i-1]))/dx


pylab.plot(x_d, d_exact*(x_d - x_0) + f_0, color="0.5", lw=2, ls=":", label="exact")
pylab.plot(x_d, d_l*(x_d - x_0) + f_0, label="left-sided")
pylab.plot(x_d, d_r*(x_d - x_0) + f_0, label="right-sided")
pylab.plot(x_d, d_c*(x_d - x_0) + f_0, label="centered")

ax = pylab.gca()

ax.spines['left'].set_position('zero')
ax.spines['right'].set_color('none')
ax.spines['bottom'].set_position('zero')
ax.spines['top'].set_color('none')
ax.spines['left'].set_smart_bounds(True)
ax.spines['bottom'].set_smart_bounds(True)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')

pylab.legend(frameon=False, loc="best")

pylab.ylim(-0.1, 1.1)

pylab.tight_layout()

pylab.savefig("derivs.eps")



