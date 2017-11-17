import numpy
from matplotlib import pyplot

ng = 2

def stencil_5_pts(q):
    """
    Do reconstruction to q_{i+1/2} using five point stencil.
    
    Parameters
    ----------
    
    q : numpy array
        Scalar data to reconstruct
        
    Returns
    -------
    
    qL : numpy array
        Reconstructed data - boundary points are zero
    """
    qL = numpy.zeros_like(q)
    np = len(q) - 2 * ng
    for i in range(ng, np+ng):
        qL[i] = ( 2 * q[i-2] - 13 * q[i-1] + 47 * q[i] + 
                  27 * q[i+1] - 3 * q[i+2] ) / 60
    return qL

def stencils_3_pts(q):
    """
    Do reconstruction to q_{i+1/2} using three point stencils.
    
    Parameters
    ----------
    
    q : numpy array
        Scalar data to reconstruct
        
    Returns
    -------
    
    qL : numpy array
        Reconstructed data - boundary points are zero
    """
    qL = numpy.zeros((3, len(q)))
    np = len(q) - 2 * ng
    for i in range(ng, np+ng):
        qL[2, i] = ( 2 * q[i-2] - 7 * q[i-1] + 11 * q[i  ]) / 6
        qL[1, i] = (-1 * q[i-1] + 5 * q[i  ] +  2 * q[i+1]) / 6
        qL[0, i] = ( 2 * q[i  ] + 5 * q[i+1] -  1 * q[i+2]) / 6
    return qL

def weno3(q):
    """
    Do WENO reconstruction, r=3
    
    Parameters
    ----------
    
    q : numpy array
        Scalar data to reconstruct
        
    Returns
    -------
    
    qL : numpy array
        Reconstructed data - boundary points are zero
    """
    order = 3
    C = numpy.array([1, 6, 3]) / 10
    a = numpy.array([[11, -7, 2], [2, 5, -1], [-1, 5, 2]]) / 6
    sigma = numpy.array([[[40, 0, 0],
                            [-124, 100, 0],
                            [44, -76, 16] ],
                           [[16, 0, 0],
                            [-52, 52, 0],
                            [20, -52, 16] ],
                           [[16, 0, 0],
                            [-76, 100, 0],
                            [44, -124, 40] ] ]) / 12

    qL = numpy.zeros_like(q)
    beta = numpy.zeros((3, len(q)))
    w = numpy.zeros_like(beta)
    np = len(q) - 2 * ng
    epsilon = 1e-16
    for i in range(ng, np+ng):
        q_stencils = numpy.zeros(order)
        alpha = numpy.zeros(order)
        for k in range(order):
            for l in range(order):
                for m in range(l+1):
                    beta[k, i] += sigma[k, l, m] * q[i+k-l] * q[i+k-m]
            alpha[k] = C[k] / (epsilon + beta[k, i]**2)
            for l in range(order):
                q_stencils[k] += a[k, l] * q[i+k-l]
        w[:, i] = alpha / numpy.sum(alpha)
        qL[i] = numpy.dot(w[:, i], q_stencils)
    
    return qL, beta, w

# First, use a smooth sine function to check convergence of the
# various reconstructions

f_smooth = lambda x: numpy.sin(2 * numpy.pi * x)
f_antiderivative = lambda x: -numpy.cos(2 * numpy.pi * x) / (2 * numpy.pi)

x_plot = numpy.linspace(0, 1, 200)
x_recon, dx = numpy.linspace(0, 1, 20, retstep=True)
x_recon_interior = x_recon[ng:-ng] + dx / 2
s_plot = f_smooth(x_plot)
s_recon = f_smooth(x_recon)
recon_exact = f_smooth(x_recon_interior)
s_5 = stencil_5_pts(s_recon)
s_3 = stencils_3_pts(s_recon)

# The fifth order case
pyplot.figure()
pyplot.plot(x_plot, s_plot, 'k--', lw=2)
pyplot.plot(x_recon, s_recon, 'bo', ms=5)
pyplot.plot(x_recon_interior, s_5[ng:-ng], 'r^', ms=5)
pyplot.xlabel(r"$x$")
pyplot.xlim(0, 1)

# The third order cases
fig, axes = pyplot.subplots(1, 3)
for r in range(3):
    axes[r].plot(x_plot, s_plot, 'k--', lw=2)
    axes[r].plot(x_recon, s_recon, 'bo', ms=5)
    axes[r].plot(x_recon_interior, s_3[r, ng:-ng], 'r^', ms=5)
    axes[r].set_xlabel(r"$x$")
    axes[r].set_xlim(0, 1)
fig.tight_layout()

# Convergence check, using the lazy fix to get error norms.
points = 20 * 2**numpy.arange(1, 7)
error_1norm_5 = numpy.zeros((len(points),))
error_1norm_3 = numpy.zeros((3, len(points)))
error_1norm_weno = numpy.zeros((len(points)))
error_2norm_5 = numpy.zeros((len(points),))
error_2norm_3 = numpy.zeros((3, len(points)))
error_2norm_weno = numpy.zeros((len(points),))
for i, pts in enumerate(points):
    x_recon, dx = numpy.linspace(0, 1, pts, retstep=True)
    x_recon_interior = x_recon[ng:-ng] + dx / 2
    s_recon = (f_antiderivative(x_recon + dx/2) -
               f_antiderivative(x_recon - dx/2)) / dx
    recon_exact = f_smooth(x_recon_interior)
    s_5 = stencil_5_pts(s_recon)
    s_3 = stencils_3_pts(s_recon)
    s_weno, beta, w = weno3(s_recon)
    error_1norm_5[i] = numpy.linalg.norm(s_5[ng:-ng] - recon_exact, 1) / pts
    error_2norm_5[i] = numpy.linalg.norm(s_5[ng:-ng] - recon_exact, 2) / numpy.sqrt(pts)
    error_1norm_weno[i] = numpy.linalg.norm(s_weno[ng:-ng] - recon_exact, 1) / pts
    error_2norm_weno[i] = numpy.linalg.norm(s_weno[ng:-ng] - recon_exact, 2) / numpy.sqrt(pts)
    for r in range(3):
        error_1norm_3[r, i] = numpy.linalg.norm(s_3[r, ng:-ng] - recon_exact, 1) / pts
        error_2norm_3[r, i] = numpy.linalg.norm(s_3[r, ng:-ng] - recon_exact, 2) / numpy.sqrt(pts)


# Plot everything
x_plot = numpy.linspace(0, 1, 200)
x_recon, dx = numpy.linspace(0, 1, 20, retstep=True)
x_recon_interior = x_recon[ng:-ng] + dx / 2
s_plot = f_smooth(x_plot)
s_recon = f_smooth(x_recon)
recon_exact = f_smooth(x_recon_interior)
s_5 = stencil_5_pts(s_recon)
s_3 = stencils_3_pts(s_recon)
s_weno, beta, w = weno3(s_recon)
fig = pyplot.figure(figsize=(10,8))
ax1 = pyplot.subplot2grid((2,2),(0,0),colspan=2)
ax2 = pyplot.subplot2grid((2,2),(1,0))
ax3 = pyplot.subplot2grid((2,2),(1,1))
ax1.plot(x_plot, s_plot, 'k--', lw=2)
ax1.plot(x_recon, s_recon, 'bo', ms=7)
ax1.plot(x_recon_interior, s_5[ng:-ng], 'r^', ms=7)
ax1.set_xlabel(r"$x$")
ax1.set_xlim(0, 1)
ax1.set_ylabel("Smooth function")
ax2.loglog(1/points, error_1norm_5, 'k^', ms=5, label='stencil 5')
ax2.loglog(1/points, error_1norm_weno, 'r+', ms=5, label='WENO')
for r in range(3):
    ax2.loglog(1/points, error_1norm_3[r,:], 'bo', ms=3, 
                  label="stencil 3 ({})".format(r))
ax2.loglog(1/points, error_1norm_3[0,0]*(points[0]/points)**3, 
              label=r"$\propto \Delta x^3$")
ax2.loglog(1/points, error_1norm_5[0]*(points[0]/points)**5, 
              label=r"$\propto \Delta x^5$")
ax2.legend()
ax2.set_xlabel(r"$\Delta x$")
ax2.set_ylabel(r"$|$Error$|_1$")

ax3.loglog(1/points, error_2norm_5, 'k^', ms=5, label='stencil 5')
ax3.loglog(1/points, error_2norm_weno, 'r+', ms=5, label='WENO')
for r in range(3):
    ax3.loglog(1/points, error_2norm_3[r,:], 'bo', ms=3, 
                  label="stencil 3 ({})".format(r))
ax3.loglog(1/points, error_2norm_3[0,0]*(points[0]/points)**3, 
              label=r"$\propto \Delta x^3$")
ax3.loglog(1/points, error_2norm_5[0]*(points[0]/points)**5, 
              label=r"$\propto \Delta x^5$")
ax3.legend()
ax3.set_xlabel(r"$\Delta x$")
ax3.set_ylabel(r"$|$Error$|_2$")
fig.tight_layout()
pyplot.savefig("weno-convergence.pdf")


# Just the convergence plots
fig, axes = pyplot.subplots(1, 2)
axes[0].loglog(1/points, error_1norm_5, 'k^', ms=5, label='stencil 5')
axes[0].loglog(1/points, error_1norm_weno, 'r+', ms=5, label='WENO')
for r in range(3):
    axes[0].loglog(1/points, error_1norm_3[r,:], 'bo', ms=3, 
                  label="stencil 3 ({})".format(r))
axes[0].loglog(1/points, error_1norm_3[0,0]*(points[0]/points)**3, 
              label=r"$\propto \Delta x^3$")
axes[0].loglog(1/points, error_1norm_5[0]*(points[0]/points)**5, 
              label=r"$\propto \Delta x^5$")
axes[0].legend()
axes[0].set_xlabel(r"$\Delta x$")
axes[0].set_ylabel(r"$|$Error$|_1$")

axes[1].loglog(1/points, error_2norm_5, 'k^', ms=5, label='stencil 5')
axes[1].loglog(1/points, error_2norm_weno, 'r+', ms=5, label='WENO')
for r in range(3):
    axes[1].loglog(1/points, error_2norm_3[r,:], 'bo', ms=3, 
                  label="stencil 3 ({})".format(r))
axes[1].loglog(1/points, error_2norm_3[0,0]*(points[0]/points)**3, 
              label=r"$\propto \Delta x^3$")
axes[1].loglog(1/points, error_2norm_5[0]*(points[0]/points)**5, 
              label=r"$\propto \Delta x^5$")
axes[1].legend()
axes[1].set_xlabel(r"$\Delta x$")
axes[1].set_ylabel(r"$|$Error$|_2$")
fig.tight_layout()    

# Now we work with a non-smooth function

f_nonsmooth = lambda x: numpy.where(x<0.5, numpy.sin(2 * numpy.pi * x),
                                           numpy.cos(2 * numpy.pi * x))

x_plot = numpy.linspace(0, 1, 200)
x_recon, dx = numpy.linspace(0, 1, 20, retstep=True)
x_recon_interior = x_recon[ng:-ng] + dx / 2
s_plot = f_nonsmooth(x_plot)
s_recon = f_nonsmooth(x_recon)
recon_exact = f_nonsmooth(x_recon_interior)
s_5 = stencil_5_pts(s_recon)
s_3 = stencils_3_pts(s_recon)
s_weno, beta, w = weno3(s_recon)

# 5 point stencil
pyplot.figure()
pyplot.plot(x_plot, s_plot, 'k--', lw=2)
pyplot.plot(x_recon, s_recon, 'bo', ms=5)
pyplot.plot(x_recon_interior, s_5[ng:-ng], 'r^', ms=5)
pyplot.xlabel(r"$x$")
pyplot.xlim(0, 1)

# 3 point stencils
fig, axes = pyplot.subplots(1, 3)
for r in range(3):
    axes[r].plot(x_plot, s_plot, 'k--', lw=2)
    axes[r].plot(x_recon, s_recon, 'bo', ms=5)
    axes[r].plot(x_recon_interior, s_3[r, ng:-ng], 'r^', ms=5)
    axes[r].set_xlabel(r"$x$")
    axes[r].set_xlim(0, 1)
fig.tight_layout()

# WENO reconstruction
pyplot.figure()
pyplot.plot(x_plot, s_plot, 'k--', lw=2)
pyplot.plot(x_recon, s_recon, 'bo', ms=5)
pyplot.plot(x_recon_interior, s_weno[ng:-ng], 'r^', ms=5)
pyplot.xlabel(r"$x$")
pyplot.xlim(0, 1)

# The smoothness indicators and the weights
colors=['r','b','g']
markers=['o','^','s']
optimal_C = numpy.array([1, 6, 3]) / 10
fig, axes = pyplot.subplots(2, 1, sharex=True, figsize=(6,10))
fig.subplots_adjust(hspace=0)
for r, (color, marker) in enumerate(zip(colors, markers)):
    axes[0].scatter(x_recon_interior, beta[r, ng:-ng], 
        marker=marker, color=color,
        label=r"Jiang-Shu $\beta_{}$".format(r))
    axes[1].scatter(x_recon_interior, w[r, ng:-ng], 
        marker=marker, color=color,
        label=r"WENO weight $\omega_{}$".format(r))
    axes[1].plot(x_recon, optimal_C[r] * numpy.ones_like(x_recon),
        linestyle='--', color=color, lw=2,
        label=r"Optimal weight $C_{}$".format(r))
axes[1].set_xlabel(r"$x$")
axes[0].set_ylabel(r"$\beta$")
axes[0].set_xlim(0,1)
axes[1].set_xlim(0,1)
axes[0].legend()
axes[1].legend()
fig.tight_layout()


# Plot everything
fig, axes = pyplot.subplots(3, 2, sharex = True, figsize=(12,14))
fig.subplots_adjust(hspace=0)
# 5 point stencil
axes[0,0].plot(x_plot, s_plot, 'k--', lw=2)
axes[0,0].plot(x_recon, s_recon, 'kx', ms=7)
axes[0,0].plot(x_recon_interior, s_5[ng:-ng], 'r^', ms=7)
axes[0,0].set_xlim(0, 1)
axes[0,0].set_ylabel('5 point stencil')

# 3 point stencils
colors=['r','b','g']
markers=['o','^','s']
axes[1,0].plot(x_plot, s_plot, 'k--', lw=2)
axes[1,0].plot(x_recon, s_recon, 'kx', ms=7)
for r, (color, marker) in enumerate(zip(colors, markers)):
    axes[1,0].plot(x_recon_interior, s_3[r, ng:-ng], 
                   '{}{}'.format(color, marker), ms=7)
axes[1,0].set_xlim(0, 1)
axes[1,0].set_ylabel('3 point stencils')

# WENO reconstruction
axes[2,0].plot(x_plot, s_plot, 'k--', lw=2)
axes[2,0].plot(x_recon, s_recon, 'kx', ms=7)
axes[2,0].plot(x_recon_interior, s_weno[ng:-ng], 'r^', ms=7)
axes[2,0].set_xlabel(r"$x$")
axes[2,0].set_xlim(0, 1)
axes[2,0].set_ylabel('WENO stencils')

# Shift the right axes
for r in range(3):
    axes[r,1].yaxis.tick_right()
    axes[r,1].yaxis.set_label_position("right")

# The function again
axes[0,1].plot(x_plot, s_plot, 'k-', lw=2)
axes[0,1].set_ylabel('Non-smooth function')

# The smoothness indicators and the weights
colors=['r','b','g']
markers=['o','^','s']
optimal_C = numpy.array([1, 6, 3]) / 10
for r, (color, marker) in enumerate(zip(colors, markers)):
    axes[1,1].scatter(x_recon_interior, beta[r, ng:-ng], 
        marker=marker, color=color,
        label=r"Jiang-Shu $\beta_{}$".format(r))
    axes[2,1].scatter(x_recon_interior, w[r, ng:-ng], 
        marker=marker, color=color,
        label=r"WENO weight $\omega_{}$".format(r))
    axes[2,1].plot(x_recon, optimal_C[r] * numpy.ones_like(x_recon),
        linestyle='--', color=color, lw=2,
        label=r"Optimal weight $C_{}$".format(r))
axes[2,1].set_xlabel(r"$x$")
axes[1,1].set_ylabel(r"$\beta$")
axes[2,1].set_ylabel(r"$C$, $\omega$")
axes[1,1].set_xlim(0,1)
axes[2,1].set_xlim(0,1)
axes[1,1].legend()
axes[2,1].legend()
fig.tight_layout()
pyplot.savefig("weno-weights.pdf")