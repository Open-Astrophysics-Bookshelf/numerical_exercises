import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['mathtext.rm'] = 'serif'


def fun(x):
    state = np.zeros_like(x)
    state[x < 0.5] = 1.0
    return state

def funs(x):
    state = 0.5 + 0.25*np.sin(2.0*np.pi*x)
    return state

def funr(x):
    state = np.ones_like(x)
    state[x < 0.5] = 0.2
    return state

def make_plot(icfun=None, label=""):

    npts_plot = 1000
    xplot = np.linspace(0.0, 1.0, npts_plot)

    nchar = 20
    xchar = np.linspace(0.0, 1.0, nchar)

    state = icfun(xplot)
    print(type(state))

    plt.clf()

    plt.subplot(211)

    plt.plot(xplot, state, lw=2)


    ax = plt.gca()
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    plt.ylim(-0.1, 1.1)

    plt.xlabel(r"$x$", fontsize="large")
    plt.ylabel(r"$u$", fontsize="large")



    plt.subplot(212)

    uchar = icfun(xchar)
    t = np.linspace(0.0, 1.0, 100)
    for n in range(nchar):
        xc = xchar[n] + uchar[n]*t
        plt.plot(xc, t, color="0.5")

    ax = plt.gca()
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    plt.xlabel(r"$x$", fontsize="large")
    plt.ylabel(r"$t$", fontsize="large")

    plt.xlim(0.0, 1.0)
    plt.ylim(0.0, 1.1)

    f = plt.gcf()
    f.set_size_inches(6.0, 8.0)

    plt.tight_layout()

    plt.savefig(f"burgers-characteristics-{label}.png")

if __name__ == "__main__":

    make_plot(funr, "rare")
    make_plot(funs, "shock")
