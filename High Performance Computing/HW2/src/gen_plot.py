"""
Reads in *.asc files produced by simple_wave and generates *.png files
"""
import matplotlib
import matplotlib.pyplot as plt
import numpy as np


def gen_figure_at_iteration(iteration):
    """
    takes iteration number (iteration) as input, loads the appropriately named file and
    generates a plot
    """
    name_prefix = "{0:07d}".format(iteration)
    dname = name_prefix + ".asc"
    pname = name_prefix + ".png"

    data = np.loadtxt(dname)

    fig, ax = plt.subplots()

    xdata = data.T[0]
    udata = data.T[1]
    vdata = data.T[2]

    u_fig = ax.plot(xdata, udata, 'r-x', label=r'$U$')
    v_fig = ax.plot(xdata, vdata, 'b-o', label=r'$V$')

    whole_fig = u_fig + v_fig
    labels = [f.get_label() for f in whole_fig]

    plt.legend(whole_fig, labels, fancybox=True, framealpha=.5, borderpad=.5, shadow=True, loc=0)
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$U$ and $V$')

    fig.savefig(pname)

    plt.close(fig)

    print("Finished processing {0}".format(dname))

num_files = 500

# generate figures from the first num_files .asc files
for i in range(num_files):
    gen_figure_at_iteration(i)
