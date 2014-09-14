#plot formatting scripts and tools
#plot settings
import pylab as PL

plparams = {'backend': 'pdf',
          'axes.labelsize': 14,
          'text.fontsize': 14,
          'title.fontsize': 22,
          'legend.fontsize': 13,
          'xtick.labelsize': 14,
          'ytick.labelsize': 14,
          'text.usetex': False}
PL.rcParams.update(plparams)


def simpleaxis(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()