import matplotlib as mpl

def set_graph_specifications():
    rc_dict = {'font.family': 'Univers 57 Condensed',
               'axes.labelsize': 9,
               'axes.titlesize': 9,
               'axes.linewidth': 0.5,
               'xtick.labelsize': 8,
               'xtick.top': True,
               'xtick.bottom': True,
               'xtick.major.size': 7.2,
               'xtick.minor.size': 3.6,
               'xtick.major.width': 0.5,
               'xtick.minor.width': 0.5,
               'xtick.direction': 'in',
               'ytick.labelsize': 8,
               'ytick.left': True,
               'ytick.right': True,
               'ytick.major.size': 7.2,
               'ytick.minor.size': 3.6,
               'ytick.major.width': 0.5,
               'ytick.minor.width': 0.5,
               'ytick.direction': 'in',
               'pdf.fonttype': 42,
               'savefig.dpi': 300,
               'savefig.transparent': True,
               'legend.fontsize': 9,
               'legend.frameon': False,
               'legend.markerscale': 1.
               }
    mpl.rcParams.update(rc_dict)
