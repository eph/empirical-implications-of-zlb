from contextlib import contextmanager
#import matplotlib

import matplotlib.pyplot as plt
from matplotlib import rc
import seaborn as sns

preamble = r"""
\usepackage[helvet]{sfmath}
\usepackage[scaled]{helvet}
\renewcommand\familydefault{\sfdefault}
\usepackage[T1]{fontenc}
"""

def figure_defaults(n=4):
    rc('text', usetex=True)
    rc('text.latex', preamble=preamble)
    
    rc('lines', linewidth=3)
    sns.set_palette(sns.color_palette("cubehelix", n))
    #sns.set_palette('gray')
    sns.set_style('white', {'axes.grid' : False})
    sns.set_context(rc={'xtick.labelsize': 14, 'ytick.labelsize': 14, 'text.usetex': True, 
                        'font.family':'serif', 'font.sans-serif':[u'Helvetica']})

    return sns.color_palette()

@contextmanager
def saved_figure(fname, **kwds):
    """
    Saves a figure in `fname`.
    """
    fig, ax = plt.subplots(**kwds)
    try:
        yield (fig, ax)
    finally:
        fig.savefig(fname, bbox_inches='tight')
        plt.close(fig)

despine = sns.despine
