import ast
import os
from pathlib import Path as P
from pathlib import PurePath as PP
from readDataFile import read

import PIL
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
import pandas as pd

plt.style.use('science')
rc('text.latex', preamble=r'\usepackage{cmbright}')
plt.rcParams['font.family'] = 'sans-serif'


FOLDER = '/Volumes/GoogleDrive/My Drive/Research/protein-dynamics/data/AsLOV DEER'
def plot(folder):
    fs = [ii for ii in P(folder).iterdir() if ii.stem.endswith('DEER') or ii.stem.endswith('DIST') or 'ERROR' in ii.stem]
    expts = [ii for ii in fs if ii.stem.endswith('DEER')]    

    fig, ax = plt.subplots()
    axin = ax.inset_axes([0.4, 0.65, 0.55, 0.3], transform=ax.transAxes)
    lw=1.25
    for exp in expts:
        dist = [ii for ii in fs if ii.stem.startswith(exp.stem.split('DEER')[0]) and ii.stem.endswith('DIST')][0]
        err = [ii for ii in fs if ii.stem.startswith(exp.stem.split('DEER')[0]) and 'ERROR' in ii.stem][0]
        exp_db = pd.read_csv(exp, delimiter='\t', encoding='mac_greek')
        dist_db = pd.read_csv(dist, delimiter='\t')
        err_db = pd.read_csv(err, delimiter='\t')
        if 'light' in exp.stem:
            color = '#00A7CA'
            label = 'Lit'
            style = '--'
        else:
            color = 'black'
            label = 'Dark'
            style = '-'
        n = np.trapz(dist_db['Probability'])
        ax.plot(0.1*dist_db['Distance(A)'], dist_db['Probability']/n, color=color, lw=lw, label=label, ls=style)
        ax.fill_between(0.1*err_db['Distance(A)'], y1=(dist_db['Probability']-2*err_db['Standard  Deviation'])/n, y2=(dist_db['Probability']+2*err_db['Standard  Deviation'])/n, color=color, interpolate=True, alpha=0.5)
        axin.plot(exp_db['time (Βs)'], exp_db['Data'], color=color, lw=lw, label=label, ls=style)
        # axin.plot(exp_db['time (Βs)'], exp_db['Fit'], color='#00A7CA', lw=lw, label='Fit', alpha=0.75)

    ax.set_ylabel('$P(r_{ee})$')
    ax.set_xlabel('Distance $r_{ee}$ (nm)')
    axin.set_xlabel('Time $t$ ($\mu$s)')
    axin.set_ylabel('$F(t)$')
    axin.set_yticklabels([])
    ax.legend(loc=(0.575, 0.125))
    # axin.legend()

    fig.tight_layout()
    fig.savefig(exp.parent.joinpath('comparison_figure.png'), dpi=300)

    # plt.show()

if __name__ == "__main__":
    plot(FOLDER)
