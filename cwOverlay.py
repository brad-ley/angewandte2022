import ast
import os
import PIL
from pathlib import Path as P
from pathlib import PurePath as PP

import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
from scipy.integrate import cumtrapz

from readDataFile import read
from makeAbsDisp import make

plt.style.use(['science'])
rc('text.latex', preamble=r'\usepackage{cmbright}')
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 14
plt.rcParams['axes.linewidth'] = 1
plt.rcParams['xtick.major.size'] = 5
plt.rcParams['xtick.major.width'] = 1
plt.rcParams['xtick.minor.size'] = 2
plt.rcParams['xtick.minor.width'] = 1
plt.rcParams['ytick.major.size'] = 5
plt.rcParams['ytick.major.width'] = 1
plt.rcParams['ytick.minor.size'] = 2
plt.rcParams['ytick.minor.width'] = 1


def overlay(targ, low=-1, high=-1):
    """overlay.
    :param targ: directory with 'absorption' files
    :param low: low field for plotting
    :param high: high field for plotting
    """
    if P(targ).is_file():
        targ = str(P(targ).parent)
    fs = [ii for ii in P(targ).iterdir() if ii.name.startswith(
        'absorption') and ii.name.endswith('_exp.txt')]
        # 'dispersion') and ii.name.endswith('_exp.txt')]
    fig, ax = plt.subplots()
    figtrapz, axtrapz = plt.subplots()
    fs.sort()
    fs.insert(2, fs.pop(1))
    for i, f in enumerate(fs):
        print(f)
        h, d = read(str(f))
        d[:, 1] = d[:, 1] - np.average(d[:, 1]) # remove baseline
        if low == -1:
            if 'plotlow' in locals():
                if d[0, 0] < plotlow:
                    plotlow = d[0, 0]
            else:
                plotlow = d[0, 0]
        else:
            plotlow = low
        if high == -1:
            if 'plothigh' in locals():
                if d[-1, 0] > plothigh:
                    plothigh = d[-1, 0]
            else:
                plothigh = d[-1, 0]
        else:
            plothigh = high

        ### for generic plotting ###
        # legend = f"{f.name.split('_')[1]}"
        ############################
        ### for manuscript plotting ###
        if "DL" in f.name.split('_')[1]:
            legend = 'DL'
            color = 'blue'
            style = ':'
        elif "537" in f.name.split('_')[1]:
            legend = 'SL E537C'
            color ='black'
            style = '-'
        elif "406" in f.name.split('_')[1]:
            legend = 'SL T406C'
            color ='red'
            style = '--'
        # if False:
        #     pass
        else:
            legend = f"{f.name.split('_')[1]}"
            # color='black'
        ############################
       
        lw=2
        try:
            scale = np.trapz(cumtrapz(d[np.where(plotlow < d[:, 0])[
                             0][0]:np.where(plothigh < d[:, 0])[0][0], 1]))
            ax.plot(d[:, 0], d[:, 1] / scale, label=legend, lw=lw, color=color, linestyle=style)
            axtrapz.plot(d[1:, 0], cumtrapz(d[:, 1]) / scale, label=legend)
        except IndexError:
            ax.plot(d[:, 0], d[:, 1] / np.max(d[:, 1]), label=legend, lw=lw, color=color, linestyle=style)
            axtrapz.plot(d[1:, 0], cumtrapz(d[:, 1]) / np.max(cumtrapz(d[: 1])), label=legend)
        except (UnboundLocalError, NameError):
            ax.plot(d[:, 0], d[:, 1] / np.max(d[:, 1]), label=legend, lw=lw)
            axtrapz.plot(d[1:, 0], cumtrapz(d[:, 1]) / np.max(cumtrapz(d[:, 1])), label=legend)
    

    ax.text(0.05, 0.11, '$T=87$ K',
            horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)
    try:
        for a in [ax, axtrapz]:
            a.set_xlim(right=8.6325)
            # a.set_yticklabels([])
            a.set_yticks([-1,0,1])
            a.set_ylabel('cwEPR signal (arb. u)')
            a.set_xlabel('Field (T)')
            # a.spines['right'].set_visible(False)
            # a.spines['top'].set_visible(False)
            a.legend(loc='upper right',markerfirst=False,handlelength=1,handletextpad=0.4,labelspacing=0.2)

        plt.tight_layout()
        fig.savefig(P(targ).joinpath('figure_comp.png'), dpi=300)
        fig.savefig(P(targ).joinpath('figure_comp.tif'), dpi=300)

    except UnboundLocalError:
        print("Need to create absorption files with makeAbsDisp.py first!")


def main():
    # make(targ, keyw='sweep', field=8.6165, numerical_keyw=False)
    overlay(targ)


if __name__ == "__main__":
    targ = '/Volumes/GoogleDrive/My Drive/Research/Data/2022/4/comparing SLDL'
    main()
    # plt.show()
