from pathlib import Path as P
from pathlib import PurePath as PP

import PIL
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np

plt.style.use('science')
rc('text.latex', preamble=r'\usepackage{cmbright}')
plt.rcParams['font.family'] = 'sans-serif'

"""
filename need to be full path to folder:
-- on Mac you can get this by right-clicking and then
holding option and selecting 'Copy as Pathname'
-- on Windows you can get this by holding shift before rightclicking and then
selecting 'Copy as path'
this code will assume your data is a .txt, .dat, or .csv file with the x-axis in the first column and all y-axis data in columns 2, 3, etc.
"""
"""
CHANGE STUFF BELOW
"""
FOLDER = '/Volumes/GoogleDrive/My Drive/Research/Data/2022/1/compare light on off' 
MUTANT = ''
TEMP = 287
FIELDS = [8.623, 8.625, 8.627, 8.628] # leave this list empty if you don't want field markers
legend_names = ['line 1'] # add as many legend names as you want to be plotted on same x axis
colors = ['red'] # as many color names as legend names
# for light on I use color '#00A7CA' and style '--'
styles = [':'] # same number here
x_Label = 'sample x label' # x axis label
y_Label = 'sample y label' # y axis label
x_tickLabels = True # set to false to remove x axis tick labels
x_ticks = True # make False to get rid of x axis tick labels (good for arbitrary unit data)
y_tickLabels = True # set to false to remove y axis tick labels
y_ticks = False # False will leave ticks but they won't be manual, they will be [-1, 0, 1]
legend = True # make False if you don't want a legend to show up
savename = 'sample save name' # how to name file that will be saved
skiprows = 0
delimiter = ','
lower_limit = 0 # mT
upper_limit = 50 # mT
Q_BAND = False
Qx = 1 # column for x axis data (python starts at 0)
Qy = 2 # column for y axis data (python starts at 0)
"""
CHANGE STUFF ABOVE
"""

def main(folder):
    fig, ax = plt.subplots()
    lw = 1.25

    if P(folder).is_file():
        folder = P(folder).parent
    data_suffixes = ['.txt', '.dat', '.asc', '.csv']
    files = [ii for ii in P(folder).iterdir() if ii.suffix in data_suffixes]
    idx = 0

    if not y_ticks:
        scale = 0
        for j, f in enumerate(files):
            data = np.loadtxt(f, skiprows=skiprows, delimiter=delimiter)
            data = data[np.logical_and(data[:, 0] >= lower_limit, data[:, 0] < upper_limit)]
            if np.max(np.abs(data[:, 1])) > scale:
                scale = np.max(np.abs(data[:, 1]))
    else:
        scale = 1
    
    for j, f in enumerate(files):
        data = np.loadtxt(f, skiprows=skiprows, delimiter=delimiter)
        data = data[np.logical_and(data[:, 0] >= lower_limit, data[:, 0] < upper_limit)]
        cols = np.shape(data)[1]

        try:
            print(f"Plot number {j + 1} is {f.stem}, label={legend_names[idx]}, color={colors[idx]}, style={styles[idx]}")
        except:
            print(f"Plot number {j + 1} is {f.stem}, styles may be undefined")

        if not Q_BAND:
            for i in range(1, cols):
                try:
                    ax.plot(data[:, 0], data[:, i]/scale, label=legend_names[idx], color=colors[idx], linestyle=styles[idx], lw=lw)
                except IndexError:
                    ax.plot(data[:, 0], data[:, i]/scale, lw=lw, label='no label')
                idx += 1
        else:
            try:
                ax.plot(data[:, Qx], data[:, Qy]/scale, label=legend_names[idx], color=colors[idx], linestyle=styles[idx], lw=lw)
            except IndexError:
                ax.plot(data[:, Qx], data[:, Qy]/scale, lw=lw, label='no label')


    if FIELDS:
        alphabet = 'abcdefghijklmnopqrstuvwxyz'
        alphas = np.linspace(0.25,0.75,len(FIELDS))
        FIELDS.sort()
        ticks = ax.get_ylim()
        n = 1.175
        fact = 0.99
        if np.min(ticks) < 0:
            ax.set_ylim(bottom=np.min(ticks)*n)
        else:
            ax.set_ylim(bottom=np.min(ticks)/n)
            fact = 1/fact
        t = ax.get_ylim()
        for i, f in enumerate(FIELDS):
            ax.axvline(x=f, c='gray',
                       alpha=0.5, lw=lw) # , label=f'$B_0^{{({alphabet[i]})}}$')
            if alphabet[i] == 'a':
                ax.text(f, fact*np.min(t), f'$B_0^{{({alphabet[i]})}}$', horizontalalignment='right', verticalalignment='bottom')        
            elif alphabet[i] == 'b':
                ax.text(f, fact*np.min(t), f'$B_0^{{({alphabet[i]})}}$', horizontalalignment='left', verticalalignment='bottom')        
            else:
                ax.text(f, fact*np.min(t), f'$B_0^{{({alphabet[i]})}}$', horizontalalignment='center', verticalalignment='bottom')        

    ax.text(0.68, 0.68, f'$T={TEMP}$ K\n{MUTANT}',
            horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)
    ax.set_xlabel(x_Label)
    ax.set_ylabel(y_Label)

    if not x_ticks:
        ax.set_xticks([])

    if not x_tickLabels:
        ax.set_xticklabels([])

    if not y_ticks:
        ax.set_yticks([-1, 0, 1])

    if not y_tickLabels:
        ax.set_yticklabels([])

    if legend:
        ax.legend()

    plt.savefig(P(FOLDER).joinpath(f"{savename}.tif"),dpi=300)
    plt.savefig(P(FOLDER).joinpath(f"{savename}.png"),dpi=300)


if __name__ == "__main__":
    main(FOLDER)
    plt.show()
