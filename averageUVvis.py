from pathlib import Path as P
from scipy.optimize import curve_fit

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

FOLDER = '/Volumes/GoogleDrive/My Drive/Research/protein-dynamics/manuscripts/methods 2022/data/uv-vis/Q513A figure'
p0 = [1, -1, 300] # initial guesses

def exp(x, a, b, c):
    return a + b*np.exp(-x / c)

def average(folder):
    if P(folder).is_file():
        folder = P(folder).parent
    excel_sheets = [ii for ii in P(folder).iterdir() if ii.suffix == '.xlsx' and not ii.name.startswith('~')]
    data = {}
    for e in excel_sheets:
        d = pd.read_excel(e) 
        data[e.stem] = d

    fig, ax = plt.subplots()
    
    for d in data:
        t_col = [ii for ii in list(data[d].columns.values) if 'time' in ii.lower()][0]
        others = [ii for ii in list(data[d].columns.values) if ii != t_col and 'Unnamed' not in ii]
        t = data[d][t_col].to_numpy()
        data_dict = {}
        for i, o in enumerate(others):
            data[d][o] = data[d][o] / np.max(data[d][o])
            data[d][o] = data[d][o].to_numpy()
            time = t[~np.isnan(data[d][o])]
            y = data[d][o][~np.isnan(data[d][o])]
            data_dict[o] = [time, y] 
            ax.plot(time, y, alpha=0.5, color='gray')
            popt, pcov = curve_fit(exp, time, y, p0=[1, -1, 100])
            print(f'{o}: {popt[-1]:.2f} s')

        ## averaging with different lengths of time acq ##
        avg = {}
        for key in data_dict:
            for i, val in enumerate(data_dict[key][0]):
                if val in avg:
                    avg[val][0] += 1
                    avg[val][1] += data_dict[key][1][i]
                else:
                    avg[val] = [1, data_dict[key][1][i]]
    
        pdat = np.zeros((len(avg.keys()), 2))
        for v, i in enumerate(avg):
            pdat[v, 0] = v
            pdat[v, 1] = avg[v][1] / avg[v][0]

        # avg = np.mean(mat, axis=1)
        popt, pcov = curve_fit(exp, pdat[:, 0], pdat[:, 1])
        print(f"average of {popt[-1]:.1f} s")
        ax.plot(pdat[:, 0], pdat[:, 1])
        outstr = ''
        for i, row in enumerate(pdat):
            outstr += f"{row[0]},{row[1]}\n"
        P(folder).joinpath('Averaged UV-Vis.txt').write_text(outstr)
        plt.savefig(P(folder).joinpath('averaged.png'), dpi=300)

if __name__ == "__main__":
    average(FOLDER)
    # plt.show()
