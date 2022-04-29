import ast
import os
from pathlib import Path as P
from pathlib import PurePath as PP

import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
import PIL
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter

from readDataFile import read

plt.style.use('science')
rc('text.latex', preamble=r'\usepackage{cmbright}')
plt.rcParams['font.family'] = 'sans-serif'


def main(filename):
    process(
        filename,
        on=5,
        off=175,
        window_frac=0,
        order=2,
        bi=False
    )


def biexponential(x, A, a1, t1, a2, t2):
    return A + a1 * np.exp(- x / t1) + a2 * np.exp(-x / t2)


def exponential(x, A, a1, t1):
    return A + a1 * np.exp(- x / t1)


def process(filename, on, off, window_frac=10, order=2, bi=True):
    """process.

    :param filename: input filename
    :param on: on-time for laser
    :param off: off-time for laser
    :window_frac: reciprocal of data length for window choice
    :order: savgol_filter polyfit order
    :bi: biexponential fit
    """
    header, data = read(filename)
    t = data[:, 0]
    r = data[:, 2] + 1j * data[:, 3]
    i = data[:, 4] + 1j * data[:, 5]
    sig = np.abs(r) + 1j * np.abs(i)

    loops = 0
    avgsig = []
    avgr = []
    avgi = []

    prevlen = 0
    smoothlen = sum(t < on + off)
    # smoothlen = 5000

    print(f"{t[-1] / (on + off):.2f} loops detected")

    while loops < t[-1] // (on + off):
        # try:
        lo = np.where(t >= loops * (on + off))[0][0]
        hi = np.where(t < (loops + 1) * (on + off))[0][-1]

        s = sig[lo:hi]
        ar = np.real(r[lo:hi]) + 1j * np.imag(r[lo:hi])
        ai = np.real(i[lo:hi]) + 1j * np.imag(i[lo:hi])

        tt = t[lo:hi]

        spacet = np.linspace(t[lo], t[hi - 1], smoothlen)

        fs = interp1d(tt, s)
        far = interp1d(tt, ar)
        fai = interp1d(tt, ai)

        s = fs(spacet)
        ar = far(spacet)
        ai = fai(spacet)

        if window_frac == 0:  # use 0 to not do savgol filtering
            avgsig.append(s)
            avgr.append(ar)
            avgi.append(ai)
        else:
            avgsig.append(savgol_filter(
                s, (2 * (np.abs(hi - lo) // window_frac) + 1), order))
            avgr.append(savgol_filter(
                np.real(ar), (2 * (np.abs(hi - lo) // window_frac) + 1), order) + 1j * savgol_filter(
                np.imag(ar), (2 * (np.abs(hi - lo) // window_frac) + 1), order))
            avgi.append(savgol_filter(
                np.real(ai), (2 * (np.abs(hi - lo) // window_frac) + 1), order) + 1j * savgol_filter(
                np.imag(ai), (2 * (np.abs(hi - lo) // window_frac) + 1), order))
        # except:
        #     print('excepted')
        #     pass

        loops += 1

    smootht = np.linspace(0, on + off, smoothlen)

    avgsig1 = [np.abs(np.real(ii)) for ii in avgsig]
    avgsig2 = [np.abs(np.imag(ii)) for ii in avgsig]

    # plots = {"Ch1 mag": avgsig1, "Ch2 mag": avgsig2, "Ch1 real": avgr, "Ch1 imag":
    #          avgr, "Ch2 real": avgi, "Ch2 imag": avgi}
    plots = {"Ch1 mag": avgsig1, "Ch2 mag": avgsig2}

    lw = 1.25

    for plot, dat in plots.items():
        if 'real' in plot:
            dat = np.real(dat)
        elif 'imag' in plot:
            dat = np.imag(dat)

        plt.figure(plot)

        full = np.mean(dat, axis=0)

        bl_guess = 1e-5  # starting point of decay
        amp_guess = 1e-5  # amplitude of exponential decay
        tau_guess = 50  # time constant guess

        # if bi:
        #     popt, pcov = curve_fit(biexponential, smootht[smootht > on], full[smootht > on], maxfev=10000000, p0=[
        #                            bl_guess, amp_guess, tau_guess, amp_guess, tau_guess])
        #     plt.plot(smootht[smootht > on], biexponential(smootht, *popt), color="red", linestyle="--", lw=lw,
        #              label=r"fit: $\tau_1$=" + f"{popt[2]:.1f}" + r" s$^{-1}$;$\tau_2$=" + f"{popt[4]:.1f}" + r" s$^{-1}$")
        # else:
        popt, pcov = curve_fit(
            exponential, smootht[smootht > on], full[smootht > on], p0=[bl_guess, amp_guess, tau_guess])
        perr = np.sqrt(np.diag(pcov))
        sd2 = 2*perr[2]
        if sd2 == np.inf:
            sd2 = 0
        if popt[1] < 0:
            m = -1
            c = 2*popt[0]
            full -= c
            full *= m
            popt, pcov = curve_fit(
                exponential, smootht[smootht > on], full[smootht > on], p0=[bl_guess, amp_guess, tau_guess])
            perr = np.sqrt(np.diag(pcov))
        else:
            m = 1
            c = 0

        for i, row in enumerate(dat):
            row -= c
            row *= m

            if i == 0:
                plt.plot(smootht, row / np.max(full),
                         alpha=0.3, color='lightgray', label="Single scans", lw=lw)
            else:
                plt.plot(smootht, row / np.max(full),
                         alpha=0.3, color='lightgray', lw=lw)

        plt.plot(smootht, full/np.max(full), label="Average", c='k', lw=lw)
        fit = exponential(smootht[smootht > on], *popt)
        plt.plot(smootht[smootht > on], fit/np.max(full), color="red", linestyle="--", label=fr"$\tau={popt[2]:.1f}\pm{sd2:.1f}$ s", lw=lw)

        
        plt.axvspan(0, on, facecolor='#00A7CA', label='Laser on')

        handles, labels = plt.gca().get_legend_handles_labels()

        # specify order of items in legend
        order = [0, 1, 2, 3]

        # add legend to plot
        plt.legend([handles[idx] for idx in order], [labels[idx]
                                                     for idx in order],
                                                     loc='upper right')

        # plt.legend()
        rang = np.max(full) - np.min(full)
        # plt.annotate('Laser\npulse', (on, np.min(full) - rang/6),
        #              color='black', horizontalalignment='left')

        # plt.ylim(np.min(full) - rang / 4, np.max(full) + rang / 4)

        # if plot != "Ch2 mag":
        #     plt.title(plot)

        # if plot == "Ch2 imag":
        #     plt.title("Time-resolved EPR absoption vs. time")

        # out = dict(name='name here', time=list(smootht), avg=list(full), fit=list(fit), tc=popt[2], te=perr[2])
        outstr = ''
        for i, v in enumerate(full):
            outstr += f"{smootht[i]}, {v}\n"

        P(filename).parent.joinpath(f"output {plot}.txt").write_text(outstr)

        plt.ylabel('Intensity (arb. u)')
        plt.xlabel('Time (s)')

        plt.savefig(P(filename).parent.joinpath(f"{plot}.png"), dpi=300)
        plt.savefig(P(filename).parent.joinpath(f"{plot}.tif"), dpi=300)


if __name__ == "__main__":
    f = '/Volumes/GoogleDrive/My Drive/Research/Data/2022/4/26/try 2/M14_pulsing_rephased.dat'
    main(f)
    plt.show()
