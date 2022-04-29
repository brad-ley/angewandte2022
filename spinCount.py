import ast
import os
from pathlib import Path as P
from pathlib import PurePath as PP
from readDataFile import read

import PIL
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import trapz, cumtrapz

def main(folder):
    """
    Used to computed the baseline-corrected double integral of files in a folder
    """
    if P(folder).is_file():
        folder = P(folder).parent
    fs = [ii for ii in P(folder).iterdir() if ii.suffix == '.txt']
    for i, f in enumerate(fs):
        data = np.loadtxt(f, delimiter=', ')
        y = data[:, 1] - trapz(data[:, 1], data[:, 0]) 
        spins = trapz(cumtrapz(data[:, 1], data[:, 0]), data[1:, 0])
        print(f"{f.stem:<50} count: {spins:.2e}")


if __name__ == "__main__":
    FOLDER = '/Volumes/GoogleDrive/My Drive/Research/Data/2022/2/wildtype compare'
    main(FOLDER)
