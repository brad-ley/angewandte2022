from pathlib import Path as P

import numpy as np


def isNumber(num):
    """
    :param s:
    :return:
    stackexchange function to check if user input is a feasible number
    """
    try:
        float(num)

        return True
    except ValueError:
        return False


def read(filename, delimiter=',', flipX=True):
    """
    Takes file from EPR computer, removes header, returns header, numpy array
    :args: filename
    :kwargs: delimiter
    :return: header, data
    """
    file = P(filename).read_text().split("\n")

    header = ''
    found = False

    for line in file:
        if line.startswith('[Data]'):
            data_idx = file.index(line)
            header += line + "\n"
            header += file[file.index(line) + 1] + "\n"
            skiprows = data_idx + 2
            found = True
        elif all([isNumber(ii) for ii in line.split(delimiter)]):
            data_idx = file.index(line)
            skiprows = data_idx
            found = True

        if found:
            break
        header += line + "\n"

    idx_list = []
    datatypes = []
    data = np.loadtxt(filename, delimiter=delimiter, skiprows=skiprows)

    for line in header.split('\n')[::-1]:
        if line != '':
            datatypes = line
            break
    
    idx = 0
    if datatypes:
        idx_dict = {'field': ['field' in ii.strip().lower()
                              for ii in datatypes.split(delimiter)],
                              'time': ['time' in ii.strip().lower()
                                  for ii in datatypes.split(delimiter)]}

        if True in idx_dict['field']:
            idx = idx_dict['field'].index(1)
        elif True in idx_dict['time']:
            idx = idx_dict['time'].index(1)

    if flipX:
        if data[0, idx] > data[-1, idx]:
            data = np.flipud(data)

    return header, data
