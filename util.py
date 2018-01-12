import os

import numpy as np
import pandas as pd


def guassian_curve(x, a, b, c):
    """
    Builds a Guassian curve over the range specified by x.
    :param x: Range over which to build Gaussian Curve
    :param a: Amplitude
    :param b: Controls where center is located
    :param c: Controls how wide the curve is
    :return: the y points
    """

    y = a * np.exp(-(x-b)**2 / (2*c**2))

    return y


def read_reference_data(filename, basedir=None):
    """
    Reads in data from the reference txt file using the pandas library
    :param filename: The filename of the reference file
    :param basedir: The base directory of the filename. If `None`, filename is
        assumed to be the full path to the file.
    :return: optical_delay: An array of optical delay values (time array)
             ref_amp: The amplitude of the reference signal at a given time
    """

    if basedir is not None:
        filename = os.path.join(basedir, filename)
    
    # Read in the reference waveform and separate out the optical delay (time)
    # and the reference amplitude
    reference_data = pd.read_csv(filename, delimiter='\t')
    optical_delay = reference_data['Optical Delay/ps'].values
    ref_amp = reference_data['Raw_Data/a.u.'].values

    return optical_delay, ref_amp


def b_scan_slicer(data):
    """
    Creates an image of each B-Scan slice and saves them to the directory specified by save_loc.
    The direction of the B-Scan is specified by the data.b_scan_dir; either horizontal or vertical.
    data: An instance of a THzData class
    """

    if data.b_scan_dir == 'horizontal':
        b_scan_layers = np.zeros((data.y_step, data.wave_length, data.x_step))
        for i in range(data.y_step):
            data.make_b_scan(i, 0)
            b_scan_layers[i] = data.b_scan

    elif data.b_scan_dir == 'vertical':
        b_scan_layers = np.zeros((data.x_step, data.wave_length, data.y_step))
        for i in range(data.x_step):
            data.make_b_scan(0, i)
            b_scan_layers[i] = data.b_scan

    else:
        raise ValueError("B-Scan direction must be either 'horizontal' or 'vertical'")

    return b_scan_layers
