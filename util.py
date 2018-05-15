import pdb
import numpy as np


def b_scan_slicer(data):
    """
    Creates an image of each B-Scan slice and saves them to the directory 
    specified by save_loc. The direction of the B-Scan is specified by the 
    data.b_scan_dir; either horizontal or vertical.
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
