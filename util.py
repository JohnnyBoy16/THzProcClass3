import numpy as np


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