import pdb

import matplotlib.pyplot as plt

from InterpolatedCScanFrame import InterpolatedCScanFrame


class ZoomedCScanFrame(InterpolatedCScanFrame):
    def __init__(self, holder, data):

        super().__init__(holder, data, title='Zoomed C-Scan Frame')

    # Override
    def plot(self):
        """
        Plots a zoomed C-Scan image so the colorbar can focus on the sample itself instead of
        background noise
        """
