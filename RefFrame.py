import pdb

import wx
import matplotlib.pyplot as plt
import numpy as np

from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.lines import Line2D
from util import read_reference_data

from ParentFrame import ParentFrame


class ReferenceFrame(ParentFrame):

    def __init__(self, filename, basedir='None', title=None):

        if title is None:
            title = 'Reference Frame'

        super().__init__(title)

        self.time_axis = None
        self.freq_axis = None
        self.line_held = None
        self.gate_held = None

        self.time, self.time_amp = read_reference_data(filename, basedir)

        self.time -= self.time[0]

        self.dt = (self.time[-1] / (len(self.time)-1))

        self.freq = np.linspace(0, 1/(self.dt*2), len(self.time)+1)

        self.freq_amp = np.fft.rfft(self.time_amp) * self.dt

        self.gate = [0, len(self.time)]

        self.plot_time_waveform()
        self.plot_freq_waveform()

    # override from ParentFrame
    def initialize_figure(self):
        """
        Creates a figure, adds two axes to that figure (one for time domain
        information and one for frequency information). Also initializes a
        figure canvas.
        """
        # override the method in ParentFrame because we want to add two axes
        # to this plot

        self.figure = plt.figure()
        self.time_axis = self.figure.add_subplot(211)
        self.freq_axis = self.figure.add_subplot(212)

        self.freq_axis.set_xlim(0, 3)

        self.figure.suptitle('Reference Waveform')

        self.time_axis.set_xlabel('Time (ps)')
        self.time_axis.set_ylabel('Amplitude')

        self.freq_axis.set_xlabel('Frequency (THz)')
        self.freq_axis.set_ylabel('Amplitude')

        self.figure_canvas = FigureCanvas(self, -1, self.figure)

    def plot_time_waveform(self):
        """
        Plots the time domain waveform and two movable gates
        """

        lead_gate_time = self.time[-1] * 0.1
        back_gate_time = self.time[-1] * 0.5

        self.time_axis.plot(self.time, self.time_amp, 'r')

        self.lead_gate = self.time_axis.axvline(lead_gate_time, color='k', linestyle='--', linewidth=1.0, picker=2)
        self.back_gate = self.time_axis.axvline(back_gate_time, color='k', linestyle='--', linewidth=1.0, picker=2)

        self.figure_canvas.draw()

    def plot_freq_waveform(self):
        """
        Plots the frequency domain waveform
        """
        self.freq_axis.plot(self.freq, self.freq_amp, 'r')
        self.figure_canvas.draw()

    def motion_handler(self, event):
        """
        Prints the current (x,y) values that the mouse is over to the status bar along with the
        amplitude of the waveform at that time.
        :param event: matplotlib motion notify event
        """
        xid = event.xdata
        yid = event.ydata

        # make sure that the mouse is actually inside of the axis, if it is not
        # xid, and yid will be None
        if xid is not None and yid is not None:
            status_string = '(%.4f, %.4f)' % (xid, yid)
            self.status_bar.SetStatusText(status_string)
        else:
            self.status_bar.SetStatusText('')

    def grab_gate_handler(self, event):
        """
        Stores that line that is clicked on and makes it bold for visibility
        :param event: matplotlib mouse click event
        """
        # if the user does not click on one of the gate lines
        if not isinstance(event.artist, Line2D):
            return

        # store the line that is grabbed and make it bold for visibility
        self.line_held = event.artist
        self.line_held.set_linewidth(2.0)

        xdata = event.xdata

        # convert the time location to index
        index = int(round(xdata*len(self.time) / self.time[-1], 0))

        self.gate_held = np.argmin(np.abs(self.gate - index))

    def release_gate_handler(self, event):
        """
        Updates the gate based on where the mouse was released
        :param event: matplotlib mouse release event
        """

        # if the user has not yet clicked on a line, do nothing
        if self.line_held is None:
            return

        # if the user does not release in the image, exit
        if event.xdata is None or event.ydata is None:
            return

        # make normal width again
        self.line_held.set_linewidth(1.0)

        # reset line handler
        self.line_held = None

        xdata = event.xdata
        index = int(round(event.xdata))
        self.gate[self.gate_held] = index

        self.freq_amp = np.fft.rfft(self.time_amp[self.gate[0]:self.gate[1]]) * self.dt

    def slide_gate_handler(self, event):
        """
        Moves the gate with the mouse after the user has grabbed it
        :param event: matplotlib motion notify event
        """

        if self.line_held is None:
            return

        x = event.xdata

        self.line_held.set_xdata([x, x])
        self.figure_canvas.draw()
