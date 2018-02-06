import pdb
import time

import matplotlib.pyplot as plt
import numpy as np

from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.lines import Line2D
from util import read_reference_data

from ParentFrame import ParentFrame


class ReferenceFrame(ParentFrame):

    def __init__(self, filename, basedir=None, title=None):

        if title is None:
            title = 'Reference Frame'

        # self.axis will come back as an array with each subplot axis
        super().__init__(title, (2, 1))

        # the time domain axis
        self.time_axis = self.axis[0]
        self.time_axis.grid(True)

        # frequency domain axis
        self.freq_axis = self.axis[1]
        # set axis limits to useful frequency range
        self.freq_axis.set_xlim(0, 3.5)
        self.freq_axis.grid(True)

        del self.axis  # clear up memory

        # the line on the plot that is being held when the user clicks on it
        self.line_held = None

        # the gate that the line being held represents
        self.gate_held = None

        self.time, self.time_amp = read_reference_data(filename, basedir)

        # adjust time vector so it starts at zero, the THz has a global optical
        # delay value that it uses when a waveform is saved as a txt file
        self.time -= self.time[0]

        self.wavelength = len(self.time)

        self.dt = (self.time[-1] / (len(self.time)-1))

        self.freq = np.linspace(0, 1/(self.dt*2), len(self.time)//2+1)

        self.freq_amp = np.fft.rfft(self.time_amp) * self.dt

        self.gate = [0, len(self.time)-1]

        self.plot_time_waveform()
        self.plot_freq_waveform()
        self.connect_events()

    def plot_time_waveform(self):
        """
        Plots the time domain waveform and two movable gates
        """

        lead_gate_time = self.time[self.gate[0]]
        back_gate_time = self.time[self.gate[1]]

        self.time_axis.plot(self.time, self.time_amp, 'r')

        self.lead_gate = self.time_axis.axvline(lead_gate_time, color='k', linestyle='--',
                                                linewidth=1.0, picker=2)
        self.back_gate = self.time_axis.axvline(back_gate_time, color='k', linestyle='--',
                                                linewidth=1.0, picker=2)

        self.time_axis.set_title('Gate: [%d, %d]' % (self.gate[0], self.gate[1]))

        self.figure_canvas.draw()

    def plot_freq_waveform(self):
        """
        Plots the frequency domain waveform
        """
        self.freq_axis.plot(self.freq, np.abs(self.freq_amp), 'r')
        self.figure_canvas.draw()

    def connect_events(self):
        self.figure_canvas.mpl_connect('pick_event', self.grab_gate_handler)
        self.figure_canvas.mpl_connect('button_release_event', self.release_gate_handler)
        self.figure_canvas.mpl_connect('motion_notify_event', self.motion_handler)
        self.figure_canvas.mpl_connect('motion_notify_event', self.slide_gate_handler)

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

        xdata = self.line_held.get_xdata()[0]

        # convert the time location to index
        index = int(round(xdata*len(self.time) / self.time[-1], 0))

        self.gate_held = np.argmin(np.abs(np.asarray(self.gate) - index))

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

        index = int(round(event.xdata / self.dt))

        # handle what happens if the line is released outside of the given time
        if index < 0:
            index = 0
        elif index > self.wavelength - 1:
            index = self.wavelength - 1

        self.gate[self.gate_held] = index

        # create a new waveform to ensure FFT is take of time domain waveform
        # with the correct number of data points
        waveform = np.zeros_like(self.time_amp)
        waveform[self.gate[0]:self.gate[1]] = self.time_amp[self.gate[0]:self.gate[1]]

        self.freq_amp = np.fft.rfft(waveform) * self.dt

        # remove the current line from the frequency axis
        old_line = self.freq_axis.lines.pop(0)
        del old_line

        # set the title on the time axis to show gate position
        self.time_axis.set_title('Gate: [%d, %d]' % (self.gate[0], self.gate[1]))

        self.plot_freq_waveform()

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


if __name__ == '__main__':
    # if program is run as main, open a file dialog so the user can select the
    # reference file that they wish to open

    import sys

    import wx

    # must pass through self parameter to wx.FileDialog, but there is no class
    # for self. Setting it to None seems to work though.
    self = None

    app = wx.App(False)

    with wx.FileDialog(self, 'Open Reference', wildcard='txt files (*.txt)|*.txt',
                       style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST) as dlg:

        if dlg.ShowModal() == wx.ID_CANCEL:  # user clicks cancel
            sys.exit(0)

        full_path = dlg.GetPath()

    frame = ReferenceFrame(full_path)

    app.MainLoop()
