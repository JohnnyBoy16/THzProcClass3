import pdb

import numpy as np
from matplotlib.lines import Line2D

from ParentFrame import ParentFrame


class AScanFrame(ParentFrame):
    """
    Displays the A-Scan of the last point clicked in the gray scale C-Scan image. Also contains
    the movable gates that update the C-Scan.
    """
    def __init__(self, holder, data, title=None):
        # TODO add functionality for frequency domain plots

        if title is None:
            title = 'A-Scan Frame'

        super().__init__(title)  # call parent class

        self.holder = holder

        # The THz Data
        self.data = data

        # the last i (row) and j (column) coordinates to be clicked on
        self.i_index = None
        self.j_index = None

        # tells itself when it has plotted the first time,
        # otherwise motion_handler() throws errors if nothing has been plotted yet
        self.is_initialized = False

        self.connect_events()

    def connect_events(self):
        self.figure_canvas.mpl_connect('pick_event', self.grab_gate_handler)
        self.figure_canvas.mpl_connect('button_release_event', self.release_gate_handler)
        self.figure_canvas.mpl_connect('motion_notify_event', self.motion_handler)

    def plot(self, i, j):
        """
        Plot a waveform from data matrix location (i, j, :)
        :param i: The row from which to plot the waveform
        :param j: The column from which to plot the waveform
        """

        if not self.is_initialized:  # set value, so motion handler method can function
            self.is_initialized = True

        # store the (i, j) that were clicked on
        self.i_index = i
        self.j_index = j

        title_string = 'Point Clicked: i=%d, j=%d' % (self.i_index, self.j_index)

        self.axis.cla()
        self.axis.plot(self.data.time, self.data.waveform[i, j, :], 'r')
        self.axis.axvline(self.data.time[self.data.gate[0][0]], color='k', linestyle='--',
                          linewidth=1.0, picker=2)
        self.axis.axvline(self.data.time[self.data.gate[0][1]], color='k', linestyle='--',
                          linewidth=1.0, picker=2)

        # if follow gate is active; plot them
        if self.data.follow_gate_on:
            follow1_idx = self.data.peak_bin[3, 1, i, j]
            follow2_idx = self.data.peak_bin[4, 1, i, j]
            self.axis.axvline(self.data.time[follow1_idx], color='b', linewidth=1.0, picker=2)
            self.axis.axvline(self.data.time[follow2_idx], color='g', linewidth=1.0, picker=2)

        self.axis.set_xlabel('Time (ps)')
        self.axis.set_ylabel('Amplitude')
        self.axis.set_title(title_string, fontsize=14)
        self.axis.grid()
        self.figure_canvas.draw()

    def motion_handler(self, event):
        """
        Prints the current (x,y) values that the mouse is over to the status bar along with the
        amplitude of the waveform at that time.
        """

        # if no C-Scan point has been clicked on yet, do nothing
        if not self.is_initialized:
            return

        xid = event.xdata
        yid = event.ydata
        if xid is not None and yid is not None:
            if xid < 0:
                t_index = 0
            elif xid > self.data.time_length:
                t_index = self.data.wave_length - 1
            else:
                t_index = int(round(xid / self.data.delta_t, 0))

            point_amp = self.data.waveform[self.i_index, self.j_index, t_index]

            status_string = '(%.4f, %.4f), [%.4f]' % (xid, yid, point_amp)
            self.status_bar.SetStatusText(status_string)
        else:
            self.status_bar.SetStatusText('')

    def grab_gate_handler(self, event):
        """
        Determines that gate that is clicked on
        """
        # if no point on the C-Scan has been clicked on yet, do nothing and return
        if not self.is_initialized:
            return

        # if the user does not click on one of the gate lines
        if not isinstance(event.artist, Line2D):
            return

        line = event.artist
        self.gate_held = event.artist
        xdata = line.get_xdata()[0]
        print(len(line.get_xdata()))

        line.set_linewidth(2.0)  # make the line clicked on bold
        self.figure_canvas.draw()

    def release_gate_handler(self, event):
        """
        Determines where the new gate will be when mouse is released.
        Updates A-Scan with new gate line
        Updates C-Scan with images between new gate values
        """
        # if no point on the C-Scan has been clicked yet, do nothing and return
        if not self.is_initialized:
            return

        # if the user does not click in the image, exit
        if event.xdata is None or event.ydata is None:
            return

        self.gate_held.set_xdata([event.xdata, event.xdata])
        self.gate_held.set_linewidth(1.0)  # return to normal width
        self.figure_canvas.draw()

        print('You released at', event.xdata, event.ydata)
        # convert to index
        index = int(round(event.xdata * self.data.wave_length / self.data.time_length, 0))
        print('Index released at =', index)
        # ensure that the gate is inside of the bounds
        if index < 0:
            index = 0
        elif index > self.data.wave_length:
            index = self.data.wave_length
        print(index)
        self.data.gate[0][self.data.gate_index] = index

        print('New gate value =', self.data.gate)

        if self.data.follow_gate_on:
            self.data.find_peaks()  # call find peaks to update peak_bin if follow gate is on

        # Generate the new C-Scan data based on the new gate location
        self.data.make_c_scan(self.data.signal_type)

        self.plot(self.i_index, self.j_index)  # update A-Scan to show new gate
        self.holder.raw_c_scan_frame.update()
        self.holder.interpolated_c_scan_frame.update()
