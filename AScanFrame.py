import pdb
import copy

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
        self.i_index = int
        self.j_index = int

        # value that allows the gate that is grabbed to be stored so it can be moved
        # 0 = front lead gate
        # 1 = back lead gate
        # 2 = front follow gate
        # 3 = back follow gate
        self.gate_held = None

        # stores the line object that is clicked on when moving the gates
        self.line_held = None

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

        xloc = self.data.x[self.j_index]
        yloc = self.data.y[self.i_index]

        title_string = 'Location: x=%0.2f, y=%0.2f' % (xloc, yloc)

        self.axis.cla()
        self.axis.plot(self.data.time, self.data.waveform[i, j, :], 'r')
        self.axis.set_xlabel('Time (ps)')
        self.axis.set_ylabel('Amplitude')
        self.axis.set_title(title_string, fontsize=14)
        self.axis.grid()

        # if the follow gate is off and signal type is 1, program uses pk to pk values across the
        # entire waveform, so plotting gates is not necessary
        if not self.data.follow_gate_on and self.data.signal_type == 1:
            self.figure_canvas.draw()
            return

        # plot the lead gates
        self.axis.axvline(self.data.time[self.data.gate[0][0]], color='k', linestyle='--',
                          linewidth=1.0, picker=2)
        self.axis.axvline(self.data.time[self.data.gate[0][1]], color='k', linestyle='--',
                          linewidth=1.0, picker=2)

        # if follow gate is on and being used with current signal type; plot them
        # signal type = 0 is to use the pk to pk voltage within the lead gates
        if self.data.follow_gate_on and self.data.signal_type != 0:
            followL_idx = self.data.peak_bin[3, 1, i, j]
            followR_idx = self.data.peak_bin[4, 1, i, j]
            self.axis.axvline(self.data.time[followL_idx], color='b', linewidth=1.0, picker=2)
            self.axis.axvline(self.data.time[followR_idx], color='g', linewidth=1.0, picker=2)

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
        self.line_held = event.artist
        xdata = line.get_xdata()[0]

        # convert xdata of point clicked on to  nearest time index
        index = int(round(xdata * self.data.wave_length / self.data.time_length, 0))

        # determine which gate was clicked on
        diff = np.abs(index - self.data.gate[0][0])  # check front lead gate
        diff1 = np.abs(index - self.data.gate[0][1])  # check back lead gate

        if diff1 > diff:
            self.gate_held = 0  # clicked on the lead lead gate
        else:
            self.gate_held = 1  # clicked on the back front gate
            diff = diff1

        # if follow gate is on we need to check the follow gates also
        # use peak_bin to check with the follow gate location for this particular (i, j) pair
        if self.data.follow_gate_on:
            # difference between point clicked and lead follow gate
            diff2 = np.abs(index - self.data.peak_bin[3, 1, self.i_index, self.j_index])

            # back follow gate
            diff3 = np.abs(index - self.data.peak_bin[4, 1, self.i_index, self.j_index])

            # determine if either is closer than the front gate
            if diff2 < diff:  # point clicked is closest to lead follow gate
                self.gate_held = 2
                diff = diff2
            if diff3 < diff:  # point clicked is closest to back follow gate
                self.gate_held = 3

        print()
        print('Gate Grabbed =', self.gate_held)
        print()

        line.set_linewidth(2.0)  # make the line clicked on bold
        self.figure_canvas.draw()  # update plot

    def release_gate_handler(self, event):
        """
        Determines where the new gate will be when mouse is released.
        Updates A-Scan with new gate line
        Updates C-Scan with images between new gate values
        """
        # if no point on the C-Scan has been clicked yet, do nothing and return
        if not self.is_initialized:
            return

        # if the user has not yet clicked on a line, do nothing
        if self.line_held is None:
            return

        # if the user does not release in the image, exit
        if event.xdata is None or event.ydata is None:
            return

        # convert time value (xdata) to index
        index = int(round(event.xdata * self.data.wave_length / self.data.time_length, 0))

        # ensure that the gate is inside of the bounds
        if index < 0:
            index = 0
        elif index > self.data.wave_length:
            index = self.data.wave_length

        new_gate = copy.deepcopy(self.data.gate)

        # if we are adjusting the lead gates, we just use the index as is
        if self.gate_held == 0:  # front lead gate
            new_gate[0][0] = index
        elif self.gate_held == 1:  # back lead gate
            new_gate[0][1] = index
        # to adjust the follow gate, we need to use new index in relation to the old one found in
        # peak_bin
        elif self.gate_held == 2:  # front follow gate
            new_gate[1][0] = self.data.gate[1][0] - \
                (self.data.peak_bin[3, 1, self.i_index, self.j_index] - index)
        else:  # gate_held = 3; back follow gate
            new_gate[1][1] = self.data.gate[1][1] - \
                (self.data.peak_bin[4, 1, self.i_index, self.j_index] - index)

        # if follow gate is on, use the change gate method, which updates bin_range, and calls
        # find_peaks function to update peak_bin
        # if a value error is encountered it usually means that the left follow gate is greater
        # than the last time length index
        if self.data.follow_gate_on:
            try:
                self.data.change_gate(new_gate)
                self.line_held.set_xdata([event.xdata, event.xdata])
                self.line_held.set_linewidth(1.0)  # return to normal width
                self.figure_canvas.draw()
            except ValueError:
                # return to normal width only, don't update position
                print('\nNeed to adjust follow gates first!\n')
                self.line_held.set_linewidth(1.0)
                self.figure_canvas.draw()
                return
        else:  # otherwise change the gate manually, this is faster
            self.data.gate = copy.deepcopy(new_gate)

        # Generate the new C-Scan data based on the new gate location
        self.data.make_c_scan(self.data.signal_type)

        self.plot(self.i_index, self.j_index)  # update A-Scan to show new gate
        self.holder.raw_c_scan_frame.update()
        self.holder.interpolated_c_scan_frame.update()

        # set line held to None, so user must first click on a line to move it
        # updated to be not None in grab_gate_handler()
        self.line_held = None
