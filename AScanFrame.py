import pdb

import numpy as np

from ParentFrame import ParentFrame


class AScanFrame(ParentFrame):
    def __init__(self, holder, data):
        ParentFrame.__init__(self, 'Raw A-Scan Test')

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
        self.figure_canvas.mpl_connect('button_press_event', self.grab_gate_handler)
        self.figure_canvas.mpl_connect('button_release_event', self.release_gate_handler)
        self.figure_canvas.mpl_connect('motion_notify_event', self.motion_handler)

    def plot(self, i, j):
        """
        Plot a waveform from data matrix location (i, j, :)
        :param i: The row from which to plot the waveform
        :param j: The column from which to plot the waveform
        """
        if not self.is_initialized:
            self.is_initialized = True

        # store the (i, j) that were clicked on
        self.i_index = i
        self.j_index = j

        title_string = 'Point Clicked: i=%d, j=%d' % (self.i_index, self.j_index)

        self.axis.cla()
        self.axis.plot(self.data.time, self.data.waveform[i, j, :], 'r')
        self.axis.axvline(self.data.time[self.data.gate[0][0]], color='k', linestyle='--', picker=2)
        self.axis.axvline(self.data.time[self.data.gate[0][1]], color='k', linestyle='--', picker=2)
        self.axis.set_xlabel('Time (ps)', fontsize=14)
        self.axis.set_ylabel('Amplitude', fontsize=14)
        # self.axis.set_title(title_string, fontsize=14)
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

            status_string = '(%.6f, %.6f), [%.6f]' % (xid, yid, point_amp)
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

        if event.xdata is not None and event.ydata is not None:
            print('You clicked', event.xdata, event.ydata)
            index = event.xdata * self.data.wave_length / self.data.time_length  # convert to index

            print('Index grabbed =', index)
            # this is the gate to replace
            self.data.gate_index = np.argmin(np.abs(np.asarray(self.data.gate) - index))

            print('Gate Index grabbed =', self.data.gate_index)

    def release_gate_handler(self, event):
        """
        Determines where the new gate will be when mouse is released.
        Updates A-Scan with new gate line
        Updates C-Scan with images between new gate values
        """
        # if no point on the C-Scan has been clicked yet, do nothing and return
        if not self.is_initialized:
            return

        if event.xdata is not None and event.ydata is not None:
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
            self.data.make_c_scan()  # Generate the new C-Scan data based on gate location

            self.plot(self.i_index, self.j_index)  # update A-Scan to show new gate
            self.holder.raw_c_scan_frame.update()
            self.holder.interpolated_c_scan_frame.update()
