import pdb
import copy

import wx
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from THzProc.ParentFrame import ParentFrame


class AScanFrame(ParentFrame):
    """
    Displays the A-Scan of the last point clicked in the gray scale C-Scan
    image. Also contains the movable gates that update the C-Scan.
    """

    def __init__(self, holder, data, a_scan_only=True, title=None):
        """
        Constructor method
        :param holder: an instance of FrameHolder, used to link actions in
                one frame to actions in another frame
        :param data: an instance of the THzData class
        :param title: the title for the frame
        """

        if title is None:
            title = 'A-Scan Frame'

        if a_scan_only:
            subplot_grid = (1, 1)  # only need one axis
        else:
            subplot_grid = (2, 1)  # need more than one axis

        # call inherited class constructor, specify AScanFrame and self as
        # arguments so this runs in python 2
        super(AScanFrame, self).__init__(title, subplot_grid)

        self.holder = holder

        # The THz Data
        self.data = data

        if a_scan_only:
            self.time_axis = self.axis
            self.freq_axis = None
        else:
            # declare and leave around in case user wants to look at frequency
            # information later
            self.time_axis = self.axis[0]
            self.freq_axis = self.axis[1]

        del self.axis  # free up memory for this class instance

        # the last i (row) and j (column) coordinates to be clicked on
        self.i_index = int
        self.j_index = int

        # controls whether or not the title is displayed above the graph
        self.show_title = True

        # value that allows the gate that is grabbed to be stored so it can be
        # moved
        # 0 = front lead gate
        # 1 = back lead gate
        # 2 = front follow gate
        # 3 = back follow gate
        self.gate_held = None

        self.gate0 = None
        self.gate1 = None

        # wx menu item that allows user to switch the A-Scan view to include
        # frequency information or not
        self.a_scan_view_menu = None

        # wx menu item that allows the user to switch between using the (i, j)
        # or (x, y) coordinate system
        self.location_view_menu = None

        # wx menu item that allows the user to switch between having the title
        # visible or not
        self.toggle_title_menu = wx.MenuItem

        # stores the line object that is clicked on when moving the gates
        self.line_held = None

        # tells itself when it has plotted the first time, otherwise
        # motion_handler() throws errors if nothing has been plotted yet
        self.is_initialized = False

        # store whether or not frequency domain plots are visible
        self.a_scan_only = a_scan_only

        # store whether or not we are showing (i, j) or (x, y) indexing
        # default to False because we generally want to see (x, y) coords to
        # match with physical location instead of index
        self.ij_indexing = False

        self.freq_gate = np.zeros(2, dtype=int)
        self.freq_gate[0] = np.argmin(np.abs(self.data.freq - 0.05))
        self.freq_gate[1] = np.argmin(np.abs(self.data.freq - 3.0))

        self.modify_menu()
        self.connect_events()

        # this will prevent a Matplotlib figure from popping up if the user
        # creates other figures and calls plt.show() after initializing the
        # application
        plt.close(self.figure)

    def connect_events(self):
        """
        Binds (connects) events to their appropriate method handlers
        """
        # matplotlib events
        self.figure_canvas.mpl_connect('pick_event', self.grab_gate_handler)
        self.figure_canvas.mpl_connect('button_release_event',
                                       self.release_gate_handler)
        self.figure_canvas.mpl_connect('motion_notify_event',
                                       self.motion_handler)
        self.figure_canvas.mpl_connect('motion_notify_event',
                                       self.gate_slider)

        # wxPython events
        self.Bind(wx.EVT_MENU, self.switch_a_scan_view, self.a_scan_view_menu)
        self.Bind(wx.EVT_MENU, self.change_index_system, self.index_view_menu)
        self.Bind(wx.EVT_MENU, self.on_change_gate_button,
                  self.change_gate_menu)
        self.Bind(wx.EVT_MENU, self.on_toggle_title_button,
                  self.toggle_title_menu)

    def modify_menu(self):
        """
        Modifies the menu bar to include an options drop down menu with options
        that are specific to the A-Scan
        """
        options_menu = wx.Menu()

        title = 'Toggle A-Scan Only'
        description = 'Turn A-Scan on/off'
        self.a_scan_view_menu = wx.MenuItem(options_menu, wx.ID_ANY, title,
                                            description)
        options_menu.Append(self.a_scan_view_menu)

        title = 'Indexing Option'
        description = 'Change indexing from (x, y) to (i, j)'
        self.index_view_menu = wx.MenuItem(options_menu, wx.ID_ANY, title,
                                           description)
        options_menu.Append(self.index_view_menu)

        self.change_gate_menu = wx.MenuItem(options_menu, wx.ID_ANY,
                                            'Change Gate', 'Change Gate')
        options_menu.Append(self.change_gate_menu)

        title = 'Toggle Title Visibility'
        description = 'Toggle Title Visibility on/off'
        self.toggle_title_menu = wx.MenuItem(options_menu, wx.ID_ANY,
                                             title, description)

        options_menu.Append(self.toggle_title_menu)

        self.menu_bar.Append(options_menu, '&Options')

        self.SetMenuBar(self.menu_bar)
        self.Fit()

    def plot(self, i=None, j=None):
        """
        Plot a waveform from data matrix location (i, j, :)
        :param i: The row from which to plot the waveform
        :param j: The column from which to plot the waveform
        """

        # set value, so motion handler method can function
        if not self.is_initialized:
            self.is_initialized = True

        if i is None and j is None:
            i = self.i_index
            j = self.j_index
        else:
            # store the (i, j) that were clicked on, this way the same line
            # can be replotted when the gate is changed
            self.i_index = i
            self.j_index = j

        self.time_plot(i, j)

        if not self.a_scan_only:
            self.freq_plot(i, j)

        # have call to draw() here as opposed to being in time_plot() or
        # freq_plot() so it will always be called regardless of whether or not
        # frequency domain plotting is on
        self.figure_canvas.draw()

    def time_plot(self, i, j):

        xloc = self.data.x[j]
        yloc = self.data.y[i]

        # can have the title give the (x, y) location or (i, j) location
        # this could also be updated to show both if one wanted to
        if self.ij_indexing:
            title_string = 'Location: i=%d, j=%d' % (i, j)
        else:
            title_string = 'Location: x=%0.2f, y=%0.2f' % (xloc, yloc)

        title_string += ', Follow Gate: [%d, %d]' % (self.data.gate[1][0],
                                                     self.data.gate[1][1])

        self.time_axis.cla()
        self.time_axis.plot(self.data.time, self.data.waveform[i, j, :], 'r')
        self.time_axis.set_xlabel('Time (ps)')
        self.time_axis.set_ylabel('Amplitude')
        self.time_axis.grid()

        if self.show_title:
            self.time_axis.set_title(title_string)
        else:
            self.time_axis.set_title('')

        # if the follow gate is off and signal type is 1, program uses pk to pk
        # values across the entire waveform, so plotting gates is not necessary
        if not self.data.follow_gate_on and self.data.signal_type == 1:
            return

        # plot the lead gates
        self.time_axis.axvline(self.data.time[self.data.gate[0][0]], color='k',
                               linestyle='--', linewidth=1.0, picker=2)
        self.time_axis.axvline(self.data.time[self.data.gate[0][1]], color='k',
                               linestyle='--', linewidth=1.0, picker=2)

        # if follow gate is on and being used with current signal type; plot
        # them signal type = 0 is to use the pk to pk voltage within the lead
        # gates
        if self.data.follow_gate_on and self.data.signal_type != 0:
            followL_idx = self.data.peak_bin[3, 1, i, j]
            followR_idx = self.data.peak_bin[4, 1, i, j]
            self.time_axis.axvline(self.data.time[followL_idx], color='b',
                                   linewidth=1.0, picker=2)
            self.time_axis.axvline(self.data.time[followR_idx], color='g',
                                   linewidth=1.0, picker=2)

        # if follow gate is on and using peak to peak voltage within follow
        # gates plot the peak locations
        if self.data.follow_gate_on and self.data.signal_type == 1:
            pos_peak = self.data.peak_bin[0, 1, i, j]
            neg_peak = self.data.peak_bin[1, 1, i, j]
            self.time_axis.plot(self.data.time[pos_peak],
                                self.data.waveform[i, j, pos_peak], 'b+')
            self.time_axis.plot(self.data.time[neg_peak],
                                self.data.waveform[i, j, neg_peak], 'gx')

            # now plot the peaks in the front gate
            pos_peak = self.data.peak_bin[0, 0, i, j]
            neg_peak = self.data.peak_bin[1, 0, i, j]
            self.time_axis.plot(self.data.time[pos_peak],
                                self.data.waveform[i, j, pos_peak], 'k*')
            self.time_axis.plot(self.data.time[neg_peak],
                                self.data.waveform[i, j, neg_peak], 'k*')

    def freq_plot(self, i, j):
        """
        Plots the frequency domain information between the active time domain
        gates
        :param i: The row
        :param j: The column
        """

        # if follow gate is on, we want to use that dimension of peak bin
        if self.data.follow_gate_on:
            idx = 1
        else:
            idx = 0

        if self.data.signal_type > 0:
            left_gate = self.data.peak_bin[3, idx, i, j]
            right_gate = self.data.peak_bin[4, idx, i, j]
        else:
            left_gate = self.data.gate[0][0]
            right_gate = self.data.gate[0][1]

        # create time domain waveform as zeros so it has the correct number
        # of points
        time_waveform = np.zeros(self.data.wave_length)
        time_waveform[left_gate:right_gate] = \
            self.data.waveform[i, j, left_gate:right_gate]

        freq_waveform = np.fft.rfft(time_waveform) * self.data.dt

        max_freq = 3.5
        max_f_idx = np.argmin(np.abs(self.data.freq - max_freq))

        freq = self.data.freq[:max_f_idx]
        freq_waveform = freq_waveform[:max_f_idx]

        self.freq_axis.cla()
        # self.freq_axis.semilogy(freq, np.abs(freq_waveform), 'r')
        self.freq_axis.plot(freq, np.abs(freq_waveform), 'r')
        self.freq_axis.set_xlabel('Frequency (THz)')
        self.freq_axis.set_ylabel('Amplitude')
        self.freq_axis.set_xlim(0, 3.5)
        self.freq_axis.grid()

        if self.data.in_freq_mode:
            # self.gate0 is the left frequency domain gate
            gate0_x = self.data.freq[self.freq_gate[0]]
            self.gate0 = self.freq_axis.axvline(gate0_x, color='b', picker=2)
            # gate1 is the right frequency domain gate
            gate1_x = self.data.freq[self.freq_gate[1]]
            self.gate1 = self.freq_axis.axvline(gate1_x, color='g', picker=2)

    def switch_a_scan_view(self, event):
        """
        Switches the A-Scan view between being A-Scan only or A-Scan and
        frequency domain information
        """

        # switch boolean value
        self.a_scan_only = not self.a_scan_only

        # clear figure of all information, including axes
        self.figure.clf()

        if self.a_scan_only:
            self.time_axis = self.figure.add_subplot(111)
            self.freq_axis = None
        else:
            self.time_axis = self.figure.add_subplot(211)
            self.freq_axis = self.figure.add_subplot(212)

        if self.is_initialized:
            # if a point has been clicked on in the C-Scan show the previous
            # waveform
            self.plot(self.i_index, self.j_index)
        else:
            # if not, update so that new axes view is visible
            self.figure_canvas.draw()

    def change_index_system(self, event):
        """
        Controls whether the title in the A-Scan frame displays (i, j) or (x, y)
        location.
        :param event: a wxPython event
        """
        # only thing we have to do is change the state of the value and then
        # replot
        self.ij_indexing = not self.ij_indexing

        if self.is_initialized:
            self.plot(self.i_index, self.j_index)

    def on_change_gate_button(self, event):
        """
        Opens a dialog box that allows the user to type in a new gate and set
        the gate at a specific location
        """

        ChangeGateFrame(self, self.holder, self.data)

    def on_toggle_title_button(self, event):
        """
        Changes whether the title is visible or not. This can be used in case
        we save a figure and do not want to to have the title displayed.
        """
        self.show_title = not self.show_title
        self.plot()

    def motion_handler(self, event):
        """
        Prints the current (x,y) values that the mouse is over to the status bar
        along with the amplitude of the waveform at that time.
        :param event: a matplotlib event
        """

        # if no C-Scan point has been clicked on yet, do nothing
        if not self.is_initialized:
            return

        xid = event.xdata
        yid = event.ydata

        # if user clicks in the window but outside of the plot itself, xid &
        # will be None
        if xid is not None and yid is not None:
            if xid < 0:
                t_index = 0
            elif xid > self.data.time_length:
                t_index = self.data.wave_length - 1
            else:
                t_index = int(round(xid / self.data.dt, 0))

            point_amp = self.data.waveform[self.i_index, self.j_index, t_index]

            status_string = '(%.4f, %.4f), [%.4f]' % (xid, yid, point_amp)
            self.status_bar.SetStatusText(status_string)
        else:
            self.status_bar.SetStatusText('')

    def grab_gate_handler(self, event):
        """
        Determines that gate that is clicked on
        """
        # if no point on the C-Scan has been clicked on yet, do nothing and
        # return
        if not self.is_initialized:
            return

        # if the user does not click on one of the gate lines
        if not isinstance(event.artist, Line2D):
            return

        # pass the event through so the handler that is called has access to
        # the event location and stuff
        if event.artist.axes == self.time_axis:
            self.time_gate_handler(event)
        else:  # on the frequency axis
            self.freq_gate_handler(event)

    def time_gate_handler(self, event):
        line = event.artist
        self.line_held = event.artist
        xdata = line.get_xdata()[0]

        # convert xdata of point clicked on to  nearest time index
        index = int(round(xdata / self.data.dt, 0))

        # determine which gate was clicked on
        diff = np.abs(index - self.data.gate[0][0])  # check front lead gate
        diff1 = np.abs(index - self.data.gate[0][1])  # check back lead gate

        if diff1 > diff:
            self.gate_held = 0  # clicked on the lead lead gate
        else:
            self.gate_held = 1  # clicked on the back front gate
            diff = diff1

        # if follow gate is on we need to check the follow gates also use
        # peak_bin to check with the follow gate location for this particular
        # (i, j) pair
        if self.data.follow_gate_on:
            # difference between point clicked and lead follow gate
            diff2 = np.abs(index - self.data.peak_bin[3, 1, self.i_index,
                                                      self.j_index])

            # back follow gate
            diff3 = np.abs(index - self.data.peak_bin[4, 1, self.i_index,
                                                      self.j_index])

            # determine if either is closer than the front gate
            if diff2 < diff:  # point clicked is closest to lead follow gate
                self.gate_held = 2
                diff = diff2
            if diff3 < diff:  # point clicked is closest to back follow gate
                self.gate_held = 3

        line.set_linewidth(2.0)  # make the line clicked on bold
        self.figure_canvas.draw()  # update plot

    def freq_gate_handler(self, event):
        """
        """
        # store the line object that was clicked on in a variable so we can
        # move it around with gate_slider() function
        self.line_held = event.artist

        # make the line that was clicked on bold
        event.artist.set_linewidth(2.0)
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

        # if the user has not yet clicked on a line, do nothing
        if self.line_held is None:
            return

        # if the user does not release in the image, exit
        if event.xdata is None or event.ydata is None:
            return

        # if line_held is one of the frequency domain gates, that means we are
        # doing frequency domain plotting
        if self.line_held == self.gate0 or self.line_held == self.gate1:
            self.line_held.set_linewidth(1.0)  # return to normal line width

            f = self.line_held.get_xdata()[0]
            if self.line_held == self.gate0:
                self.freq_gate[0] = np.argmin(np.abs(self.data.freq - f))
            else:  # holding self.gate1
                self.freq_gate[1] = np.argmin(np.abs(self.data.freq - f))

            self.data.make_freq_c_scan(0, self.freq_gate[0], self.freq_gate[1])
            # update A-Scan, Raw C-Scan and Interpolated C-Scan plots with the
            # current information
            self.plot()
            self.holder.raw_c_scan_frame.update()
            self.holder.interpolated_c_scan_frame.update()

            self.line_held = None

            return

        # convert time value (xdata) to index
        index = int(round(event.xdata * (self.data.wave_length-1) /
                          self.data.time_length, 0))

        # ensure that the gate is inside of the bounds
        if index < 0:
            index = 0
        elif index > self.data.wave_length - 1:
            index = self.data.wave_length - 1

        new_gate = copy.deepcopy(self.data.gate)

        # if we are adjusting the lead gates, we just use the index as is
        if self.gate_held == 0:  # front lead gate
            new_gate[0][0] = index
        elif self.gate_held == 1:  # back lead gate
            new_gate[0][1] = index
        # to adjust the follow gate, we need to use new index in relation to
        # the old one found in peak_bin
        elif self.gate_held == 2:  # front follow gate
            new_gate[1][0] = self.data.gate[1][0] - \
                (self.data.peak_bin[3, 1, self.i_index, self.j_index] - index)
        else:  # gate_held = 3; back follow gate
            new_gate[1][1] = self.data.gate[1][1] - \
                (self.data.peak_bin[4, 1, self.i_index, self.j_index] - index)

        # if follow gate is on, use the change gate method, which updates
        # bin_range, and calls find_peaks function to update peak_bin
        # if a value error is encountered it usually means that the left
        # follow gate is greater than the last time length index
        if self.data.follow_gate_on:
            self.data.change_gate(new_gate)
            self.line_held.set_xdata([event.xdata, event.xdata])
            self.line_held.set_linewidth(1.0)  # return to normal width
            self.figure_canvas.draw()
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

    def gate_slider(self, event):
        """
        Moves and redraws the gate that has been grabbed
        """
        # if no point on the C-Scan has been clicked yet, do nothing and return
        if not self.is_initialized:
            return

        # if the user has not yet clicked on a line, do nothing
        if self.line_held is None:
            return

        # if the mouse leaves the image
        if event.xdata is None or event.ydata is None:
            return

        # slide the line around with the mouse
        self.line_held.set_xdata([event.xdata, event.xdata])

        # update figure
        self.figure_canvas.draw()


class ChangeGateFrame(wx.Frame):
    """
    Frame to hold the change gate panel
    """

    def __init__(self, parent, holder, data):
        """
        Constructor Method
        :param parent: The parent frame that called this frame
        :param holder: The holder class that is binding everything together
        :param data: An instance of the THzData class
        """
        # this frame exists to hold the ChangeGatePanel that's about it

        super(ChangeGateFrame, self).__init__(parent, title='Change Gates ...',
                                              size=(450, 200))

        self.holder = holder

        self.panel = _ChangeGatePanel(self, data)

        self.Show(True)


# Private class
class _ChangeGatePanel(wx.Panel):
    """
    Private class that contains type boxes for the user to manually type
    in values to change the gate to
    """

    def __init__(self, parent, data):
        """
        Constructor method
        :param parent: The parent frame that this panel is being called by
        :param data: An instance of the THzData class
        """
        super(_ChangeGatePanel, self).__init__(parent)

        # store the parent Frame that created this panel
        self.parent = parent

        # The THz Data
        self.data = data

        self.main_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.gbs = wx.GridBagSizer(vgap=10, hgap=10)

        self.add_controls()
        self.bind_controls()

        # set the main sizer to the panel
        self.SetSizer(self.main_sizer)

    def add_controls(self):
        """
        Add the Text Control boxes to the main frame
        """
        # text boxes for the front and back lead gates
        self.front1_text_control = wx.TextCtrl(self, -1)
        self.front2_text_control = wx.TextCtrl(self, -1)

        # text boxes for the front and back follow gates
        self.follow1_text_control = wx.TextCtrl(self, -1)
        self.follow2_text_control = wx.TextCtrl(self, -1)

        # OK and cancel buttons
        self.ok_button = wx.Button(self, label='OK')
        self.cancel_button = wx.Button(self, label='Cancel')

        # set the text for the front and follow gate text control boxes to
        # contain the position of the gates when the panel was created
        self.front1_text_control.SetValue(str(self.data.gate[0][0]))
        self.front2_text_control.SetValue(str(self.data.gate[0][1]))
        self.follow1_text_control.SetValue(str(self.data.gate[1][0]))
        self.follow2_text_control.SetValue(str(self.data.gate[1][1]))

        self.gbs.Add(wx.StaticText(self, -1, 'Front Gate 1'), (0, 0))
        self.gbs.Add(self.front1_text_control, (0, 1))

        self.gbs.Add(wx.StaticText(self, -1, 'Front Gate 2'), (0, 2))
        self.gbs.Add(self.front2_text_control, (0, 3))

        self.gbs.Add(wx.StaticText(self, -1, 'Follow Gate 1'), (1, 0))
        self.gbs.Add(self.follow1_text_control, (1, 1))

        self.gbs.Add(wx.StaticText(self, -1, 'Follow Gate 2'), (1, 2))
        self.gbs.Add(self.follow2_text_control, (1, 3))

        self.gbs.Add(self.ok_button, (2, 0))
        self.gbs.Add(self.cancel_button, (2, 1))

        self.main_sizer.Add(self.gbs, 1, wx.ALL | wx.EXPAND, 10)

    def bind_controls(self):
        """
        Binds the controls to their appropriate method handler
        """
        self.ok_button.Bind(wx.EVT_BUTTON, self.on_ok_button)
        self.cancel_button.Bind(wx.EVT_BUTTON, self.on_cancel_button)

    def on_cancel_button(self, event):
        """
        Closes the dialog and does nothing but exit
        :param event: a wxPython event
        """
        # call the parent destroy method to close the window
        self.parent.Destroy()

    def on_ok_button(self, event):
        """
        Attempt to get the new gate values and call THzData's change_gate method
        to change the gate, then exit the Frame
        """

        # pre-allocate gate as a list
        gate = [[None, None], [None, None]]

        gate[0][0] = int(self.front1_text_control.GetValue())
        gate[0][1] = int(self.front2_text_control.GetValue())

        gate[1][0] = int(self.follow1_text_control.GetValue())
        gate[1][1] = int(self.follow2_text_control.GetValue())

        # change the gate
        self.data.change_gate(gate)

        # create the new C-Scan
        self.data.make_c_scan()

        # update the A-Scan, Raw C-Scan, and Interpolated C-Scan frames
        if self.parent.holder.a_scan_frame.is_initialized:
            self.parent.holder.a_scan_frame.plot()
        self.parent.holder.raw_c_scan_frame.update()
        self.parent.holder.interpolated_c_scan_frame.update()

        # call the parent frame's destroy method to close the window
        self.parent.Destroy()
