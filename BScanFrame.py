import pdb

import wx
import numpy as np
import matplotlib.pyplot as plt

from THzProc.ParentFrame import ParentFrame


class BScanFrame(ParentFrame):
    """
    Frame to display the B-Scan at the last location clicked on in the gray
    scale C-Scan
    """

    def __init__(self, holder, data, title=None):
        """
        Constructor method
        :param holder: an instance of FrameHolder, used to link actions in
                one frame to actions in another frame
        :param data: an instance of the THzData class
        :param title: the title for the frame
        """

        if title is None:
            title = 'B-Scan Frame'

        # call the inherited ParentFrame constructor using the arguments so
        # that the code works in python 2
        super(BScanFrame, self).__init__(title)

        # the class instance holding stuff together
        self.holder = holder

        # the THz data
        self.data = data

        # whether or not a point in the C-Scan has been clicked on yet
        self.is_initialized = False

        # contains the B-Scan picture
        self.image = None

        # the last (i, j) coordinate that was clicked on
        self.i_index = None
        self.j_index = None

        # menu item to toggle the B-Scan orientation between vertical and
        # horizontal
        self.toggle_orientation_menu = None

        self.modify_menu()
        self.connect_events()

        plt.close(self.figure)

    def modify_menu(self):
        """
        Modifies the menu to add an options menu that includes an item to switch
        the B-Scan orientation
        """
        options_menu = wx.Menu()

        status_string = 'Toggles orientation between horizontal and vertical'
        self.toggle_orientation_menu = wx.MenuItem(options_menu, wx.ID_ANY,
                                                   'Switch Orientation', status_string)

        options_menu.Append(self.toggle_orientation_menu)

        self.menu_bar.Append(options_menu, '&Options')
        self.SetMenuBar(self.menu_bar)
        self.Fit()

    def connect_events(self):
        """
        Connects events to their appropriate method
        """
        # matplotlib events
        self.figure_canvas.mpl_connect('motion_notify_event', self.motion_handler)

        # wx events
        self.Bind(wx.EVT_MENU, self.switch_orientation, self.toggle_orientation_menu)

    def plot(self, i, j):
        """
        Plot the B-Scan cross section
        :param i: The column index to look at
        :param j: The row index to look at
        """
        # show B-Scan for given (i, j) index, the THzData class handles the
        # direction (either horizontal or vertical)

        if not self.is_initialized:
            self.is_initialized = True

        # store the last points that were clicked on in the Raw C-Scan Frame
        self.i_index = i
        self.j_index = j

        if self.data.b_scan_dir == 'vertical':
            title_string = 'Line at x = %0.2f' % self.data.x[j]
            xlabel_string = 'Y Scan Location (mm)'
        else:  # b_scan_dir is horizontal
            title_string = 'Line at y = %0.2f' % self.data.y[i]
            xlabel_string = 'X Scan Location (mm)'

        self.axis.cla()
        self.data.make_b_scan(i, j)  # make the B-Scan image
        self.image = self.axis.imshow(self.data.b_scan, interpolation='none',
                                      cmap='seismic', aspect='auto',
                                      extent=self.data.b_scan_extent)
        self.axis.set_xlabel(xlabel_string)
        self.axis.set_ylabel('Time (ps)')
        self.axis.set_title(title_string)
        self.axis.grid()
        self.figure_canvas.draw()

    def switch_orientation(self, event):
        """
        Switches the B-Scan orientation
        """
        if self.data.b_scan_dir == 'horizontal':
            self.data.b_scan_dir = 'vertical'
        else:  # data.b_scan_dir should be vertical, so set it to horizontal
            self.data.b_scan_dir = 'horizontal'

        # if a point in the C-Scan has not been clicked on yet, we don't want
        # to do anything other than change the direction string above
        if not self.is_initialized:
            return

        # make the new B-Scan in THzProc at the location specified
        self.data.make_b_scan(self.i_index, self.j_index)

        # plot the new B-Scan
        self.plot(self.i_index, self.j_index)

        self.holder.raw_c_scan_frame.flash_b_scan_line()

    def motion_handler(self, event):
        """
        Prints the current (x,y) values that the mouse is over to the status
        bar along with the pixel value at that (x,y) location.
        """
        # Since we are looking at the B-Scan, the y-axis now represents time
        # and the x-axis on the plot represents either the x or y axis on the
        # sample depending on the B-Scan orientation

        # if now point in the C-Scan has been clicked on yet there is no
        # information to display on the status bar
        if not self.is_initialized:
            return

        # if the event is not in the axis
        if not event.inaxes:
            self.status_bar.SetStatusText('')
            return

        xid = event.xdata
        yid = event.ydata

        if self.data.b_scan_dir == 'horizontal':
            x_index = np.argmin(np.abs(self.data.x - xid))
        elif self.data.b_scan_dir == 'vertical':
            x_index = np.argmin(np.abs(self.data.y - xid))
        else:
            raise ValueError('B-Scan direction must be either "horizontal" or '
                             '"vertical"!')

        # need to flip the time array to index correctly. In the B-Scan the
        # bottom of the y-axis corresponds to t[0], but in the B-Scan image
        # this is the bottom row and accessed [-1].
        t_index = np.argmin(np.abs(np.flipud(self.data.time) - yid))
        point_amp = self.data.b_scan[t_index, x_index]
        status_string = '(%.4f, %.4f), [%.4f]' % (xid, yid, point_amp)
        self.status_bar.SetStatusText(status_string)
