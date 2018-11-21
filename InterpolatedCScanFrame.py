import pdb

import matplotlib.pyplot as plt
import wx

from THzProc.ParentFrame import ParentFrame


class InterpolatedCScanFrame(ParentFrame):
    """
    Frame for the Interpolated C-Scan.
    """
    # this plot is not interactive. Clicking in the image will not change the
    # A-Scan. However, the image will change when the gates are updated.

    def __init__(self, holder, data, title=None):
        """
        Constructor method
        :param holder: an instance of FrameHolder, used to link actions in
                one frame to actions in another frame
        :param data: an instance of the THzData class
        :param title: the title for the frame
        """

        if title is None:
            title = 'Interpolated C-Scan Test'

        # call inherited ParentFrame constructor using the arguments so it runs
        # in python 2
        super(InterpolatedCScanFrame, self).__init__(title)

        # the instance of holder class for this frame
        self.holder = holder

        # the THz data
        self.data = data

        # the C-Scan image
        self.image = None

        # the options menu for the this frame
        self.options_menu = None

        # menu options to change the orientation of the colorbar
        self.change_colorbar_dir_button = None

        # the orientation of the colorbar for this specific frame
        self.colorbar_dir = 'horizontal'

        # the colorbar for the figure
        self.colorbar = None

        self.add_options_menu()
        self.connect_events()

        # make sure to call plot so the interpolated image is plotted when
        # starting
        self.plot()

        # close the figure, if this is not done, a normal matplotlib figure for
        # this image will pop up if plt.show() is called after creating the GUI
        plt.close(self.figure)

    def add_options_menu(self):
        """
        Adds an options menu to the frame with on options to change to the
        orientation of the colorbar
        """
        self.options_menu = wx.Menu()

        title = 'Change Colorbar Orientation'
        description = 'Changes the colorbar orientation between vertical and horizontal'
        self.change_colorbar_dir_button = wx.MenuItem(self.options_menu, wx.ID_ANY,
                                                      title, description)

        self.options_menu.Append(self.change_colorbar_dir_button)

        self.menu_bar.Append(self.options_menu, 'Options')
        self.SetMenuBar(self.menu_bar)
        self.Fit()

    def connect_events(self):
        """
        Connects events to their appropriate method
        """
        # matplotlib events
        self.figure_canvas.mpl_connect('motion_notify_event', self.motion_handler)

        # wx events
        self.Bind(wx.EVT_MENU, self.on_change_colorbar_dir, self.change_colorbar_dir_button)

    def plot(self):
        """
        Plot the interpolated C-Scan for the current gate location
        """
        self.image = self.axis.imshow(self.data.c_scan, interpolation='bilinear', cmap='jet',
                                      extent=self.data.c_scan_extent)
        self.axis.set_xlabel('X Scan Location (mm)')
        self.axis.set_ylabel('Y Scan Location (mm)')
        self.colorbar = plt.colorbar(self.image, orientation=self.colorbar_dir)
        self.axis.grid()
        self.figure_canvas.draw()

    def update(self):
        """
        Updates the figure, by clearing the axis and remaking stuff, including
        the colorbar
        """
        # At the moment it is necessary to have a separate method from plot(),
        # because it appears that there must be a call to
        # colorbar.update_bruteforce() to have the colorbar actually update.
        # Using self.axis.cla() does not actually remove the colorbar
        self.axis.cla()
        self.image = self.axis.imshow(self.data.c_scan, interpolation='bilinear', cmap='jet',
                                      extent=self.data.c_scan_extent)
        self.axis.set_xlabel('X Scan Location (mm)')
        self.axis.set_ylabel('Y Scan Location (mm)')
        self.colorbar.update_bruteforce(self.image)
        self.axis.grid()
        self.figure_canvas.draw()

    def motion_handler(self, event):
        """
        Prints the current (x,y) values that the mouse is over to the status
        bar along with the pixel value at that (x,y) location.
        """
        xid = event.xdata
        yid = event.ydata
        if xid is not None and yid is not None:
            x_index = int((xid - self.data.x_min) / self.data.dx)
            y_index = int((yid - self.data.y_min) / self.data.dy)
            point_amp = self.data.c_scan[y_index, x_index]
            status_string = '(%.6f, %.6f), [%.6f]' % (xid, yid, point_amp)
            self.status_bar.SetStatusText(status_string)
        else:
            self.status_bar.SetStatusText('')

    def on_change_colorbar_dir(self, event):
        """
        Changes the orientation of the colorbar
        """

        if self.colorbar_dir == 'horizontal':
            self.colorbar_dir = 'vertical'
        else:
            self.colorbar_dir = 'horizontal'

        self.colorbar.remove()

        # need to specify which axis the colorbar will be drawn on. Otherwise
        # it will attach itself to the most recent figure object, which may not
        # be the figure from this frame if other plots have been created
        self.colorbar = plt.colorbar(self.image, ax=self.axis,
                                     orientation=self.colorbar_dir)

        self.figure_canvas.draw()
