import pdb

import matplotlib.pyplot as plt
import wx
import numpy as np

from THzProc.ParentFrame import ParentFrame


class RawCScanFrame(ParentFrame):
    """
    Frame for the gray-scale interactive C-Scan.
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
            title = 'Raw C-Scan Frame'

        super().__init__(title)  # call parent class

        self.holder = holder

        # the THzData
        self.data = data

        # the orientation of the colorbar for this specific frame
        self.colorbar_dir = 'vertical'

        # the C-Scan image
        self.image = None

        # stores the colorbar
        self.colorbar = None

        # the last (i,j) [y,x] indices that were clicked
        self.i_index = None
        self.j_index = None

        # the menu item that rescales the colorbar based on the part of the
        # C-Scan that is currently in view
        self.rescale_colorbar_menu = None

        # menu item that allows the user to change the colorbar orientation
        self.colorbar_dir_menu_button = None

        self.modify_menu()
        self.connect_events()
        self.plot()  # make sure to plot the C-Scan to start out with

        plt.close(self.figure)

    def modify_menu(self):
        """
        Adds an Options menu to the menu system. The options menu will include
        and options to rescale the colorbar
        """
        self.options_menu = wx.Menu()

        self.rescale_colorbar_menu = wx.MenuItem(self.options_menu, wx.ID_ANY, 'Rescale Colorbar',
                                                 'Rescale colorbar with current image')

        title = 'Change Colorbar Orientation'
        description = 'Changes the colorbar orientation between horizontal and vertical'
        self.colorbar_dir_menu_button = wx.MenuItem(self.options_menu, wx.ID_ANY, title,
                                                    description)

        self.options_menu.Append(self.rescale_colorbar_menu)
        self.options_menu.Append(self.colorbar_dir_menu_button)

        self.menu_bar.Append(self.options_menu, '&Options')
        self.SetMenuBar(self.menu_bar)
        self.Fit()

    def plot(self):
        """
        Plots the Raw C-Scan initially.
        """
        self.image = self.axis.imshow(self.data.c_scan, interpolation='none', cmap='gray',
                                      extent=self.data.c_scan_extent, picker=True, origin='upper')
        self.axis.set_xlabel('X Scan Location (mm)')
        self.axis.set_ylabel('Y Scan Location (mm)')
        self.colorbar = plt.colorbar(self.image, ax=self.axis, orientation=self.colorbar_dir)
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
        # Using self.axis.cla() does not actually remove the colorbar.

        # Other methods of updating the plot involve changing the data of
        # the image object with (image.set_data(new_data)), and then manually
        # updating the colorbar range with colorbar.set_clim(vmin, vmax), but
        # this seems to result in the colorbar losing its connection to the
        # the axis on which it is attached for some reason, the axis that is
        # associated with the colorbar seems to become None type. Then when
        # change_colorbar_dir is called, it will crash.

        self.axis.cla()
        self.image = self.axis.imshow(self.data.c_scan, interpolation='none', cmap='gray',
                                      extent=self.data.c_scan_extent, picker=True)
        self.axis.set_xlabel('X Scan Location (mm)')
        self.axis.set_ylabel('Y Scan Location (mm)')
        self.colorbar.update_bruteforce(self.image)
        self.axis.grid()
        self.figure_canvas.draw()

    def connect_events(self):
        """
        Binds the matplotlib and wx events to their method handlers
        """
        # matplotlib events
        self.figure_canvas.mpl_connect('motion_notify_event', self.motion_handler)
        self.figure_canvas.mpl_connect('button_press_event', self.select_point)

        # wx events
        self.Bind(wx.EVT_MENU, self.on_rescale_click, self.rescale_colorbar_menu)
        self.Bind(wx.EVT_MENU, self.change_colorbar_dir, self.colorbar_dir_menu_button)

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
            status_string = '(%.4f, %.4f), [%.4f]' % (xid, yid, point_amp)
            self.status_bar.SetStatusText(status_string)
        else:
            self.status_bar.SetStatusText('')

    def select_point(self, event):
        """
        Allows the user to click on the C-Scan and update the A-Scan
        """
        # if the user does not click within the axis of the plot just return
        # otherwise lots of error will be thrown
        if not event.inaxes:
            return

        # if the user has either the zoom or pan button selected just return
        # because we don't want to change A-Scan when the user is zooming or
        # panning around
        if self.toolbar._active == 'ZOOM' or self.toolbar._active == 'PAN':
            return

        x_data = event.xdata
        y_data = event.ydata

        self.i_index = int((y_data - self.data.y_min) / self.data.dy)
        self.j_index = int((x_data - self.data.x_min) / self.data.dx)

        self.holder.a_scan_frame.plot(self.i_index, self.j_index)
        self.holder.b_scan_frame.plot(self.i_index, self.j_index)

    def flash_b_scan_line(self):
        """
        Flashes a horizontal or vertical line (depending or the orientation
        that is currently in use) to correspond with the location of the B-Scan
        """
        # using a hex color for the color because using 'g' was not bright
        # enough in my opinion. Hex color '#81ff6b' corresponds to bright
        # green

        n_flashes = 3

        if self.data.b_scan_dir == 'horizontal':
            for i in range(n_flashes):
                line = self.axis.axhline(self.data.y[self.i_index],
                                         color='#81ff6b', linestyle='--')
                self.figure_canvas.draw()
                plt.pause(0.2)
                line.remove()
                self.figure_canvas.draw()
                plt.pause(0.2)

        if self.data.b_scan_dir == 'vertical':
            for i in range(n_flashes):
                line = self.axis.axvline(self.data.x[self.j_index],
                                         color='#81ff6b', linestyle='--')
                self.figure_canvas.draw()
                plt.pause(0.2)
                line.remove()
                self.figure_canvas.draw()
                plt.pause(0.2)

    def on_rescale_click(self, event):
        """
        Rescales the colorbar with the maximum and minimum values that are
        currently in the visible plot boundary.
        """
        # get the current x and y boundaries of the image
        # convert to array so they can be edited if necessary
        x_bounds = np.array(self.axis.get_xlim())
        y_bounds = np.array(self.axis.get_ylim())

        # since the y_bounds are normally inverted to show the data how we see
        # it in the lab rearrange them in low -> high order
        if y_bounds[0] > y_bounds[1]:
            temp = y_bounds[0]
            y_bounds[0] = y_bounds[1]
            y_bounds[1] = temp

        # it is possible that x_bounds may be inverted also
        if x_bounds[0] > x_bounds[1]:
            temp = x_bounds[0]
            x_bounds[0] = x_bounds[1]
            x_bounds[1] = temp

        j0 = np.argmin(np.abs(self.data.x - x_bounds[0]))
        j1 = np.argmin(np.abs(self.data.x - x_bounds[1]))

        i0 = np.argmin(np.abs(self.data.y - y_bounds[0]))
        i1 = np.argmin(np.abs(self.data.y - y_bounds[1]))

        # want to include the bounds in calculation, so add 1 to be inclusive
        area = self.data.c_scan[i0:i1+1, j0:j1+1]

        # resent vmin and vmax to be the min and max of area inside of plot
        # bounds
        self.image.set_clim(vmin=area.min(), vmax=area.max())

        # update colorbar with new vmin and vmax values
        self.colorbar.update_normal(self.image)

        self.figure_canvas.draw()  # redraw image

    def change_colorbar_dir(self, event):
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
        # be the figure from this frame
        self.colorbar = plt.colorbar(self.image, ax=self.axis,
                                     orientation=self.colorbar_dir)

        self.figure_canvas.draw()
