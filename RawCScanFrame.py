import pdb

import matplotlib.pyplot as plt
import wx
import numpy as np

from ParentFrame import ParentFrame


class RawCScanFrame(ParentFrame):
    """
    Frame for the gray-scale interactive C-Scan.
    """

    def __init__(self, holder, data, title=None):

        if title is None:
            title = 'Raw C-Scan Frame'

        super().__init__(title)  # call parent class

        self.holder = holder

        # the THzData
        self.data = data

        # the C-Scan image
        self.image = None

        # stores the colorbar
        self.colorbar = None

        # the last (i,j) [y,x] indices that were clicked
        self.i_index = None
        self.j_index = None

        self.rescale_button = wx.Button(self, -1, label='Rescale Colorbar')

        self.modify_sizer()
        self.bind_controls()
        self.plot()  # make sure to plot the C-Scan to start out with
        self.connect_events()

    def modify_sizer(self):
        """
        Adds a button to rescale the colorbar with the maximum and minimum values in the C-Scan
        image after it has been zoomed in or out using the magnify tool from the matplotlib
        toolbar.
        """
        # TODO: in the future this should be moved to the start GUI if I get around to making that

        # create a new sizer to display the button
        button_sizer = wx.BoxSizer(wx.HORIZONTAL)

        button_sizer.Add(self.rescale_button)

        self.sizer.Add(button_sizer)
        self.SetSizer(self.sizer)
        self.Fit()

    def plot(self):
        """
        Plots the Raw C-Scan in the object figure windows after clearing what was in the image
        previously
        """
        self.image = self.axis.imshow(self.data.c_scan, interpolation='none', cmap='gray',
                                      extent=self.data.c_scan_extent, picker=True, origin='upper')
        # self.axis.set_title('Raw C-Scan Test', fontsize=16)
        self.axis.set_xlabel('X Scan Location (mm)')
        self.axis.set_ylabel('Y Scan Location (mm)')
        self.colorbar = plt.colorbar(self.image, orientation=self.data.colorbar_dir)
        self.axis.grid()
        self.figure_canvas.draw()

    def update(self):
        """
        Updates the figure, by clearing the axis and remaking stuff, including the colorbar
        """
        # At the moment it is necessary to have a separate method from plot(), because it appears
        # that there must be a call to colorbar.update_bruteforce() to have the colorbar actually
        # update. Using self.axis.cla() does not actually remove the colorbar

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
        Binds the events to methods
        """
        self.figure_canvas.mpl_connect('motion_notify_event', self.motion_handler)
        self.figure_canvas.mpl_connect('button_press_event', self.select_point)

    def bind_controls(self):
        """
        Binds controls with their method handlers
        """
        self.rescale_button.Bind(wx.EVT_BUTTON, self.on_rescale_click)

    def motion_handler(self, event):
        """
        Prints the current (x,y) values that the mouse is over to the status bar along with the
        pixel value at that (x,y) location.
        """
        xid = event.xdata
        yid = event.ydata
        if xid is not None and yid is not None:
            x_index = int((xid - self.data.x_min) / self.data.delta_x)
            y_index = int((yid - self.data.y_min) / self.data.delta_y)
            point_amp = self.data.c_scan[y_index, x_index]
            status_string = '(%.4f, %.4f), [%.4f]' % (xid, yid, point_amp)
            self.status_bar.SetStatusText(status_string)
        else:
            self.status_bar.SetStatusText('')

    def select_point(self, event):
        """
        Allows the user to click on the C-Scan and update the A-Scan
        """
        x_data = event.xdata
        y_data = event.ydata

        # apparent method to check if matplotlib toolbar is clicked
        # doesn't seem to work with the wx backend yet
        # print(self.figure_canvas.manager.toolmanager.active_toggle())

        if x_data and y_data is not None:  # make sure that the user clicks inside of the plot
            self.i_index = int((y_data - self.data.y_min) / self.data.delta_y)
            self.j_index = int((x_data - self.data.x_min) / self.data.delta_x)

            self.holder.a_scan_frame.plot(self.i_index, self.j_index)
            self.holder.b_scan_frame.plot(self.i_index, self.j_index)

    def on_rescale_click(self, event):
        """
        Rescales the colorbar with the maximum and minimum values that are currently in the visible
        plot bounday.
        """
        # get the current x and y boundaries of the image
        # convert to array so they can be edited if necessary
        x_bounds = np.array(self.axis.get_xlim())
        y_bounds = np.array(self.axis.get_ylim())

        # since the y_bounds are normally inverted to show the data how we see it in the lab
        # rearrange them in low -> high order
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

        area = self.data.c_scan[i0:i1+1, j0:j1+1]  # want to include the bounds in calculation

        # resent vmin and vmax to be the min and max of area inside of plot bounds
        self.image.set_clim(vmin=area.min(), vmax=area.max())
        self.colorbar.update_normal(self.image)  # update colorbar with new vmin and vmax values
        self.figure_canvas.draw()  # redraw image
