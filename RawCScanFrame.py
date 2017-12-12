import pdb

import matplotlib.pyplot as plt

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

        self.plot()  # make sure to plot the C-Scan to start out with
        self.connect_events()

    def plot(self):
        """
        Plots the Raw C-Scan in the object figure windows after clearing what was in the image
        previously
        """
        self.image = self.axis.imshow(self.data.c_scan, interpolation='none', cmap='gray',
                                      extent=self.data.c_scan_extent, picker=True)
        # self.axis.set_title('Raw C-Scan Test', fontsize=16)
        self.axis.set_xlabel('X Scan Location (mm)', fontsize=14)
        self.axis.set_ylabel('Y Scan Location (mm)', fontsize=14)
        self.colorbar = plt.colorbar(self.image, orientation='vertical')
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
            status_string = '(%.6f, %.6f), [%.6f]' % (xid, yid, point_amp)
            self.status_bar.SetStatusText(status_string)
        else:
            self.status_bar.SetStatusText('')

    def select_point(self, event):
        """
        Allows the user to click on the C-Scan and update the A-Scan
        """
        x_data = event.xdata
        y_data = event.ydata

        if x_data and y_data is not None:  # make sure that the user clicks inside of the plot
            self.i_index = int((y_data - self.data.y_min) / self.data.delta_y)
            self.j_index = int((x_data - self.data.x_min) / self.data.delta_x)

            self.holder.a_scan_frame.plot(self.i_index, self.j_index)
            self.holder.b_scan_frame.plot(self.i_index, self.j_index)

            print(self.data.waveform[self.i_index, self.j_index, :].max())
            print(self.data.time[self.data.waveform[self.i_index, self.j_index, :].argmax()])
            print()
