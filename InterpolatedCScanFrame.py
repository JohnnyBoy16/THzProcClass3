import pdb

import matplotlib.pyplot as plt

from THzProc.ParentFrame import ParentFrame


class InterpolatedCScanFrame(ParentFrame):
    """
    Frame for the Interpolated C-Scan.
    """
    # this plot is not interactive. Clicking in the image will not change the A-Scan. However,
    # the image will change when the gates are updated.

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

        super().__init__(title)  # call parent class

        # the instance of holder class for this frame
        self.holder = holder

        # the THz data
        self.data = data

        # the C-Scan image
        self.image = None

        # the colorbar for the figure
        self.colorbar = None

        self.connect_events()
        self.plot()  # make sure to call plot so the interpolated image is plotted when starting

    def connect_events(self):
        """
        Connects events to their appropriate method
        """
        self.figure_canvas.mpl_connect('motion_notify_event', self.motion_handler)

    def plot(self):
        """
        Plot the interpolated C-Scan for the current gate location
        """
        self.image = self.axis.imshow(self.data.c_scan, interpolation='bilinear', cmap='jet',
                                      extent=self.data.c_scan_extent)
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
        self.image = self.axis.imshow(self.data.c_scan, interpolation='bilinear', cmap='jet',
                                      extent=self.data.c_scan_extent)
        self.axis.set_xlabel('X Scan Location (mm)')
        self.axis.set_ylabel('Y Scan Location (mm)')
        self.colorbar.update_bruteforce(self.image)
        self.axis.grid()
        self.figure_canvas.draw()

    def motion_handler(self, event):
        """
        Prints the current (x,y) values that the mouse is over to the status bar along with the
        pixel value at that (x,y) location.
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
