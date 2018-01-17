import pdb

import numpy as np

from ParentFrame import ParentFrame


class BScanFrame(ParentFrame):
    """
    Frame to display the B-Scan at the last location clicked on in the gray scale C-Scan
    """
    def __init__(self, holder, data, title=None):

        if title is None:
            title = 'B-Scan Frame'

        super().__init__(title)

        # the class instance holding stuff together
        self.holder = holder

        # the THz data
        self.data = data

        # contains the B-Scan picture
        self.image = None

        # the last (i, j) coordinate that was clicked on
        self.i_index = None
        self.j_index = None

    def connect_events(self):
        """
        Connects events to their appropriate method
        """
        self.figure_canvas.mpl_connect('motion_notify_event', self.motion_handler)

    def plot(self, i, j):
        """
        Plot the B-Scan cross section
        :param i: The column index to look at
        :param j: The row index to look at
        """
        # show B-Scan for given (i, j) index, the THzData class handles the direction
        # (either horizontal or vertical)

        if self.data.b_scan_dir == 'vertical':
            title_string = 'Line at x = %0.2f' % self.data.x[j]
        else:  # b_scan_dir is horizontal
            title_string = 'Line at y = %0.2f' % self.data.y[i]

        self.axis.cla()
        self.data.make_b_scan(i, j)  # make the B-Scan image
        self.image = self.axis.imshow(self.data.b_scan, interpolation='none', cmap='seismic',
                                      extent=self.data.b_scan_extent)
        self.axis.set_xlabel('X Scan Location')
        self.axis.set_ylabel('Time (ps)')
        self.axis.set_title(title_string)
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
