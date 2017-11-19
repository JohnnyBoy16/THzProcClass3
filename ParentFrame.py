import pdb

import matplotlib.pyplot as plt
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_wx import NavigationToolbar2Wx
import wx


class ParentFrame(wx.Frame):
    """
    Parent Frame for all other frames to inherit
    """
    # This organizes the code that all frames have in common into one place

    def __init__(self, title):
        super().__init__(None, -1, title)

        # holds the figure object
        self.figure = None

        # holds the axis object
        self.axis = None

        # holds the FigureCanvas Object
        self.figure_canvas = None

        # holds the status bar object
        self.status_bar = None

        # holds the matplotlib toolbar that allows the user to
        # save, zoom, ect.
        self.toolbar = None

        # the main sizer for the frame
        self.sizer = None

        self.initialize_figure()
        self.initialize_toolbar()
        self.initialize_sizer()

        self.Show(True)

    def initialize_figure(self):
        """
        Create a figure, adds an axis to that figure and initializes a figure canvas
        """
        self.figure = plt.figure()
        self.axis = self.figure.add_subplot(111)

        self.figure_canvas = FigureCanvas(self, -1, self.figure)

    def initialize_toolbar(self):
        """
        Creates a matplotlib toolbar and adds it to the status bar
        """
        # toolbar which allows users to save, zoom, pan, ect...
        self.status_bar = wx.StatusBar(self, -1)
        self.status_bar.SetFieldsCount(1)
        self.SetStatusBar(self.status_bar)

        # add toolbar to allow user to zoom, save, ect
        self.toolbar = NavigationToolbar2Wx(self.figure_canvas)

    def initialize_sizer(self):
        """
        Creates the frame sizer and adds the toolbar and figure canvas to itself
        """
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.toolbar, 0, wx.LEFT | wx.EXPAND)
        self.sizer.Add(self.figure_canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
        self.SetSizer(self.sizer)
        self.Fit()
        self.toolbar.Show()
