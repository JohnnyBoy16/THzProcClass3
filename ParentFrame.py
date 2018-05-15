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

    def __init__(self, title, subplot_grid=(1, 1)):
        """
        Constructor method.
        :param title: The title for the Frame. This is the name that goes in
            the very top Window's navigation bar.
        :param subplot_grid: The axes grid layout the figure is to have.
            Expected as (nrows, ncols)
        """
        super().__init__(None, -1, title)

        # holds the figure object
        self.figure = None

        # holds the axis object(s). If initialize_figure() is called with either
        # nrows or ncols greater than 1. This will be an array with each axis
        # instance in order
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

        # the menu bar the contains the file menu. Make this an instance
        # attribute so it can be modified by child classes
        self.menu_bar = None

        # menu button to open a reference txt file and look at it
        self.open_ref_menu_button = None

        # the exit menu for the frames
        self.exit_menu = None

        # Start calling the methods
        self.create_menu()
        self.connect_menu()
        self.initialize_figure(subplot_grid[0], subplot_grid[1])
        self.initialize_toolbar()
        self.initialize_sizer()

        # show the frame so the use doesn't have to call this manually
        # in their driving script
        self.Show(True)

    def initialize_figure(self, nrows=1, ncols=1):
        """
        Create a figure, adds an axis or axes to that figure and initializes a figure canvas. If
        either nrows or ncols is greater than 1, axis attribute will be an array containing axis
        instances.
        :param nrows: The number of rows for the subplot. Default: 1
        :param ncols: The number of columns for the subplot. Default: 1
        """

        self.figure, self.axis = plt.subplots(nrows, ncols)

        self.figure_canvas = FigureCanvas(self, -1, self.figure)

    def initialize_toolbar(self):
        """
        Creates a matplotlib toolbar and adds it to the status bar
        """
        # status bar which shows information about current mouse location
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

    def create_menu(self):
        """
        Creates a menu that includes a exit option to terminate the program
        """
        file_menu = wx.Menu()

        self.open_ref_menu_button = wx.MenuItem(file_menu, wx.ID_ANY, 'Open Ref',
                                                'Open a Reference File')
        file_menu.Append(self.open_ref_menu_button)

        self.exit_menu = wx.MenuItem(file_menu, wx.ID_EXIT, 'E&xit', 'Terminate the Program')
        file_menu.Append(self.exit_menu)

        self.menu_bar = wx.MenuBar()
        self.menu_bar.Append(file_menu, '&File')

        self.SetMenuBar(self.menu_bar)

    def connect_menu(self):
        """
        Binds the menu options with their method handlers
        """
        self.Bind(wx.EVT_MENU, self.on_exit, self.exit_menu)
        self.Bind(wx.EVT_MENU, self.on_open_ref, self.open_ref_menu_button)

    def on_open_ref(self, event):
        """
        Opens a new frame that shows the time and frequency information of a
        reference txt file
        """
        # need to have this import statement here, otherwise there will be
        # circular imports between RefFrame and ParentFrame
        from THzProc.RefFrame import ReferenceFrame

        dlg = wx.FileDialog(self, 'Open Reference', 
                            wildcard='txt files (*.txt)|*.txt',
                            style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST)

        # set the default directory to the reference directory
        # this is the directory that the dialog will open in
        dlg.SetDirectory('C:\\Work\\Refs')

        if dlg.ShowModal() == wx.ID_CANCEL:
            return

        full_path = dlg.GetPath()

        ReferenceFrame(full_path)

    @staticmethod
    def on_exit(event):
        """
        Terminates the program when the user clicks the exit option from the 
        file menu
        """
        for window in wx.GetTopLevelWindows():
            window.Destroy()
