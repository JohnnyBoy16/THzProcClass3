import os

import wx
import matplotlib.pyplot as plt


class StartFrame(wx.Frame):
    """
    Frame to hold information that starts the THzProc Data Class
    """
    def __init__(self, title):
        super().__init__(None, -1, title, size=(-1, -1))  # size=(-1, -1) gives default size

        # Variable declarations
        # define the menus that are used
        self.clear_menu = None
        self.open_menu = None
        self.save_menu = None
        self.load_menu = None
        self.exit_menu = None
        self.about_menu = None

        self.create_menus()
        self.connect_menus()

        self.Show(True)  # always show the frame

    def create_menus(self):
        """
        Create menu system at the top of the window
        """
        # create a status bar in the bottom of the window
        self.CreateStatusBar()

        file_menu = wx.Menu()

        self.clear_menu = wx.MenuItem(file_menu, wx.ID_NEW, '&Clear', 'Clear THz Plots')
        file_menu.Append(self.clear_menu)

        self.open_menu = wx.MenuItem(file_menu, wx.ID_OPEN, '&Open', 'Open a .tvl file')
        file_menu.Append(self.open_menu)

        self.save_menu = wx.MenuItem(file_menu, wx.ID_SAVE, '&Save', 'Save the state of the GUI')
        file_menu.Append(self.save_menu)

        self.load_menu = wx.MenuItem(file_menu, wx.ID_ANY, '&Load',
                                     'Load a previously saved GUI state')
        file_menu.Append(self.load_menu)

        self.exit_menu = wx.MenuItem(file_menu, wx.ID_EXIT, 'E&xit', 'Terminate the Program')
        file_menu.Append(self.exit_menu)

        info_menu = wx.Menu()

        self.about_menu = wx.MenuItem(info_menu, wx.ID_ABOUT, '&About', 'Information')
        info_menu.Append(self.about_menu)

        # create  the menu bar and add the menus
        menu_bar = wx.MenuBar()
        menu_bar.Append(file_menu, '&File')
        menu_bar.Append(info_menu, '&Info')

        # add the menu bar to the frame
        self.SetMenuBar(menu_bar)

    def connect_menus(self):
        """
        Binds the menus with their appropriate handler
        """
        self.Bind(wx.EVT_MENU, self.on_open, self.open_menu)
        self.Bind(wx.EVT_MENU, self.exit, self.exit_menu)
        self.Bind(wx.EVT_MENU, self.about, self.about_menu)
        self.Bind(wx.EVT_MENU, self.clear, self.clear_menu)
        self.Bind(wx.EVT_MENU, self.save, self.save_menu)
        self.Bind(wx.EVT_MENU, self.load, self.load_menu)

    def on_open(self, event):
        """
        Opens a dialog box so the user can select a .tvl file to load
        """

        print('You clicked the open button!')

        # give dialog to find the file
        dlg = wx.FileDialog(self, 'Open file...', os.getcwd(), style=wx.FD_OPEN,
                            wildcard='TVL Files (*.tvl)|*.tvl')

        # if a user selects okay in the dialog
        if dlg.ShowModal() == wx.ID_OK:
            # save the pathname from the dialog
            full_path = dlg.GetPath()
            print(full_path)

        dlg.Destroy()

    def save(self, event):
        """
        Saves the settings in the GUI to a csv file
        """
        print('You clicked the save button!')
        # dlg = wx.FileDialog(self, 'Save GUI State...', os.getcwd(), style=wx.SAVE | wx.OVERWRITE_PROMPT,
        #                     wildcard='csv Files (*.csv)|*.csv')
        #
        # # if the user selects okay withing the dialog
        # if dlg.ShowModal() == wx.ID_OK:
        #     full_path = dlg.GetPath()
        #     self.store_gui_values(full_path)

    def load(self, event):
        """
        Allows the user to load settings from a previously saved GUI state
        """
        print('You clicked the load button!')

        # dlg = wx.FileDialog(self, 'Load GUI State...', os.getcwd(), style=wx.OPEN, wildcard='csv Files (*.csv)|*.csv')
        # if dlg.ShowModal() == wx.ID_OK:
        #     full_path = dlg.GetPath()
        #     self.panel.read_stored_gui_data(full_path)

    def about(self, event):
        """
        Provides information about the GUI
        """
        print('You clicked the about button!')

    def clear(self, event):
        """
        Clears the figures that have been generated by THzProc. Allows the program to be run again.
        """
        print('You clicked the clear button!')

        # close the pyplot windows (basically terminate THzProc)
        # plt.close('all')

    def exit(self, event):
        """
        Stores the values in the GUI to a csv file and exits the GUI
        """
        print('You clicked the exit button!')

        self.Destroy()


class StartPanel(wx.Panel):
    """
    Panel to control information that is needed to run THzProc Data Class
    """
    def __init__(self):
        super().__init__()  # call parent constructor

        # Initialize filename, basedir, and gate
        self.filename = ''
        self.basedir = ''
        self.gate = None

        # create the main sizer
        main_sizer = wx.BoxSizer(wx.VERTICAL)

        # add the sizer to the panel
        self.SetSizer(main_sizer)

    def add_controls(self):
        """
        Add controls to the panel
        """


if __name__ == '__main__':
    from FrameHolder import FrameHolder
    from THzData import THzData

    # create a new app and don't redirect stdout/stderr to a window
    # not sure what this does exactly
    app = wx.App(False)

    # instantiate the class
    frame = StartFrame(title='THzProc Start GUI')

    # holder = FrameHolder(data)

    # start the event loop
    app.MainLoop()
