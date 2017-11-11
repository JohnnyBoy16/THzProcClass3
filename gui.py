import os
import csv
import sys

import wx


class THzFrame(wx.Frame):
    def __init__(self, parent, title):
        # call the parent initializer
        wx.Frame.__init__(self, parent, title=title, size=(600, 400))

        # create the menu items and connect them to handlers
        self.create_menus()
        self.connect_menus()

        # add the panel to the frame
        self.panel = THzPanel(self)

        # show the GUI so user does not have to call this manually
        self.Show(True)

    def create_menus(self):
        """
        Creates a menu system
        """

        # a status bar in the bottom of the window
        self.CreateStatusBar()

        # create the file menu
        file_menu = wx.Menu()

        # create items that we want in the file menu
        self.clear_menu = wx.MenuItem(file_menu, wx.ID_NEW, '&Clear', 'Clear the THz Plots')
        file_menu.Append(self.clear_menu)

        self.open_menu = wx.MenuItem(file_menu, wx.ID_OPEN, '&Open', 'Open a .tvl File')
        file_menu.Append(self.open_menu)

        self.save_menu = wx.MenuItem(file_menu, wx.ID_SAVE, '&Save', 'Save the state of the GUI')
        file_menu.Append(self.save_menu)

        string = 'Load a previously save GUI state'
        self.load_menu = wx.MenuItem(file_menu, wx.ID_ANY, '&Load', string)
        file_menu.Append(self.load_menu)

        # create the menu bar add the menus
        menu_bar = wx.MenuBar()
        menu_bar.Append(file_menu, '&File')

        # Add the MenuBar to the Frame
        self.SetMenuBar(menu_bar)

    def connect_menus(self):
        """
        Binds the menu item with their appropriate handler
        """
        pass

    def on_open(self, event):
        """
        Opens the file handler
        """
        # create dialog that allows user to look for file
        dlg = wx.FileDialog(self, 'Open file...', os.getcwd(), style=wx.OPEN,
                            wildcard='TVL Files (*tvl)|*.tvl')

        # if a user selects ok in the dialog
        if dlg.ShowModl() == wx.ID_OK:
            # save the pathname from the dialog
            full_path = dlg.GetPath()

            # add code to change the path in the GUI screen

    def on_save(self, event):
        """
        Saves the setting in the GUI to a csv file
        """
        pass

    def on_load(self, event):
        """
        Loads the settings from a previously saved GUI state
        """
        pass

    def on_clear(self, event):
        """
        Clears the figures that have been generated. Allows the program to be run again.
        """
        pass

    def on_exit(self, event):
        """
        Stores the values in the GUI to a csv file and exits the GUI
        """
        pass


class THzPanel(wx.Panel):

    def __init__(self, parent):
        wx.Panel.__init__(self, parent=parent)

        # Initialize filename and basedir
        self.filename = ''
        self.basedir = ''
        self.gate = [[None, None], [None, None]]

        # create a sizer
        self.main_sizer = wx.BoxSizer(wx.VERTICAL)

        # self.add_controls()
        # self.set_initial_values()
        # self.b_scan_dir = None

        # # Set the sizer to the panel
        # self.SetSizer(self.main_sizer)

        # # Connect the control event to an event handler
        # self.bind_controls()


if __name__ == '__main__':
    app = wx.App(False)

    frame = THzFrame(None, 'Basic THz GUI')

    app.MainLoop()
