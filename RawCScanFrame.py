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

        # the menu item that calculates the signal to noise ratio from the part
        # of the C-Scan that is currently in view
        self.calculate_sn_ratio_menu_button = None

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
        options_menu = wx.Menu()

        self.rescale_colorbar_menu = wx.MenuItem(options_menu, wx.ID_ANY, 'Rescale Colorbar',
                                                 'Rescale colorbar with current image')

        description = 'Calculate C-Scan S/N ratio based on current view'
        title = 'Calculate C-Scan S/N Ratio'
        self.calculate_cscan_sn_ratio_button = wx.MenuItem(options_menu, wx.ID_ANY, title,
                                                           description)

        description = 'Calculate A-Scan S/N ratio based on current C-Scan View'
        title = 'Calculate A-Scan S/N Ratio'
        self.calculate_sn_ratio_menu_button = wx.MenuItem(options_menu, wx.ID_ANY, title,
                                                          description)

        title = 'Calculate A-Scan S/N Ratio from Defect only'
        description = 'Calculate S/N ratio from defect waveforms only'
        self.calculate_sn_ratio_menu_button2 = wx.MenuItem(options_menu, wx.ID_ANY, title,
                                                           description)

        title = 'Change Colorbar Orientation'
        description = 'Changes the colorbar orientation between horizontal and vertical'
        self.colorbar_dir_menu_button = wx.MenuItem(options_menu, wx.ID_ANY, title, description)

        options_menu.Append(self.rescale_colorbar_menu)
        options_menu.Append(self.calculate_cscan_sn_ratio_button)
        options_menu.Append(self.calculate_sn_ratio_menu_button)
        options_menu.Append(self.calculate_sn_ratio_menu_button2)
        options_menu.Append(self.colorbar_dir_menu_button)

        self.menu_bar.Append(options_menu, '&Options')
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
        self.colorbar = plt.colorbar(self.image, ax=self.axis, orientation=self.data.colorbar_dir)
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
        Binds the matplotlib and wx events to their method handlers
        """
        # matplotlib events
        self.figure_canvas.mpl_connect('motion_notify_event', self.motion_handler)
        self.figure_canvas.mpl_connect('button_press_event', self.select_point)

        # wx events
        self.Bind(wx.EVT_MENU, self.on_rescale_click, self.rescale_colorbar_menu)
        self.Bind(wx.EVT_MENU, self.on_cscan_signal_noise_click,
                  self.calculate_cscan_sn_ratio_button)
        self.Bind(wx.EVT_MENU, self.on_ascan_signal_noise_click,
                  self.calculate_sn_ratio_menu_button)
        self.Bind(wx.EVT_MENU, self.on_ascan_signal_noise_click_defect_only,
                  self.calculate_sn_ratio_menu_button2)
        self.Bind(wx.EVT_MENU, self.change_colorbar_dir, self.colorbar_dir_menu_button)

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
            self.i_index = int((y_data - self.data.y_min) / self.data.dy)
            self.j_index = int((x_data - self.data.x_min) / self.data.dx)

            self.holder.a_scan_frame.plot(self.i_index, self.j_index)
            self.holder.b_scan_frame.plot(self.i_index, self.j_index)

    def on_rescale_click(self, event):
        """
        Rescales the colorbar with the maximum and minimum values that are
        currently in the visible plot boundary.
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

        # want to include the bounds in calculation, so add 1 to be inclusive
        area = self.data.c_scan[i0:i1+1, j0:j1+1]

        # resent vmin and vmax to be the min and max of area inside of plot bounds
        self.image.set_clim(vmin=area.min(), vmax=area.max())
        self.colorbar.update_normal(self.image)  # update colorbar with new vmin and vmax values
        self.figure_canvas.draw()  # redraw image

    def on_cscan_signal_noise_click(self, event):
        """
        Calculate the signal to noise ratio based on the area that is currently
        in view in the C-Scan. This method is based only on the pixel values
        that are in the C-Scan image.

        Calculate the signal to noise ratio defined as ...
            (Peak Defect - Avg. Noise) / (Peak Noise - Avg. Noise)
        """
        import skimage.filters

        # this is the signal to noise ratio metric that was used during the
        # engine titanium consortium

        x_bounds = np.array(self.axis.get_xlim())
        y_bounds = np.array(self.axis.get_ylim())

        # the y bounds are normally inverted on the C-Scan to show the data as
        # how it is viewed in lab (-y corresponds to the top of the image).
        # Rearrange them in low -> high order
        if y_bounds[0] > y_bounds[1]:
            temp = y_bounds[0]
            y_bounds[0] = y_bounds[1]
            y_bounds[1] = temp

        # it could be possible that x_bounds may be inverted also
        if x_bounds[0] > x_bounds[1]:
            temp = x_bounds[0]
            x_bounds[0] = x_bounds[1]
            x_bounds[1] = temp

        j0 = np.argmin(np.abs(self.data.x - x_bounds[0]))
        j1 = np.argmin(np.abs(self.data.x - x_bounds[1]))

        i0 = np.argmin(np.abs(self.data.y - y_bounds[0]))
        i1 = np.argmin(np.abs(self.data.y - y_bounds[1]))

        # we want to area to be inclusive of the end points so add +1 to them
        area = self.data.c_scan[i0:i1+1, j0:j1+1]

        thresh = skimage.filters.threshold_otsu(area)

        # create a new figure that shows the the thresholded image
        plt.figure('C-Scan After Threshold')
        plt.imshow(area > thresh, cmap='gray', extent=self.axis.axis())
        plt.xlabel('X Scan Location (mm)')
        plt.ylabel('Y Scan Location (mm)')
        plt.grid()
        plt.show()

        peak_defect = area[np.where(area > thresh)].max()
        avg_noise = area[np.where(area < thresh)].mean()
        peak_noise = area[np.where(area < thresh)].max()

        sn_ratio = (peak_defect - avg_noise) / (peak_noise - avg_noise)

        # open a message dialog that displays the signal to noise ratio
        # calculation
        mssg_string = 'Signal to Noise Ratio = %0.4f' % sn_ratio
        title_string = 'Signal to Noise Ratio'
        dlg = wx.MessageDialog(self, mssg_string, title_string, wx.OK |
                               wx.ICON_INFORMATION)
        dlg.ShowModal()

    def on_ascan_signal_noise_click(self, event):
        """
        Calculates the signal to noise ratio based on the signal and noise in
        the A-Scan between the follow gates. This method uses the defect and
        the surrounding waveforms to calculate the signal to noise ratio.
        """
        import skimage.filters

        # this method calculates the signal to noise ratio using the largest
        # peak to peak value from a defect waveform (based on thresholding) as
        # the peak signal value. It then searches all of the non-defect
        # waveforms for the largest peak to peak value and calls that the peak
        # noise value. While searching it rectifies all noise waveforms and
        # averages each waveform to find the average noise value.

        x_bounds = np.array(self.axis.get_xlim())
        y_bounds = np.array(self.axis.get_ylim())

        # the y bounds are normally inverted on the C-Scan to show the data as
        # how it is viewed in lab (-y corresponds to the top of the image).
        # Rearrange them in low -> high order
        if y_bounds[0] > y_bounds[1]:
            temp = y_bounds[0]
            y_bounds[0] = y_bounds[1]
            y_bounds[1] = temp

        # it could be possible that x_bounds may be inverted also
        if x_bounds[0] > x_bounds[1]:
            temp = x_bounds[0]
            x_bounds[0] = x_bounds[1]
            x_bounds[1] = temp

        j0 = np.argmin(np.abs(self.data.x - x_bounds[0]))
        j1 = np.argmin(np.abs(self.data.x - x_bounds[1]))

        i0 = np.argmin(np.abs(self.data.y - y_bounds[0]))
        i1 = np.argmin(np.abs(self.data.y - y_bounds[1]))

        # we want to area to be inclusive of the end points so add +1 to them
        area = self.data.c_scan[i0:i1+1, j0:j1+1]

        thresh = skimage.filters.threshold_otsu(area)

        binary_image = area > thresh

        max_defect = 0
        max_noise = 0
        avg_noise = 0
        count = 0
        for i in range(binary_image.shape[0]):
            ii = i+i0
            for j in range(binary_image.shape[1]):
                jj = j+j0
                # these are the left and right gates for the waveform
                left = self.data.peak_bin[3, 1, ii, jj]
                right = self.data.peak_bin[4, 1, ii, jj]

                # extract the part of the wave that is within the gates
                wave = self.data.waveform[ii, jj, left:right]

                # due to find peaks looking for peaks within pulse length
                # the maximum or minimum of the waveform may not actually be
                # one of the peaks, so use peak bin to get index
                max_pos = self.data.peak_bin[0, 1, ii, jj]
                min_pos = self.data.peak_bin[1, 1, ii, jj]

                max_val = self.data.waveform[ii, jj, max_pos]
                min_val = self.data.waveform[ii, jj, min_pos]
                vpp = max_val - min_val

                # check for max defect signal in waveforms
                if binary_image[i, j] == 1:
                    if vpp > max_defect:
                        max_defect = vpp

                # get average noise and max noise signal
                else:
                    if vpp > max_noise:
                        max_noise = vpp
                    avg_noise += np.mean(np.abs(wave))
                    count += 1

        avg_noise /= count
        sn_ratio = (max_defect-avg_noise) / (max_noise-avg_noise)

        # create a new figure that shows the the thresholded image
        plt.figure('C-Scan After Threshold')
        plt.imshow(area > thresh, cmap='gray', extent=self.axis.axis())
        plt.xlabel('X Scan Location (mm)')
        plt.ylabel('Y Scan Location (mm)')
        plt.grid()
        plt.show()

        # open a message dialog that displays the signal to noise ratio
        # calculation
        mssg_string = ('Signal to Noise Ratio = %0.4f\n'
                       'Max Peak = %0.4f\n'
                       'Max Noise = %0.4f\n'
                       'Avg. Noise = %0.4f'
                       % (sn_ratio, max_defect, max_noise, avg_noise))
        title_string = 'Signal to Noise Ratio'

        dlg = wx.MessageDialog(self, mssg_string, title_string, wx.OK |
                               wx.ICON_INFORMATION)
        dlg.ShowModal()

    def on_ascan_signal_noise_click_defect_only(self, event):
        """
        Calculates the signal to noise ratio based on the signal and noise in
        the A-Scan between the follow gates. This method uses the waveforms
        from the defect only to calculate the signal to noise ratio.
        """
        import skimage.filters

        # this method calculates the signal to noise ratio using the largest
        # peak to peak value from a defect waveform (based on thresholding) as
        # the peak signal value. It then searches all of the non-defect
        # waveforms for the largest peak to peak value and calls that the peak
        # noise value. While searching it rectifies all noise waveforms and
        # averages each waveform to find the average noise value.

        x_bounds = np.array(self.axis.get_xlim())
        y_bounds = np.array(self.axis.get_ylim())

        # the y bounds are normally inverted on the C-Scan to show the data as
        # how it is viewed in lab (-y corresponds to the top of the image).
        # Rearrange them in low -> high order
        if y_bounds[0] > y_bounds[1]:
            temp = y_bounds[0]
            y_bounds[0] = y_bounds[1]
            y_bounds[1] = temp

        # it could be possible that x_bounds may be inverted also
        if x_bounds[0] > x_bounds[1]:
            temp = x_bounds[0]
            x_bounds[0] = x_bounds[1]
            x_bounds[1] = temp

        j0 = np.argmin(np.abs(self.data.x - x_bounds[0]))
        j1 = np.argmin(np.abs(self.data.x - x_bounds[1]))

        i0 = np.argmin(np.abs(self.data.y - y_bounds[0]))
        i1 = np.argmin(np.abs(self.data.y - y_bounds[1]))

        # we want to area to be inclusive of the end points so add +1 to them
        area = self.data.c_scan[i0:i1+1, j0:j1+1]

        thresh = skimage.filters.threshold_otsu(area)

        binary_image = area > thresh

        max_defect = 0
        max_noise = 0
        avg_noise = 0
        count = 0
        for i in range(binary_image.shape[0]):
            ii = i+i0
            print(ii)
            for j in range(binary_image.shape[1]):
                jj = j+j0

                # due to find peaks looking for peaks within pulse length
                # the maximum or minimum of the waveform may not actually be
                # one of the peaks, so use peak bin to get index
                max_pos = self.data.peak_bin[0, 1, ii, jj]
                min_pos = self.data.peak_bin[1, 1, ii, jj]

                max_val = self.data.waveform[ii, jj, max_pos]
                min_val = self.data.waveform[ii, jj, min_pos]
                vpp = max_val - min_val

                if binary_image[i, j] == 1:  # on defect
                    count += 1
                    if vpp > max_defect:
                        max_defect = vpp
                    # return max_noise and avg_noise from this (i, j) waveform
                    max_noise_wave, avg_noise_wave = self._noise_helper(max_pos, min_pos, ii, jj)

                    avg_noise += avg_noise_wave
                    if max_noise_wave > max_noise:
                        max_noise = max_noise_wave

        avg_noise /= count
        sn_ratio = (max_defect-avg_noise) / (max_noise-avg_noise)

        # create a new figure that shows the the thresholded image
        plt.figure('C-Scan After Threshold')
        plt.imshow(area > thresh, cmap='gray', extent=self.axis.axis())
        plt.xlabel('X Scan Location (mm)')
        plt.ylabel('Y Scan Location (mm)')
        plt.grid()
        plt.show()

        # open a message dialog that displays the signal to noise ratio
        # calculation
        mssg_string = ('Signal to Noise Ratio = %0.4f\n'
                       'Max Peak = %0.4f\n'
                       'Max Noise = %0.4f\n'
                       'Avg. Noise = %0.4f'
                       % (sn_ratio, max_defect, max_noise, avg_noise))
        title_string = 'Signal to Noise Ratio'

        dlg = wx.MessageDialog(self, mssg_string, title_string, wx.OK |
                               wx.ICON_INFORMATION)
        dlg.ShowModal()

    def _noise_helper(self, max_pos, min_pos, i, j):
        """
        Helper function to find the max noise and avg noise inside of an A-Scan
        waveform that is located on a defect
        :param max_pos: the index of the maximum peak in the part of the
            waveform that we are interested in. Should be peak_bin[0, 1, i, j]
        :param min_pos: the index of the minimum peak in the part of the
            waveform that we are interested in. Should be peak_bin[1, 1, i, j]
        :param i: The row the waveform is from
        :param j: The column the waveform is from
        """
        # determining where the signal from the defect ends and noise starts
        # could be tricky. I am thinking that I will look for two zero
        # crossings and call that the end of the defect signal. The rest of the
        # A-Scan that is inside of the gates will then be called noise
        cross_counter = 0

        wave = self.data.waveform[i, j, :]

        idx = min(max_pos, min_pos)
        while cross_counter < 2:
            idx -= 1
            if wave[idx] * wave[idx+1] < 0:
                cross_counter += 1
        # this is the farthest left we can go before running into defect signal
        back_left = idx

        idx = max(max_pos, min_pos)
        cross_counter = 0
        while cross_counter < 2:
            idx += 1
            if wave[idx] * wave[idx-1] < 0:
                cross_counter += 1
        # this is the farthest right we can go before running into defect signal
        front_right = idx

        # these are the follow gates that are stored in peak_bin
        left_gate = self.data.peak_bin[3, 1, i, j]
        right_gate = self.data.peak_bin[4, 1, i, j]

        # make sure that we are looking inside of the gates that are set in
        # peak_bin
        if back_left < left_gate:
            back_left = left_gate
        if front_right > right_gate:
            front_right = right_gate

        if front_right == right_gate or back_left == left_gate:
            pdb.set_trace()

        try:
            left_max = wave[left_gate:back_left].max()
            left_min = wave[left_gate:back_left].min()
        except ValueError:
            string = ('There is not enough room before the defect signal\n'
                      'You should move the left gate forwards!')
            raise ValueError(string)
        left_wave = wave[left_gate:back_left]

        try:
            right_max = wave[front_right:right_gate].max()
            right_min = wave[front_right:right_gate].min()
        except ValueError:
            string = ('There is not enough room after the defect signal\n'
                      'You should move the right gate backwards!')
            raise ValueError(string)
        right_wave = wave[front_right:right_gate]

        noise_wave = np.concatenate((left_wave, right_wave))

        # this is the max noise signal as vpp on either side of the defect
        # signal in the A-Scan
        max_noise = max(left_max, right_max) - min(left_min, right_min)

        avg_noise = np.mean(np.abs(noise_wave))

        return max_noise, avg_noise

    def change_colorbar_dir(self, event):
        """
        Changes the orientation of the colorbar
        """
        # TODO should remove colorbar_dir as an attribute to THzData and make
        # TODO an attribute of self instead. That way, Raw C-Scan and
        # TODO interpolated C-Scan can have their own colorbar direction
        if self.data.colorbar_dir == 'horizontal':
            self.data.colorbar_dir = 'vertical'
        else:
            self.data.colorbar_dir = 'horizontal'

        self.colorbar.remove()

        # need to specify which axis the colorbar will be drawn on. Otherwise
        # it will attach itself to the most recent figure object, which may not
        # be the figure from this frame
        self.colorbar = plt.colorbar(self.image, ax=self.axis, orientation=self.data.colorbar_dir)

        self.figure_canvas.draw()
