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
        self.calculate_c_scan_sn_ratio_button = None

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

        description = 'Calculate C-Scan S/N ratio based on current view'
        title = 'Calculate C-Scan S/N Ratio'
        self.calculate_c_scan_sn_ratio_button = wx.MenuItem(self.options_menu, wx.ID_ANY, title,
                                                            description)

        description = 'Calculate A-Scan S/N ratio based on current C-Scan View'
        title = 'Calculate A-Scan S/N Ratio'
        self.calculate_sn_ratio_menu_button = wx.MenuItem(self.options_menu, wx.ID_ANY, title,
                                                          description)

        title = 'Calculate A-Scan S/N Ratio from Defect only'
        description = 'Calculate S/N ratio from defect waveforms only'
        self.calculate_sn_ratio_menu_button2 = wx.MenuItem(self.options_menu, wx.ID_ANY, title,
                                                           description)

        title = 'Change Colorbar Orientation'
        description = 'Changes the colorbar orientation between horizontal and vertical'
        self.colorbar_dir_menu_button = wx.MenuItem(self.options_menu, wx.ID_ANY, title,
                                                    description)

        self.options_menu.Append(self.rescale_colorbar_menu)
        self.options_menu.Append(self.calculate_c_scan_sn_ratio_button)
        self.options_menu.Append(self.calculate_sn_ratio_menu_button)
        self.options_menu.Append(self.calculate_sn_ratio_menu_button2)
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
                  self.calculate_c_scan_sn_ratio_button)
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

        from base_util.base_util import clear_small_defects

        area, thresh, bound_coords = self._sn_define_area_helper()

        binary_image = np.zeros(area.shape)
        binary_image[np.where(area > thresh)] = 1

        binary_image = clear_small_defects(binary_image, 4)

        # create a new figure that shows the the thresholded image
        plt.figure('C-Scan After Threshold')
        plt.imshow(binary_image, cmap='gray', extent=self.axis.axis())
        plt.xlabel('X Scan Location (mm)')
        plt.ylabel('Y Scan Location (mm)')
        plt.grid()
        plt.show()

        peak_defect = area[np.where(binary_image == 1)].max()
        avg_noise = area[np.where(binary_image == 0)].mean()
        peak_noise = area[np.where(binary_image == 0)].max()

        sn_ratio = (peak_defect - avg_noise) / (peak_noise - avg_noise)

        # open a message dialog that displays the signal to noise ratio
        # calculation
        mssg_string = ('Signal to Noise Ratio = %0.4f\n'
                       'Max Peak = %0.4f\n'
                       'Max Noise = %0.4f\n'
                       'Avg. Noise = %0.4f'
                       % (sn_ratio, peak_defect, peak_noise, avg_noise))
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
        # this method calculates the signal to noise ratio using the largest
        # peak to peak value from a defect waveform (based on thresholding) as
        # the peak signal value. It then searches all of the non-defect
        # waveforms for the largest peak to peak value and calls that the peak
        # noise value. While searching it rectifies all noise waveforms and
        # averages each waveform to find the average noise value.

        import matplotlib.patches as mpatches
        from skimage.measure import label, regionprops
        from skimage.morphology import closing, square
        from skimage.segmentation import clear_border
        from base_util.base_util import clear_small_defects, combine_close_defects

        area, thresh, bound_coords = self._sn_define_area_helper()
        i0, i1, j0, j1 = bound_coords

        binary_image = np.zeros(area.shape)
        binary_image[np.where(area > thresh)] = 1
        binary_image = clear_small_defects(binary_image, 5)

        # use skimage closing function to close up the defects that are close
        # to each other. This forces the label function to treat them as one
        # defect. The line below with size 3 seems to close up defects that
        # happen to become separated by a single line of black.
        binary_image = closing(binary_image, square(3))
        binary_image = clear_border(binary_image)

        labeled_image = label(binary_image)

        plt.figure('Labled Image')
        plt.imshow(labeled_image, cmap='gray', extent=self.axis.axis())
        plt.xlabel('X Scan Location (mm)')
        plt.ylabel('Y Scan Location (mm)')
        plt.colorbar()
        plt.grid()

        # expand the bbox for each of the regions to contain a larger area
        bbox_list = list()
        region_list = regionprops(labeled_image)
        for region in region_list:
            bbox = np.asarray(region.bbox)
            bbox_area = (bbox[2]-bbox[0]) * (bbox[3]-bbox[1])
            while bbox_area < region.area * 10:
                bbox[0] -= 1
                bbox[1] -= 1
                bbox[2] += 1
                bbox[3] += 1
                bbox_area = (bbox[2]-bbox[0]) * (bbox[3]-bbox[1])
            bbox_list.append(bbox)

        fig = plt.figure('Raw C-Scan with Bounding Boxes')
        axis = fig.add_subplot(111)
        axis.imshow(area, cmap='gray', extent=self.axis.axis())
        for i in range(len(bbox_list)):
            min_row, min_col, max_row, max_col = bbox_list[i]
            min_x = self.data.x[j0 + bbox_list[i][1]]
            max_x = self.data.x[j0 + bbox_list[i][3]]
            min_y = self.data.y[i0 + bbox_list[i][0]]
            max_y = self.data.y[i0 + bbox_list[i][2]]
            height = max_y - min_y
            width = max_x - min_x
            bbox_rectangle = mpatches.Rectangle((min_x, min_y), width, height,
                                                fill=False, edgecolor='red')
            axis.add_patch(bbox_rectangle)
        plt.xlabel('X Scan Location (mm)')
        plt.ylabel('Y Scan Location (mm)')
        plt.grid()

        # create a new figure that shows the the thresholded image
        fig = plt.figure('C-Scan After Threshold')
        axis = fig.add_subplot(111)
        axis.imshow(binary_image, cmap='gray', extent=self.axis.axis())
        for i in range(len(bbox_list)):
            min_row, min_col, max_row, max_col = bbox_list[i]
            min_x = self.data.x[j0 + bbox_list[i][1]]
            max_x = self.data.x[j0 + bbox_list[i][3]]
            min_y = self.data.y[i0 + bbox_list[i][0]]
            max_y = self.data.y[i0 + bbox_list[i][2]]
            height = max_y - min_y
            width = max_x - min_x
            bbox_rectangle = mpatches.Rectangle((min_x, min_y), width, height,
                                                fill=False, edgecolor='red')
            axis.add_patch(bbox_rectangle)
        plt.xlabel('X Scan Location (mm)')
        plt.ylabel('Y Scan Location (mm)')
        plt.grid()

        plt.show(block=False)

        defect_coords = combine_close_defects(region_list, bbox_list)
        print(len(defect_coords))
        print(len(region_list))
        test_image = np.zeros(area.shape)
        for i, defect in enumerate(defect_coords):
            for loc in defect:
                test_image[loc[0], loc[1]] = i + 1

        fig = plt.figure('Test Image')
        axis = fig.add_subplot(111)
        im = axis.imshow(test_image, cmap='gray', extent=self.axis.axis())
        for i in range(len(bbox_list)):
            min_row, min_col, max_row, max_col = bbox_list[i]
            min_x = self.data.x[j0 + bbox_list[i][1]]
            max_x = self.data.x[j0 + bbox_list[i][3]]
            min_y = self.data.y[i0 + bbox_list[i][0]]
            max_y = self.data.y[i0 + bbox_list[i][2]]
            height = max_y - min_y
            width = max_x - min_x
            bbox_rectangle = mpatches.Rectangle((min_x, min_y), width, height,
                                                fill=False, edgecolor='red')
            axis.add_patch(bbox_rectangle)
        plt.xlabel('X Scan Location (mm)')
        plt.ylabel('Y Scan Location (mm)')
        plt.grid()
        plt.colorbar(im)

        plt.show(block=False)

        sn_list = list()
        max_defect_list = list()
        max_noise_list = list()
        avg_noise_list = list()
        for i, region in enumerate(region_list):
            bbox = bbox_list[i]
            max_defect = 0
            max_noise = 0
            avg_noise = 0
            count = 0
            for i in range(bbox[0], bbox[2]):
                ii = i+i0
                for j in range(bbox[1], bbox[3]):
                    jj = j+j0
                    # these are the left and right gates for the waveform
                    try:
                        left = self.data.peak_bin[3, 1, ii, jj]
                        right = self.data.peak_bin[4, 1, ii, jj]
                    except IndexError:
                        pdb.set_trace()

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

            sn_list.append(sn_ratio)
            max_defect_list.append(max_defect)
            max_noise_list.append(max_noise)
            avg_noise_list.append(avg_noise)

        for i in range((len(sn_list))):
            sn_ratio = sn_list[i]
            max_defect = max_defect_list[i]
            max_noise = max_noise_list[i]
            avg_noise = avg_noise_list[i]
            print('\nDefect %d' % (i+1))
            print('S/N Ratio %0.4f' % sn_ratio)
            print('Max Defect:  %0.4f' % max_defect)
            print('Max Noise: %0.4f' % max_noise)
            print('Avg. Noise: %0.4f' % avg_noise)

    def on_ascan_signal_noise_click_defect_only(self, event):
        """
        Calculates the signal to noise ratio based on the signal and noise in
        the A-Scan between the follow gates. This method uses the waveforms
        from the defect only to calculate the signal to noise ratio.
        """
        from skimage.measure import label, regionprops
        from skimage.morphology import closing, square
        from skimage.segmentation import clear_border
        from base_util.base_util import clear_small_defects

        area, thresh, bound_coords = self._sn_define_area_helper()
        i0, i1, j0, j1 = bound_coords

        binary_image = np.zeros(area.shape)
        binary_image[np.where(area > thresh)] = 1
        binary_image = clear_small_defects(binary_image, 5)

        # use skimage closing function to close up the defects that are close
        # to each other. This forces the label function to treat them as one
        # defect and makes
        binary_image = closing(binary_image, square(3))
        binary_image = clear_border(binary_image)

        labeled_image = label(binary_image)

        plt.figure('Labled Image')
        plt.imshow(labeled_image, cmap='gray', extent=self.axis.axis())
        plt.xlabel('X Scan Location (mm)')
        plt.ylabel('Y Scan Location (mm)')
        plt.colorbar()
        plt.grid()

        # create a new figure that shows the the thresholded image
        plt.figure('C-Scan After Threshold')
        plt.imshow(binary_image, cmap='gray', extent=self.axis.axis())
        plt.xlabel('X Scan Location (mm)')
        plt.ylabel('Y Scan Location (mm)')
        plt.grid()
        plt.show()

        sn_list = list()
        max_defect_list = list()
        max_noise_list = list()
        avg_noise_list = list()

        # we don't care about bbox for this method because we are calculating
        # noise from A-Scan on defect only
        for region in regionprops(labeled_image):
            max_defect = 0
            max_noise = 0
            avg_noise = 0
            for loc in region.coords:
                i = loc[0]
                j = loc[1]
                ii = i+i0
                jj = j+j0

                # due to find peaks looking for peaks within pulse length
                # the maximum or minimum of the waveform may not actually be
                # one of the peaks, so use peak bin to get index
                max_pos = self.data.peak_bin[0, 1, ii, jj]
                min_pos = self.data.peak_bin[1, 1, ii, jj]

                max_val = self.data.waveform[ii, jj, max_pos]
                min_val = self.data.waveform[ii, jj, min_pos]
                vpp = max_val - min_val

                if vpp > max_defect:
                    max_defect = vpp

                max_noise_wave, avg_noise_wave = self._noise_helper(max_pos, min_pos, ii, jj)

                avg_noise += avg_noise_wave
                if max_noise_wave > max_noise:
                    max_noise = max_noise_wave

            avg_noise /= region.area
            sn_ratio = (max_defect-avg_noise) / (max_noise-avg_noise)

            sn_list.append(sn_ratio)
            max_defect_list.append(max_defect)
            max_noise_list.append(max_noise)
            avg_noise_list.append(avg_noise)

        for i in range((len(sn_list))):
            sn_ratio = sn_list[i]
            max_defect = max_defect_list[i]
            max_noise = max_noise_list[i]
            avg_noise = avg_noise_list[i]
            print('\nDefect %d' % (i+1))
            print('S/N Ratio %0.4f' % sn_ratio)
            print('Max Defect:  %0.4f' % max_defect)
            print('Max Noise: %0.4f' % max_noise)
            print('Avg. Noise: %0.4f' % avg_noise)

        plt.figure('Labled Image')
        plt.imshow(labeled_image, cmap='gray', extent=self.axis.axis())
        plt.xlabel('X Scan Location (mm)')
        plt.ylabel('Y Scan Location (mm)')
        plt.colorbar()
        plt.grid()

        # create a new figure that shows the the thresholded image
        plt.figure('C-Scan After Threshold')
        plt.imshow(binary_image, cmap='gray', extent=self.axis.axis())
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
        :return max_noise: The maximum noise signal in that waveform
        :return avg_noise: The average noise in that waveform
        """
        # determining where the signal from the defect ends and noise starts
        # could be tricky. I am thinking that I will look for two zero
        # crossings and call that the end of the defect signal. The rest of the
        # A-Scan that is inside of the gates will then be called noise

        wave = self.data.waveform[i, j, :]

        idx = min(max_pos, min_pos)
        cross_counter = 0
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

        right_only = False
        try:
            left_max = wave[left_gate:back_left].max()
            left_min = wave[left_gate:back_left].min()
        except ValueError:
            string = ('There is not enough room before the defect signal\n'
                      'You should move the left gate forwards!')
            # print(string)
            left_max = 0
            left_min = 0
            right_only = True
        left_wave = wave[left_gate:back_left]

        left_only = False
        try:
            right_max = wave[front_right:right_gate].max()
            right_min = wave[front_right:right_gate].min()
        except ValueError:
            string = ('There is not enough room after the defect signal\n'
                      'You should move the right gate backwards!')
            # print(string)
            right_max = 0
            right_min = 0
            left_only = True
        right_wave = wave[front_right:right_gate]

        if left_only and right_only:
            raise ValueError('Something went very wrong!')
        elif left_only:
            noise_wave = left_wave
        elif right_only:
            noise_wave = right_wave
        else:
            noise_wave = np.concatenate((left_wave, right_wave))

        # this is the max noise signal as vpp on either side of the defect
        # signal in the A-Scan
        max_noise = max(left_max, right_max) - min(left_min, right_min)

        avg_noise = np.mean(np.abs(noise_wave))

        return max_noise, avg_noise

    def _sn_define_area_helper(self):
        """
        Method to determine the area in the image that was visible when a
        signal to noise calculation button was hit. The keeps common code in
        one place that will be called for any signal to noise ratio calculation
        :return area: The area that is inside of the x & y limits in the C-Scan
        :return thresh: The threshold for the area to create a binary image.
        :return coord_bounds: A tuple that has the i & j indices that bound the
            area of interest. Returned as (i0, i1, j0, j1).
        """
        import skimage.filters

        # using the get_bound method instead of get_lim method returns the
        # lower and upper bounds of the image in increasing order regardless of
        # which is left or right; or top or bottom. This will allow us to
        # access the data x & y array without rearranging the boundary values
        x_bounds = np.array(self.axis.get_xbound())
        y_bounds = np.array(self.axis.get_ybound())

        j0 = np.argmin(np.abs(self.data.x - x_bounds[0]))
        j1 = np.argmin(np.abs(self.data.x - x_bounds[1]))

        i0 = np.argmin(np.abs(self.data.y - y_bounds[0]))
        i1 = np.argmin(np.abs(self.data.y - y_bounds[1]))

        # we want to area to be inclusive of the end points so add +1 to them
        area = self.data.c_scan[i0:i1+1, j0:j1+1]

        # get the y bounds again with get_ylim() so they are in bottom-top
        # order to compare to c_scan_extent
        y_bounds = np.array(self.axis.get_ylim())

        # the threshold that we want to use depends on whether or not the full
        # sample is in view. Yen's Threshold seems to do better on the full
        # sample, while Otsu's Threshold does better on a smaller area with one
        # or two defects visible
        # if np.array_equal(bounds, self.data.c_scan_extent):
        #     thresh = skimage.filters.threshold_yen(area)
        # else:
        #     thresh = skimage.filters.threshold_otsu(area)

        thresh = skimage.filters.threshold_triangle(area)
        skimage.filters.try_all_threshold(area)

        # put the coords in a tuple to return them
        coord_bounds = (i0, i1, j0, j1)

        return area, thresh, coord_bounds

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
