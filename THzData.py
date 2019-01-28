from __future__ import print_function, division
import os
import copy
import struct
import sys
import pdb

import numpy as np
import skimage.filters

from THzProc.thz_functions import AmpCor300, FindPeaks, ReMap

# THINGS THAT STILL HAVE TO BE IMPLEMENTED
# Depth Map
# trend off != 0


# inherit object class so this will still be a new-style class even if run
# using python 2
class THzData(object):
    """
    Class representation of THz data from a tvl file. Can be used in other
    scripts by instantiating the class like so ...
    data = THzData(filename, basedir, etc.)
    """

    def __init__(self, filename, basedir=None, gate=None, follow_gate_on=True,
                 signal_type=1, center=False, print_on=True, trend_off=0,
                 pulse_length=3, x_cor_tol=0.1):
        """
        Constructor method
        :param filename: Either the filename or full path to the tvl file. If
            only the filename is given, basedir must be provided
        :param basedir: The base directory for the tvl file. If left as None
            (default), it will be assumed that filename contains the full path
        :param gate: The 2x2 gate that is used to bracket the front surface and
            the portion of the inside of the sample that we are interested in.
        :param follow_gate_on: If True (default), follow gate is on. If False,
            follow gate is off.
        :param center: If True, will adjust x and y axis values so that (0, 0)
            is at the center of the C-Scan image (default False)
        :param print_on: If True (default), will print information about the
            scan to the console, such as number of X & Y steps, resolution, etc.
        :param trend_off: Method to choose to remove the baseline trend in near-
            field imaging. Leave as 0 (default) to do not baseline removal.
        :param pulse_length: The length in picoseconds within which to search
            when determining peak to peak amplitude.
        """

        # handle what happens if they don't pass through gate
        if gate is None:
            self.gate = [[100, 1500], [-200, 200]]

        # if the user passed through a value for gate, make sure that it is 2x2
        elif np.shape(gate) != (2, 2):
            raise ValueError('gate must by a 2x2 list or numpy array!')
        else:
            # make self.gate a separate copy of the gate that is passed into
            # the constructor. Otherwise if self.gate is changed by this program
            # the variable that was passed in will change as well
            self.gate = copy.deepcopy(gate)

        # make sure that follow_gate_on is either True or False
        if not isinstance(follow_gate_on, bool):
            raise ValueError('follow_gate_on must be of a boolean!')
        else:
            self.follow_gate_on = follow_gate_on

        # Open the file using TeraView's DataFile class
        dat = DataFile(filename, basedir=basedir)

        # set signal_type to pass to make_c_scan()
        self.signal_type = signal_type

        self.x = dat.data[3:]['x']
        self.y = dat.data[3:]['y']
        self.waveform = dat.data['waveform'][3:]
        self.wave_length = len(self.waveform[0])

        # B Scan direction is usually horizontal by default
        self.b_scan_dir = 'horizontal'

        # whether or not to correct for the excessive amplification on the
        # edges of the 300 ps waveforms
        self.amp_correction300_on = True

        # initialize all of the variables that are wanted from the header,
        # these are calculated in the method header_info()
        self.time_length = None  # the time length of the scan in ps
        self.x_res = None  # the attempted spacing between x points in the scan in mm
        self.y_res = None  # the attempted spacing between y points in the scan in mm
        self.x_step = None  # the number of steps per row on major axis (usually x-dimension)
        self.y_step = None  # the number of steps per column (usually y-dimension)
        self.x_min = None  # the smallest x value (mm)
        self.y_min = None  # the smallest y value (mm)
        self.x_max = None  # the largest x value (mm)
        self.y_max = None  # the largest y value (mm)
        self.scan_type = None  # the type of scan performed usually '2D Image Scan with Encoder'

        # the first axis, usually x, but can be turntable if rotational scan in performed.
        self.axis = None

        # the frequency domain representation of the data between either the
        # follow gates, if follow_gate_on is True or the front gates
        self.freq_waveform = np.ndarray

        self.dt = None  # the spacing between time values
        self.time = None  # array of time values that is used for plotting
        self.freq = None  # array of frequency values that is used for plotting
        self.df = None  # spacing between frequency values

        # the number of data points in a half pulse
        self.n_half_pulse = int

        self.true_x_res = None  # the spacing between x points after Remap is called
        self.true_y_res = None  # the spacing between y points after Remap is called
        self.dx = None  # the average difference between x data points
        self.dy = None  # the average difference between y data points
        self.b_scan = None  # the B-Scan values
        self.tof_c_scan = None  # the time of flight of the front surface

        # if True there will be little crosses displayed on the peak locations
        # in the lead gates and the follow
        self.flag_peak_on = True

        # incoming angle of the THz Beam (17.5 degrees) converted to radians
        self.theta0 = 17.5 * np.pi / 180

        # whether or not a call to resize() has been issued. This attribute is
        # for other functions that base their decisions on whether of not the
        # resize method has been called
        self.has_been_resized = False

        # whether or not the C-Scan will be generated using frequency domain
        # information. Default to False because we are usually interested in
        # time domain info.
        self.in_freq_mode = False

        # CONSTANTS ------------------------------------------------------------

        # 1 for the positive peak, 2 for the negative peak, prefer 1
        self.PICK_PEAK = 1

        # provides algorithm for removing baseline trend in near-field imaging.
        self.TREND_OFF = trend_off

        # VERY constant: do not change unless you have a good reason to
        self.FSE_TOLERANCE = -0.17

        # parameters to pass to the AmpCor300 function
        self.AMP_CORRECTION_300_PAR = [0., 35., 5.0, 1., 240., 300., 1., 4.5,
                                       4.]

        # controls how much (as a multiplier * resolution) how much the lines
        # shift during remapping
        self.X_CORRECTION_TOLERANCE = x_cor_tol
        # self.X_CORRECTION_TOLERANCE = 0.1
        # self.X_CORRECTION_TOLERANCE = 0.75  # for 50 um resolution
        # self.X_CORRECTION_TOLERANCE = 2.0  # for 25 um resolution
        # self.X_CORRECTION_TOLERANCE = 2.5  # for 15 um resolution

        # TODO: (John): The X_CORRECTION_TOLERANCE seems to be a constant value
        # TODO: instead of a constant multiplied by the scan resolution, this
        # TODO: value seems to be around 0.35 ~ 0.40 mm, might actually be 50

        # half pulse width for bracketing the front surface echo this value was
        # originally 2 in Thomas's THzProc
        self.HALF_PULSE = 2

        # the number of half pulses within which to search for the peak to
        # peak value within the follow gate. For example if the maximum peak
        # in the follow gates is found, then the minimum peak that defines the
        # Vpp with can only a maximum of PULSE_LENGTH * HALF_PULSE ps away. The
        # length of the half pulse is defined by HALF_PULSE
        self.PULSE_LENGTH = pulse_length

        # difference that a point is allowed to deviate (as a ratio with
        # respect to resolution) from its desired coordinate
        self.X_DIFFERENCE_TOLERANCE = 0.4
        self.Y_DIFFERENCE_TOLERANCE = 0.3

        # tolerance ratio for number of scan points in a X line
        # this is nlim from base THzProc
        self.N_LIMIT = 0.8

        # threshold for lead gate signal
        self.FSE_THRESHOLD = 0.5

        # tolerance ratio for number of actual scan points in a X line compared
        # to how many are supposed to be in that line
        self.X_RATIO_LIMIT = 0.8

        # a small value that is used to see if floating points are close
        # together.
        self.TINY = 1e-4

        # END OF CONSTANTS -------------------------------------------------------------------------

        # Begin driving the methods
        self.header_info(dat.header)
        self.delta_calculator()

        if print_on:
            self.printer()

        # call remap to correct for the backlash of the system
        # the 0 parameter is for XChkTol, which is not used in the function
        self.waveform, self.c_scan, self.x, self.y, self.pos, self.x_step, self.y_step = \
            ReMap(self.waveform, self.x, self.y, self.x_max, self.x_min, self.y_min,
                  self.true_x_res, self.true_y_res, self.x_step, self.y_step, self.scan_type,
                  self.axis, self.wave_length, self.time_length, self.X_CORRECTION_TOLERANCE,
                  self.X_DIFFERENCE_TOLERANCE, self.Y_DIFFERENCE_TOLERANCE, 0, self.TINY,
                  self.TREND_OFF, self.N_LIMIT)

        # have this after ReMap because x may be adjusted during call to ReMap
        if center:
            self.center_coordinates()

        # run a check to make sure that the time length is actually 300
        if self.amp_correction300_on and abs(300 - self.wave_length) < self.TINY:
            AmpCor300(self.AMP_CORRECTION_300_PAR, self.wave_length, self.dt, self.x_step,
                      self.y_step, self.waveform)

        if self.follow_gate_on:
            self.bin_range = copy.deepcopy(self.gate)
            self.peak_bin = np.zeros((5, 2, self.x_step, self.y_step))
        elif not self.follow_gate_on:
            self.bin_range = [[0, self.wave_length]]
            self.peak_bin = np.zeros((5, 1, self.x_step, self.y_step))
        else:
            raise ValueError('follow_gate_on must be set to either True or False')

        # Return peak_bin to keep track of index values for various points on the waveform
        # peak_bin is 5 dimensional
        # index 0 = positive peak
        # index 1 = negative peak
        # index 2 = half way between positive and negative peak
        # index 3 = left gate index
        # index 4 = right gate index
        self.find_peaks()

        # do not use x_min and y_min for extent because they do not represent
        # the true x and y values after remap has been called.
        # Add half of a pixel length to all sides of the extent boundary
        # because matplotlib uses the coordinate system as the center of pixels.
        # On an image with (i, j) index coordinates, the center of the top left
        # pixel is treated as (0, 0), see documentation on plt.imshow() here
        # https://matplotlib.org/api/_as_gen/matplotlib.pyplot.imshow.html
        # especially look at the description of `extent` parameter

        self.c_scan_extent = (self.x[0] - self.true_x_res/2,
                              self.x[-1] + self.true_x_res/2,
                              self.y[-1] + self.true_y_res/2,
                              self.y[0] - self.true_y_res/2)

        self.b_scan_extent = (self.x[0] - self.true_x_res/2,
                              self.x[-1] + self.true_x_res/2,
                              0 - self.dt/2,
                              self.time_length + self.dt/2)

        # Generate the C-Scan
        self.make_c_scan(self.signal_type)

    def make_c_scan(self, signal_type=None):
        """
        Generates the C-Scan based on the gate location and signal type
        At the moment this just does peak to peak voltage response
        :param: signal_type - determines how the C-Scan is calculated
                    choices:
                        0: (default) Use Peak to Peak voltage with the front
                           gates regardless of whether follow gate is on or not.
                        1: Use Peak to Peak voltage with the follow gates if on.
                           If follow gate is not on then use peak to peak
                           voltage across entire waveform
                        2: avg mag (mag = abs amp)
                        3: median mag
                        4: min mag
                        5: max mag
                        6: avg amp ~ integrated gate
                        7: median amp
                        8: min amp
                        9: max amp
                        10: ratio between FSE/BSE.
        """

        # if user does not pass through a value default to the one that is set
        # by the class instance
        if signal_type is None:
            signal_type = self.signal_type

        if self.follow_gate_on:
            idx = 1
        else:
            idx = 0

        # TODO: work on vectorizing code to avoid the for loops
        # it seems that if using peak_bin we can't vectorize, look into solution for this

        # use Vpp within the front gates regardless of whether follow gate is on or not
        # the gates are ABSOLUTE INDEX POSITION and DO NOT account for differences in height of the
        # front surface
        if signal_type == 0:
            max_amp = np.amax(self.waveform[:, :, self.gate[0][0]:self.gate[0][1]], axis=2)
            min_amp = np.amin(self.waveform[:, :, self.gate[0][0]:self.gate[0][1]], axis=2)
            self.c_scan = max_amp - min_amp

        # use Vpp within the follow gates if on, else use Vpp across entire A-Scan
        # It appears that using peak bin does not allow for vectorization
        elif signal_type == 1:
            for i in range(self.y_step):
                for j in range(self.x_step):
                    max_amp = self.waveform[i, j, self.peak_bin[0, idx, i, j]]
                    min_amp = self.waveform[i, j, self.peak_bin[1, idx, i, j]]
                    self.c_scan[i, j] = max_amp - min_amp

        # use avg mag within the gate (abs amplitude within the gate, sum up then avg)
        elif signal_type == 2:
            for i in range(self.y_step):
                for j in range(self.x_step):
                    L = self.peak_bin[3, self.follow_gate_on, i, j]
                    R = self.peak_bin[4, self.follow_gate_on, i, j]
                    self.c_scan[i, j] = np.sum(np.abs(self.waveform[i, j, L:R])) / (R - L)

        # median magnitude within the gate
        elif signal_type == 3:
            for i in range(self.y_step):
                for j in range(self.x_step):
                    L = self.peak_bin[3, self.follow_gate_on, i, j]
                    R = self.peak_bin[4, self.follow_gate_on, i, j]
                    self.c_scan[i, j] = np.median(np.abs(self.waveform[i, j, L:R]))

        # min magnitude within the gate
        elif signal_type == 4:
            for i in range(self.y_step):
                for j in range(self.x_step):
                    L = self.peak_bin[3, self.follow_gate_on, i, j]
                    R = self.peak_bin[4, self.follow_gate_on, i, j]
                    self.c_scan[i, j] = np.min(np.abs(self.waveform[i, j, L:R]))

        # max magnitude within the gate
        elif signal_type == 5:
            for i in range(self.y_step):
                for j in range(self.x_step):
                    L = self.peak_bin[3, self.follow_gate_on, i, j]
                    R = self.peak_bin[4, self.follow_gate_on, i, j]
                    self.c_scan[i, j] = np.amax(np.abs(self.waveform[i, j, L:R]))

        # avg amp within the gate ~ integrated gate
        elif signal_type == 6:
            for i in range(self.y_step):
                for j in range(self.x_step):
                    L = self.peak_bin[3, self.follow_gate_on, i, j]
                    R = self.peak_bin[4, self.follow_gate_on, i, j]
                    self.c_scan[i, j] = np.sum(self.waveform[i, j, L:R]) / (R - L)

        # median amp within the gate
        elif signal_type == 7:
            for i in range(self.y_step):
                for j in range(self.x_step):
                    L = self.peak_bin[3, self.follow_gate_on, i, j]
                    R = self.peak_bin[4, self.follow_gate_on, i, j]
                    self.c_scan[i, j] = np.median(self.waveform[i, j, L:R])

        # min amp within the gate
        elif signal_type == 8:
            for i in range(self.y_step):
                for j in range(self.x_step):
                    L = self.peak_bin[3, self.follow_gate_on, i, j]
                    R = self.peak_bin[4, self.follow_gate_on, i, j]
                    self.c_scan[i, j] = np.min(self.waveform[i, j, L:R])

        # max amp within the gate
        elif signal_type == 9:
            for i in range(self.y_step):
                for j in range(self.x_step):
                    L = self.peak_bin[3, self.follow_gate_on, i, j]
                    R = self.peak_bin[4, self.follow_gate_on, i, j]
                    self.c_scan[i, j] = np.amax(self.waveform[i, j, L:R])

        # use the peak to peak amplitude of the current follow-gate and then
        # find the max and minimum value of the wave that is left over after
        # the rear follow gate
        elif signal_type == 10:
            for i in range(self.y_step):
                for j in range(self.x_step):
                    # get Vpp within the follow gate
                    L = self.peak_bin[3, idx, i, j]
                    R = self.peak_bin[4, idx, i, j]
                    fse_max = np.amax(self.waveform[i, j, L:R])
                    fse_min = np.amin(self.waveform[i, j, L:R])

                    # assume that other Vpp occurs after this gate
                    bse_max = np.amax(self.waveform[i, j, R:])
                    bse_min = np.amin(self.waveform[i, j, R:])

                    # make ration of FSE / BSE
                    self.c_scan[i, j] = (fse_max-fse_min) / (bse_max-bse_min)

        # 27JAN2015 FSE handle: if FSE too small, throw away this location
        if signal_type != 0 and self.FSE_TOLERANCE > 0:
            for i in range(self.y_step):
                for j in range(self.x_step):
                    if (self.waveform[i, j, self.peak_bin[0, 0, i, j]] -
                            self.waveform[i, j, self.peak_bin[1, 0, i, j]] < self.FSE_TOLERANCE):
                        self.c_scan[i, j] = 0.

    def make_freq_c_scan(self, signal_type, gate0, gate1=None, is_index=True):
        """
        Generates a frequency domain C-Scan
        :param signal_type: The type of frequency C-Scan is generated,
                = 0: integrate between gate0 and gate1
                = 1: value at gate0 only, single index slice
        :param gate0: The 1st gate
        :param gate1: The 2nd gate
        :param indexing: Whether the values in gate0 and gate1 are the index to
            the frequency array or the frequency value itself. Should be either
            'index' or 'frequency'
        """

        if not is_index:
            gate0 = np.argmin(np.abs(self.freq - gate0))
            gate1 = np.argmin(np.abs(self.freq - gate1))

        if signal_type == 0 and gate1 is not None:  # integrate between gates
            magnitude = np.abs(self.freq_waveform[:, :, gate0:gate1+1])
            self.c_scan = np.sum(magnitude, axis=2) * self.df

        if signal_type == 1:
            self.c_scan = np.abs(self.data.freq_waveform[:, :, gate0])

    def make_time_of_flight_c_scan(self):
        """
        Generates a time of flight C-Scan that shows the time of flight in
        picoseconds of each pixel on the front surface.
        """
        self.tof_c_scan = np.zeros((self.y_step, self.x_step))

        # if self.follow_gate_on:
        #     tof_index = self.peak_bin[0, 1, :, :]
        # else:
        #     tof_index = self.peak_bin[0, 0, :, :]

        tof_index = self.peak_bin[0, 0, :, :]

        for i in range(self.y_step):
            for j in range(self.x_step):
                self.tof_c_scan[i, j] = self.time[tof_index[i, j]]

    def gate_and_take_fft(self, gate=None):
        """
        Generates a frequency domain representation of the data by taking the
        FFT of the part of the time domain waveforms that is gated by the
        follow gates.
        """

        # set up gated waveform so the area of interest from waveform (defined
        # by the gates in peak_bin) still has the proper number of data points
        # which is (wave_length//2 + 1) eg. 4096//2 + 1 = 2049
        gated_waveform = np.zeros(self.waveform.shape)

        if self.follow_gate_on:
            idx = 1
        else:
            idx = 0

        # extract the part of the waveform that we are interested in and put it
        # in gated_waveform so it is properly zero-padded. This also allows us
        # to take the FFT in one shot.
        for i in range(self.y_step):
            for j in range(self.x_step):
                left_gate = self.peak_bin[3, idx, i, j]
                right_gate = self.peak_bin[4, idx, i, j]
                gated_waveform[i, j, left_gate:right_gate] = \
                    self.waveform[i, j, left_gate:right_gate]

        self.freq_waveform = np.fft.rfft(gated_waveform, axis=2)

    def resize(self, indices, indexing='xy', return_indices=False):
        """
        Resizes the data in the bounds between x0, x1, y0, and y1. Should be
        used to remove the edges from the data if it was over scanned. This
        method creates attributes waveform_small, c_scan_small, and if time of
        flight c-scan has already been calculated it also creates
        tof_c_scan_small.
        :param indices: Expected to be an iterable containing the values from
            which to bound the data. If indexing=('xy') (Default), the iterable
            is expected to contain (x_min, x_max, y_min, y_max). If
            indexing='ij', it is expected to contain
            (i_min, i_max, j_min, j_max).

            Note that normally in the C-Scan that is output by THzProc the -y
            values are on top of the C-Scan image.

        :param return_indecices: If passed as True will return the indices that
            were used to generate the small C-Scan as (i0, i1, j0, j1). Where
            i0 is the top most index. i1 is bottom most index. j0 is the left
            most index, and j1 is the right most.
        """

        self.has_been_resized = True

        if indexing == 'xy':
            # if the 'xy' indexing is used, need to convert to (i, j). Do this
            # by getting the indices of the closest value in the x and y array
            # to ones that the user passed through
            x0 = indices[0]
            x1 = indices[1]
            y0 = indices[2]
            y1 = indices[3]

            j0 = np.argmin(np.abs(self.x - x0))
            j1 = np.argmin(np.abs(self.x - x1))
            i0 = np.argmin(np.abs(self.y - y0))
            i1 = np.argmin(np.abs(self.y - y1))

        elif indexing == 'ij':
            i0 = indices[0]
            i1 = indices[1]
            j0 = indices[2]
            j1 = indices[3]

        else:
            raise ValueError('Parameter indexing must be either "xy" or "ij"')

        # Make sure to add +1 to all of the numpy slicing because we want to
        # include the index that is as close as possible to the parameters that
        # are passed in. Normal numpy slicing is exclusive of the endpoint.
        self.x_small = self.x[j0:j1+1]
        self.x_step_small = len(self.x_small)

        self.y_small = self.y[i0:i1+1]
        self.y_step_small = len(self.y_small)

        self.waveform_small = self.waveform[i0:i1+1, j0:j1+1, :]

        self.c_scan_small = self.c_scan[i0:i1+1, j0:j1+1]

        self.small_extent = (self.x_small[0] - self.true_x_res/2,
                             self.x_small[-1] + self.true_x_res/2,
                             self.y_small[-1] + self.true_y_res/2,
                             self.y_small[0] - self.true_y_res/2)

        if self.tof_c_scan is not None:
            self.tof_c_scan_small = self.tof_c_scan[i0:i1+1, j0:j1+1]

        if return_indices:
            return (i0, i1, j0, j1)

    def center_coordinates(self):
        """
        Centers the x and y coordinates such that (0, 0) is in the center of
        the image.
        """
        self.x -= ((self.x[0] + self.x[-1]) / 2)
        self.y -= ((self.y[0] + self.y[-1]) / 2)

        # update max and min of x and y
        self.x_min = self.x[0]
        self.x_max = self.x[-1]
        self.y_min = self.y[0]
        self.y_max = self.y[-1]

        # update extent
        self.c_scan_extent = (self.x[0], self.x[-1], self.y[-1], self.y[0])

    def adjust_coordinates(self, i, j):
        """
        Adjusts the coordinate system so that the point at index (i, j) is at
        (0, 0).
        :param i: The index of the new center row
        :param j: The index of the new center column
        """
        self.x -= self.x[j]
        self.y -= self.y[i]

        # reset min and max values for each coordinate
        self.x_min = self.x[0]
        self.x_max = self.x[-1]
        self.y_min = self.y[0]
        self.y_max = self.y[-1]

        # reset extent
        self.c_scan_extent = (self.x[0], self.x[-1], self.y[-1], self.y[0])

    def make_b_scan(self, yid, xid):
        """
        Generates the B-Scan based on the last cursor location clicked and
        whether b_scan_dir is horizontal or vertical
        :param yid: the row from which to generate B-Scan
        :param xid: the column clicked from which to generate B-Scan
        """
        # in order to vectorize this method, the x or y coordinates are put in the rows, but to plot
        # we want them in the columns. Thus the call to transpose after arranging the data
        if self.b_scan_dir == 'horizontal':
            self.b_scan = self.waveform[yid, :, :]
            self.b_scan_extent = (self.x[0], self.x[-1], 0, self.time[-1])
        else:  # b_scan_dir == 'vertical'
            self.b_scan = self.waveform[:, xid, :]
            self.b_scan_extent = (self.y[0], self.y[-1], 0, self.time[-1])

        # call to transpose is necessary to flip data axis, so x or y location is on the bottom and
        # time is along the y-axis on the plot
        self.b_scan = np.transpose(self.b_scan)
        self.b_scan = np.flipud(self.b_scan)  # flip so top of sample is at bottom of image

    def set_follow_gate(self, given_boolean):
        """
        Sets the follow gate to be either on or off. If follow gate is changed,
        code will update bin_range and peak_bin.
        :param given_boolean: True: follow gate is on. False: follow gate is
            off.
        """
        # if user tries to set it to the value that it already is, do nothing
        if self.follow_gate_on == given_boolean:
            return

        self.follow_gate_on = given_boolean

        # if follow gate changes reset bin_range and peak_bin appropriately
        if self.follow_gate_on:
            self.bin_range = copy.deepcopy(self.gate)
            self.peak_bin = np.zeros((5, 2, self.x_step, self.y_step))
        elif not self.follow_gate_on:
            self.bin_range = [[0, self.wave_length]]
            self.peak_bin = np.zeros((5, 1, self.x_step, self.y_step))
        else:
            raise ValueError('follow_gate_on must be set to either True or False')

        self.find_peaks()

    def change_gate(self, incoming_gate):
        """
        Change the gate locations by providing an a new gate. Then updates
        bin_range and peak_bin accordingly. This method does NOT create a new
        C-Scan image based on the new gate.
        :param incoming_gate: The value of gate that is to be changed to, must
            be a 2x2 list or numpy array
        """
        import time

        # incoming_gate is the same as what gate currently is; do nothing
        if np.array_equal(self.gate, incoming_gate):
            return
        # make sure that incoming_gate is a 2x2 array
        elif np.shape(np.asarray(incoming_gate)) != (2, 2):
            raise ValueError('gate must by a 2x2 list or numpy array!')

        t0 = time.time()

        # if the front gates haven't changed we can avoid the call to
        # find_peaks and speed things up I think.
        if np.array_equal(self.gate[0], incoming_gate[0]):
            self.gate[1] = incoming_gate[1]
            self.bin_range[1] = incoming_gate[1]

            # adjust the left and right follow gates of the 2nd layer based on
            # front surface peaks and the incoming gate
            left_gate = self.peak_bin[0, 0, :, :] + incoming_gate[1][0]
            right_gate = self.peak_bin[0, 0, :, :] + incoming_gate[1][1]
            left_gate[np.where(left_gate < 0)] = 0
            right_gate[np.where(right_gate > self.wave_length-1)] = self.wave_length - 1

            self.peak_bin[3, 1, :, :] = left_gate
            self.peak_bin[4, 1, :, :] = right_gate

            self.fast_find_peaks()
        else:
            # update bin_range and call find peaks with new bin_range
            self.bin_range = copy.deepcopy(incoming_gate)
            self.find_peaks()

        # Update gate after everything else in case an error is thrown in
        # find_peaks(). If an error is thrown, gate will not update. This is
        # what we want
        self.gate = copy.deepcopy(incoming_gate)

        t1 = time.time()
        print('Time to change gates: %0.2f seconds' % (t1-t0))

    def fast_find_peaks(self):
        """
        Faster version of find peaks that is called if only the follow gates
        are changed.
        """
        print('In Fast Find Peaks!')
        # TODO: in the future need to add functions for finding the front peak
        #       and the peaks in the follow gate, as separate functions.

        # if PULSE_LENGTH < 0, the maximum and minimum peaks will be found with
        # the follow gate regardless of how far apart they are from each other.
        # If PULSE_LENGTH > 0, the maximum and minimum peaks will only be that
        # far apart from each other and will also be inside of the follow gate

        # if only the follow gates are changed we don't have to worry about
        # finding the FSE within the lead gates because it is already set
        left_gate = self.peak_bin[3, 1, :, :]
        right_gate = self.peak_bin[4, 1, :, :]

        # determine the width of the gate in which to search for the vpp
        # PULSE_LENGTH is the number of half pulses that are to be used in the
        # search
        pulse_width = int(self.PULSE_LENGTH * self.n_half_pulse)

        # IMPORTANT...
        # need to add left_gate or left pulse to all positions found. Numpy
        # returns the index for the position found within waveform that it is
        # passed

        for i in range(self.y_step):
            for j in range(self.x_step):
                max_pos = self.waveform[i, j, left_gate[i, j]:right_gate[i, j]].argmax()
                max_pos += left_gate[i, j]
                if self.PULSE_LENGTH < 0:
                    min_pos = self.waveform[i, j, left_gate[i, j]:right_gate[i, j]].argmin()
                    min_pos += left_gate[i, j]
                else:
                    left_pulse = max_pos - pulse_width
                    right_pulse = max_pos + pulse_width
                    # make sure pulse gates are inside of follow gates
                    if left_pulse < left_gate[i, j]:
                        left_pulse = left_gate[i, j]
                    if right_pulse > right_gate[i, j]:
                        right_pulse = right_gate[i, j]
                    # take argmin() from within pulse_gate
                    min_pos = self.waveform[i, j, left_pulse:right_pulse].argmin()
                    min_pos += left_pulse
                    vpp = self.waveform[i, j, max_pos] - self.waveform[i, j, min_pos]

                    # now get argmin() from within entire follow gate
                    # min_pos2, max_pos2, & vpp2 are from using pulse gate
                    # around argmin instead of argmax()
                    min_pos2 = self.waveform[i, j, left_gate[i, j]:right_gate[i, j]].argmin()
                    min_pos2 += left_gate[i, j]

                    # if min_pos == min_pos2, then we have already found max peak to peak value
                    # by starting search with argmax()
                    if min_pos == min_pos2:
                        pass
                    else:
                        left_pulse = min_pos2 - pulse_width
                        right_pulse = min_pos2 + pulse_width
                        if left_pulse < left_gate[i, j]:
                            left_pulse = left_gate[i, j]
                        if right_pulse > right_gate[i, j]:
                            right_pulse = right_gate[i, j]
                        max_pos2 = self.waveform[i, j, left_pulse:right_pulse].argmax()
                        max_pos2 += left_pulse
                        vpp2 = self.waveform[i, j, max_pos2] - self.waveform[i, j, min_pos2]

                        if vpp2 > vpp:
                            max_pos = max_pos2
                            min_pos = min_pos2

                self.peak_bin[0, 1, i, j] = max_pos
                self.peak_bin[1, 1, i, j] = min_pos
                self.peak_bin[2, 1, i, j] = (max_pos + min_pos) // 2

    def center_image(self):
        """
        Will center the x & y axis arrays over the middle of the sample such
        that the center of the sample is coordinate (0, 0). This method assumes
        that the sample is rectangular
        """

        binary = np.zeros(self.c_scan.shape)

        # otsu threshold is generally good at differentiating between the
        # sample as whole and the background
        threshold = skimage.filters.threshold_otsu(self.c_scan)

        # assume that the peak to peak value of the waveforms on the sample
        # will be larger than on the overscanned areas, so let values larger
        # than threshold be 1
        binary[np.where(self.c_scan > threshold)] = 1

        sample = np.where(binary)
        first_i = sample[0].min()
        last_i = sample[0].max()
        first_j = sample[1].min()
        last_j = sample[1].max()

        center_i = (last_i + first_i) // 2
        center_j = (last_j + first_j) // 2

        self.adjust_coordinates(center_i, center_j)

    def find_peaks(self):
        """
        Put FindPeaks in method so it can be called in another program.
        Initializes peak_bin, which contains information about the location of
        the front surface echo, and if follow gate is on, the location of the
        2nd peak.
        """
        self.peak_bin = FindPeaks(self.waveform, self.x_step, self.y_step, self.wave_length,
                                  self.n_half_pulse, self.FSE_THRESHOLD, self.bin_range,
                                  self.PULSE_LENGTH, self.follow_gate_on)

    def printer(self):
        """
        Prints information about the scan to the console when a THz data class
        in instantiated
        """
        print()
        print(' asn wave length =', self.wave_length, ' asn Time Length =', self.time_length,
              ' delta_t =', self.dt, ' delta_f =', self.df, ' scan type =',
              self.scan_type)

        print('X min =', self.x_min, ' max =', self.x_max, 'Y min =', self.y_min, ' max=',
              self.y_max, ' scan step  X =', self.x_step, ' Y =', self.y_step, ' res X =',
              self.x_res, ' Y =', self.y_res)

        print('True resolution:  X =', self.true_x_res, ' Y =', self.true_y_res)
        print()

    def header_info(self, header):
        """
        Retrieves information about the scan from the header file that is
        created by the DataFile class.
        :param header: the header of the data file class
        """
        # get key parameters from the header
        # change from Python2: add 'b' to the front of all strings, this tells Python that it is
        # looking at a list of bytes

        for item in header:
            b = item.strip()
            a = b.split()

            if a[0] == b'RSDLAMP:':
                self.time_length = float(a[1])
            elif a[0] == b'RESOLUTION:':
                self.x_res = float(a[1])
            elif a[0] == b'Y_RESOLUTION:':
                self.y_res = float(a[1])
            elif a[0] == b'XSTEPS:':
                self.x_step = int(float(a[1]))  # can't directly int(a[1]) of string a[1] ?!
            elif a[0] == b'YSTEPS:':
                self.y_step = int(float(a[1]))
            elif a[0] == b'XMIN:':
                self.x_min = float(a[1])
            elif a[0] == b'XMAX:':
                self.x_max = float(a[1])
            elif a[0] == b'YMIN:':
                self.y_min = float(a[1])
            elif a[0] == b'YMAX:':
                self.y_max = float(a[1])
            elif a[0] == b'SCAN_NAME:':
                self.scan_type = b[10:].strip().decode('utf-8')
            elif a[0] == b'AXIS1:':
                self.axis = b[7:].strip().decode('utf-8')

    def delta_calculator(self):
        """
        Calculates dt, df, n_half_pulse, and the time and frequency array that
        are used for plotting.
        """
        # Note that len(freq) is wave_length/2 + 1. This is so it can be used with
        # numpy's rfft function. Thomas's THzProc code uses len(freq) as wave_length/2.

        self.dt = self.time_length / (self.wave_length - 1)
        self.df = 1. / (self.wave_length * self.dt)
        self.time = np.linspace(0., self.time_length, self.wave_length)
        self.freq = np.linspace(0., (self.wave_length/2) * self.df, self.wave_length//2+1)
        self.n_half_pulse = int(self.HALF_PULSE / self.dt)
        self.true_x_res = (self.x_max - self.x_min) / float(self.x_step - 1)
        self.true_y_res = (self.y_max - self.y_min) / float(self.y_step - 1)
        self.dx = (self.x_max - self.x_min) / float(self.x_step)
        self.dy = (self.y_max - self.y_min) / float(self.y_step)


class DataFile(object):
    # Thomas says this class is from Teraview
    def __init__(self, filename, basedir=None):
        if basedir is not None:  # allows the user to pass through the full filename as 1 argument
            filename = os.path.join(basedir, filename)
        with open(filename, 'rb') as fobj:
            # Edit from THzProc in Python2, need to add 'b' char to beginning of string in Python3
            # I think that this lets the interpreter know it is an instance of bytes instead of a
            # string
            self.header = list(iter(fobj.readline, b"ENDOFHEADER\r\n"))
            # the first 4-byte word after the header gives the length of each row
            first_word = fobj.read(4)
            # convert the 4-byte string to a float
            # change from THzProc in Python2: convert col_size to type int after it has been
            # unpacked as a float. This change is necessary because Python3 no longer accepts
            # floats as array indices eg. array[4.0] is invalid. Used in bin_dtype argument to
            # np.fromfile() below
            # As far as I can tell we have to unpack it as float first, then convert it to int
            col_size = int(struct.unpack(">f", first_word)[0])
            # move the read point back so we can read the first row in its entirety
            fobj.seek(-4, 1)
            # define a compound data type for the data: the two coordinate values
            # followed by the THz waveform
            bin_dtype = [("x", ">f"), ("y", ">f"), ("waveform", ">f", col_size-2)]

            # read the data into an array
            self.data = np.fromfile(fobj, bin_dtype)


# explicitly inherit object so it will be a new-style class even if run using
# python 2
class RefData(object):
    """
    Class representation of the reference data. User has free access to the
    time and waveform information as in THzData class. Important instance
    variables are.
    time = the time array
    waveform = the amplitude array in time domain
    freq = the frequency array
    freq_waveform = the complex amplitude array in the frequency domain
    dt = the spacing between data points in time domain
    df = the spacing between data points in the frequency domain
    n_points = the number of points in the time domain waveform
    """

    def __init__(self, filename, basedir=None, zero=True, gate=[0, None],
                 fix_time=True):
        """
        Constructor method. Everything is done here.
        :param filename: Either the base filename or the full path to the file
                if basedir is left as None
        :param basedir: The path to the directory that contains the file.
                (Default: None) If left as None, filename is expected to
                contain the full path to the file.
        :param zero: If true adjust time values so time[0] is 0. If False, leave
                as the raw values from optical delay (Default: True)
        :param gate: The index slices to remove the front "blip" and the noise
                from the water vapor after the reflection.
        :param fix_time: Whether or not to fix the time domain values not
                having the correct dt. If True (default), the correct final time
                value will be inferred from the filename. If False, nothing will
                be done to the raw data. If passed as a numerical value the last
                time vale will be taken as that value.
        """

        # the array of time values and waveform amplitudes
        self.time, self.waveform = read_reference_data(filename, basedir, zero)

        if type(fix_time) is bool and fix_time:
            # use the filename to infer what the time length is supposed to be
            # reference files are usually named as '50ps waveform.txt'

            # set the max time value to be based on start time in case user
            # does not want to zero time
            if basedir is None:
                basedir, filename = os.path.split(filename)

            # if user includes a directory in the filename string attempt to
            # root it out
            if filename.find('\\') != -1:
                filename = os.path.split(filename)[1]

            # get the integer value of the waveform length, information on time
            # length of the reference signal usually starts with 'ps' or a
            # space.
            max_pos = min(filename.index('p'), filename.index(' '))
            max_time = self.time[0] + int(filename[:max_pos])
            self.time = np.linspace(self.time[0], max_time, len(self.time))

        elif type(fix_time) is int or type(fix_time) is float:
            # allow the user to set the last time value

            max_time = self.time[0] + fix_time
            self.time = np.linspace(self.time[0], max_time, len(self.time))

        # rearrange time array so it starts at zero. When a waveform is saved
        # in the TeraView software it saves the time from some absolute
        # reference point based on mirror position (I think). So the first value
        # in the time array probably won't be zero.
        if zero:
            self.time -= self.time[0]

        # the number of points in the waveform
        self.n_points = len(self.time)

        # calculate dt based on start time, end time, and number of points
        # this should work regardless of whether or not it was zeroed above
        self.dt = (self.time[-1] - self.time[0]) / (self.n_points - 1)

        # the spacing in the frequency domain
        self.df = 1.0 / (self.n_points * self.dt)

        self.freq = np.linspace(0, (self.n_points/2)*self.df,
                                self.n_points//2+1)

        # the frequency amplitude waveform
        # create a deepcopy here so we keep the origin waveform with "blip"
        self.freq_waveform = copy.deepcopy(self.waveform)

        # remove the front blip from the frequency domain waveform if they
        # passed through a value
        if gate[0] != 0:
            self.freq_waveform[:gate[0]] = 0

        # remove the trailing noise if they passed through a rear gate value
        # this cuts down on noise from water vapor
        if gate[1] is not None:
            self.freq_waveform[gate[1]:] = 0

        self.freq_waveform = np.fft.rfft(self.freq_waveform) * self.dt


def read_reference_data(filename, basedir=None, zero=True):
    """
    Reads in a reference file and returns the time values along with amplitudes

    :param filename: The name of the file to be loaded, or the full path if
        basedir is None
    :param basedir: The base directory that contains the file. If None,
        filename is expected to be the full path
    :param zero: If True (default), the reference signal time array is returned
        with zero as the first value. If False, no modification is made.

    :return: optical_delay: Array of time values
    :return: ref_amp: Array of amplitude values
    """
    import pandas as pd

    # allow the user to enter full path in filename parameter
    if basedir is not None:
        filename = os.path.join(basedir, filename)

    # Read in the reference waveform and separate out the optical delay (time)
    # and the reference amplitude
    reference_data = pd.read_csv(filename, delimiter='\t')
    optical_delay = reference_data['Optical Delay/ps'].values
    ref_amp = reference_data['Raw_Data/a.u.'].values

    if zero:  # adjust optical delay so it starts at zero
        optical_delay -= optical_delay[0]

    return optical_delay, ref_amp
