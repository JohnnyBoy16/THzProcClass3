import numpy as np
import os
import copy
import struct
import pdb
from thz_functions import ReMap, AmpCor300, FindPeaks

# THINGS THAT STILL HAVE TO BE IMPLEMENTED
# Depth Map
# trend off != 0


class THzData:
    # user can provide the full path the file in filename, or provide the filename and base
    # directory separately
    def __init__(self, filename, basedir=None, gate=None, follow_gate_on=False, signal_type=0,
                 center=False, print_on=True):

        dat = DataFile(filename, basedir=basedir)  # first thing to do is open the file

        # handle what happens if they don't pass through gate
        if gate is None:
            self.gate = [[100, 1500], [700, 900]]

        # if the user passed through a value for gate, make sure that it is 2x2
        elif np.shape(gate) != (2, 2):
            raise ValueError('gate must by a 2x2 list or numpy array!')
        else:
            self.gate = gate

        # make sure that follow_gate_on is either True or False
        if not isinstance(follow_gate_on, bool):
            raise ValueError('follow_gate_on must be of a boolean!')
        else:
            self.follow_gate_on = follow_gate_on

        # set signal_type to pass to make_c_scan()
        self.signal_type = signal_type

        self.x = dat.data[3:]['x']
        self.y = dat.data[3:]['y']
        self.waveform = dat.data['waveform'][3:]
        self.wave_length = len(self.waveform[0])

        self.a_scan_only = True
        self.b_scan_on = True
        self.b_scan_dir = 'horizontal'  # B Scan direction is usually horizontal by default

        # whether or not to correct for the excessive amplification on the edges of the
        # 300 ps waveforms
        self.amp_correction300_on = True

        # initialize all of the variables that are wanted from the header, these are calculated in
        # the method header_info
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

        self.delta_t = None  # the spacing between time values
        self.time = None  # array of time values that is used for plotting
        self.freq = None  # array of frequency values that is used for ploting
        self.delta_f = None  # spacing between frequency values
        self.n_half_pulse = None  # the number of half pulses in the scan
        self.true_x_res = None  # the spacing between x points after Remap is called
        self.true_y_res = None  # the spacing between y points after Remap is called
        self.delta_x = None  # the average difference between x data points
        self.delta_y = None  # the average difference between y data points
        self.b_scan = None  # the B-Scan values
        self.tof_c_scan = None  # the time of flight of the front surface

        # the largest frequency that you would like the frequency plots to go up to
        self.work_freq = 3.0

        # flags peaks (the FSE always and peaks in follow gate if follow gate is on)
        self.flag_peak_on = True

        # CONSTANTS --------------------------------------------------------------------------------

        # 1 for the positive peak, 2 for the negative peak, prefer 1
        self.PICK_PEAK = 1

        # provides algorithm for removing baseline trend in near-field imaging. leave 0 for now
        self.TREND_OFF = 0

        # VERY constant: do not change unless you have a good reason to
        self.FSE_TOLERANCE = -0.17

        # if true, removes excess amplitude on the ends of each waveform for a 300ps scan
        self.AMP_CORRECTION_300 = True
        # parameters to pass to the AmpCor300 function
        self.AMP_CORRECTION_300_PAR = [0., 35., 5.0, 1., 240., 300., 1., 4.5, 4.]

        self.X_CORRECTION_TOLERANCE = 0.1
        self.PULSE_LENGTH = 3.0

        # difference that a point is allowed to deviate (as a ratio with respect to resolution)
        # from its desired coordinate
        self.X_DIFFERENCE_TOLERANCE = 0.4
        self.Y_DIFFERENCE_TOLERANCE = 0.3

        # tolerance ratio for number of scan points in a X line
        # this is nlim from base THzProc
        self.N_LIMIT = 0.8

        # half pulse width for bracketing the follower signal (usually the front surface echo)
        self.HALF_PULSE = 2  # don't change unless you have a good reason to

        # threshold for lead gate signal
        self.FSE_THRESHOLD = 0.5

        # tolerance ratio for number of actual scan points in a X line compared to how many are
        # supposed to be in that line
        self.X_RATIO_LIMIT = 0.8

        # the height of the figure in inches
        self.FIGURE_HEIGHT = 8

        # a small value that is used to see if floating points are close together.
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
            AmpCor300(self.AMP_CORRECTION_300_PAR, self.wave_length, self.delta_t, self.x_step,
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

        self.c_scan_extent = (self.x_min, self.x_max, self.y_max, self.y_min)
        self.b_scan_extent = (self.x_min, self.x_max, 0, self.time_length)

        # Generate the C-Scan
        self.make_c_scan(self.signal_type)

    def make_c_scan(self, signal_type):
        """
        Generates the C-Scan based on the gate location and signal type
        At the moment this just does peak to peak voltage response
        :param: signal_type - determines how the C-Scan is calculated
                    choices:
                        0: (default) Use Peak to Peak voltage with the front gates regardless of
                           whether follow gate is on or not.
                        1: Use Peak to Peak voltage with the follow gates if on. If follow gate
                           is not on then use peak to peak voltage across entire waveform
                        2: avg mag (mag = abs amp)
                        3: median mag
                        4: min mag
                        5: max mag
                        6: avg amp ~ integrated gate
                        7: median amp
                        8: min amp
                        9: max amp
        """
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
            if self.follow_gate_on:
                idx = 1
            else:
                idx = 0
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

        # 27JAN2015 FSE handle: if FSE too small, throw away this location
        if signal_type != 0 and self.FSE_TOLERANCE > 0:
            for i in range(self.y_step):
                for j in range(self.x_step):
                    if (self.waveform[i, j, self.peak_bin[0, 0, i, j]] -
                            self.waveform[i, j, self.peak_bin[1, 0, i, j]] < self.FSE_TOLERANCE):
                        self.c_scan[i, j] = 0.

    def make_time_of_flight_c_scan(self):
        """
        Generates a time of flight C-Scan that shows the time of flight in picoseconds of each
        pixel on the front surface.
        """
        tof_index = self.waveform.argmax(axis=2)
        self.tof_c_scan = np.zeros((self.y_step, self.x_step))

        for i in range(self.y_step):
            for j in range(self.x_step):
                self.tof_c_scan[i, j] = self.time[tof_index[i, j]]

    def correct_amplitude(self, filename, basedir=None, focus=None):
        """
        Attempts to correct the amplitude for the time of flight in accordance with a gaussian
        beam profile
        :param
        """
        from scipy.interpolate import interp1d
        from scipy.optimize import curve_fit
        import util

        c = 0.2998  # speed of light in mm/ps
        theta0 = 17.5 * np.pi / 180  # incoming angle of THz beam

        # need to the time of flight info to make adjustments, so if tof_c_scan has not been made
        # made yet, make it
        if self.tof_c_scan is None:
            self.make_time_of_flight_c_scan()

        if basedir is not None:
            filename = os.path.join(basedir, filename)

        distance, amplitude = np.loadtxt(filename, skiprows=1, unpack=True)

        amplitude /= amplitude.max()  # normalize amplitude

        # only interested in the values close to maximum as ceramic is pretty flat
        # should give a better fit
        max_index = amplitude.argmax()
        amplitude = amplitude[max_index-5:max_index+6]
        distance = distance[max_index-5:max_index+6]

        # initial parameter guesses for curve_fit to Gaussian
        a = 1  # amplitude
        b = 0  # center
        c = 0.25  # width
        p0 = (a, b, c)

        # only want best fit parameters back
        p = curve_fit(util.guassian_curve, distance, amplitude, p0)[0]
        print(p)
        # p[1] = 0  # set b to be 0, center the curve

        x = np.linspace(-0.25, 0.25, 100)
        y = util.guassian_curve(x, *p)

        import matplotlib.pyplot as plt
        plt.figure('Gaussian Curve Fit Test')
        plt.plot(x, y, 'b', label='Best Fit')
        plt.plot(distance, amplitude, 'ro', label='True Points')
        plt.xlabel('Distance from Focus (mm)')
        plt.ylabel('Normalized Amplitude')
        plt.grid()

        if focus is None:
            # assume that the pixel at (0, 0) is the focus
            # all other TOF values will be in reference to this
            j_idx = np.argmin(np.abs(self.x - 0))
            i_idx = np.argmin(np.abs(self.y - 0))
        else:
            raise ValueError('Focus has to be (0,0) for now!')

        tof_temp = copy.deepcopy(self.tof_c_scan)

        tof_temp -= tof_temp[i_idx, j_idx]  # normalize TOF to focus point

        i0 = np.where(self.y == self.y_small[0])[0][0]
        i1 = np.where(self.y == self.y_small[-1])[0][0]
        j0 = np.where(self.x == self.x_small[0])[0][0]
        j1 = np.where(self.x == self.x_small[-1])[0][0]

        # relative height of the sample
        # NOTE!!!!!!!!: the height array is actually the negative of what it should be.
        # however we want this because the THz gantry considers -z to be up
        height = tof_temp[i0:i1+1, j0:j1+1] * c * np.cos(theta0) / 2
        print('Max height = %0.4f' % height.min())  # therefore min() is actually highest point
        print('Min Height = %0.4f' % height.max())  # and max() is lowest

        amplitude_correction = util.guassian_curve(height, *p)

        plt.figure('Amplitude Correction for Sample')
        plt.imshow(amplitude_correction, interpolation='none', cmap='gray',
                   extent=self.small_extent)
        plt.xlabel('X Scan Location (mm)')
        plt.ylabel('Y Scan Location (mm)')
        plt.title('Amplitude Correction Factor')
        plt.colorbar()
        plt.grid()

        c_scan_corrected = self.c_scan_small / amplitude_correction

        return c_scan_corrected, height

    def resize(self, x0, x1, y0, y1):
        """
        Resizes the data in the bounds between x0, x1, y0, and y1. Should be used to remove the
        edges from the data if it was over scanned
        :param x0: The smallest x value in the new image
        :param x1: The largest x value in the new image
        :param y0: The smallest y value in the new image
        :param y1: The largest y value in the new image
        """

        j0 = np.argmin(np.abs(self.x - x0))
        j1 = np.argmin(np.abs(self.x - x1))
        self.x_small = self.x[j0:j1]
        self.x_step_small = len(self.x_small)

        i0 = np.argmin(np.abs(self.y - y0))
        i1 = np.argmin(np.abs(self.y - y1))
        self.y_small = self.y[i0:i1]
        self.y_step_small = len(self.y_small)

        self.waveform_small = self.waveform[i0:i1, j0:j1, :]

        max_amp = np.amax(self.waveform_small[:, :, self.gate[0][0]:self.gate[0][1]], axis=2)
        min_amp = np.amin(self.waveform_small[:, :, self.gate[0][0]:self.gate[0][1]], axis=2)
        self.c_scan_small = max_amp - min_amp

        self.small_extent = (self.x_small[0], self.x_small[-1], self.y_small[0],
                             self.y_small[-1])

        if self.tof_c_scan is not None:
            self.tof_c_scan_small = self.tof_c_scan[i0:i1, j0:j1]

    def center_coordinates(self):
        """
        Centers the x and y coordinates such that (0, 0) is in the center of the image
        """
        self.x -= ((self.x[0]-self.x[-1]) / 2)
        self.y -= ((self.y[0]-self.y[-1]) / 2)

    def make_b_scan(self, yid, xid):
        """
        Generates the B-Scan based on the last cursor location clicked and whether b_scan_dir is
        horizontal or vertical
        :param yid: the row from which to generate B-Scan
        :param xid: the column clicked from which to generate B-Scan
        """
        # in order to vectorize this method, the x or y coordinates are put in the rows, but to plot
        # we want them in the columns. Thus the call to transpose after arranging the data
        if self.b_scan_dir == 'horizontal':
            self.b_scan = self.waveform[yid, :, :]
        else:  # b_scan_dir == 'vertical'
            self.b_scan = self.waveform[:, xid, :]

        # call to transpose is necessary to flip data axis, so x or y location is on the bottom and
        # time is along the y-axis on the plot
        self.b_scan = np.transpose(self.b_scan)
        self.b_scan = np.flipud(self.b_scan)  # flip so top of sample is at bottom of image

    def set_follow_gate(self, given_boolean):
        """
        Sets the follow gate to be either on or off. If follow gate is changed, code will update
        bin_range and peak_bin.
        :param given_boolean: True: follow gate is on. False: follow gate is off.
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
        Change the gate locations by providing an a new gate. Then updates bin_range and peak_bin
        accordingly
        :param incoming_gate: The value of gate that is to be changed to, must be a 2x2 list or
                    numpy array
        """

        # incoming_gate is the same as what gate currently is; do nothing
        if np.array_equal(self.gate, incoming_gate):
            return
        # make sure that incoming_gate is a 2x2 array
        elif np.shape(np.asarray(incoming_gate)) != (2, 2):
            raise ValueError('gate must by a 2x2 list or numpy array!')

        # update bin_range and call find peaks with new bin_range
        self.bin_range = copy.deepcopy(incoming_gate)
        self.find_peaks()

        # Update gate after everything else in case an error is thrown in find_peaks().
        # If an error is thrown, gate will not update. This is what we want
        self.gate = copy.deepcopy(incoming_gate)

    def set_b_scan_on(self, given_boolean):
        pass
        # TODO: write code so setting B-Scan on when it was previously off,
        # TODO: calls make B-Scan from the last point clicked

    def set_b_scan_direction(self, given_direction):
        pass
        # TODO: write code so changing B-Scan direction will call make B-Scan from last point clicked

    def find_peaks(self):
        """
        Put FindPeaks in method so it can be called in another program. Initializes peak_bin,
        which contains information about the location of the front surface echo, and if follow
        gate is on, the location of the 2nd peak.
        """
        self.peak_bin = FindPeaks(self.waveform, self.x_step, self.y_step, self.wave_length,
                                  self.n_half_pulse, self.FSE_THRESHOLD, self.bin_range,
                                  self.PULSE_LENGTH, self.follow_gate_on)

    def printer(self):
        """
        Prints information about the scan to the console when a THz data class in instantiated
        """
        print()
        print(' asn wave length =', self.wave_length, ' asn Time Length =', self.time_length,
              ' delta_t =', self.delta_t, ' delta_f =', self.delta_f, ' scan type =',
              self.scan_type)

        print('X min =', self.x_min, ' max =', self.x_max, 'Y min =', self.y_min, ' max=',
              self.y_max, ' scan step  X =', self.x_step, ' Y =', self.y_step, ' res X =',
              self.x_res, ' Y =', self.y_res)

        print('True resolution:  X =', self.true_x_res, ' Y =', self.true_y_res)
        print()

    def header_info(self, header):
        """
        Retrieves information about the scan from the header file that is created by the DataFile
        class.
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
        Calculates delta_t, delta_f, n_half_pulse, and the time and frequency array that are used
        for plotting.
        """
        # Note that len(freq) is wave_length/2 + 1. This is so it can be used with
        # numpy's rfft function. Thomas's THzProc code uses len(freq) as wave_length/2.

        self.delta_t = self.time_length / (self.wave_length - 1)
        self.delta_f = 1. / (self.wave_length * self.delta_t)
        self.time = np.linspace(0., self.time_length, self.wave_length)
        self.freq = np.linspace(0., (self.wave_length/2) * self.delta_f, self.wave_length//2+1)
        self.n_half_pulse = int(self.HALF_PULSE / self.delta_t)
        self.true_x_res = (self.x_max - self.x_min) / float(self.x_step - 1)
        self.true_y_res = (self.y_max - self.y_min) / float(self.y_step - 1)
        self.delta_x = (self.x_max - self.x_min) / float(self.x_step)
        self.delta_y = (self.y_max - self.y_min) / float(self.y_step)


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
