"""
Script for testing THzData class
"""
import wx
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import numpy as np

from THzData import THzData
from FrameHolder import FrameHolder

basedir = ''
filename = ''

data = THzData(filename, basedir, follow_gate_on=True)

app = wx.App(False)

holder = FrameHolder(data)

data.make_time_of_flight_c_scan()

data.resize(-12, 11, -11.5, 11.2)

calibration_file = 'CalibrationCurve\\CalibrationCurve06DEC2017.txt'
corrected_amplitude, thz_height = data.correct_amplitude(calibration_file)

distance, amplitude = np.loadtxt(calibration_file, skiprows=1, unpack=True)

amplitude /= amplitude.max()  # normalize amplitude

# create the interpolation function
interp_f = interp1d(distance, amplitude, kind='quadratic')

x = np.linspace(distance[0], distance[-1], 1000)
y = interp_f(x)

plt.figure('Normalized & Interpolated Calibration Curve')
plt.plot(distance, amplitude, 'bo')
plt.plot(x, y, 'r')
plt.xlabel('Distance From Focus (mm)')
plt.ylabel('Normalized Amplitude')
plt.title('Calibration Curve')
plt.grid()

plt.figure('Zoomed C-Scan')
plt.imshow(data.c_scan_small, interpolation='none', cmap='gray', extent=data.small_extent)
plt.title('THz Amplitude')
plt.xlabel('X Scan Location (mm)')
plt.ylabel('Y Scan Location (mm)')
plt.colorbar()
plt.grid()

plt.figure('Corrected C-Scan')
plt.imshow(corrected_amplitude, interpolation='none', cmap='gray', extent=data.small_extent)
plt.xlabel('X Scan Location (mm)')
plt.ylabel('Y Scan Location (mm)')
plt.title('Amplitude Corrected C-Scan')
plt.colorbar()
plt.grid()

plt.figure('Time of flight C-Scan')
plt.imshow(data.tof_c_scan, interpolation='none', cmap='gray_r', extent=data.c_scan_extent,
           vmax=11)
plt.title('Time of Flight C-Scan (ps)')
plt.xlabel('X Scan Location (mm)')
plt.ylabel('Y Scan Location (mm)')
plt.colorbar()
plt.grid()

# we want to make the amplitude at TOF C-Scan Small be zero
center = np.zeros(2, dtype=int)

center[0] = np.argmin(np.abs(0 - data.y_small))
center[1] = np.argmin(np.abs(0 - data.x_small))

center_value = data.tof_c_scan_small[center[0], center[1]]

data.tof_c_scan_small = data.tof_c_scan_small - center_value

plt.figure('Zoomed TOF C-Scan')
plt.imshow(data.tof_c_scan_small, interpolation='none', cmap='gray_r', extent=data.small_extent)
plt.xlabel('X Scan Location')
plt.ylabel('Y Scan Location')
plt.title('Time of Flight: Normalized (ps)')
plt.colorbar()
plt.grid()

plt.figure('THz Height')
plt.imshow(thz_height*1e3, interpolation='none', cmap='gray_r', extent=data.small_extent)
plt.xlabel('X Scan Location')
plt.ylabel('Y Scan Location')
plt.title(r'THz Height: Normalized ($\mu$ m)')
plt.colorbar()
plt.grid()

plt.show()

app.MainLoop()
