"""
Script for testing THzData class
"""
import wx
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import numpy as np

from THzData import THzData
from FrameHolder import FrameHolder

basedir = 'C:\\Work\\Faint Defect Testing\\Yellow Composite'
filename = 'Scan with Two Tape Defects F@FS (res=0.5mm).tvl'

data = THzData(filename, basedir, follow_gate_on=True)

app = wx.App(False)

holder = FrameHolder(data)

# data.resize(-12, 11, -11.5, 11.2)
#
# calibration_file = 'CalibrationCurve\\CalibrationCurve04DEC2015.txt'
# corrected_amplitude = data.correct_amplitude(calibration_file)
#
# plt.figure('Zoomed C-Scan')
# plt.imshow(data.c_scan_small, interpolation='none', cmap='jet', extent=data.small_extent)
# plt.title('THz Amplitude')
# plt.xlabel('X Scan Location (mm)')
# plt.ylabel('Y Scan Location (mm)')
# plt.colorbar()
# plt.grid()
#
# plt.figure('Corrected C-Scan')
# plt.imshow(corrected_amplitude, interpolation='none', cmap='jet', extent=data.small_extent)
# plt.xlabel('X Scan Location (mm)')
# plt.ylabel('Y Scan Location (mm)')
# plt.colorbar()
# plt.grid()
#
# plt.figure('Time of flight C-Scan')
# plt.imshow(data.tof_c_scan, interpolation='none', cmap='jet_r', extent=data.c_scan_extent,
#            vmax=11)
# plt.title('Time of Flight C-Scan (ps)')
# plt.xlabel('X Scan Location (mm)')
# plt.ylabel('Y Scan Location (mm)')
# plt.colorbar()
# plt.grid()
#
# plt.show()

app.MainLoop()
