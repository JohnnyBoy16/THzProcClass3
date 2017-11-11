"""
Script for testing THzData class
"""
import wx
import matplotlib.pyplot as plt

from THzData import THzData
from FrameHolder import FrameHolder

basedir = 'C:\\Work\\Faint Defect Testing\\Yellow Composite'
filename = 'Scan with Two Tape Defects F@FS Take 2 (res=0.5mm).tvl'

data = THzData(filename, basedir, follow_gate_on=True)

app = wx.App(False)

holder = FrameHolder(data)

# data.resize(-12, 11, -11.5, 11.2)
#
# plt.figure('Zoomed C-Scan')
# plt.imshow(data.c_scan_zoomed, interpolation='bilinear', cmap='jet', extent=data.zoomed_extent)
# plt.xlabel('X Scan Location (mm)')
# plt.ylabel('Y Scan Location (mm)')
# plt.colorbar()
# plt.grid()
#
# plt.show()

app.MainLoop()
