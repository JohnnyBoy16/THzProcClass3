"""
Script for testing THzData class
"""
import wx
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import numpy as np

from THzData import THzData
from FrameHolder import FrameHolder

basedir = 'D:\\RR 2015\\John Scans'
filename = 'Sample 1 (res=0.2mm 50ps).tvl'

data = THzData(filename, basedir, follow_gate_on=True)

app = wx.App(False)

holder = FrameHolder(data)

app.MainLoop()
