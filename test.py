"""
Script for testing THzData class
"""
import wx
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import numpy as np

from THzData import THzData
from FrameHolder import FrameHolder

basedir = 'F:\\RR 2016\\THz Data\\4-1'
filename = '4-1 10-19-2016 (100ps res=0.25mm).tvl'

data = THzData(filename, basedir, follow_gate_on=True)

app = wx.App(False)

holder = FrameHolder(data)

app.MainLoop()
