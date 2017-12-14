"""
Script for testing THzData class
"""
import wx
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import numpy as np

from THzData import THzData
from FrameHolder import FrameHolder

basedir = 'F:\\RR 2016\\THz Data\\Grinding Trial Sample\\1st Grind'
filename = 'Sample 4-1 After 1st Polish res=0.25mm.tvl'

data = THzData(filename, basedir, follow_gate_on=True)

app = wx.App(False)

holder = FrameHolder(data)

app.MainLoop()
