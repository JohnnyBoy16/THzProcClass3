"""
Script for testing THzData class
"""
import wx

from THzData import THzData
from FrameHolder import FrameHolder

basedir = 'F:\\RR 2017\\THz Data\\Black CMC'
filename = 'Sample 1 (res=0.25mm OD=60ps).tvl'

gate = [[100, 1000], [2250, 2550]]

data = THzData(filename, basedir, follow_gate_on=True, signal_type=1, gate=gate)

app = wx.App(False)

holder = FrameHolder(data)

app.MainLoop()
