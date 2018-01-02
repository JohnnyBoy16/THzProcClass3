"""
Script for testing THzData class
"""
import wx

from THzData import THzData
from FrameHolder import FrameHolder

basedir = 'D:\\RR 2017\\THz Data\\Coated Rods'
filename = 'Fully Coated Sample (res=0.5mm OD=40ps).tvl'

data = THzData(filename, basedir, follow_gate_on=True, signal_type=1)

app = wx.App(False)

holder = FrameHolder(data)

app.MainLoop()
