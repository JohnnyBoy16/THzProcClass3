import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Cursor

from THzData import THzData

basedir = 'C:\\Work\\Signal Modelling\\THz Data\\HDPE Lens'
filename = 'Smooth Side (res=0.5mm, OD=60ps).tvl'
gate = [[100, 1000], [700, 900]]


def c_scan_click_handler(event):
    global xid
    global yid
    xid = int((event.mouseevent.xdata - data.x_min) / data.delta_x)
    yid = int((event.mouseevent.ydata - data.y_min) / data.delta_y)

    # grab the A-Scan from the clicked point and plot it
    a_axis.cla()
    a_axis.set_title('A-Scan Test', fontsize=20)
    a_axis.set_xlabel('Time (ps)')
    a_axis.set_ylabel('Amplitude (a.u.)')
    a_axis.plot(data.time, data.waveform[yid, xid, :], 'r')
    a_axis.axvline(data.time[data.gate[0][0]], color='k', linestyle='--', picker=5)
    a_axis.axvline(data.time[data.gate[0][1]], color='k', linestyle='--', picker=5)
    a_axis.grid()
    a_fig.canvas.draw()

    a_fig.canvas.mpl_connect('button_press_event', grab_gate_handler)
    a_fig.canvas.mpl_connect('button_release_event', release_gate_handler)

    b_axis.cla()
    data.make_b_scan(yid, xid)  # make the B-Scan
    b_axis.set_title('B-Scan Test', fontsize=20)
    b_axis.set_xlabel('X Scan Location (mm)')
    b_axis.set_ylabel('Time (ps)')
    b_axis.imshow(data.b_scan, interpolation='bilinear', cmap='seismic', extent=data.b_scan_extent)
    b_axis.grid()
    b_fig.canvas.draw()


def grab_gate_handler(event):
    print('You clicked', event.xdata, event.ydata)
    index = event.xdata * data.wave_length / data.time_length  # convert to index

    print(index)
    data.gate_index = np.argmin(abs(np.asarray(data.gate) - index))  # this is the gate to replace
    print(data.gate_index)


def release_gate_handler(event):
    print('You released at', event.xdata, event.ydata)

    # convert the (x,y) coordinates from the image to index
    gate[0][data.gate_index] = int(round(event.xdata * data.wave_length / data.time_length, 0))

    data.make_c_scan()  # Generate the new C-Scan data based on gate location

    a_scan_updater()  # update A-Scan with new gate location
    c_scan_updater()  # call the image handler


def a_scan_updater():
    global xid
    global yid

    a_axis.cla()
    a_axis.set_title('A-Scan Test', fontsize=20)
    a_axis.set_xlabel('Time (ps)')
    a_axis.set_ylabel('Amplitude (a.u.)')
    a_axis.plot(data.time, data.waveform[yid, xid, :], 'r')
    a_axis.axvline(data.time[data.gate[0][0]], color='k', linestyle='--', picker=5)
    a_axis.axvline(data.time[data.gate[0][1]], color='k', linestyle='--', picker=5)
    a_axis.grid()
    a_fig.canvas.draw()


def c_scan_updater():
    raw_axis.cla()  # clear the old Raw C-Scan
    interp_axis.cla()  # clear the old Interpolated C-Scan

    plt.figure('Raw C-Scan Test')  # select the raw c-scan figure
    raw_axis.set_title('Raw C-Scan Test', fontsize=20)
    raw_axis.set_xlabel('X Scan Location (mm)')
    raw_axis.set_ylabel('Y Scan Location (mm)')
    raw_image = plt.imshow(data.c_scan, cmap='gray', interpolation='none',
                           extent=data.c_scan_extent, origin='upper', picker=True)
    raw_axis.grid()
    raw_cb.update_bruteforce(raw_image)  # update the colorbar for the new image

    # update the new Raw C-Scan image
    raw_fig.canvas.draw()

    # add the cursor to the raw c-scan image
    # curser seems to disappear after gates have been moved in A-Scan plot
    cursor = Cursor(raw_axis, color='red')

    plt.figure('Interpolated C-Scan Test')  # select the interpolated C-Scan figure
    interp_axis.set_title('Interpolated C-Scan Test', fontsize=20)
    interp_axis.set_xlabel('X Scan Location (mm)')
    interp_axis.set_ylabel('Y Scan Location (mm)')
    interp_image = plt.imshow(data.c_scan, cmap='jet', interpolation='bilinear',
                              extent=data.c_scan_extent, origin='upper')
    interp_axis.grid()
    interp_cb.update_bruteforce(interp_image)  # update the colorbar for the new image
    interp_fig.canvas.draw()


data = THzData(filename, basedir)
print('time length is', data.time_length)
print('wave length is', data.wave_length)

# create and define the Raw C-Scan
raw_fig = plt.figure('Raw C-Scan Test')
raw_axis = raw_fig.add_subplot(111)
# connect the figure to the handler
raw_fig.canvas.mpl_connect('pick_event', c_scan_click_handler)

# create the interpolated C-Scan Figure
interp_fig = plt.figure('Interpolated C-Scan Test')
interp_axis = interp_fig.add_subplot(111)

# create the A-Scan figure
a_fig = plt.figure('A-Scan Test')
a_axis = a_fig.add_subplot(111)

# create the B-Scan figure
b_fig = plt.figure('B-Scan Test')
b_axis = b_fig.add_subplot(111)

plt.figure('Raw C-Scan Test')  # select the raw c-scan figure
raw_axis.set_title('Raw C-Scan Test', fontsize=20)
raw_axis.set_xlabel('X Scan Location (mm)')
raw_axis.set_ylabel('Y Scan Location (mm)')
raw_image = plt.imshow(data.c_scan, cmap='gray', interpolation='none', extent=data.c_scan_extent,
                       origin='upper', picker=True)
raw_axis.grid()
# instantiate a colorbar object for the raw C-Scan, this is necessary so it can be updated
raw_cb = plt.colorbar(raw_image, orientation='horizontal')

# add the cursor to the raw c-scan image
cursor = Cursor(raw_axis, color='red')

plt.figure('Interpolated C-Scan Test')  # select the interpolated C-Scan figure
interp_axis.set_title('Interpolated C-Scan Test', fontsize=20)
interp_axis.set_xlabel('X Scan Location (mm)')
interp_axis.set_ylabel('Y Scan Location (mm)')
interp_image = plt.imshow(data.c_scan, cmap='jet', interpolation='bilinear',
                          extent=data.c_scan_extent, origin='upper')
interp_axis.grid()
# instantiate a colorbar object for interpolated C-Scan, this is necessary so it can be updated
interp_cb = plt.colorbar(interp_image, orientation='horizontal')

plt.show()
