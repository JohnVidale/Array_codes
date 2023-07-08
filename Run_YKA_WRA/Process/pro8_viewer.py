#!/usr/bin/env python3
# makes a little animation, controlled by mouse, of snapshots of tdiff as a function of slownesses
import os
os.environ['PATH'] += os.pathsep + '/usr/local/bin'
import numpy as np
import matplotlib.pyplot as plt
from obspy import UTCDateTime
from obspy import Stream, Trace
from obspy import read
import time

start_time_wc = time.time()

#%%
class IndexTracker(object):
    def __init__(self, ax, X):
        self.ax = ax
        ax.set_title('use scroll wheel to navigate images')

        self.X = X
        rows, cols, self.slices = X.shape
        self.ind = self.slices//2

        self.im = ax.imshow(self.X[:, :, self.ind])
        self.update()

    def onscroll(self, event):
        print("%s %s" % (event.button, event.step))
        if event.button == 'up':
            self.ind = (self.ind + 1) % self.slices
        else:
            self.ind = (self.ind - 1) % self.slices
        self.update()

    def update(self):
        self.im.set_data(self.X[:, :, self.ind])
        self.ax.set_ylabel('slice %s' % self.ind)
        self.im.axes.figure.canvas.draw()

#%% define common parameters
ARRAY = 1
eq_file1 = 'event1.txt'
eq_file2 = 'event2.txt'
start_buff = 20
end_buff   = 200

slowR_lo   = -0.05
slowR_hi   =  0.05
slowT_lo   = -0.05
slowT_hi   =  0.05
slow_delta =  0.005
decimate_fac   =  5

plot_scale_fac = 0.003
fig_index = 301
plot_dyn_range = 100

os.chdir('/Users/vidale/Documents/PyCode/LASA')

#%% Input parameters
# Get saved event info, also used to name files
# date_label = '2018-04-02' # date for filename
file = open('EvLocs/' + eq_file1, 'r')
lines=file.readlines()
split_line = lines[0].split()
t1           = UTCDateTime(split_line[1])
date_label1  = split_line[1][0:10]

file = open('EvLocs/' + eq_file2, 'r')
lines=file.readlines()
split_line = lines[0].split()
t2           = UTCDateTime(split_line[1])
date_label2  = split_line[1][0:10]

fname1  = 'HD' + date_label1 + '_' + date_label2 + '_tshift.mseed'
tshift  = Stream()
tshift  = read(fname1)
fname1  = 'HD' + date_label1 + '_' + date_label2 + '_amp_ratio.mseed'
amp_ratio  = Stream()
amp_ratio  = read(fname1)
fname2  = 'HD' + date_label1 + '_' + date_label2 + '_amp_ave.mseed'
amp_ave = Stream()
amp_ave = read(fname2)

tshift_full = tshift.copy()
tshift.decimate(decimate_fac, no_filter=True)
amp_ratio.decimate(decimate_fac, no_filter=True)
amp_ave.decimate(decimate_fac, no_filter=True)

#print(f'len(tshift): ' {len(tshift):4d} 'len(tshift[0].data): ' {len(tshift[0].data:4d})
print(f'len(tshift): {len(tshift):4d}')
print(f'len(amp_ave): {len(amp_ave):4d}')
print(f'len(amp_ratio): {len(amp_ratio):4d}')

nt = len(tshift[0].data)
dt = tshift[0].stats.delta
print(f'nt: {nt:6d}')

#%% Make grid of slownesses
slowR_n = int(1 + (slowR_hi - slowR_lo)/slow_delta)  # number of slownesses
slowT_n = int(1 + (slowT_hi - slowT_lo)/slow_delta)  # number of slownesses
# In English, stack_slows = range(slow_n) * slow_delta - slow_lo
a1R = range(slowR_n)
a1T = range(slowT_n)
stack_Rslows = [(x * slow_delta + slowR_lo) for x in a1R]
stack_Tslows = [(x * slow_delta + slowT_lo) for x in a1T]
print(str(slowT_n) + ' trans slownesses,  hi and lo are ' + str(slowT_hi) + '  ' + str(slowT_lo))
print(str(slowR_n) + ' radial slownesses, hi and lo are ' + str(slowR_hi) + '  ' + str(slowR_lo))

total_slows = slowR_n * slowT_n
global_max = 0
for slow_i in range(total_slows): # find envelope, phase, tshift, and global max
	local_max = max(abs(amp_ave[slow_i].data))
	if local_max > global_max:
		global_max = local_max

for slow_i in range(total_slows): # ignore less robust points
	for it in range(nt):
		if ((amp_ratio[slow_i].data[it] < 0.7) or (amp_ratio[slow_i].data[it] > 1.5) or (amp_ave[slow_i].data[it] < (0.1 * global_max))):
			tshift[slow_i].data[it] = np.nan

X = np.random.rand(slowT_n, slowR_n, nt)

for slowT_i in range(slowT_n):
	for slowR_i in range(slowR_n):
		for it in range(nt):
			index = slowR_i*slowT_n + slowT_i
			X[slowR_i,slowT_i,it] = tshift[index].data[it]

for slowT_i in range(slowT_n):  # Set scale in upper left corner
	for slowR_i in range(slowR_n):
		for it in range(nt):
			X[0,0,it] = 1
			X[0,1,it] = -1
#%%
#	total_slows = slowR_n * slowT_n

plt.title('T-R stack ' + fname1[2:12] + ' ' + fname2[2:12])

fig, ax = plt.subplots(1, 1)

plt.xlabel('T Slowness (s/km)')
plt.ylabel('R Slowness (s/km)')
#plt.xlim(stack_Rslows[0],stack_Rslows[-1])
#plt.ylim(stack_Tslows[0],stack_Tslows[-1])

tracker = IndexTracker(ax, X)
fig.canvas.mpl_connect('scroll_event', tracker.onscroll)
plt.show()

elapsed_time_wc = time.time() - start_time_wc
print('This job took ' + str(elapsed_time_wc) + ' seconds')
os.system('say "Done"')