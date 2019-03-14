#!/usr/bin/env python
# Slant stack
# Input is set of hinet traces
# traces have already been aligned and corrected for near-vertical statics
#   to have specified phase start at the earthquake origin time
# This programs deals with a single event.
# John Vidale 2/2019

def pro6plotstack(eq_file, plot_scale_fac = 0.05, slow_delta = 0.0005,
			  slowR_lo = -0.1, slowR_hi = 0.1, slowT_lo = -0.1, slowT_hi = 0.1,
			  start_buff = 50, end_buff = 50, snaptime = 0, snaps = 1,
			  plot_dyn_range = 1000, fig_index = 401, skip_T = 1, skip_R = 0):

	import obspy
	import obspy.signal
	from obspy import UTCDateTime
	from obspy import Stream, Trace
	from obspy import read
	from obspy.geodetics import gps2dist_azimuth
	import numpy as np
	import os
	from obspy.taup import TauPyModel
	import obspy.signal as sign
	import matplotlib.pyplot as plt
	model = TauPyModel(model='iasp91')
	from scipy.signal import hilbert
	import math
	import time

	start_time_wc = time.time()

	file = open(eq_file, 'r')
	lines=file.readlines()
	split_line = lines[0].split()
#			ids.append(split_line[0])  ignore label for now
	t           = UTCDateTime(split_line[1])
	date_label  = split_line[1][0:10]

	#%% Input parameters
	# #%% Get saved event info, also used to name files
	# date_label = '2018-04-02' # date for filename
	fname = 'HD' + date_label + '_2dstack.mseed'
	st = Stream()
	st = read(fname)
	print('Read in: ' + str(len(st)) + ' traces')
	nt = len(st[0].data)
	dt = st[0].stats.delta
	print('First trace has : ' + str(nt) + ' time pts, time sampling of '
		  + str(dt) + ' and thus duration of ' + str((nt-1)*dt))

	#%% Make grid of slownesses
	slowR_n = int(1 + (slowR_hi - slowR_lo)/slow_delta)  # number of slownesses
	slowT_n = int(1 + (slowT_hi - slowT_lo)/slow_delta)  # number of slownesses
	stack_nt = int(1 + ((end_buff + start_buff)/dt))  # number of time points
	print(str(slowT_n) + ' trans slownesses, hi and lo are ' + str(slowT_hi) + '  ' + str(slowT_lo))
	# In English, stack_slows = range(slow_n) * slow_delta - slow_lo
	a1R = range(slowR_n)
	a1T = range(slowT_n)
	stack_Rslows = [(x * slow_delta + slowR_lo) for x in a1R]
	stack_Tslows = [(x * slow_delta + slowT_lo) for x in a1T]
	print(str(slowR_n) + ' radial slownesses, ' + str(slowT_n) + ' trans slownesses, ')

	#%% Find transverse slowness nearest zero
	if skip_R != 1:
		lowest_Tslow = 1000000
		for slow_i in range(slowT_n):
			if abs(stack_Tslows[slow_i]) < lowest_Tslow:
				lowest_Tindex = slow_i
				lowest_Tslow = abs(stack_Tslows[slow_i])

		print(str(slowT_n) + ' T slownesses, ' + str(lowest_Tindex) + ' min T slow, min is ' + str(lowest_Tslow))

		# Select only stacks with that slowness for Radial plot
		centralR_st = Stream()
		for slowR_i in range(slowR_n):
			centralR_st += st[slowR_i*slowT_n + lowest_Tindex]

	#%% If desired, find radial slowness nearest zero
	if skip_T != 1:
		lowest_Rslow = 1000000
		for slow_i in range(slowR_n):
			if abs(stack_Rslows[slow_i]) < lowest_Rslow:
				lowest_Rindex = slow_i
				lowest_Rslow = abs(stack_Rslows[slow_i])

		print(str(slowR_n) + ' R slownesses, ' + str(lowest_Rindex) + ' min R slow, min is ' + str(lowest_Rslow))

		# Select only stacks with that slowness for Radial plot
		centralT_st = Stream()
		for slowT_i in range(slowT_n):
			centralT_st += st[lowest_Rindex*slowT_n + slowT_i]
#%%
	total_slows = slowR_n * slowT_n
	global_max = 0
	for slow_i in range(total_slows): # find global max, and if requested, take envelope
		if len(st[slow_i].data) == 0:
				print('%d data has zero length ' % (slow_i))
		local_max = max(abs(st[slow_i].data))
		if local_max > global_max:
			global_max = local_max

	#%% compute timing time series
	ttt = (np.arange(len(st[0].data)) * st[0].stats.delta - start_buff) # in units of seconds

#%%  Plotting
# Regular radial-time stack
	if skip_R != 1:
		stack_array = np.zeros((slowR_n,stack_nt))

		min_allowed = global_max/plot_dyn_range
		for it in range(stack_nt):  # check points one at a time
			for slowR_i in range(slowR_n):  # for this station, loop over slownesses
				num_val = centralR_st[slowR_i].data[it]
				if num_val < min_allowed:
					num_val = min_allowed
				stack_array[slowR_i, it] = math.log10(num_val)
		stack_array[0,0] = math.log10(global_max)
		stack_array[0,1] = math.log10(min_allowed)

		y, x = np.mgrid[slice(stack_Rslows[0], stack_Rslows[-1] + slow_delta, slow_delta),
					 slice(ttt[0], ttt[-1] + dt, dt)]

		fig, ax = plt.subplots(1)
		print('Figure is set to ' + str(fig))
		c = ax.pcolormesh(x, y, stack_array, cmap=plt.cm.gist_rainbow_r)
		ax.axis([x.min(), x.max(), y.min(), y.max()])
		fig.colorbar(c, ax=ax)
		plt.xlabel('Time (s)')
		plt.ylabel('Slowness (s/km)')
		plt.title('Radial stack at 0 T slow, ' + fname[2:12])
		plt.show()

#%%  Transverse-time stack
	if skip_T != 1:
		stack_array = np.zeros((slowT_n,stack_nt))

		min_allowed = global_max/plot_dyn_range
		for it in range(stack_nt):  # check points one at a time
			for slowT_i in range(slowT_n):  # for this station, loop over slownesses
				num_val = centralT_st[slowT_i].data[it]
				if num_val < min_allowed:
					num_val = min_allowed
				stack_array[slowT_i, it] = math.log10(num_val)
		stack_array[0,0] = math.log10(global_max)
		stack_array[0,1] = math.log10(min_allowed)

		y, x = np.mgrid[slice(stack_Tslows[0], stack_Tslows[-1] + slow_delta, slow_delta),
					 slice(ttt[0], ttt[-1] + dt, dt)]

		fig, ax = plt.subplots(1)
		print('Figure is set to ' + str(fig))
		c = ax.pcolormesh(x, y, stack_array, cmap=plt.cm.gist_rainbow_r)
		ax.axis([x.min(), x.max(), y.min(), y.max()])
		fig.colorbar(c, ax=ax)
		plt.xlabel('Time (s)')
		plt.ylabel('Slowness (s/km)')
		plt.title('Transverse stack at 0 R slow, ' + fname[2:12])
		plt.show()

#%% R-T stack
	stack_array1 = np.zeros((slowR_n,slowT_n))
	for snap_num in range(snaps):
		fig_index += 1
		it = int((snaptime + start_buff)/dt) + snap_num
		for slowR_i in range(slowR_n):  # loop over radial slownesses
			for slowT_i in range(slowT_n):  # loop over transverse slownesses
				index = slowR_i*slowT_n + slowT_i
				num_val = st[index].data[it]
				if num_val < min_allowed:
					num_val = min_allowed
				stack_array1[slowR_i, slowT_i] = math.log10(num_val)
		stack_array1[0,0] = math.log10(global_max)
		stack_array1[0,1] = math.log10(min_allowed)

		y1, x1 = np.mgrid[slice(stack_Rslows[0], stack_Rslows[-1] + slow_delta, slow_delta),
					 slice(stack_Tslows[0], stack_Tslows[-1] + slow_delta, slow_delta)]

		fig, ax = plt.subplots(1)
		c = ax.pcolormesh(x1, y1, stack_array1, cmap=plt.cm.gist_rainbow_r)
		ax.axis([x1.min(), x1.max(), y1.min(), y1.max()])
		fig.colorbar(c, ax=ax)
		plt.xlabel('T Slowness (s/km)')
		plt.ylabel('R Slowness (s/km)')
		plt.title('T-R stack at rel time ' + str(snaptime + snap_num*dt) + '  ' + fname[2:12])
		plt.show()

	#  Save processed files
#	fname = 'HD' + date_label + '_slice.mseed'
#	stack.write(fname,format = 'MSEED')

	elapsed_time_wc = time.time() - start_time_wc
	print('This job took ' + str(elapsed_time_wc) + ' seconds')
	os.system('say "Done"')