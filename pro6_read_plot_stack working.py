#!/usr/bin/env python
# Slant stack
# Input is set of hinet traces
# traces have already been aligned and corrected for near-vertical statics
#   to have specified phase start at the earthquake origin time
# This programs deals with a single event.
# John Vidale 2/2019

def pro6plotstack(eq_file, plot_scale_fac = 0.05, slow_delta = 0.0005,
			  slowR_lo = -0.1, slowR_hi = 0.1, slowT_lo = -0.1, slowT_hi = 0.1,
			  start_buff = 50, end_buff = 50,
			  plot_dyn_range = 1000, fig_index = 401):

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
	# In English, stack_slows = range(slow_n) * slow_delta - slow_lo
	a1R = range(slowR_n)
	a1T = range(slowT_n)
	stack_Rslows = [(x * slow_delta + slowR_lo) for x in a1R]
	stack_Tslows = [(x * slow_delta + slowT_lo) for x in a1T]
	print(str(slowR_n) + ' radial slownesses, ' + str(slowT_n) + ' trans slownesses, ')

	#%% Find transverse slowness nearest zero
	lowest_Tslow = 1000000
	for slow_i in range(slowT_n):
		if abs(stack_Tslows[slow_i]) < lowest_Tslow:
			lowest_Tindex = slow_i
			lowest_Tslow = abs(stack_Tslows[slow_i])

	print(str(slowT_n) + ' T slownesses, ' + str(lowest_Tindex) + ' min T slow, min is ' + str(lowest_Tslow))

	#%% Select only stacks with that slowness to plot
	central_st = Stream()
	for slowR_i in range(slowR_n):
		central_st += st[slowR_i*slowT_n + lowest_Tindex]

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
	stack_array = np.zeros((slowR_n,stack_nt))

	min_allowed = global_max/plot_dyn_range
	for it in range(stack_nt):  # check points one at a time
		for slow_i in range(slowR_n):  # for this station, loop over slownesses
			num_val = central_st[slow_i].data[it]
			if num_val < min_allowed:
				num_val = min_allowed
			stack_array[slow_i, it] = math.log10(num_val)

	y, x = np.mgrid[slice(0, slowR_n, 1),
				 slice(ttt[0], ttt[-1] + dt, dt)]
#	y, x = np.mgrid[slice(stack_slows[0], stack_slows[-1] + slow_delta, slow_delta),
#				 slice(ttt[0], ttt[-1] + dt, dt)]
#	y, x = np.mgrid[ stack_slows , time ]  # make underlying x-y grid for plot
	plt.close(fig_index)
	plt.figure(fig_index,figsize=(10,10))

	fig, ax = plt.subplots(1)
	c = ax.pcolormesh(x, y, stack_array, cmap=plt.cm.gist_rainbow_r)
	ax.axis([x.min(), x.max(), y.min(), y.max()])
	fig.colorbar(c, ax=ax)
	plt.close(fig_index)
	plt.xlabel('Time (s)')
	plt.ylabel('Slowness (s/km)')
	plt.title(fname[2:12])
	plt.show()

	#  Save processed files
#	fname = 'HD' + date_label + '_slice.mseed'
#	stack.write(fname,format = 'MSEED')

	elapsed_time_wc = time.time() - start_time_wc
	print('This job took ' + str(elapsed_time_wc) + ' seconds')
	os.system('say "Done"')