#!/usr/bin/env python
# Slant stack
# Input is set of hinet traces
# traces have already been aligned and corrected for near-vertical statics
#   to have specified phase start at the earthquake origin time
# This programs deals with a single event.
# John Vidale 2/2019

def pro7stacked_seis(eq_file1, eq_file2, plot_scale_fac = 0.05, slow_delta = 0.0005,
			  slowR_lo = -0.1, slowR_hi = 0.1, slowT_lo = -0.1, slowT_hi = 0.1,
			  start_buff = 50, end_buff = 50, snaptime = 0, snaps = 1, norm = 0,
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

	file = open(eq_file1, 'r')
	lines=file.readlines()
	split_line = lines[0].split()
	t1           = UTCDateTime(split_line[1])
	date_label1  = split_line[1][0:10]

	file = open(eq_file2, 'r')
	lines=file.readlines()
	split_line = lines[0].split()
	t2           = UTCDateTime(split_line[1])
	date_label2  = split_line[1][0:10]

	#%% Input parameters
	# #%% Get saved event info, also used to name files
	# date_label = '2018-04-02' # date for filename
	fname1 = 'HD' + date_label1 + '_2dstack.mseed'
	fname2 = 'HD' + date_label2 + '_2dstack.mseed'
	st1 = Stream()
	st2 = Stream()
	st1 = read(fname1)
	st2 = read(fname2)
	print('Read in: event 1 ' + str(len(st1)) + ' event 2 ' + str(len(st2)) + ' traces')
	nt1 = len(st1[0].data)
	nt2 = len(st2[0].data)
	dt1 = st1[0].stats.delta
	dt2 = st2[0].stats.delta
	print('Event 1 - First trace has ' + str(nt1) + ' time pts, time sampling of '
		  + str(dt1) + ' and thus duration of ' + str((nt1-1)*dt1))
	print('Event 2 - First trace has ' + str(nt2) + ' time pts, time sampling of '
		  + str(dt2) + ' and thus duration of ' + str((nt2-1)*dt2))
	if nt1 != nt2 or dt1 != dt2:
		print('nt or dt not does not match')
		exit(-1)

	#%% Make grid of slownesses
	slowR_n = int(1 + (slowR_hi - slowR_lo)/slow_delta)  # number of slownesses
	slowT_n = int(1 + (slowT_hi - slowT_lo)/slow_delta)  # number of slownesses
	stack_nt = int(1 + ((end_buff + start_buff)/dt1))  # number of time points
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

		# Select only stacks with that slowness for radial plot
		centralR_st1 = Stream()
		centralR_st2 = Stream()
		for slowR_i in range(slowR_n):
			centralR_st1 += st1[slowR_i*slowT_n + lowest_Tindex]
			centralR_st2 += st2[slowR_i*slowT_n + lowest_Tindex]

	#%% If desired, find radial slowness nearest zero
	if skip_T != 1:
		lowest_Rslow = 1000000
		for slow_i in range(slowR_n):
			if abs(stack_Rslows[slow_i]) < lowest_Rslow:
				lowest_Rindex = slow_i
				lowest_Rslow = abs(stack_Rslows[slow_i])

		print(str(slowR_n) + ' R slownesses, ' + str(lowest_Rindex) + ' min R slow, min is ' + str(lowest_Rslow))

		# Select only stacks with that slowness for transverse plot
		centralT_st1 = Stream()
		centralT_st2 = Stream()
		for slowT_i in range(slowT_n):
			centralT_st1 += st1[lowest_Rindex*slowT_n + slowT_i]
			centralT_st2 += st2[lowest_Rindex*slowT_n + slowT_i]
#%%
	total_slows = slowR_n * slowT_n
	global_max = 0
	for slow_i in range(total_slows): # find global max, and if requested, take envelope
		if len(st1[slow_i].data) == 0:
				print('%d data has zero length ' % (slow_i))
		local_max = max(abs(st1[slow_i].data))
		if local_max > global_max:
			global_max = local_max

	#%% compute timing time series
	ttt = (np.arange(len(st1[0].data)) * st1[0].stats.delta - start_buff) # in units of seconds

	#%%
	# Plot radial slowness plots
	fig_index = 3
	plt.close(fig_index)
	plt.figure(fig_index,figsize=(10,10))
	plt.xlim(-start_buff,end_buff)
	plt.ylim(stack_Rslows[0], stack_Rslows[-1])
#	for tr in st1good:
	for slowR_i in range(slowR_n):  # for this station, loop over slownesses
		dist_offset = stack_Rslows[slowR_i] # trying for approx degrees
		ttt = (np.arange(len(centralR_st1[slowR_i].data)) * centralR_st1[slowR_i].stats.delta
		 + (centralR_st1[slowR_i].stats.starttime - t1))
		if norm==1:
			plt.plot(ttt, (centralR_st1[slowR_i].data - np.median(centralR_st1[slowR_i].data))*plot_scale_fac /(centralR_st1[slowR_i].data.max()
				- centralR_st1[slowR_i].data.min()) + dist_offset, color = 'green')
			plt.plot(ttt, (centralR_st2[slowR_i].data - np.median(centralR_st2[slowR_i].data))*plot_scale_fac /(centralR_st2[slowR_i].data.max()
				- centralR_st2[slowR_i].data.min()) + dist_offset, color = 'red')
		else:
			plt.plot(ttt, (centralR_st1[slowR_i].data - np.median(centralR_st1[slowR_i].data))*plot_scale_fac /global_max + dist_offset, color = 'green')
			plt.plot(ttt, (centralR_st2[slowR_i].data - np.median(centralR_st2[slowR_i].data))*plot_scale_fac /global_max + dist_offset, color = 'red')
	plt.title('Radial stacks at 0 T slowness')

	# Plot transverse slowness plots
	fig_index = 4
	plt.close(fig_index)
	plt.figure(fig_index,figsize=(10,10))
	plt.xlim(-start_buff,end_buff)
	plt.ylim(stack_Tslows[0], stack_Tslows[-1])
#	for tr in st1good:
	for slowT_i in range(slowT_n):  # for this station, loop over slownesses
		dist_offset = stack_Tslows[slowT_i] # trying for approx degrees
		ttt = (np.arange(len(centralT_st1[slowT_i].data)) * centralT_st1[slowT_i].stats.delta
		 + (centralT_st1[slowT_i].stats.starttime - t1))
		if norm==1:
			plt.plot(ttt, (centralT_st1[slowT_i].data - np.median(centralT_st1[slowT_i].data))*plot_scale_fac /(centralT_st1[slowT_i].data.max()
				- centralT_st1[slowT_i].data.min()) + dist_offset, color = 'green')
			plt.plot(ttt, (centralT_st2[slowT_i].data - np.median(centralT_st2[slowT_i].data))*plot_scale_fac /(centralT_st2[slowT_i].data.max()
				- centralT_st2[slowT_i].data.min()) + dist_offset, color = 'red')
		else:
			plt.plot(ttt, (centralT_st1[slowT_i].data - np.median(centralT_st1[slowT_i].data))*plot_scale_fac /global_max + dist_offset, color = 'green')
			plt.plot(ttt, (centralT_st2[slowT_i].data - np.median(centralT_st2[slowT_i].data))*plot_scale_fac /global_max + dist_offset, color = 'red')
	plt.title('Transverse stacks at 0 R slowness')

	elapsed_time_wc = time.time() - start_time_wc
	print('This job took ' + str(elapsed_time_wc) + ' seconds')
	os.system('say "Done"')
