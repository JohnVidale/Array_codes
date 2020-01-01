#!/usr/bin/env python
# Read in 2D stacks for two events
# Compute tdiff, ave_amp, amp_ratio
# Plot radial and transverse cuts through stack, plus beam sum
# Write out tdiff, ave_amp, amp_ratio results
# John Vidale 3/2019

def pro6stacked_seis(eq_file1, eq_file2, plot_scale_fac = 0.03, slow_delta = 0.0005,
			  slowR_lo = -0.1, slowR_hi = 0.1, slowT_lo = -0.1, slowT_hi = 0.1,
			  start_buff = 50, end_buff = 50, norm = 0, freq_corr = 1.0,
			  plot_dyn_range = 1000, fig_index = 401, get_stf = 0, ref_phase = 'blank',
			  ARRAY = 0, min_rat = 0.6, max_rat = 1.8, min_amp = 0.3):

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
	import statistics

#%% Get info
	#%% get locations
	print('Running pro6_plot_stacked_seis')
	start_time_wc = time.time()

	if ARRAY == 0:
		file = open(eq_file1, 'r')
	elif ARRAY == 1:
		file = open('EvLocs/' + eq_file1, 'r')
	lines=file.readlines()
	split_line = lines[0].split()
	t1           = UTCDateTime(split_line[1])
	date_label1  = split_line[1][0:10]

	if ARRAY == 0:
		file = open(eq_file2, 'r')
	elif ARRAY == 1:
		file = open('EvLocs/' + eq_file2, 'r')
	lines=file.readlines()
	split_line = lines[0].split()
	t2           = UTCDateTime(split_line[1])
	date_label2  = split_line[1][0:10]

	#%% read files
	# #%% Get saved event info, also used to name files
	# date_label = '2018-04-02' # date for filename
	if ARRAY == 0:
		fname1 = 'HD' + date_label1 + '_2dstack.mseed'
		fname2 = 'HD' + date_label2 + '_2dstack.mseed'
	elif ARRAY == 1:
		fname1 = 'Pro_Files/HD' + date_label1 + '_2dstack.mseed'
		fname2 = 'Pro_Files/HD' + date_label2 + '_2dstack.mseed'
	st1 = Stream()
	st2 = Stream()
	st1 = read(fname1)
	st2 = read(fname2)

	tshift    = st1.copy()  # make array for time shift
	amp_ratio = st1.copy()  # make array for relative amplitude
	amp_ave   = st1.copy()  # make array for relative amplitude

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
	print(str(slowT_n) + ' trans slownesses, hi and lo are ' + str(slowT_hi) + '  ' + str(slowT_lo))
	# In English, stack_slows = range(slow_n) * slow_delta - slow_lo
	a1R = range(slowR_n)
	a1T = range(slowT_n)
	stack_Rslows = [(x * slow_delta + slowR_lo) for x in a1R]
	stack_Tslows = [(x * slow_delta + slowT_lo) for x in a1T]
	print(str(slowR_n) + ' radial slownesses, ' + str(slowT_n) + ' trans slownesses, ')

#%%  Loop over slowness
	total_slows = slowR_n * slowT_n
	global_max = 0
	for slow_i in range(total_slows): # find envelope, phase, tshift, and global max
		if slow_i % 200 == 0:
			print('At line 101, ' +str(slow_i) + ' slowness out of ' + str(total_slows))
		if len(st1[slow_i].data) == 0: # test for zero-length traces
				print('%d data has zero length ' % (slow_i))

		seismogram1 = hilbert(st1[slow_i].data)  # make analytic seismograms
		seismogram2 = hilbert(st2[slow_i].data)

		env1 = np.abs(seismogram1) # amplitude
		env2 = np.abs(seismogram2)
		amp_ave[slow_i].data    = 0.5 * (env1 + env2)
		amp_ratio[slow_i].data  = env1/env2

		angle1 = np.angle(seismogram1) # time shift
		angle2 = np.angle(seismogram2)
		phase1 = np.unwrap(angle1)
		phase2 = np.unwrap(angle2)
		dphase = (angle1 - angle2)
#		dphase = phase1 - phase2
		for it in range(nt1):
			if dphase[it] > math.pi:
				dphase[it] -= 2 * math.pi
			elif dphase[it] < -1 * math.pi:
				dphase[it] += 2 * math.pi
			if dphase[it] > math.pi or dphase[it] < -math.pi:
				print(f'Bad dphase value {dphase[it]:.2f}  {it:4d}')
		freq1 = np.diff(phase1) #freq in radians/sec
		freq2 = np.diff(phase2)
		ave_freq = 0.5*(freq1 + freq2)
		ave_freq_plus = np.append(ave_freq,[1]) # ave_freq one element too short
#		tshift[slow_i].data     = dphase / ave_freq_plus # 2*pi top and bottom cancels
		tshift[slow_i].data     = dphase/(2*math.pi*freq_corr)

		local_max = max(abs(amp_ave[slow_i].data))
		if local_max > global_max:
			global_max = local_max
#%% Extract slices
	tshift_full = tshift.copy()  # make array for time shift
	for slow_i in range(total_slows): # ignore less robust points
		if slow_i % 200 == 0:
			print('At line 140, ' +str(slow_i) + ' slowness out of ' + str(total_slows))
		for it in range(nt1):
			if ((amp_ratio[slow_i].data[it] < min_rat) or (amp_ratio[slow_i].data[it] > max_rat) or (amp_ave[slow_i].data[it] < (min_amp * global_max))):
				tshift[slow_i].data[it] = np.nan
	#%% Find transverse slowness nearest zero
	lowest_Tslow = 1000000
	for slow_i in range(slowT_n):
		if abs(stack_Tslows[slow_i]) < lowest_Tslow:
			lowest_Tindex = slow_i
			lowest_Tslow = abs(stack_Tslows[slow_i])

	print(str(slowT_n) + ' T slownesses, ' + str(lowest_Tindex) + ' min T slow, min is ' + str(lowest_Tslow))

	# Select only stacks with that slowness for radial plot
	centralR_st1 = Stream()
	centralR_st2 = Stream()
	centralR_amp   = Stream()
	centralR_ampr  = Stream()
	centralR_tdiff = Stream()
	for slowR_i in range(slowR_n):
		ii = slowR_i*slowT_n + lowest_Tindex
		centralR_st1 += st1[ii]
		centralR_st2 += st2[ii]
		centralR_amp   += amp_ave[ii]
		centralR_ampr  += amp_ratio[ii]
		centralR_tdiff += tshift[ii]

	#%% If desired, find radial slowness nearest zero
	lowest_Rslow = 1000000
	for slow_i in range(slowR_n):
		if abs(stack_Rslows[slow_i]) < lowest_Rslow:
			lowest_Rindex = slow_i
			lowest_Rslow = abs(stack_Rslows[slow_i])

	print(str(slowR_n) + ' R slownesses, ' + str(lowest_Rindex) + ' min R slow, min is ' + str(lowest_Rslow))

	# Select only stacks with that slowness for transverse plot
	centralT_st1 = Stream()
	centralT_st2 = Stream()
	centralT_amp   = Stream()
	centralT_ampr  = Stream()
	centralT_tdiff = Stream()

	#%% to extract stacked time functions
	event1_sample = Stream()
	event2_sample = Stream()

	for slowT_i in range(slowT_n):
		ii = lowest_Rindex*slowT_n + slowT_i
		centralT_st1 += st1[ii]
		centralT_st2 += st2[ii]
		centralT_amp   += amp_ave[ii]
		centralT_ampr  += amp_ratio[ii]
		centralT_tdiff += tshift[ii]

	#%% compute timing time series
	ttt = (np.arange(len(st1[0].data)) * st1[0].stats.delta - start_buff) # in units of seconds

#%% Plot radial amp and tdiff vs time plots
	fig_index = 6
#	plt.close(fig_index)
	plt.figure(fig_index,figsize=(30,10))
	plt.xlim(-start_buff,end_buff)
	plt.ylim(stack_Rslows[0], stack_Rslows[-1])
	for slowR_i in range(slowR_n):  # loop over radial slownesses
		dist_offset = stack_Rslows[slowR_i] # trying for approx degrees
		ttt = (np.arange(len(centralR_st1[slowR_i].data)) * centralR_st1[slowR_i].stats.delta
		 + (centralR_st1[slowR_i].stats.starttime - t1))
		plt.plot(ttt, (centralR_st1[slowR_i].data - np.median(centralR_st1[slowR_i].data))*plot_scale_fac /global_max + dist_offset, color = 'green')
		plt.plot(ttt, (centralR_st2[slowR_i].data - np.median(centralR_st2[slowR_i].data))*plot_scale_fac /global_max + dist_offset, color = 'red')
		# extract stacked time functions
		if get_stf != 0:
			if np.abs(stack_Rslows[slowR_i]- 0.005) < 0.000001: # kludge, not exactly zero when desired
				event1_sample = centralR_st1[slowR_i].copy()
				event2_sample = centralR_st2[slowR_i].copy()
#		plt.plot(ttt, (centralR_amp[slowR_i].data)  *plot_scale_fac/global_max + dist_offset, color = 'purple')
		plt.plot(ttt, (centralR_tdiff[slowR_i].data)*plot_scale_fac/1 + dist_offset, color = 'black')
		plt.plot(ttt, (centralR_amp[slowR_i].data)*0.0 + dist_offset, color = 'lightgray') # reference lines
	plt.title('Seismograms and tdiff at 0 T slowness')
	# Plot transverse amp and tdiff vs time plots
	fig_index = 7
#	plt.close(fig_index)
	plt.figure(fig_index,figsize=(30,10))
	plt.xlim(-start_buff,end_buff)
	plt.ylim(stack_Tslows[0], stack_Tslows[-1])

	for slowT_i in range(slowT_n):  # loop over transverse slownesses
		dist_offset = stack_Tslows[slowT_i] # trying for approx degrees
		ttt = (np.arange(len(centralT_st1[slowT_i].data)) * centralT_st1[slowT_i].stats.delta
		 + (centralT_st1[slowT_i].stats.starttime - t1))
		plt.plot(ttt, (centralT_st1[slowT_i].data - np.median(centralT_st1[slowT_i].data))*plot_scale_fac /global_max + dist_offset, color = 'green')
		plt.plot(ttt, (centralT_st2[slowT_i].data - np.median(centralT_st2[slowT_i].data))*plot_scale_fac /global_max + dist_offset, color = 'red')
#		plt.plot(ttt, (centralT_amp[slowT_i].data)  *plot_scale_fac/global_max + dist_offset, color = 'purple')
		plt.plot(ttt, (centralT_tdiff[slowT_i].data)*plot_scale_fac/1 + dist_offset, color = 'black')
		plt.plot(ttt, (centralT_amp[slowT_i].data)*0.0 + dist_offset, color = 'lightgray') # reference lines
	plt.title(ref_phase + ' seismograms and tdiff at 0 R slowness')
#%% R-T tshift averaged over time window
	fig_index = 8
	stack_slice = np.zeros((slowR_n,slowT_n))
	for slowR_i in range(slowR_n):  # loop over radial slownesses
		for slowT_i in range(slowT_n):  # loop over transverse slownesses
			index = slowR_i*slowT_n + slowT_i
			num_val = np.nanmedian(tshift[index].data)
#			num_val = statistics.median(tshift_full[index].data)
			stack_slice[slowR_i, slowT_i] = num_val # adjust for dominant frequency of 1.2 Hz, not 1 Hz
#	stack_slice[0,0] = -0.25
#	stack_slice[0,1] =  0.25
#	tdiff_clip = 0.4/1.2
	tdiff_clip_max =  0.1  # DO NOT LEAVE COMMENTED OUT!!
	tdiff_clip_min = -0.2

	y1, x1 = np.mgrid[slice(stack_Rslows[0], stack_Rslows[-1] + slow_delta, slow_delta),
				 slice(stack_Tslows[0], stack_Tslows[-1] + slow_delta, slow_delta)]

	fig, ax = plt.subplots(1, figsize=(4.8,6))
#		fig, ax = plt.subplots(1, figsize=(9,2))
#		fig.subplots_adjust(bottom=0.3)
#	c = ax.pcolormesh(x1, y1, stack_slice, cmap=plt.cm.coolwarm, vmin = -tdiff_clip, vmax = tdiff_clip)
	c = ax.pcolormesh(x1, y1, stack_slice, cmap=plt.cm.bwr, vmin = tdiff_clip_min, vmax = tdiff_clip_max)
	ax.axis([x1.min(), x1.max(), y1.min(), y1.max()])
	circle1 = plt.Circle((0, 0), 0.019, color='black', fill=False)
	ax.add_artist(circle1)
	circle2 = plt.Circle((0, 0), 0.040, color='black', fill=False)
	ax.add_artist(circle2)  #outer core limit
	fig.colorbar(c, ax=ax)
	plt.ylabel('R Slowness (s/km)')
	plt.title(ref_phase + ' time shift')
#	plt.title('T-R average time shift ' + date_label1 + ' ' + date_label2)
	plt.show()

#%% R-T amplitude averaged over time window
	fig_index = 9
	stack_slice = np.zeros((slowR_n,slowT_n))
	smax = 0
	for slowR_i in range(slowR_n):  # loop over radial slownesses
		for slowT_i in range(slowT_n):  # loop over transverse slownesses
			index = slowR_i*slowT_n + slowT_i
			num_val = np.nanmedian(amp_ave[index].data)
			stack_slice[slowR_i, slowT_i] = num_val
			if num_val > smax:
				smax = num_val
#	stack_slice[0,0] = 0

	y1, x1 = np.mgrid[slice(stack_Rslows[0], stack_Rslows[-1] + slow_delta, slow_delta),
				 slice(stack_Tslows[0], stack_Tslows[-1] + slow_delta, slow_delta)]

#	fig, ax = plt.subplots(1)
	fig, ax = plt.subplots(1, figsize=(7,6))
#	c = ax.pcolormesh(x1, y1, stack_slice/smax, cmap=plt.cm.gist_yarg, vmin = 0.5)
	c = ax.pcolormesh(x1, y1, stack_slice, cmap=plt.cm.gist_rainbow_r, vmin = 0)
	ax.axis([x1.min(), x1.max(), y1.min(), y1.max()])
	circle1 = plt.Circle((0, 0), 0.019, color='black', fill=False)
	ax.add_artist(circle1)  #inner core limit
	circle2 = plt.Circle((0, 0), 0.040, color='black', fill=False)
	ax.add_artist(circle2)  #outer core limit
	fig.colorbar(c, ax=ax)
	plt.xlabel('Transverse Slowness (s/km)')
	plt.ylabel('Radial Slowness (s/km)')
	plt.title(ref_phase + ' beam amplitude')
#	plt.title('Beam amplitude ' + date_label1 + ' ' + date_label2)
	plt.show()

#%%  Save processed files
	if ARRAY == 0:
		fname = 'HD' + date_label1 + '_' + date_label2 + '_tshift.mseed'
	elif ARRAY == 1:
		fname = 'Pro_Files/HD' + date_label1 + '_' + date_label2 + '_tshift.mseed'
	tshift_full.write(fname,format = 'MSEED')
	if ARRAY == 0:
		fname = 'HD' + date_label1 + '_' + date_label2 + '_amp_ave.mseed'
	elif ARRAY == 1:
		fname = 'Pro_Files/HD' + date_label1 + '_' + date_label2 + '_amp_ave.mseed'
	amp_ave.write(fname,format = 'MSEED')
	if ARRAY == 0:
		fname = 'HD' + date_label1 + '_' + date_label2 + '_amp_ratio.mseed'
	elif ARRAY == 1:
		fname = 'Pro_Files/HD' + date_label1 + '_' + date_label2 + '_amp_ratio.mseed'
	amp_ratio.write(fname,format = 'MSEED')

#%% Option to write out stf
	if get_stf != 0:
		event1_sample.taper(0.1)
		event2_sample.taper(0.1)
		fname = 'Pro_Files/HD' + date_label1 + '_stf.mseed'
		event1_sample.write(fname,format = 'MSEED')
		fname = 'Pro_Files/HD' + date_label2 + '_stf.mseed'
		event2_sample.write(fname,format = 'MSEED')

	elapsed_time_wc = time.time() - start_time_wc
	print('This job took ' + str(elapsed_time_wc) + ' seconds')
	os.system('say "Done"')
