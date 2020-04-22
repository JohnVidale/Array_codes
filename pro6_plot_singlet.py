#!/usr/bin/env python
# Read in 2D stacks for two events
# Compute tdiff, ave_amp, amp_ratio
# Plot radial and transverse cuts through stack, plus beam sum
# Write out tdiff, ave_amp, amp_ratio results
# John Vidale 3/2019

def pro6stacked_singlet(eq_file, plot_scale_fac = 0.03, slow_delta = 0.0005,
			  slowR_lo = -0.1, slowR_hi = 0.1, slowT_lo = -0.1, slowT_hi = 0.1,
			  start_buff = -50, end_buff = 50, start_beam = 0, end_beam = 0,
			  norm = 0, dphase = 'PKiKP',
			  plot_dyn_range = 1000, fig_index = 401, get_stf = 0, ref_phase = 'blank',
			  ARRAY = 0, R_slow_plot = 0, T_slow_plot = 0, event_no = 0,
			  ref_loc = 0, ref_lat = 36.3, ref_lon = 138.5, NS = 1):

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

	start_time_wc = time.time()
	plot_trace_sums = 0
	print('Running pro6_plot_singlet')
#%% Get info
	#%% get big info file with event location and more

	sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/events_good.txt'
	with open(sta_file, 'r') as file:
		lines = file.readlines()
	event_count = len(lines)

	print(str(event_count) + ' lines read from ' + sta_file)
	station_index = range(event_count)
	event_names        = []

	event_index = np.zeros(event_count)
	event_year  = np.zeros(event_count)
	event_mo    = np.zeros(event_count)
	event_day   = np.zeros(event_count)
	event_hr    = np.zeros(event_count)
	event_min   = np.zeros(event_count)
	event_sec   = np.zeros(event_count)
	event_lat   = np.zeros(event_count)
	event_lon   = np.zeros(event_count)
	event_dep   = np.zeros(event_count)
	event_mb    = np.zeros(event_count)
	event_ms    = np.zeros(event_count)
	event_tstart       = np.zeros(event_count)
	event_tend         = np.zeros(event_count)
	event_gcdist       = np.zeros(event_count)
	event_dist         = np.zeros(event_count)
	event_baz          = np.zeros(event_count)
	event_SNR          = np.zeros(event_count)
	event_Sflag        = np.zeros(event_count)
	event_PKiKPflag    = np.zeros(event_count)
	event_ICSflag      = np.zeros(event_count)
	event_PKiKP_radslo = np.zeros(event_count)
	event_PKiKP_traslo = np.zeros(event_count)
	event_PKiKP_qual   = np.zeros(event_count)
	event_ICS_qual     = np.zeros(event_count)

	iii = 0
	for ii in station_index:   # read file
		line = lines[ii]
		split_line = line.split()

		event_index[ii]  = float(split_line[0])
		event_names.append(split_line[1])
		event_year[ii]   = float(split_line[2])
		event_mo[ii]     = float(split_line[3])
		event_day[ii]    = float(split_line[4])
		event_hr[ii]     = float(split_line[5])
		event_min[ii]    = float(split_line[6])
		event_sec[ii]    = float(split_line[7])
		event_lat[ii]    = float(split_line[8])
		event_lon[ii]    = float(split_line[9])
		event_dep[ii]    = float(split_line[10])
		event_mb[ii]     = float(split_line[11])
		event_ms[ii]     = float(split_line[12])
		event_tstart[ii] = float(split_line[13])
		event_tend[ii]   = float(split_line[14])
		event_gcdist[ii] = float(split_line[15])
		event_dist[ii]   = float(split_line[16])
		event_baz[ii]    = float(split_line[17])
		event_SNR[ii]    = float(split_line[18])
		event_Sflag[ii]  = float(split_line[19])
		event_PKiKPflag[ii]     = float(split_line[20])
		event_ICSflag[ii]       = float(split_line[21])
		event_PKiKP_radslo[ii]  = float(split_line[22])
		event_PKiKP_traslo[ii]  = float(split_line[23])
		event_PKiKP_qual[ii]    = float(split_line[24])
		event_ICS_qual[ii]      = float(split_line[25])
		if event_index[ii] == event_no:
			iii = ii

	if iii == 0:
		print('Event ' + str(event_no) + ' not found')
	else:
		print('Event ' + str(event_no) + ' is ' + str(iii))

	if ref_loc == 0:
		if ARRAY == 0:
			ref_lat = 36.3  # °N, around middle of Japan
			ref_lon = 138.5 # °E
		elif ARRAY == 1:
			ref_lat = 46.7  # °N keep only inner rings A-D
			ref_lon = -106.22   # °E
		elif ARRAY == 2: # China set and center
			ref_lat = 38      # °N
			ref_lon = 104.5   # °E
	ref_dist_az = gps2dist_azimuth(event_lat[iii],event_lon[iii],ref_lat,ref_lon)
	ref_back_az = ref_dist_az[2]
	ref_dist    = ref_dist_az[0]/(1000*111)

	#  find predicted slowness
	arrivals1 = model.get_travel_times(source_depth_in_km=event_dep[iii],distance_in_degree=ref_dist-0.5,phase_list=[dphase])
	arrivals2 = model.get_travel_times(source_depth_in_km=event_dep[iii],distance_in_degree=ref_dist+0.5,phase_list=[dphase])
	dtime = arrivals2[0].time - arrivals1[0].time
	event_pred_slo  = dtime/111.  # s/km

	# convert to pred rslo and tslo
	sin_baz = np.sin(ref_back_az * np.pi /180)
	cos_baz = np.cos(ref_back_az * np.pi /180)
	if NS == 1:
		pred_Nslo = event_pred_slo * cos_baz
		pred_Eslo = event_pred_slo * sin_baz
	else:
		pred_Nslo = event_pred_slo
		pred_Eslo = 0


	#  rotate observed slowness to N and E
	obs_Nslo = (event_PKiKP_radslo[iii] * cos_baz) - (event_PKiKP_traslo[iii] * sin_baz)
	obs_Eslo = (event_PKiKP_radslo[iii] * sin_baz) + (event_PKiKP_traslo[iii] * cos_baz)
	#  find observed back-azimuth
#	bazi_rad = np.arctan(event_PKiKP_traslo[ii]/event_PKiKP_radslo[ii])
#	event_obs_bazi  = event_baz[ii] + (bazi_rad * 180 / np.pi)
	print('PR '+ str(pred_Nslo) + ' PT ' + str(pred_Eslo) + ' OR ' + str(obs_Nslo) + ' OT ' + str(obs_Eslo))

	#%% get date and time
	goto = '/Users/vidale/Documents/PyCode/EvLocs'
	os.chdir(goto)

	file = open(eq_file, 'r')
	lines=file.readlines()
	split_line = lines[0].split()
	t           = UTCDateTime(split_line[1])
	date_label  = split_line[1][0:10]

	#%% read stack files
	# date_label = '2018-04-02' # date for filename
	goto = '/Users/vidale/Documents/PyCode/Pro_files'
	os.chdir(goto)
	fname = 'HD' + date_label + '_2dstack.mseed'
	st = Stream()
	st = read(fname)

	amp_ave   = st.copy()  # make array for amplitude

	print('Read in: event ' + str(len(st)) + ' traces')
	nt = len(st[0].data)
	dt = st[0].stats.delta
	print('First trace has ' + str(nt) + ' time pts, time sampling of ' + str(dt) + ' and thus duration of ' + str((nt-1)*dt))

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
	for slow_i in range(total_slows): # find envelope and global max
		if slow_i % 200 == 0:
			print('At line 101, ' +str(slow_i) + ' slowness out of ' + str(total_slows))
		if len(st[slow_i].data) == 0: # test for zero-length traces
				print('%d data has zero length ' % (slow_i))

		seismogram = hilbert(st[slow_i].data)  # make analytic seismograms

		env = np.abs(seismogram) # amplitude envelope
		amp_ave[slow_i].data    = env

		local_max = max(abs(amp_ave[slow_i].data))
		if local_max > global_max:
			global_max = local_max
	#%% Make plots summing across all R slownesses at a given T slowness, then the converse
	#find transverse slowness nearest T_slow_plot
	if plot_trace_sums == 1:
		lowest_Tslow = 1000000
		for slow_i in range(slowT_n):
			if abs(stack_Tslows[slow_i] - T_slow_plot) < lowest_Tslow:
				lowest_Tindex = slow_i
				lowest_Tslow = abs(stack_Tslows[slow_i] - T_slow_plot)

		print(str(slowT_n) + ' T slownesses, index ' + str(lowest_Tindex) + ' is closest to input parameter ' + str(T_slow_plot) + ', slowness diff there is ' + str(lowest_Tslow) + ' and slowness is ' + str(stack_Tslows[lowest_Tindex]))
		# Select only stacks with that slowness for radial plot
		centralR_st  = Stream()
		centralR_amp   = Stream()
		for slowR_i in range(slowR_n):
			ii = slowR_i*slowT_n + lowest_Tindex
			centralR_st += st[ii]
			centralR_amp   += amp_ave[ii]

		#%% find radial slowness nearest R_slow_plot
		lowest_Rslow = 1000000
		for slow_i in range(slowR_n):
			if abs(stack_Rslows[slow_i] - R_slow_plot) < lowest_Rslow:
				lowest_Rindex = slow_i
				lowest_Rslow = abs(stack_Rslows[slow_i] - R_slow_plot)

		print(str(slowR_n) + ' R slownesses, index ' + str(lowest_Rindex) + ' is closest to input parameter ' + str(R_slow_plot) + ', slowness diff there is ' + str(lowest_Rslow) + ' and slowness is ' + str(stack_Rslows[lowest_Rindex]))

		# Select only stacks with that slowness for transverse plot
		centralT_st = Stream()
		centralT_amp   = Stream()

		for slowT_i in range(slowT_n):
			ii = lowest_Rindex*slowT_n + slowT_i
			centralT_st += st[ii]
			centralT_amp   += amp_ave[ii]

	#%% compute timing time series
	ttt = (np.arange(len(st[0].data)) * st[0].stats.delta + start_buff) # in units of seconds

#%% Plot radial amp and tdiff vs time plots
	if plot_trace_sums == 1:
		fig_index = 6
		plt.figure(fig_index,figsize=(30,10))
		plt.xlim(start_buff,end_buff)
		plt.ylim(stack_Rslows[0], stack_Rslows[-1])
		for slowR_i in range(slowR_n):  # loop over radial slownesses
			dist_offset = stack_Rslows[slowR_i] # trying for approx degrees
			ttt = (np.arange(len(centralR_st[slowR_i].data)) * centralR_st[slowR_i].stats.delta
			 + (centralR_st[slowR_i].stats.starttime - t))
			plt.plot(ttt, (centralR_st[slowR_i].data - np.median(centralR_st[slowR_i].data))*plot_scale_fac /global_max + dist_offset, color = 'green')
		plt.xlabel('Time (s)')
		plt.ylabel('R Slowness (s/km)')
		plt.title(str(event_no) + '  ' + date_label + '  ' + dphase + ' seismograms at ' + str(T_slow_plot) + ' T slowness')
		# Plot transverse amp and tdiff vs time plots
		fig_index = 7
		plt.figure(fig_index,figsize=(30,10))
		plt.xlim(start_buff,end_buff)
		plt.ylim(stack_Tslows[0], stack_Tslows[-1])

		for slowT_i in range(slowT_n):  # loop over transverse slownesses
			dist_offset = stack_Tslows[slowT_i] # trying for approx degrees
			ttt = (np.arange(len(centralT_st[slowT_i].data)) * centralT_st[slowT_i].stats.delta
			 + (centralT_st[slowT_i].stats.starttime - t))
			plt.plot(ttt, (centralT_st[slowT_i].data - np.median(centralT_st[slowT_i].data))*plot_scale_fac /global_max + dist_offset, color = 'green')
		plt.xlabel('Time (s)')
		plt.ylabel('T Slowness (s/km)')
		plt.title(str(event_no) + '  ' + date_label + '  ' + dphase + ' seismograms at ' + str(R_slow_plot) + ' R slowness')
	#	os.chdir('/Users/vidale/Documents/PyCode/Plots')
	#	plt.savefig(date_label1 + '_' + str(start_buff) + '_' + str(end_buff) + '_stack.png')

#%% R-T amplitude averaged over time window
	fig_index = 9
	stack_slice = np.zeros((slowR_n,slowT_n))
	smax = 0
	if start_beam == 0 and end_beam == 0:
		full_beam = 1
	else:  # beam just part of stack volume
		full_beam = 0
		start_index = int((start_beam - start_buff) / dt)
		end_index   = int((end_beam   - start_buff) / dt)
		print('beam is ' + str(start_beam) + ' to ' + str(end_beam) + 's, out of ' + str(start_buff)
			+ ' to ' + str(end_buff) + 's, dt is ' + str(dt)  + 's, and indices are '+ str(start_index) + ' ' + str(end_index))
	for slowR_i in range(slowR_n):  # loop over radial slownesses
		for slowT_i in range(slowT_n):  # loop over transverse slownesses
			index = slowR_i*slowT_n + slowT_i
			if full_beam == 1:
				num_val = np.nanmedian(amp_ave[index].data)
			else:
				num_val = np.nanmedian(amp_ave[index].data[start_index:end_index])
			stack_slice[slowR_i, slowT_i] = num_val
			if num_val > smax:
				smax = num_val

	y1, x1 = np.mgrid[slice(stack_Rslows[0], stack_Rslows[-1] + slow_delta, slow_delta),
				 slice(stack_Tslows[0], stack_Tslows[-1] + slow_delta, slow_delta)]

	fig, ax = plt.subplots(1, figsize=(7,6))
	c = ax.pcolormesh(x1, y1, stack_slice/smax, cmap=plt.cm.gist_rainbow_r, vmin = 0)
	fig.colorbar(c, ax=ax, label='linear amplitude')
	ax.axis([x1.min(), x1.max(), y1.min(), y1.max()])
	circle1 = plt.Circle((0, 0), 0.019, color='black', fill=False)
	ax.add_artist(circle1)  # inner core limit
	circle2 = plt.Circle((0, 0), 0.040, color='black', fill=False)
	ax.add_artist(circle2)  # outer core limit

	c = ax.scatter(pred_Eslo, pred_Nslo, color='black'  , s=50, alpha=0.75)
#	c = ax.scatter( obs_Eslo,  obs_Nslo, color='purple', s=100, alpha=0.75)
	c = ax.scatter(        0,         0, color='black' , s=50,  alpha=0.75)

	if NS == 1:
		plt.xlabel('East Slowness (s/km)')
		plt.ylabel('North Slowness (s/km)')
	else:
		plt.xlabel('Transverse Slowness (s/km)')
		plt.ylabel('Radial Slowness (s/km)')
	if start_beam == 0 and end_beam == 0:
		plt.title(str(event_no) + '  ' + date_label + '  ' + dphase + ' beam amplitude' + ' ' + str(start_buff) + '-' + str(end_buff) + 's')
	else:
		plt.title(str(event_no) + '  ' + date_label + '  ' + dphase + ' selected beam amplitude' + ' ' + str(start_beam) + '-' + str(end_beam) + 's')
	os.chdir('/Users/vidale/Documents/PyCode/Plots')
	plt.savefig(date_label + '_' + str(start_buff) + '_' + str(end_buff) + '_' + str(start_beam) + '-' + str(end_beam) + '_beam.png')
	plt.show()

#%%  Save processed files
#	goto = '/Users/vidale/Documents/PyCode/Pro_Files'
#	os.chdir(goto)
#
#	fname = 'HD' + date_label + '_amp_ave.mseed'
#	amp_ave.write(fname,format = 'MSEED')

#%% Option to write out stf
	elapsed_time_wc = time.time() - start_time_wc
	print('This job took ' + str(elapsed_time_wc) + ' seconds')
	os.system('say "Done"')
