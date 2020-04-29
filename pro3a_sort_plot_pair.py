#!/usr/bin/env python
# input is pair of repeated events from hinet, LASA, or NORSAR
# this program tapers, filters, selects range and SNR
# plots seismograms with traveltime curves, either raw or reduced against traveltimes
# outputs selected traces, "*sel.mseed"
# John Vidale 2/2019

def pro3pair(eq_file1, eq_file2, stat_corr = 1, simple_taper = 0, skip_SNR = 0,
			dphase = 'PKIKP', dphase2 = 'PKiKP', dphase3 = 'PKIKP', dphase4 = 'PKiKP',
			rel_time = 1, start_buff = -200, end_buff = 500,
			plot_scale_fac = 0.05, qual_threshold = 0, corr_threshold = 0.5,
			freq_min = 1, freq_max = 3, min_dist = 0, max_dist = 180, auto_dist = 1,
			alt_statics = 0, statics_file = 'nothing', ARRAY = 0,
			ref_loc = 0, ref_rad = 0.4, ref_lat = 36.3, ref_lon = 138.5):


#%% Import functions
	from obspy import UTCDateTime
	from obspy import Stream
	from obspy import read
	from obspy.geodetics import gps2dist_azimuth
	import numpy as np
	import os
	from obspy.taup import TauPyModel
	import matplotlib.pyplot as plt
	import time
	model = TauPyModel(model='iasp91')

	import sys # don't show any warnings
	import warnings

	if not sys.warnoptions:
	    warnings.simplefilter("ignore")

	print('Running pro3a_sort_plot_pair')
	start_time_wc = time.time()

#%%  Set some parameters
	verbose = 0           # more output
#	rel_time = 1          # timing is relative to a chosen phase, otherwise relative to OT
	taper_frac = .05      #Fraction of window tapered on both ends
	signal_dur = 5.       # signal length used in SNR calculation
	plot_tt = 1           # plot the traveltimes?
	do_decimate = 0         # 0 if no decimation desired
	#ref_loc = 0  # 1 if selecting stations within ref_rad of ref_lat and ref_lon
	             # 0 if selecting stations by distance from earthquake
	if ref_loc == 0:
		if ARRAY == 0:
			ref_lat = 36.3  # °N, around middle of Japan
			ref_lon = 138.5 # °E
		if ARRAY == 1:
			ref_lat = 46.7      # °N keep only inner rings A-D if radius is 0.4°
			ref_lon = -106.22   # °E
		if ARRAY == 2:
			ref_lat = 38      # °N
			ref_lon = 104.5   # °E
#		ref_rad = 0.4    # ° radius (°) set by input or at top

	if rel_time == 0: # SNR requirement not implemented for unaligned traces
		qual_threshold = 0

	# Plot with reduced velocity?
	red_plot = 0
	red_dist = 55
	red_time = 300
	red_slow = 7.2 # seconds per degree
#%% Get saved event info, also used to name files
	#  event 2016-05-28T09:47:00.000 -56.241 -26.935 78
	fname1 = '/Users/vidale/Documents/PyCode/EvLocs/' + eq_file1
	print('Opening ' + eq_file1)
	file = open(fname1, 'r')
	lines=file.readlines()
	split_line = lines[0].split()
#			ids.append(split_line[0])  ignore label for now
	t1           = UTCDateTime(split_line[1])
	date_label1  = split_line[1][0:10]
	year1        = split_line[1][0:4]
	ev_lat1      = float(      split_line[2])
	ev_lon1      = float(      split_line[3])
	ev_depth1    = float(      split_line[4])
	print('1st event: date_label ' + date_label1 + ' time ' + str(t1) + ' lat '
	   + str(ev_lat1) + ' lon ' + str( ev_lon1) + ' depth ' + str(ev_depth1))

	fname2 = '/Users/vidale/Documents/PyCode/EvLocs/' + eq_file2
	print('Opening ' + eq_file2)
	file = open(fname2, 'r')
	lines=file.readlines()
	split_line = lines[0].split()
#			ids.append(split_line[0])  ignore label for now
	t2           = UTCDateTime(split_line[1])
	date_label2  = split_line[1][0:10]
	year2        = split_line[1][0:4]
	ev_lat2      = float(      split_line[2])
	ev_lon2      = float(      split_line[3])
	ev_depth2    = float(      split_line[4])
	print('2nd event: date_label ' + date_label2 + ' time ' + str(t2) + ' lat '
	   + str(ev_lat2) + ' lon ' + str( ev_lon2) + ' depth ' + str(ev_depth2))

#%% Get station location file
	if stat_corr == 1:  # load static terms, only applies to Hinet and LASA
		if ARRAY == 0:
			if alt_statics == 0: # standard set
				sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/sta_statics_hinet.txt'
			else: # custom set made by this event for this event
				sta_file = ('/Users/vidale/Documents/PyCode/Hinet/Array_codes/Files/' + 'HA' +
				   date_label1[:10] + 'pro4_' + dphase + '.statics')
		elif ARRAY == 1:
			sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/sta_statics_LASA.txt'
		elif ARRAY == 2:
			sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/sta_statics_ch.txt'
		with open(sta_file, 'r') as file:
			lines = file.readlines()
		print(str(len(lines)) + ' stations read from ' + sta_file)
		# Load station coords into arrays
		station_index = range(len(lines))
		st_names = []
		st_dist  = []
		st_lats  = []
		st_lons  = []
		st_shift = []
		st_corr  = []
		for ii in station_index:
			line = lines[ii]
			split_line = line.split()
			st_names.append(split_line[0])
			st_dist.append(split_line[1])
			st_lats.append( split_line[2])
			st_lons.append( split_line[3])
			st_shift.append(split_line[4])
			st_corr.append(split_line[5])
	else: # no static terms, always true for LASA or NORSAR
		if ARRAY == 0: # Hinet set
			sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/sta_hinet.txt'
		elif ARRAY == 1: #         LASA set
			sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/sta_LASA.txt'
		elif ARRAY == 2: #         LASA set
			sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/sta_ch.txt'
		else: #         NORSAR set
			sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/sta_NORSAR.txt'
		with open(sta_file, 'r') as file:
			lines = file.readlines()
		print(str(len(lines)) + ' stations read from ' + sta_file)
		# Load station coords into arrays
		station_index = range(len(lines))
		st_names = []
		st_lats  = []
		st_lons  = []
		for ii in station_index:
			line = lines[ii]
			split_line = line.split()
			st_names.append(split_line[0])
			st_lats.append( split_line[1])
			st_lons.append( split_line[2])
	if ARRAY == 0:  # shorten and make upper case Hi-net station names to match station list
		for ii in station_index:
			this_name = st_names[ii]
			this_name_truc = this_name[0:5]
			st_names[ii]  = this_name_truc.upper()

#%% Is taper too long compared to noise estimation window?
	totalt = end_buff - start_buff
	noise_time_skipped = taper_frac * totalt
	if simple_taper == 0:
		if noise_time_skipped >= 0.5 * (-start_buff):
			print('Specified taper of ' + str(taper_frac * totalt) +
			   ' is not big enough compared to available noise estimation window ' +
			   str(-start_buff - noise_time_skipped) + '. May not work well.')
			old_taper_frac = taper_frac
			taper_frac = -0.5*start_buff/totalt
			print('Taper reset from ' + str(old_taper_frac * totalt) + ' to '
			   + str(taper_frac * totalt) + ' seconds.')

#%% Load waveforms and decimate to 10 sps
	st1 = Stream()
	st2 = Stream()
	fname1     = '/Users/vidale/Documents/PyCode/Mseed/HD' + date_label1 + '.mseed'
	fname2     = '/Users/vidale/Documents/PyCode/Mseed/HD' + date_label2 + '.mseed'
	st1=read(fname1)
	st2=read(fname2)

	if do_decimate != 0:
		st1.decimate(do_decimate, no_filter=True)
		st2.decimate(do_decimate, no_filter=True)

	print('1st trace has : ' + str(len(st1[0].data)) + ' time pts ')
	print('st1 has ' + str(len(st1)) + ' traces')
	print('st2 has ' + str(len(st2)) + ' traces')
	print('1st trace starts at ' + str(st1[0].stats.starttime) + ', event at ' + str(t1))
	print('2nd trace starts at ' + str(st2[0].stats.starttime) + ', event at ' + str(t2))

#%% Select by distance, window and adjust start time to align picked times
	st_pickalign1 = Stream()
	st_pickalign2 = Stream()
	tra_in_range  = 0
	tra_sta_found = 0
	nodata        = 0
	min_dist_auto = 180
	max_dist_auto = 0
	min_time_plot =  1000000
	max_time_plot = -1000000

	# not used in all cases, but printed out below
	# only used if rel_slow == 1, preserves 0 slowness, otherwise 0 is set to phase slowness
	ref_distance = gps2dist_azimuth(ref_lat,ref_lon,ev_lat1,ev_lon1)
	ref1_dist  = ref_distance[0]/(1000*111)
	dist_minus = ref1_dist - 0.5
	dist_plus  = ref1_dist + 0.5
	arrivals_ref   = model.get_travel_times(source_depth_in_km=ev_depth1,distance_in_degree=ref1_dist, phase_list=[dphase])
	arrivals_minus = model.get_travel_times(source_depth_in_km=ev_depth1,distance_in_degree=dist_minus,phase_list=[dphase])
	arrivals_plus  = model.get_travel_times(source_depth_in_km=ev_depth1,distance_in_degree=dist_plus ,phase_list=[dphase])
	atime_ref = arrivals_ref[0].time  # phase arrival time at reference distance
	ref_slow = arrivals_plus[0].time - arrivals_minus[0].time  # dt over 1 degree at ref distance

	for tr in st1: # find lat-lon from list, chop, statics, traces one by one
		if float(year1) < 1970: # fix the damn 1969 -> 2069 bug in Gibbon's LASA data
			temp_t = str(tr.stats.starttime)
			temp_tt = '19' + temp_t[2:]
			tr.stats.starttime = UTCDateTime(temp_tt)
		if tr.stats.station in st_names:  # find station in station list
			ii = st_names.index(tr.stats.station)
			tra_sta_found += 1

			if stat_corr != 1 or float(st_corr[ii]) > corr_threshold: # if using statics, reject low correlations
				stalat = float(st_lats[ii]) # look up lat & lon again to find distance
				stalon = float(st_lons[ii])

				distance = gps2dist_azimuth(stalat,stalon,ev_lat1,ev_lon1) # Get traveltimes again, hard to store
				tr.stats.distance=distance[0] # distance in km
				dist = distance[0]/(1000*111)

				in_range = 0  # flag for whether this trace goes into stack
				if ref_loc == 0:  # check whether trace is in distance range from earthquake
					if min_dist < dist and dist < max_dist:
						in_range = 1
						tra_in_range += 1
				elif ref_loc == 1:  # alternately, check whether trace is close enough to ref_location
					ref_distance = gps2dist_azimuth(ref_lat,ref_lon,stalat,stalon)
					ref2_dist = ref_distance[0]/(1000*111)
					if ref2_dist < ref_rad:
						in_range = 1
						tra_in_range += 1
				if in_range == 1:   # trace fulfills the specified criteria for being in range
					s_t = t1 + start_buff
					e_t = t1 + end_buff
					if stat_corr == 1: # apply static station corrections
						tr.stats.starttime -= float(st_shift[ii])
					if rel_time == 0:  #  don't adjust absolute time
						tr.trim(starttime=s_t,endtime = e_t)
					else:              # shift relative to a chosen phase
						arrivals_each = model.get_travel_times(source_depth_in_km=ev_depth1,distance_in_degree=dist,phase_list=[dphase])
						atime_each = arrivals_each[0].time
						if rel_time == 1: # each window has a shift proportional to ref_dist at phase slowness at ref_dist
							s_t += atime_each
							e_t += atime_each
							tr.trim( starttime=s_t, endtime = e_t)
							tr.stats.starttime -= atime_each - (dist-ref1_dist) * ref_slow
						elif rel_time == 2: # each window has a distinct shift, but offset is common to all stations
							s_t += atime_each
							e_t += atime_each
							tr.trim( starttime=s_t, endtime = e_t)
							tr.stats.starttime -= atime_ref
						elif rel_time == 3:  # each station has an individual, chosen-phase shift, phase arrival set to zero
							s_t += atime_each
							e_t += atime_each
							tr.trim( starttime=s_t, endtime = e_t)
							tr.stats.starttime -= atime_each
						elif rel_time == 4: # use same window around chosen phase for all stations, phase arrival set to zero
							s_t += atime_ref
							e_t += atime_ref
							tr.trim( starttime=s_t, endtime = e_t)
							tr.stats.starttime -= atime_ref
						else:
							print('invalid rel_time, must be integer 0 to 4')
							sys.exit()
					if len(tr.data) > 0:
						st_pickalign1 += tr
					else:
						nodata += 1
		else:
			print(tr.stats.station + ' not found in station list ')
#			sys.exit()

	print('After alignment + range and correlation selection - event: ' + str(len(st_pickalign1)) + ' traces')
	print('Traces found: ' + str(tra_sta_found) + ' Traces in range: ' + str(tra_in_range) + ' Traces with no data: ' + str(nodata))
	print(f'ref1_distance  {ref1_dist:.3f}  relative start time  {atime_ref:.3f}')
	print('ref_loc == 1, ref_lat: ' + str(ref_lat) + ' ref_lon: ' + str(ref_lon))
	print(f'last station: distance {dist:.3f}  last station lat: {stalat:.3f}   last station lon: {stalat:.3f}')

	for tr in st2: # find lat-lon from list, chop, statics, traces one by one
		if float(year2) < 1970: # fix the damn 1969 -> 2069 bug in Gibbon's LASA data
			temp_t = str(tr.stats.starttime)
			temp_tt = '19' + temp_t[2:]
			tr.stats.starttime = UTCDateTime(temp_tt)
		if tr.stats.station in st_names:  # find station in station list
			ii = st_names.index(tr.stats.station)
			tra_sta_found += 1

			if stat_corr != 1 or float(st_corr[ii]) > corr_threshold: # if using statics, reject low correlations
				stalat = float(st_lats[ii]) # look up lat & lon again to find distance
				stalon = float(st_lons[ii])

				distance = gps2dist_azimuth(stalat,stalon,ev_lat2,ev_lon2) # Get traveltimes again, hard to store
				tr.stats.distance=distance[0] # distance in km
				dist = distance[0]/(1000*111)

				in_range = 0  # flag for whether this trace goes into stack
				if ref_loc == 0:  # check whether trace is in distance range from earthquake
					if min_dist < dist and dist < max_dist:
						in_range = 1
						tra_in_range += 1
				elif ref_loc == 1:  # alternately, check whether trace is close enough to ref_location
					ref_distance = gps2dist_azimuth(ref_lat,ref_lon,stalat,stalon)
					ref2_dist = ref_distance[0]/(1000*111)
					if ref2_dist < ref_rad:
						in_range = 1
						tra_in_range += 1
				if in_range == 1:   # trace fulfills the specified criteria for being in range
					s_t = t2 + start_buff
					e_t = t2 + end_buff
					if stat_corr == 1: # apply static station corrections
						tr.stats.starttime -= float(st_shift[ii])
					if rel_time == 0:  #  don't adjust absolute time
						tr.trim(starttime=s_t,endtime = e_t)
					else:              # shift relative to a chosen phase
						arrivals_each = model.get_travel_times(source_depth_in_km=ev_depth2,distance_in_degree=dist,phase_list=[dphase])
						atime_each = arrivals_each[0].time
						if rel_time == 1: # each window has a shift proportional to ref_dist at phase slowness at ref_dist
							s_t += atime_each
							e_t += atime_each
							tr.trim( starttime=s_t, endtime = e_t)
							tr.stats.starttime -= atime_each - (dist-ref1_dist) * ref_slow
						elif rel_time == 2: # each window has a distinct shift, but offset is common to all stations
							s_t += atime_each
							e_t += atime_each
							tr.trim( starttime=s_t, endtime = e_t)
							tr.stats.starttime -= atime_ref
						elif rel_time == 3:  # each station has an individual, chosen-phase shift, phase arrival set to zero
							s_t += atime_each
							e_t += atime_each
							tr.trim( starttime=s_t, endtime = e_t)
							tr.stats.starttime -= atime_each
						elif rel_time == 4: # use same window around chosen phase for all stations, phase arrival set to zero
							s_t += atime_ref
							e_t += atime_ref
							tr.trim( starttime=s_t, endtime = e_t)
							tr.stats.starttime -= atime_ref
						else:
							print('invalid rel_time, must be integer 0 to 4')
							sys.exit()
					if len(tr.data) > 0:
						st_pickalign2 += tr
					else:
						nodata += 1
		else:
			print(tr.stats.station + ' not found in station list')

	print('After alignment + range and correlation selection - event: ' + str(len(st_pickalign2)) + ' traces')
	print('Traces found: ' + str(tra_sta_found) + ' Traces in range: ' + str(tra_in_range) + ' Traces with no data: ' + str(nodata))
	print(f'ref1_distance  {ref1_dist:.3f}  relative start time  {atime_ref:.3f}')
	print('ref_loc == 1, ref_lat: ' + str(ref_lat) + ' ref_lon: ' + str(ref_lon))
	print(f'last station: distance {dist:.3f}  last station lat: {stalat:.3f}   last station lon: {stalat:.3f}')

	print('After alignment and range selection: ' + str(len(st_pickalign1)) + ' traces')
	print('After alignment and range selection: ' + str(len(st_pickalign2)) + ' traces')

	#%%
	#print(st) # at length
	if verbose:
		print(st1.__str__(extended=True))
		print(st2.__str__(extended=True))
		if rel_time == 1:
			print(st_pickalign1.__str__(extended=True))
			print(st_pickalign2.__str__(extended=True))

#%%  Detrend, taper, filter
	st_pickalign1.detrend(type='simple')
	st_pickalign2.detrend(type='simple')
	st_pickalign1.taper(taper_frac)
	st_pickalign2.taper(taper_frac)
	st_pickalign1.filter('bandpass', freqmin=freq_min, freqmax=freq_max, corners=4, zerophase=True)
	st_pickalign2.filter('bandpass', freqmin=freq_min, freqmax=freq_max, corners=4, zerophase=True)
	st_pickalign1.taper(taper_frac)
	st_pickalign2.taper(taper_frac)

#%%  Cull further by imposing SNR threshold on both traces
	st1good = Stream()
	st2good = Stream()
	for tr1 in st_pickalign1:
		for tr2 in st_pickalign2:
			if ((tr1.stats.network  == tr2.stats.network) &
			    (tr1.stats.station  == tr2.stats.station)):
				if skip_SNR == 1:
					st1good += tr1
					st2good += tr2
				else:
					# estimate median noise
					t_noise1_start  = int(len(tr1.data) * taper_frac)
					t_noise2_start  = int(len(tr2.data) * taper_frac)
					t_noise1_end    = int(len(tr1.data) * (-start_buff)/(end_buff - start_buff))
					t_noise2_end    = int(len(tr2.data) * (-start_buff)/(end_buff - start_buff))
					noise1          = np.median(abs(tr1.data[t_noise1_start:t_noise1_end]))
					noise2          = np.median(abs(tr2.data[t_noise2_start:t_noise2_end]))
		# estimate median signal
					t_signal1_start = int(len(tr1.data) * (-start_buff)/(end_buff - start_buff))
					t_signal2_start = int(len(tr2.data) * (-start_buff)/(end_buff - start_buff))
					t_signal1_end   = t_signal1_start + int(len(tr1.data) * signal_dur/(end_buff - start_buff))
					t_signal2_end   = t_signal2_start + int(len(tr2.data) * signal_dur/(end_buff - start_buff))
					signal1         = np.median(abs(tr1.data[t_signal1_start:t_signal1_end]))
					signal2         = np.median(abs(tr2.data[t_signal2_start:t_signal2_end]))
		#			test SNR
					SNR1 = signal1/noise1;
					SNR2 = signal2/noise2;
					if (SNR1 > qual_threshold and SNR2 > qual_threshold):
						st1good += tr1
						st2good += tr2
	if skip_SNR == 1:
		print('Matches (no SNR test): ' + str(len(st1good)) + ' traces')
	else:
		print('Match and above SNR threshold: ' + str(len(st1good)) + ' traces')

	#%%  get station lat-lon, compute distance for plot
	min_dist_auto = 180
	max_dist_auto = 0
	min_time_plot =  1000000
	max_time_plot = -1000000

	for tr in st1good:
		if tr.stats.station in st_names:  # find station in station list
			ii = st_names.index(tr.stats.station)
			stalon = float(st_lons[ii]) # look up lat & lon again to find distance
			stalat = float(st_lats[ii])
			distance = gps2dist_azimuth(stalat,stalon,ev_lat1,ev_lon1)
			tr.stats.distance=distance[0]/(1000*111) # distance in km
			if tr.stats.distance < min_dist_auto:
				min_dist_auto = tr.stats.distance
			if tr.stats.distance > max_dist_auto:
				max_dist_auto = tr.stats.distance
			if tr.stats.starttime - t1 < min_time_plot:
				min_time_plot = tr.stats.starttime - t1
			if ((tr.stats.starttime - t1) + ((len(tr.data)-1) * tr.stats.delta)) > max_time_plot:
				max_time_plot =  ((tr.stats.starttime - t1) + ((len(tr.data)-1) * tr.stats.delta))
	for tr in st2good:
		if tr.stats.station in st_names:  # find station in station list
			ii = st_names.index(tr.stats.station)
			stalon = float(st_lons[ii]) # look up lat & lon again to find distance
			stalat = float(st_lats[ii])
			distance = gps2dist_azimuth(stalat,stalon,ev_lat2,ev_lon2)
			tr.stats.distance=distance[0]/(1000*111) # distance in km

	print(f'Min distance is   {min_dist_auto:.3f}   Max distance is {max_dist_auto:.3f}')
	print(f'Min time is   {min_time_plot:.2f}   Max time is {max_time_plot:.2f}')

	#%%
	# plot traces
	fig_index = 3
	plt.close(fig_index)
	plt.figure(fig_index,figsize=(8,8))
	plt.xlim(start_buff,end_buff)

	if auto_dist == 1:
		dist_diff = max_dist_auto - min_dist_auto # add space at extremes
		plt.ylim(min_dist_auto - 0.1 * dist_diff, max_dist_auto + 0.1 * dist_diff)
	else:
		plt.ylim(min_dist,max_dist)

	for tr in st1good:
		dist_offset = tr.stats.distance # trying for approx degrees
		ttt = np.arange(len(tr.data)) * tr.stats.delta + (tr.stats.starttime - t1)
		if red_plot == 1:
			shift = red_time + (dist_offset - red_dist) * red_slow
			ttt = ttt - shift
		plt.plot(ttt, (tr.data - np.median(tr.data))*plot_scale_fac /(tr.data.max()
			- tr.data.min()) + dist_offset, color = 'green')

	for tr in st2good:
		dist_offset = tr.stats.distance # trying for approx degrees
		ttt = np.arange(len(tr.data)) * tr.stats.delta + (tr.stats.starttime - t2)
		if red_plot == 1:
			shift = red_time + (dist_offset - red_dist) * red_slow
			ttt = ttt - shift
		ttt = ttt
		plt.plot(ttt, (tr.data - np.median(tr.data))*plot_scale_fac /(tr.data.max()
			- tr.data.min()) + dist_offset, color = 'red')
	print('And made it to here.')

#%% Plot traveltime curves
	if rel_time != 1:
		if plot_tt:
			# first traveltime curve
			line_pts = 50
			dist_vec  = np.arange(min_dist_auto, max_dist_auto, (max_dist_auto - min_dist_auto)/line_pts) # distance grid
			time_vec1 = np.arange(min_dist_auto, max_dist_auto, (max_dist_auto - min_dist_auto)/line_pts) # empty time grid of same length (filled with -1000)
			for i in range(0,line_pts):
				arrivals = model.get_travel_times(source_depth_in_km=ev_depth1,distance_in_degree
											=dist_vec[i],phase_list=[dphase])
				num_arrivals = len(arrivals)
				found_it = 0
				for j in range(0,num_arrivals):
					if arrivals[j].name == dphase:
						time_vec1[i] = arrivals[j].time
						found_it = 1
				if found_it == 0:
					time_vec1[i] = np.nan
			# second traveltime curve
			if dphase2 != 'no':
				time_vec2 = np.arange(min_dist_auto, max_dist_auto, (max_dist_auto - min_dist_auto)/line_pts) # empty time grid of same length (filled with -1000)
				for i in range(0,line_pts):
					arrivals = model.get_travel_times(source_depth_in_km=ev_depth1,distance_in_degree
												=dist_vec[i],phase_list=[dphase2])
					num_arrivals = len(arrivals)
					found_it = 0
					for j in range(0,num_arrivals):
						if arrivals[j].name == dphase2:
							time_vec2[i] = arrivals[j].time
							found_it = 1
					if found_it == 0:
						time_vec2[i] = np.nan
				if   rel_time == 3 or rel_time == 4:
					time_vec2 = time_vec2 - time_vec1
				elif rel_time == 2:
					time_vec2 = time_vec2 - atime_ref
				plt.plot(time_vec2,dist_vec, color = 'orange')
			# third traveltime curve
			if dphase3 != 'no':
				time_vec3 = np.arange(min_dist_auto, max_dist_auto, (max_dist_auto - min_dist_auto)/line_pts) # empty time grid of same length (filled with -1000)
				for i in range(0,line_pts):
					arrivals = model.get_travel_times(source_depth_in_km=ev_depth1,distance_in_degree
												=dist_vec[i],phase_list=[dphase3])
					num_arrivals = len(arrivals)
					found_it = 0
					for j in range(0,num_arrivals):
						if arrivals[j].name == dphase3:
							time_vec3[i] = arrivals[j].time
							found_it = 1
					if found_it == 0:
						time_vec3[i] = np.nan
				if   rel_time == 3 or rel_time == 4:
					time_vec2 = time_vec2 - time_vec1
				elif rel_time == 2:
					time_vec2 = time_vec2 - atime_ref
				plt.plot(time_vec3,dist_vec, color = 'yellow')
			# fourth traveltime curve
			if dphase4 != 'no':
				time_vec4 = np.arange(min_dist, max_dist_auto, (max_dist_auto - min_dist)/line_pts) # empty time grid of same length (filled with -1000)
				for i in range(0,line_pts):
					arrivals = model.get_travel_times(source_depth_in_km=ev_depth1,distance_in_degree
												=dist_vec[i],phase_list=[dphase4])
					num_arrivals = len(arrivals)
					found_it = 0
					for j in range(0,num_arrivals):
						if arrivals[j].name == dphase4:
							time_vec4[i] = arrivals[j].time
							found_it = 1
					if found_it == 0:
						time_vec4[i] = np.nan
				if   rel_time == 3 or rel_time == 4:
					time_vec2 = time_vec2 - time_vec1
				elif rel_time == 2:
					time_vec2 = time_vec2 - atime_ref
				plt.plot(time_vec4,dist_vec, color = 'purple')

			if   rel_time == 3 or rel_time == 4:
				time_vec1 = time_vec1 - time_vec1
			elif rel_time == 2:
				time_vec1 = time_vec1 - atime_ref
			plt.plot(time_vec1,dist_vec, color = 'blue')
			plt.show()

	plt.xlabel('Time (s)')
	plt.ylabel('Epicentral distance from event (°)')
	plt.title(dphase + ' for ' + fname1[39:49] + ' vs ' + fname2[39:49])
	plt.show()

#%%  Save processed files
	fname1 = '/Users/vidale/Documents/PyCode/Pro_Files/HD' + date_label1 + 'sel.mseed'
	fname2 = '/Users/vidale/Documents/PyCode/Pro_Files/HD' + date_label2 + 'sel.mseed'
	st1good.write(fname1,format = 'MSEED')
	st2good.write(fname2,format = 'MSEED')

	elapsed_time_wc = time.time() - start_time_wc
	print('This job took ' + str(elapsed_time_wc) + ' seconds')
	os.system('say "Done"')