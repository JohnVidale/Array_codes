#!/usr/bin/env python
# reads in raw mseed traces from a single event
# this program tapers, filters, selects range and SNR
# plots against traveltime curves, either raw or reduced against traveltimes
# writes out "*sel.mseed" file
# plot lines are blue, orange, yellow, purple for phases 1 through 4
# John Vidale 2/2019

def pro3singlet(eq_file, stat_corr = 0, rel_time = 1, simple_taper = 0, skip_SNR = 0,
			dphase = 'PKiKP', dphase2 = 'PKKP', dphase3 = 'PKIKP', dphase4 = 'PPP',
			start_buff = -10, end_buff = 30,
			plot_scale_fac = 0.05, qual_threshold = 0, corr_threshold = 0,
			freq_min = 0.25, freq_max = 1, do_filt = 1,
			min_dist = 0, max_dist = 180, auto_dist = 0, do_decimate = 0,
			alt_statics = 0, statics_file = 'nothing', ARRAY = 0, JST = 0, ref_loc = 0, ref_rad = 0.4,
			verbose = 0, fig_index = 101, event_no = 0, ref_lat = 36.3, ref_lon = 138.5):
# 0 is Hinet, 1 is LASA, 2 is NORSAR

#%% Import functions
	from obspy import UTCDateTime
	from obspy import Stream
	from obspy import read
	from obspy.geodetics import gps2dist_azimuth
	import numpy as np
	import os
	import sys
	from obspy.taup import TauPyModel
	import matplotlib.pyplot as plt
	import time
	model = TauPyModel(model='iasp91')

#	import sys # don't show any warnings
#	import warnings
#
#	if not sys.warnoptions:
#	    warnings.simplefilter("ignore")

	print('Running pro3b_sort_plot_singlet')
	start_time_wc = time.time()

#%% Get saved event info, also used to name files
	#  input event data with 1-line file of format
	#  event 2016-05-28T09:47:00.000 -56.241 -26.935 78
	file = open('/Users/vidale/Documents/PyCode/EvLocs/' + eq_file, 'r')
	lines=file.readlines()
	split_line = lines[0].split()
#			ids.append(split_line[0])  ignore label for now
	t           = UTCDateTime(split_line[1])
	date_label  = split_line[1][0:10]
	date_long   = split_line[1][0:22]
	year        = split_line[1][0:4]
	ev_lat      = float(      split_line[2])
	ev_lon      = float(      split_line[3])
	ev_depth    = float(      split_line[4])
	print('date_label ' + date_label + ' time ' + str(t) + ' lat ' + str(ev_lat) + ' lon ' + str( ev_lon) + ' depth ' + str(ev_depth))

#%% Get station location file
	if stat_corr == 1:  # load static terms, only applies to Hinet, LASA, and China
		if ARRAY == 0:
			if alt_statics == 0: # standard set
				sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/sta_statics_hinet.txt'
			else: # custom set made by this event for this event
				sta_file = ('/Users/vidale/Documents/PyCode/Hinet/Array_codes/Files/' + 'HA' +
				   date_label[:10] + 'pro4_' + dphase + '.statics')
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
			if ARRAY == 0 or ARRAY == 1:
				st_dist.append(split_line[1])
				st_lats.append( split_line[2])
				st_lons.append( split_line[3])
				st_shift.append(split_line[4])
				st_corr.append(split_line[5])
			elif ARRAY == 2:
#				st_dist.append(split_line[1])
				st_lats.append( split_line[1])
				st_lons.append( split_line[2])
				st_shift.append(split_line[3])
				st_corr.append(split_line[4]) # but really std dev

	else: # no static terms, always true for NORSAR
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

#%%  Set some parameters
#	fig_index = 101
#	stat_corr = 1 # apply station static corrections
#	rel_time = 1          # timing is relative to a chosen phase, otherwise relative to OT
#	dphase  = 'PKIKP'       # phase to be aligned
#	dphase2 = 'PKiKP'      # another phase to have traveltime plotted
#	dphase3 = 'PKP'        # another phase to have traveltime plotted
#	dphase4 = 'pP'        # another phase to have traveltime plotted
	taper_frac = .05      #Fraction of window tapered on both ends
	signal_dur = 10.       # signal length used in SNR calculation
#	plot_scale_fac = 0.5    #  Bigger numbers make each trace amplitude bigger on plot
#	qual_threshold =  2 # minimum SNR
#	corr_threshold = 0.7  # minimum correlation in measuring shift to use station
	plot_tt = 1           # plot the traveltimes?
#	ref_loc = 0  # 1 if selecting stations within ref_rad of ref_lat and ref_lon
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

#%% Is taper too long compared to noise estimation window?
	totalt = end_buff - start_buff
	noise_time_skipped = taper_frac * totalt
	if simple_taper == 0:
		if noise_time_skipped >= -0.5 * start_buff:
			print('Specified taper of ' + str(taper_frac * totalt) +
			   ' is not big enough compared to available noise estimation window ' +
			   str(-start_buff - noise_time_skipped) + '. May not work well.')
			old_taper_frac = taper_frac
			taper_frac = -0.5*start_buff/totalt
			if start_buff > 0:
					taper_frac = 0.05 # pick random minimal window if there is no leader
			print('Taper reset from ' + str(old_taper_frac * totalt) + ' to '
			   + str(taper_frac * totalt) + ' seconds.')

	if rel_time == 0: # SNR requirement not implemented for unaligned traces
		qual_threshold = 0

	# Plot with reduced velocity?
	red_plot = 0
	red_dist = 55
	red_time = 300
	red_slow = 7.2 # seconds per degree

#%% Load waveforms and decimate to 10 sps
	st = Stream()
	if ARRAY == 0 or ARRAY == 2:
		fname     = '/Users/vidale/Documents/PyCode/Mseed/HD' + date_label + '.mseed'
	if ARRAY == 1:
		fname     = '/Users/vidale/Documents/PyCode/Mseed/L' + date_long[2:4] + date_long[5:7] + date_long[8:10] + '_' + date_long[11:13] + date_long[14:16] + '.mseed'
	st=read(fname)
	if do_decimate != 0:
		st.decimate(do_decimate, no_filter=True)

	print(fname)
	print('Read in: ' + str(len(st)) + ' traces')
	print('First trace has : ' + str(len(st[0].data)) + ' time pts ')
	print('Start time : ' + str(st[0].stats.starttime) + '  event time : ' + str(t))
	print('After decimation: ' + str(len(st)) + ' traces')
	nt = len(st[0].data)
	dt = st[0].stats.delta
	print('First trace has : ' + str(nt) + ' time pts, time sampling of '
		  + str(dt) + ' and thus duration of ' + str((nt-1)*dt))
	#	print(f'Sta lat-lon {stalat:.4f}  {stalon:.4f}')

#%% Select by distance, window and adjust start time to align picked times
	st_pickalign = Stream()

	tra_in_range  = 0
	tra_sta_found = 0
	nodata        = 0
	min_dist_auto = 180
	max_dist_auto = 0
	min_time_plot =  1000000
	max_time_plot = -1000000

	# not used in all cases, but printed out below
	# only used if rel_slow == 1, preserves 0 slowness, otherwise 0 is set to phase slowness
	ref_distance = gps2dist_azimuth(ref_lat,ref_lon,ev_lat,ev_lon)
	ref1_dist  = ref_distance[0]/(1000*111)
	dist_minus = ref1_dist - 0.5
	dist_plus  = ref1_dist + 0.5
	arrivals_ref   = model.get_travel_times(source_depth_in_km=ev_depth,distance_in_degree=ref1_dist, phase_list=[dphase])
	arrivals_minus = model.get_travel_times(source_depth_in_km=ev_depth,distance_in_degree=dist_minus,phase_list=[dphase])
	arrivals_plus  = model.get_travel_times(source_depth_in_km=ev_depth,distance_in_degree=dist_plus ,phase_list=[dphase])
	atime_ref = arrivals_ref[0].time  # phase arrival time at reference distance
	ref_slow = arrivals_plus[0].time - arrivals_minus[0].time  # dt over 1 degree at ref distance

	for tr in st: # traces one by one, find lat-lon
		if float(year) < 1970: # fix the damn 1969 -> 2069 bug in Gibbon's LASA data
			temp_t = str(tr.stats.starttime)
			temp_tt = '19' + temp_t[2:]
			tr.stats.starttime = UTCDateTime(temp_tt)
		if JST == 1: # if necessary, convert JST -> UTC, time in Greenwich 9 hours earlier than Japan
			tr.stats.starttime = tr.stats.starttime - 9*60*60
#			temp_t = str(tr.stats.starttime)
#			temp_tt = '19' + temp_t[2:]
#			tr.stats.starttime = UTCDateTime(temp_tt)
#			times.append(   tr.stats.starttime.strptime(split_line[1], dtformat) - dt.timedelta(seconds=offset))
#			timesUTC.append(dt.datetime.strptime(split_line[1], dtformat)) # keep UTC version for output
		if tr.stats.station in st_names:  # find station in station list
			ii = st_names.index(tr.stats.station)
			tra_sta_found += 1

			if stat_corr != 1 or float(st_corr[ii]) > corr_threshold: # if using statics, reject low correlations
				stalat = float(st_lats[ii]) # look up lat & lon again to find distance
				stalon = float(st_lons[ii])

				distance = gps2dist_azimuth(stalat,stalon,ev_lat,ev_lon) # Get traveltimes again, hard to store
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
					s_t = t + start_buff
					e_t = t + end_buff
					if stat_corr == 1: # apply static station corrections
						tr.stats.starttime -= float(st_shift[ii])
					if rel_time == 0:  #  don't adjust absolute time
						tr.trim(starttime=s_t,endtime = e_t)
					else:              # shift relative to a chosen phase
						arrivals_each = model.get_travel_times(source_depth_in_km=ev_depth,distance_in_degree=dist,phase_list=[dphase])
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
						st_pickalign += tr
					else:
						nodata += 1
		else:
			print(tr.stats.station + ' not found in station list')
	print('After alignment + range and correlation selection - event: ' + str(len(st_pickalign)) + ' traces')
	print('Traces found: ' + str(tra_sta_found) + ' Traces in range: ' + str(tra_in_range) + ' Traces with no data: ' + str(nodata))
	print(f'ref1_distance  {ref1_dist:.3f}  relative start time  {atime_ref:.3f}')
	print('ref_loc == 1, ref_lat: ' + str(ref_lat) + ' ref_lon: ' + str(ref_lon))
	print(f'last station: distance {dist:.3f}  last station lat: {stalat:.3f}   last station lon: {stalat:.3f}')

	#print(st) # at length
	if verbose:
		print(st.__str__(extended=True))
		if rel_time == 1:
			print(st_pickalign.__str__(extended=True))

#%%  Detrend, taper, filter
	st_pickalign.detrend(type='simple')
	print('taper_frac is ' + str(taper_frac))
	st_pickalign.taper(taper_frac)
	if do_filt == 1:
		st_pickalign.filter('bandpass', freqmin=freq_min, freqmax=freq_max, corners=4, zerophase=True)
	st_pickalign.taper(taper_frac)

#%%  Cull further by imposing SNR threshold on both traces
	if skip_SNR == 1:
		stgood = st_pickalign.copy()
	else:
		stgood = Stream()
		for tr in st_pickalign:
		# estimate median noise
			t_noise_start  = int(len(tr.data) * taper_frac)
			t_noise_end    = int(len(tr.data) * start_buff/(start_buff-end_buff))
			noise          = np.median(abs(tr.data[t_noise_start:t_noise_end]))
		# estimate median signal
			t_signal_start = int(len(tr.data) * start_buff/(start_buff-end_buff))
			t_signal_end   = t_signal_start + int(len(tr.data) * signal_dur/(end_buff - start_buff))
			signal         = np.median(abs(tr.data[t_signal_start:t_signal_end]))
		#			test SNR
			SNR = signal/noise;
			if (SNR > qual_threshold):
				stgood += tr

	print('Above SNR threshold: ' + str(len(stgood)) + ' traces')
	if verbose:
		for tr in stgood:
			print('Distance is ' + str(tr.stats.distance/(1000*111)) + ' for station ' + tr.stats.station)

	#%%  get station lat-lon, compute distance and time limits for plot
	for tr in stgood:
		if tr.stats.station in st_names:  # find station in station list
			ii = st_names.index(tr.stats.station)
			stalon = float(st_lons[ii]) # look up lat & lon again to find distance
			stalat = float(st_lats[ii])
			distance = gps2dist_azimuth(stalat,stalon,ev_lat,ev_lon)
			tr.stats.distance=distance[0]/(1000*111) # distance in degrees
			if tr.stats.distance < min_dist_auto:
				min_dist_auto = tr.stats.distance
			if tr.stats.distance > max_dist_auto:
				max_dist_auto = tr.stats.distance
			if tr.stats.starttime - t < min_time_plot:
				min_time_plot = tr.stats.starttime - t
			if ((tr.stats.starttime - t) + ((len(tr.data)-1) * tr.stats.delta)) > max_time_plot:
				max_time_plot =  ((tr.stats.starttime - t) + ((len(tr.data)-1) * tr.stats.delta))
	print(f'Min distance is   {min_dist_auto:.3f}   Max distance is {max_dist_auto:.3f}')
	print(f'Min time is   {min_time_plot:.2f}   Max time is {max_time_plot:.2f}')

#%%  This section causes a crash in Spyder
	# plot traces
	plt.close(fig_index)
	plt.figure(fig_index,figsize=(10,10))
	plt.xlim(min_time_plot,max_time_plot)

	if auto_dist == 1:
		dist_diff = max_dist_auto - min_dist_auto # add space at extremes
		plt.ylim(min_dist_auto - 0.1 * dist_diff, max_dist_auto + 0.1 * dist_diff)
	else:
		plt.ylim(min_dist,max_dist)

	for tr in stgood:
		dist_offset = tr.stats.distance
		ttt = np.arange(len(tr.data)) * tr.stats.delta + (tr.stats.starttime - t)
		if red_plot == 1:
			shift = red_time + (dist_offset - red_dist) * red_slow
			ttt = ttt - shift
		if len(tr.data) > 0:
			if tr.data.max() - tr.data.min() > 0:
				plt.plot(ttt, (tr.data - np.median(tr.data))*plot_scale_fac /(tr.data.max() - tr.data.min()) + dist_offset, color = 'black')
			else:
				print('Max ' + str(tr.data.max()) + ' equals min ' + str(tr.data.min()) + ', skip plotting')
		else:
			nodata += 1
			print('Trace ' + tr.stats.station + ' has : ' + str(len(tr.data)) + ' time pts, skip plotting')
#%% Plot traveltime curves
	if rel_time != 1:
		if plot_tt:
			# first traveltime curve
			line_pts = 50
			dist_vec  = np.arange(min_dist_auto, max_dist_auto, (max_dist_auto - min_dist_auto)/line_pts) # distance grid
			time_vec1 = np.arange(min_dist_auto, max_dist_auto, (max_dist_auto - min_dist_auto)/line_pts) # empty time grid of same length (filled with -1000)
			for i in range(0,line_pts):
				arrivals = model.get_travel_times(source_depth_in_km=ev_depth,distance_in_degree
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
					arrivals = model.get_travel_times(source_depth_in_km=ev_depth,distance_in_degree
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
					arrivals = model.get_travel_times(source_depth_in_km=ev_depth,distance_in_degree
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
					arrivals = model.get_travel_times(source_depth_in_km=ev_depth,distance_in_degree
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
	plt.title(dphase + ' for ' + date_label + ' event # ' + str(event_no))
	os.chdir('/Users/vidale/Documents/PyCode/Plots')
#	plt.savefig(date_label + '_' + str(event_no) + '_raw.png')
	plt.show()

#%%  Save processed files
	fname3 = '/Users/vidale/Documents/PyCode/Pro_Files/HD' + date_label + 'sel.mseed'

	stgood.write(fname3,format = 'MSEED')

	elapsed_time_wc = time.time() - start_time_wc
	print(f'This job took   {elapsed_time_wc:.1f}   seconds')
	os.system('say "Done"')