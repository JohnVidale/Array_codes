#!/usr/bin/env python
# reads in raw mseed traces from a single event
# this program tapers, filters, selects range and SNR
# plots against traveltime curves, either raw or reduced against traveltimes
# writes out "*sel.mseed" file
# plot lines are blue, orange, yellow, purple for phases 1 through 4
# John Vidale 2/2019

def pro3singlet(eq_file, stat_corr = 0, rel_time = 1, simple_taper = 0, skip_SNR = 0,
			dphase = 'PKIKP', dphase2 = 'PKiKP', dphase3 = 'PKIKP', dphase4 = 'PKiKP',
			start_buff = -10, end_buff = 30,
			plot_scale_fac = 0.05, qual_threshold = 0, corr_threshold = 0,
			freq_min = 0.25, freq_max = 1,
			min_dist = 0, max_dist = 180, do_decimate = 0,
			alt_statics = 0, statics_file = 'nothing', ARRAY = 0, ref_loc = 0,
			verbose = 0, fig_index = 101):
# 0 is Hinet, 1 is LASA, 2 is NORSAR

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
	if ARRAY == 0:
		file = open(eq_file, 'r')
	elif ARRAY == 1:
		file = open('EvLocs/' + eq_file, 'r')
	lines=file.readlines()
	split_line = lines[0].split()
#			ids.append(split_line[0])  ignore label for now
	t           = UTCDateTime(split_line[1])
	date_label  = split_line[1][0:10]
	year        = split_line[1][0:4]
	ev_lat      = float(      split_line[2])
	ev_lon      = float(      split_line[3])
	ev_depth    = float(      split_line[4])
	print('date_label ' + date_label + ' time ' + str(t) + ' lat ' + str(ev_lat) + ' lon ' + str( ev_lon) + ' depth ' + str(ev_depth))

#%% Get station location file
	if stat_corr == 1:  # load static terms, only applies to Hinet and LASA
		if ARRAY == 0:
			if alt_statics == 0: # standard set
				sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/hinet_sta_statics.txt'
			else: # custom set made by this event for this event
				sta_file = ('/Users/vidale/Documents/PyCode/Hinet/Array_codes/Files/' + 'HA' +
				   date_label[:10] + 'pro4_' + dphase + '.statics')
		elif ARRAY == 1:
			sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/L_sta_statics.txt'
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
			sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/hinet_sta.txt'
		elif ARRAY == 1: #         LASA set
			sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/LASA_sta.txt'
		else: #         NORSAR set
			sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/NORSAR_sta.txt'
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
	if ref_loc == 1:
		if ARRAY == 0:
			ref_lat = 36.3  # °N, around middle of Japan
			ref_lon = 138.5 # °E
			ref_rad = 1.5   # ° radius (°)
		elif ARRAY == 1:
			ref_lat = 46.7  # °N keep only inner rings A-D
			ref_lon = -106.22   # °E
			ref_rad = 0.4    # ° radius (°)

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

	#%% In case one wants to manually enter data here
#	date_label = '2018-04-02' # date for filename
#	ev_lon   = -63.006
#	ev_lat   = -20.659
#	ev_depth = 559
#	t        = UTCDateTime('2018-04-02T13:40:34.840')

#%% Load waveforms and decimate to 10 sps
	st = Stream()
	if ARRAY == 0:
		fname     = 'HD' + date_label + '.mseed'
	elif ARRAY == 1:
		fname     = 'Mseed/HD' + date_label + '.mseed'
	st=read(fname)
	if do_decimate != 0:
		st.decimate(do_decimate)

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

	tra_located   = 0
	tra_in_range  = 0
	tra_sta_found = 0
#	for ii in station_index:
#		print('Station name of index ' + str(ii) + ' is ' + str(st_names[ii])) # enumerate stations
#	for tr in st: # traces one by one, find lat-lon by searching entire inventory.  Inefficient
#		print('Station name of tr ' + str(tr.stats.station)) # enumerate stations
	for tr in st: # traces one by one, find lat-lon by searching entire inventory.  Inefficient
		if float(year) < 1970: # fix the damn 1969 -> 2069 bug in Gibbon's LASA data
			temp_t = str(tr.stats.starttime)
			temp_tt = '19' + temp_t[2:]
			tr.stats.starttime = UTCDateTime(temp_tt)
		for ii in station_index:
			if ARRAY == 0:  # have to chop off last letter, always 'h'
				this_name = st_names[ii]
				this_name_truc = this_name[0:5]
				name_truc_cap  = this_name_truc.upper()
			elif ARRAY == 1:
				name_truc_cap = st_names[ii]
			if (tr.stats.station == name_truc_cap): # find station in inventory
				tra_sta_found += 1
				if stat_corr != 1 or float(st_corr[ii]) > corr_threshold: # if using statics, reject low correlations
					stalat = float(st_lats[ii])
					stalon = float(st_lons[ii]) # look up lat & lon again to find distance
					if ref_loc == 1:
						ref_distance = gps2dist_azimuth(stalat,stalon,ref_lat,ref_lon)
						ref_dist = ref_distance[0]/(1000*111)
					distance = gps2dist_azimuth(stalat,stalon,ev_lat,ev_lon) # Get traveltimes again, hard to store
					tr.stats.distance=distance[0] # distance in km
					dist = distance[0]/(1000*111)
					if ref_loc != 1:
						tra_located += 1
						if min_dist < dist and dist < max_dist: # select distance range from earthquake
							tra_in_range += 1
							try:
								arrivals = model.get_travel_times(source_depth_in_km=ev_depth,distance_in_degree=dist,phase_list=[dphase])
								atime = arrivals[0].time
								if stat_corr == 1: # apply static station corrections
									tr.stats.starttime -= float(st_shift[ii])
								if rel_time == 1:
									s_t = t + atime + start_buff
									e_t = t + atime + end_buff
								else:
									s_t = t + start_buff
									e_t = t + end_buff
								tr.trim(starttime=s_t,endtime = e_t)
								# deduct theoretical traveltime and start_buf from starttime
								if rel_time == 1:
									tr.stats.starttime -= atime
								st_pickalign += tr
							except:
								pass
					elif ref_loc == 1:
						if ref_dist < ref_rad: # alternatively, select based on distance from ref location
							try:
								arrivals = model.get_travel_times(source_depth_in_km=ev_depth,distance_in_degree=dist,phase_list=[dphase])
								atime = arrivals[0].time
								if stat_corr == 1: # apply static station corrections
									tr.stats.starttime -= float(st_shift[ii])
								if rel_time == 1:
									s_t = t + atime + start_buff
									e_t = t + atime + end_buff
								else:
									s_t = t + start_buff
									e_t = t + end_buff
								tr.trim(starttime=s_t,endtime = e_t)
								# deduct theoretical traveltime and start_buf from starttime
								if rel_time == 1:
									tr.stats.starttime -= atime
								st_pickalign += tr
							except:
								pass
	print('After alignment + range and correlation selection - event: ' + str(len(st_pickalign)) + ' traces')
	print('Traces found: ' + str(tra_sta_found) + ' traces with range examined: ' + str(tra_located) + ' traces in range: ' + str(tra_in_range))

	#print(st) # at length
	if verbose:
		print(st.__str__(extended=True))
		if rel_time == 1:
			print(st_pickalign.__str__(extended=True))

#%%  Detrend, taper, filter
	st_pickalign.detrend(type='simple')
	print('taper_frac is ' + str(taper_frac))
	st_pickalign.taper(taper_frac)
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

	#%%  get station lat-lon, compute distance for plot
	for tr in stgood:
		for ii in station_index:
			if (tr.stats.station == st_names[ii]): # find station in inventory
				stalon = float(st_lons[ii]) # look up lat & lon again to find distance
				stalat = float(st_lats[ii])
				distance = gps2dist_azimuth(stalat,stalon,ev_lat,ev_lon)
				tr.stats.distance=distance[0] # distance in km

#%%  This section causes a crash in Spyder
	# plot traces
	plt.close(fig_index)
	plt.figure(fig_index,figsize=(10,10))
	plt.xlim(start_buff,end_buff)
	plt.ylim(min_dist,max_dist)
	for tr in stgood:
		dist_offset = tr.stats.distance/(1000*111) # trying for approx degrees
		ttt = np.arange(len(tr.data)) * tr.stats.delta + (tr.stats.starttime - t)
		if red_plot == 1:
			shift = red_time + (dist_offset - red_dist) * red_slow
			ttt = ttt - shift
		plt.plot(ttt, (tr.data - np.median(tr.data))*plot_scale_fac /(tr.data.max()
			- tr.data.min()) + dist_offset, color = 'black')
#%% Plot traveltime curves
	if plot_tt:
		# first traveltime curve
		line_pts = 50
		dist_vec  = np.arange(min_dist, max_dist, (max_dist - min_dist)/line_pts) # distance grid
		time_vec1 = np.arange(min_dist, max_dist, (max_dist - min_dist)/line_pts) # empty time grid of same length (filled with -1000)
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
			time_vec2 = np.arange(min_dist, max_dist, (max_dist - min_dist)/line_pts) # empty time grid of same length (filled with -1000)
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
			if rel_time == 1:
				time_vec2 = time_vec2 - time_vec1
			plt.plot(time_vec2,dist_vec, color = 'orange')
		# third traveltime curve
		if dphase3 != 'no':
			time_vec3 = np.arange(min_dist, max_dist, (max_dist - min_dist)/line_pts) # empty time grid of same length (filled with -1000)
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
			if rel_time == 1:
				time_vec3 = time_vec3 - time_vec1
			plt.plot(time_vec3,dist_vec, color = 'yellow')
		# fourth traveltime curve
		if dphase4 != 'no':
			time_vec4 = np.arange(min_dist, max_dist, (max_dist - min_dist)/line_pts) # empty time grid of same length (filled with -1000)
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
			if rel_time == 1:
				time_vec4 = time_vec4 - time_vec1
			plt.plot(time_vec4,dist_vec, color = 'purple')

		if rel_time == 1:
			time_vec1 = time_vec1 - time_vec1
		plt.plot(time_vec1,dist_vec, color = 'blue')
		plt.show()

	plt.xlabel('Time (s)')
	plt.ylabel('Epicentral distance from event (°)')
	if ARRAY == 0:
		plt.title(dphase + ' for ' + fname)
	elif ARRAY == 1:
		plt.title(dphase + ' for ' + fname[8:18])
	plt.show()

#%%  Save processed files
	if ARRAY == 0:
		fname3 = 'HD' + date_label + 'sel.mseed'
	elif ARRAY == 1:
		fname3 = 'Pro_Files/HD' + date_label + 'sel.mseed'

	stgood.write(fname3,format = 'MSEED')

	elapsed_time_wc = time.time() - start_time_wc
	print('This job took ' + str(elapsed_time_wc) + ' seconds')
	os.system('say "Done"')