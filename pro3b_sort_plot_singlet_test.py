#!/usr/bin/env python
# Add three plane waves to test beam resolution
# This programs deals with a single event.
# John Vidale 3/2019

def pro3singlet(eq_file, stat_corr = 0,
			dphase = 'PKIKP', dphase2 = 'PKiKP', dphase3 = 'PKIKP', dphase4 = 'PKiKP',
			start_buff = 10, end_buff = 30,
			plot_scale_fac = 0.05, qual_threshold = 0, corr_threshold = 0,
			freq_min = 0.25, freq_max = 1,
			min_dist = 0, max_dist = 180, do_decimate = 0,
			alt_statics = 0, statics_file = 'nothing', ARRAY = 0, ref_loc = 0,
			verbose = 0):
# 0 is Hinet, 1 is LASA, 2 is NORSAR

	from obspy import UTCDateTime
	from obspy import Stream
	from obspy import read
	from obspy.geodetics import gps2dist_azimuth
	import numpy as np
	import os
	from obspy.taup import TauPyModel
	import matplotlib.pyplot as plt
	import time
	import math
	model = TauPyModel(model='iasp91')

	import sys # don't show any warnings
	import warnings

	if not sys.warnoptions:
	    warnings.simplefilter("ignore")

	start_time_wc = time.time()

	#%% input event data with 1-line file of format
	#  event 2016-05-28T09:47:00.000 -56.241 -26.935 78
	file = open(eq_file, 'r')
	lines=file.readlines()
	split_line = lines[0].split()
#			ids.append(split_line[0])  ignore label for now
	t           = UTCDateTime(split_line[1])
	date_label  = split_line[1][0:10]
	ev_lat      = float(      split_line[2])
	ev_lon      = float(      split_line[3])
	ev_depth    = float(      split_line[4])
	print('date_label ' + date_label + ' time ' + str(t) + ' lat ' + str(ev_lat) + ' lon ' + str( ev_lon) + ' depth ' + str(ev_depth))

	#%% Get Hinet, LASA, or NORSAR station location file
	sta_file = '/Users/vidale/Documents/GitHub/Hinet-codes/LASA_sta.txt'
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

	#%%
	fig_index = 101
	rel_time = 1          # timing is relative to a chosen phase, otherwise relative to OT
	taper_frac = .05      #Fraction of window tapered on both ends
	plot_tt = 1           # plot the traveltimes?
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
	totalt = start_buff + end_buff
	noise_time_skipped = taper_frac * totalt
	if noise_time_skipped >= 0.5 * start_buff:
		print('Specified taper of ' + str(taper_frac * totalt) +
		   ' is not big enough compared to available noise estimation window ' +
		   str(start_buff - noise_time_skipped) + '. May not work well.')
		old_taper_frac = taper_frac
		taper_frac = 0.5*start_buff/totalt
		print('Taper reset from ' + str(old_taper_frac * totalt) + ' to '
		   + str(taper_frac * totalt) + ' seconds.')

	#%% Load waveforms and decimate to 10 sps
	st = Stream()
	fname     = 'HD' + date_label + '.mseed'
	st=read(fname)
	if do_decimate != 0:
		st.decimate(do_decimate)

	print('Read in: ' + str(len(st)) + ' traces')
	print('First trace has : ' + str(len(st[0].data)) + ' time pts ')
	print('After decimation: ' + str(len(st)) + ' traces')
	nt = len(st[0].data)
	dt = st[0].stats.delta
	print('First trace has : ' + str(nt) + ' time pts, time sampling of '
		  + str(dt) + ' and thus duration of ' + str((nt-1)*dt))

	#%%
	# select by distance, window and adjust start time to align picked times
	st_pickalign = Stream()

	for tr in st: # traces one by one, find lat-lon by searching entire inventory.  Inefficient
		for ii in station_index:
			if (tr.stats.station == st_names[ii]): # find station in inventory
				stalat = float(st_lats[ii])
				stalon = float(st_lons[ii]) # look up lat & lon again to find distance
				if ref_loc == 1:
					ref_distance = gps2dist_azimuth(stalat,stalon,ref_lat,ref_lon)
					ref_dist = ref_distance[0]/(1000*111)
				distance = gps2dist_azimuth(stalat,stalon,ev_lat,ev_lon) # Get traveltimes again, hard to store
				tr.stats.distance=distance[0] # distance in km
				dist = distance[0]/(1000*111)
				if ref_loc != 1 and min_dist < dist and dist < max_dist: # select distance range from earthquake
					try:
						arrivals = model.get_travel_times(source_depth_in_km=ev_depth,distance_in_degree=dist,phase_list=[dphase])
						atime = arrivals[0].time
						s_t = t + atime - start_buff
						e_t = t + atime + end_buff
						tr.trim(starttime=s_t,endtime = e_t)
						# deduct theoretical traveltime and start_buf from starttime
						tr.stats.starttime -= atime
						st_pickalign += tr
					except:
						pass
				elif ref_loc == 1:
					if ref_dist < ref_rad: # alternatively, select based on distance from ref location
						try:
							arrivals = model.get_travel_times(source_depth_in_km=ev_depth,distance_in_degree=dist,phase_list=[dphase])
							atime = arrivals[0].time
							s_t = t + atime - start_buff
							e_t = t + atime + end_buff
							tr.trim(starttime=s_t,endtime = e_t)
							# deduct theoretical traveltime and start_buf from starttime
							tr.stats.starttime -= atime
							st_pickalign += tr
						except:
							pass

				for it in range(20):  # add three strong signals
#					tr.data[it+100] += 100000*math.sin(it*(180/10.))
#					tr.data[it+360 - int(35*(dist-60.))] += 100000*math.sin(it*(180/5.))
#					tr.data[it+260 + int(35*(dist-60.))] += 100000*math.sin(it*(180/8.))
					tr.data[it+100] += 100000*math.sin(it*(math.pi/10.))
					tr.data[it+360 - int(35*(dist-60.))] += 100000*math.sin(it*2*math.pi/10.)
					tr.data[it+260 + int(35*(dist-60.))] += 100000*math.sin(it*4*math.pi/10.)
	print('After alignment + range and correlation selection - event: ' + str(len(st_pickalign)) + ' traces')

	#%%  Cull further by imposing SNR threshold on both traces
	stgood = Stream()
	for tr in st_pickalign:
		stgood += tr
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

	#%%
	# plot traces
	plt.close(fig_index)
	plt.figure(fig_index,figsize=(10,10))
	plt.xlim(-start_buff,end_buff)
	plt.ylim(min_dist,max_dist)
	for tr in stgood:
		dist_offset = tr.stats.distance/(1000*111) # trying for approx degrees
		ttt = np.arange(len(tr.data)) * tr.stats.delta + (tr.stats.starttime - t)
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

		if rel_time == 1:
			time_vec1 = time_vec1 - time_vec1
		plt.plot(time_vec1,dist_vec, color = 'blue')
		plt.show()

	plt.xlabel('Time (s)')
	plt.ylabel('Epicentral distance from event (°)')
	plt.title(dphase + ' for ' + fname[2:12])
	plt.show()

	#  Save processed files
	fname3 = 'HD' + date_label + 'sel.mseed'
	stgood.write(fname3,format = 'MSEED')

	elapsed_time_wc = time.time() - start_time_wc
	print('This job took ' + str(elapsed_time_wc) + ' seconds')
	os.system('say "Done"')