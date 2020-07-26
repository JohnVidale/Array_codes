#!/usr/bin/env python
#  reads in "*sel.mseed" for an event
#  align traces by shifting, and record time shifts for station statics
#  plots traces before and after time shift
#  saves aligned traces (not generally used) and static corrections, used in pro3 codes
#  John Vidale, 2/2019
#   revisited vy Vidale 7/2020

def pro4statics(eq_file, use_ref_trace = 0, ref_trace = 'nothing',
				dphase = 'PcP', dphase2 = 'PKiKP', dphase3 = 'P', dphase4 = 'PP',
				start_corr_win = -1, end_corr_win = 3, plot_scale_fac = 0.05,
				start_buff = -10, end_buff = 30, qual_threshold = 0, corr_threshold = 0,
				max_time_shift = 2, min_dist = 17, max_dist = 21, ARRAY = 0):

	from obspy import UTCDateTime
	from obspy.signal.cross_correlation import xcorr_pick_correction
	from obspy import Stream
	from obspy import Trace
	from obspy import read
	from obspy.geodetics import gps2dist_azimuth
	import numpy as np
	import os
	import sys
	from obspy.taup import TauPyModel
	import matplotlib.pyplot as plt
	model = TauPyModel(model='iasp91')

# 	import warnings   # don't show any warnings
# 	if not sys.warnoptions:
# 	    warnings.simplefilter("ignore")
	print('    Running pro4')

	#%% Get station location file
	if   ARRAY == 0: # Hinet set
		sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/sta_hinet.txt'
	elif ARRAY == 1: # LASA set
		sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/sta_LASA.txt'
	elif ARRAY == 2: # China set
		sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/sta_ch.txt.txt'
	with open(sta_file, 'r') as file:
		lines = file.readlines()
	print(str(len(lines)) + ' stations read from ' + sta_file)
	# Load station coords into arrays
		# old line: station_index = range(343)
	station_index = range(len(lines))
	st_lats  = []
	st_lons  = []
	st_deps  = []
	st_names = []
	for ii in station_index:
		line = lines[ii]
		split_line = line.split()
		st_names.append(split_line[0])
		st_lats.append( split_line[1])
		st_lons.append( split_line[2])
		st_deps.append( split_line[3])

	# initialize lists of statics
	sta_names   = []
	sta_dists   = []
	sta_lats    = []
	sta_lons    = []
	sta_statics = []
	sta_corrs   = []

	#%% Parameter list
	#dphase  = 'PKIKP'       # phase to be aligned
	#dphase2 = 'PKiKP'      # another phase to have traveltime plotted
	#dphase3 = 'PKP'        # another phase to have traveltime plotted
	#dphase4 = 'pP'        # another phase to have traveltime plotted
	#ref_trace = 'N.SZW'   # trace with reference waveform
	#start_corr_win = 2       # plots start Xs before PKiKP
	#end_corr_win   = 7       # plots end Xs before PKiKP
	#max_time_shift = 2       # searches up to this time shift for alignment
	start_plot_win = start_buff       # plots start Xs before PKiKP
	end_plot_win   =   end_buff       # plots end Xs before PKiKP
	#corr_threshold = 0.  # threshold that correlation is good enough to keep trace
	#max_dist = 151
	#min_dist = 150.6
	#plot_scale_fac = 0.2    #  Bigger numbers make each trace amplitude bigger on plot
	#qual_threshold =  0 # minimum SNR
	plot_tt = True           # plot the traveltimes?
	plot_flag = False     # plot for each trace?  Watch out, can be lots, one for each station pair!!

	#%% Get saved event info, also used to name files
	#  event 2016-05-28T09:47:00.000 -56.241 -26.935 78
	file = open('/Users/vidale/Documents/PyCode/EvLocs/' + eq_file, 'r')
	lines=file.readlines()
	split_line = lines[0].split()
#			ids.append(split_line[0])  ignore label, now "event"
	t           = UTCDateTime(split_line[1])
	date_label  = split_line[1][0:10]
	ev_lat      = float(      split_line[2])
	ev_lon      = float(      split_line[3])
	ev_depth    = float(      split_line[4])

	print('Date label ' + date_label + ' lat ' + str(ev_lat) + ' lon ' + str(ev_lon))

	st = Stream()
#	fname     = 'HD' + date_label + '.mseed'
	fname     = 'HD' + date_label + 'sel.mseed'  # sel file has windowing, shift?, filtering

	print('fname ' + fname)

	os.chdir('/Users/vidale/Documents/PyCode/Pro_Files/')
	os.system('pwd')
	st=read(fname)
	print('Read in: ' + str(len(st)) + ' traces')
	print('First trace has : ' + str(len(st[0].data)) + ' time pts ')

	#%% Reference trace
	trim_start = t + start_buff
	trim_end   = t +   end_buff
	time_buff = end_buff - start_buff
	tr_ref = Trace()
	#%% Stack reference trace
	if use_ref_trace == 0:
		counter = 0
		for tr in st: # loop over seismograms to find reference trace, put it in tr_ref
			if counter == 0:  # copy first trace to stack
				tr_ref = tr.copy()
				tr_ref.stats.station = 'STACK'
				tr_ref.trim(starttime = trim_start - time_buff, endtime = trim_end + time_buff)
				nt_ref = len(tr_ref.data)
				tr_ref.normalize()
				counter = counter + 1
			else:         # add the rest of the traces to stack
				tr_add = tr.copy()
				tr_add.trim(starttime = trim_start - time_buff, endtime = trim_end + time_buff)
				nt_add = len(tr_ref.data)
				tr_add.normalize()

				for it in range(nt_ref):  # add seismogram one point at a time
					if nt_ref != nt_add: # are seismograms the same length?
						print('trying to stack seismograms of different lengths, debug!')
					tr_ref.data[it] += tr_add[it]
				counter = counter + 1
		tr_ref.data = tr_ref.data/counter

	#%% Pick reference trace
	if use_ref_trace == 1:
		for tr in st: # loop over seismograms to find reference trace, put it in tr_ref
			if (tr.stats.station == ref_trace): # found it
				tr_ref = tr.copy()
				tr_ref.trim(starttime = trim_start - time_buff, endtime = trim_end + time_buff)
				nt_ref = len(tr_ref.data)
				tr_ref.normalize()
				print('found reference station ' + tr.stats.station)

	#%% Plot reference trace
	plt.close(4)
	plt.figure(4,figsize=(10,10))
	plt.xlim(start_corr_win - time_buff, end_corr_win + time_buff)
	plt.ylim(min(tr_ref.data), max(tr_ref.data))

	time = np.arange(nt_ref) * tr_ref.stats.delta + start_corr_win - time_buff
	plt.plot(time, tr_ref.data, color = 'black')


	stgood = Stream()
	st2 = st.copy()  # hard to measure timing of traces without adjusting entire thing
	print('st2 has: ' + str(len(st)) + ' traces' + ' t (origin time) ' + str(t))
	print('t ' + str(t) + ' start_corr_win end_corr_win max_time_shift ' + str(start_corr_win) + '  ' + str(end_corr_win) + '  ' + str(max_time_shift) + '  ')

	#  get station lat-lon, compute distance for plot
	good_corr = 0
	bad_corr = 0
# 	sys.exit()
	for tr in st: # do all seismograms
		for ii in station_index: # find station in inventory
			tested_name = st_names[ii]
			actual_trace = tr.stats.station
			if ARRAY == 0: # convoluted patch for long Hinet names
				this_name_truc = tested_name[0:5]
				name_truc_cap  = this_name_truc.upper()
				this_name = name_truc_cap
				actual_trace = tr.stats.station.upper()
			if (actual_trace == this_name): # found it
# 				print('found Station ' + this_name + '  ' + actual_trace)
				stalon = float(st_lons[ii]) # look up lat & lon to find distance
				stalat = float(st_lats[ii])
				distance = gps2dist_azimuth(stalat,stalon,ev_lat,ev_lon)
				tr.stats.distance=distance[0]/(1000.*111) # distance for phase time and plotting
				arrivals = model.get_travel_times(source_depth_in_km=ev_depth,distance_in_degree
									=tr.stats.distance,phase_list=[dphase])
# 				print(tr.stats.station + '  ' + tr_ref.stats.station + ' start_corr ' +
# 				   str(start_corr_win) + ' end ' + str(end_corr_win))
				try:
					dt, coeff = xcorr_pick_correction(t, tr_ref, t, tr, start_corr_win, end_corr_win, max_time_shift, plot=plot_flag)
					if dt > max_time_shift:
						print('Hey!  Excess shift: %.3f' % dt)
						print('Station ' + tr.stats.station + ' corr is ' + str(coeff))
					if coeff > 1:
						print('Hey!  Excess coeff: %.3f' % coeff)
						print('Station ' + tr.stats.station + ' corr is ' + str(coeff))
					if coeff > corr_threshold:
						good_corr += 1
						if plot_flag == True:
							print('Time correction for pick 2: %.6f' % dt)
							print('Correlation coefficient: %.2f' % coeff)
						tr.stats.starttime -= dt
						sta_names.extend([tr.stats.station])
						sta_dists.extend([tr.stats.distance])
						sta_lats.extend([stalat])
						sta_lons.extend([stalon])
						sta_statics.extend([dt])
						sta_corrs.extend([coeff])
						stgood += tr
					else:
						bad_corr += 1
				except:
					print('No time shift for ' + tr.stats.station + ' at distance ' + str(tr.stats.distance))

	##		# store shift to write out
	##			if coeff > corr_threshold:
	#			# write out station_name, dt, coeff
	#			# record shifted waveform in stgood
	print('Good ' + str(good_corr) + ' bad ' + str(bad_corr) + ' out of ' + str(good_corr + bad_corr) + ' are greater than corr threshhold of ' + str(corr_threshold))


	#%% Make plots
	plt.close(5)
	plt.figure(5,figsize=(10,10))
	plt.xlim(start_plot_win, end_plot_win)
	plt.ylim(min_dist, max_dist)

	for tr in stgood:
		dist_offset = tr.stats.distance # trying for approx degrees
		time = np.arange(len(tr.data)) * tr.stats.delta + (tr.stats.starttime - t)
		plt.plot(time, (tr.data - np.median(tr.data))*plot_scale_fac /(tr.data.max()
			- tr.data.min()) + dist_offset, color = 'black')

	#%% Plot before shift
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
			plt.plot(time_vec4,dist_vec, color = 'purple')

		plt.plot(time_vec1,dist_vec, color = 'blue')
		plt.show()

	plt.xlabel('Time (s)')
	plt.ylabel('Epicentral distance from event (°)')
	plt.title('Post-alignment ' + dphase + ' for ' + fname[2:12])
	plt.show()

	# plot traces
	plt.close(6)
	plt.figure(6,figsize=(10,10))
	plt.xlim(start_plot_win, end_plot_win)
	plt.ylim(min_dist, max_dist)

	for tr in st2: # regenerate distances into st2 as they were loaded into st for plots
		for ii in station_index: # find station in inventory
			tested_name = st_names[ii]
			actual_trace = tr.stats.station
# 			print('actual trace  ' + actual_trace + '  tested_name ' + tested_name)
			if ARRAY == 0: # convoluted patch for long Hinet names
				this_name_truc = tested_name[0:5]
				name_truc_cap  = this_name_truc.upper()
				this_name = name_truc_cap
				actual_trace = tr.stats.station.upper()
			if (actual_trace == this_name): # found it
				stalon = float(st_lons[ii]) # look up lat & lon to find distance
				stalat = float(st_lats[ii])
				distance = gps2dist_azimuth(stalat,stalon,ev_lat,ev_lon)
				tr.stats.distance=distance[0]/(1000.*111) # distance for phase time and plotting
# 				sys.exit()

	for tr in st2: # generate plot
		dist_offset = tr.stats.distance # trying for approx degrees
		time = np.arange(len(tr.data)) * tr.stats.delta + (tr.stats.starttime - t)
		plt.plot(time, (tr.data - np.median(tr.data))*plot_scale_fac /(tr.data.max()
			- tr.data.min()) + dist_offset, color = 'black')

	#%% Plot after shift
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
			plt.plot(time_vec4,dist_vec, color = 'purple')

		plt.plot(time_vec1,dist_vec, color = 'blue')
		plt.show()

	plt.xlabel('Time (s)')
	plt.ylabel('Epicentral distance from event (°)')
	plt.title('Pre-alignment ' + dphase + ' for ' + fname[2:12])
	plt.show()

	#  Save stats
	fname_stats = '/Users/vidale/Documents/PyCode/Mseed/H.fine_statics.txt'

	#  Save station static correction files
	#fname_stats = 'Statics' + etime[:10] + dphase + ref_trace + '.txt'
	stats_file = open(fname_stats, 'w')
	len_file1 = len(sta_names)
	for j in range(0,len_file1):
		dist_str = '{:.2f}'.format(  sta_dists[j]) # 3 digits after decimal place
		lat_str  = '{:.4f}'.format(   sta_lats[j]) # 2 digits after decimal place
		lon_str  = '{:.4f}'.format(   sta_lons[j])
		stat_str = '{:.3f}'.format(sta_statics[j])
		corr_str = '{:.3f}'.format(  sta_corrs[j])
		write_line = sta_names[j] +' ' + dist_str +' ' + lat_str +' ' + lon_str +' ' + stat_str + ' ' + corr_str + '\n'
		stats_file.write(write_line)
	file.close()
	print('Correlation files have: ' + str(len_file1) + ' traces')

# 	os.system('say "Done"')