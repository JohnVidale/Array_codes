#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 29 10:44:55 2020
Makes a map of the LASA events
@author: vidale
"""
def map_plot(min_dist = 0, max_dist = 180):
	import numpy as np
	import matplotlib.pyplot as plt

	sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/Peng_event_table.txt'
	with open(sta_file, 'r') as file:
		lines = file.readlines()

	print(str(len(lines)) + ' lines read from ' + sta_file)
	# Load station coords into arrays
	event_count = 78
	station_index = range(event_count)
	event_names = []
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
	event_tstart      = np.zeros(event_count)
	event_tend        = np.zeros(event_count)
	event_gcdist      = np.zeros(event_count)
	event_dist        = np.zeros(event_count)
	event_baz         = np.zeros(event_count)
	event_SNR         = np.zeros(event_count)
	event_Sflag       = np.zeros(event_count)
	event_PKiKPflag   = np.zeros(event_count)
	event_ICSflag     = np.zeros(event_count)


	for ii in station_index:
		line = lines[ii]
		split_line = line.split()

		event_names.append(split_line[0])
		event_year[ii] = float(split_line[1])
		event_mo[ii] = float(split_line[2])
		event_day[ii] = float(split_line[3])
		event_hr[ii] = float(split_line[4])
		event_min[ii] = float(split_line[5])
		event_sec[ii] = float(split_line[6])
		event_lat[ii] = float(split_line[7])
		event_lon[ii] = float(split_line[8])
		event_dep[ii] = float(split_line[9])
		event_mb[ii] = float(split_line[10])
		event_ms[ii] = float(split_line[11])
		event_tstart[ii] = float(split_line[12])
		event_tend[ii] = float(split_line[13])
		event_gcdist[ii] = float(split_line[14])
		event_dist[ii] = float(split_line[15])
		event_baz[ii] = float(split_line[16])
		event_SNR[ii] = float(split_line[17])
		event_Sflag[ii] = float(split_line[18])
		event_PKiKPflag[ii] = float(split_line[19])
		event_ICSflag[ii] = float(split_line[20])

	print(event_names[77] + '  ' + str(event_lat[77]) + '  ' + str(event_lon[77]) + '  ' + str(event_tend[77]))


#%% plot data vs prediction
	fig_index = 1
#	plt.figure(1, figsize=(10,5))
	fig, ax = plt.subplots(1)

	#ax = fig.gca()
	ax.set_xticks(np.arange(0, 360, 30))
	ax.set_yticks(np.arange(0, 100, 30))

	# old
#	obs = [0.05,  0.10,  0.03, -0.04, -0.20,  0.10, -0.05, -0.18, -0.10,  0.20, -0.05,  -0.10,  0.05,  0.00]
#	pred = [0.3, -0.6, -0.1,  0.2,  0.5, -0.5,  0.1,  0.5,  0.5, -0.5,  0.4,  0.7, -0.5, 0.3]

#	FILL NUMBERS

#	plt.scatter(event_baz, event_dist, c='k', s=100, alpha=1, marker='.')
	plt.scatter(event_baz, event_dist, c=event_ICSflag, s=100, alpha=1, marker='.', cmap=plt.cm.autumn)
#	plt.scatter(pred,obs, c='k', s=100, alpha=1, marker='.')
#	plt.plot(np.unique(pred), np.poly1d(np.polyfit(pred, obs, 1))(np.unique(pred)))
	plt.grid()
	plt.rc('grid', linestyle="-", color='black')
	plt.xlabel('Back-azimuth')
	plt.ylabel('Distance')
	plt.title('LASA events available - baz-range')
#	plot.legend([])
	plt.colorbar()
	plt.show()

	fig_index = 2
#	plt.figure(fig_index,figsize=(10,5))
	fig, ax = plt.subplots(1)

	#ax = fig.gca()
	ax.set_xticks(np.arange(-180, 180, 30))
	ax.set_yticks(np.arange(-90, 90, 30))

	# old
#	obs = [0.05,  0.10,  0.03, -0.04, -0.20,  0.10, -0.05, -0.18, -0.10,  0.20, -0.05,  -0.10,  0.05,  0.00]
#	pred = [0.3, -0.6, -0.1,  0.2,  0.5, -0.5,  0.1,  0.5,  0.5, -0.5,  0.4,  0.7, -0.5, 0.3]

#	FILL NUMBERS

#	plt.scatter(event_baz, event_dist, c='k', s=100, alpha=1, marker='.')
	plt.scatter(event_lon, event_lat, c=event_dep, s=100, alpha=1, marker='.', cmap=plt.cm.autumn)
#	plt.scatter(pred,obs, c='k', s=100, alpha=1, marker='.')
#	plt.plot(np.unique(pred), np.poly1d(np.polyfit(pred, obs, 1))(np.unique(pred)))
	plt.grid()
	plt.rc('grid', linestyle="-", color='black')
	plt.xlabel('Longitude (°)')
	plt.ylabel('Latitude')
	plt.title('LASA events available - map')
#	plot.legend([])
	plt.colorbar()
	plt.show()