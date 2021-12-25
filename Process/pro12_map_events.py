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

#	sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/Peng_event_table.txt'
	sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/events_good.txt'
	with open(sta_file, 'r') as file:
		lines = file.readlines()
	event_count = len(lines)

	print(str(event_count) + ' lines read from ' + sta_file)
	# Load station coords into arrays
	station_index = range(event_count)
	event_names = []
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


	for ii in station_index:
		line = lines[ii]
		split_line = line.split()

		event_names.append(split_line[0])
		event_index[ii] = float(split_line[1])
		event_year[ii]  = float(split_line[2])
		event_mo[ii]    = float(split_line[3])
		event_day[ii]   = float(split_line[4])
		event_hr[ii]    = float(split_line[5])
		event_min[ii]   = float(split_line[6])
		event_sec[ii]   = float(split_line[7])
		event_lat[ii]   = float(split_line[8])
		event_lon[ii]   = float(split_line[9])
		event_dep[ii]   = float(split_line[10])
		event_mb[ii]    = float(split_line[11])
		event_ms[ii]    = float(split_line[12])
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

#	print(event_names[77] + '  ' + str(event_lat[77]) + '  ' + str(event_lon[77]) + '  ' + str(event_tend[77]))


#%% plot data vs prediction
	fig_index = 1
#	plt.figure(1, figsize=(10,5))
	fig, ax = plt.subplots(1)

	#ax = fig.gca()
	ax.set_xticks(np.arange(0, 360, 30))
	ax.set_yticks(np.arange(0, 100, 30))

#	FILL NUMBERS

#	plt.scatter(event_baz, event_dist, c='k', s=100, alpha=1, marker='.')
	plt.scatter(event_lon, event_lat, c=event_PKiKPflag, s=100, alpha=1, marker='.', cmap=plt.cm.autumn)
#	plt.scatter(pred,obs, c='k', s=100, alpha=1, marker='.')
#	plt.plot(np.unique(pred), np.poly1d(np.polyfit(pred, obs, 1))(np.unique(pred)))
	plt.grid()
	plt.rc('grid', linestyle="-", color='black')
	plt.xlabel('Longitude (°)')
	plt.ylabel('Latitude')
	plt.title('Peng PKiKP quality - map')
#	plot.legend([])
	plt.colorbar()
	plt.show()

	fig_index = 2
#	plt.figure(fig_index,figsize=(10,5))
	fig, ax = plt.subplots(1)

	#ax = fig.gca()
	ax.set_xticks(np.arange(-180, 180, 30))
	ax.set_yticks(np.arange(-90, 90, 30))

#	FILL NUMBERS

#	plt.scatter(event_baz, event_dist, c='k', s=100, alpha=1, marker='.')
	plt.scatter(event_lon, event_lat, c=event_PKiKP_qual, s=100, alpha=1, marker='.', cmap=plt.cm.autumn)
#	plt.scatter(pred,obs, c='k', s=100, alpha=1, marker='.')
#	plt.plot(np.unique(pred), np.poly1d(np.polyfit(pred, obs, 1))(np.unique(pred)))
	plt.grid()
	plt.rc('grid', linestyle="-", color='black')
	plt.xlabel('Longitude (°)')
	plt.ylabel('Latitude')
	plt.title('Our PKiKP quality - map')
#	plot.legend([])
	plt.colorbar()
	plt.show()

	fig_index = 3
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='polar')
	cc = ax.scatter(np.pi * event_baz / 180, event_dist, c=event_PKiKP_qual, s=100, cmap='brg', alpha=0.75)
#	cc = ax.scatter(np.pi * event_baz / 180, event_dist, c=event_baz, s=100, cmap='hsv', alpha=0.75)

	plt.title('PKiKP quality - polar plot')
#	plot.legend([])
#	plt.colorbar()
	ax.set_theta_zero_location("N")  # theta=0 at the top
	ax.set_theta_direction(-1)  # theta increasing clockwise
	plt.show()

	fig_index = 4
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='polar')
	cc = ax.scatter(np.pi * event_baz / 180, event_dist, c=event_ICS_qual, s=100, cmap='brg', alpha=0.75)
#	cc = ax.scatter(np.pi * event_baz / 180, event_dist, c=event_baz, s=100, cmap='hsv', alpha=0.75)

	plt.title('ICS quality - polar plot')
#	plot.legend([])
#	plt.colorbar()
	ax.set_theta_zero_location("N")  # theta=0 at the top
	ax.set_theta_direction(-1)  # theta increasing clockwise
	plt.show()