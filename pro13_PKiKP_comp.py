#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat March 21 2020
Compares predicted and observed PKiKP slownesses
@author: vidale
"""
def map_slo_plot(min_dist = 0, max_dist = 180):
	import numpy as np
	import matplotlib.pyplot as plt
	from obspy.taup import TauPyModel

	print('Starting')
	model = TauPyModel(model='iasp91')
	dphase = 'PKiKP'

	sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/event_table_good.txt'
	with open(sta_file, 'r') as file:
		lines = file.readlines()

	print(str(len(lines)) + ' lines read from ' + sta_file)
	# Load station coords into arrays
	event_count = 73
	station_index = range(event_count)
	event_names        = []
	event_year         = np.zeros(event_count)
	event_mo           = np.zeros(event_count)
	event_day          = np.zeros(event_count)
	event_hr           = np.zeros(event_count)
	event_min          = np.zeros(event_count)
	event_sec          = np.zeros(event_count)
	event_lat          = np.zeros(event_count)
	event_lon          = np.zeros(event_count)
	event_dep          = np.zeros(event_count)
	event_mb           = np.zeros(event_count)
	event_ms           = np.zeros(event_count)
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

	for ii in station_index:   # read file
		line = lines[ii]
		split_line = line.split()

		event_names.append(split_line[0])
		event_year[ii]      = float(split_line[1])
		event_mo[ii]        = float(split_line[2])
		event_day[ii]       = float(split_line[3])
		event_hr[ii]        = float(split_line[4])
		event_min[ii]       = float(split_line[5])
		event_sec[ii]       = float(split_line[6])
		event_lat[ii]       = float(split_line[7])
		event_lon[ii]       = float(split_line[8])
		event_dep[ii]       = float(split_line[9])
		event_mb[ii]        = float(split_line[10])
		event_ms[ii]        = float(split_line[11])
		event_tstart[ii]    = float(split_line[12])
		event_tend[ii]      = float(split_line[13])
		event_gcdist[ii]    = float(split_line[14])
		event_dist[ii]      = float(split_line[15])
		event_baz[ii]       = float(split_line[16])
		event_SNR[ii]       = float(split_line[17])
		event_Sflag[ii]     = float(split_line[18])
		event_PKiKPflag[ii] = float(split_line[19])
		event_ICSflag[ii]   = float(split_line[20])
		event_PKiKP_radslo[ii]  = float(split_line[21])
		event_PKiKP_traslo[ii]  = float(split_line[22])

	event_pred_bazi     = np.zeros(event_count)
	event_pred_slo     = np.zeros(event_count)
	event_obs_bazi      = np.zeros(event_count)
	event_obs_slo      = np.zeros(event_count)
	print(str(event_dep[0]))
	print(str(event_PKiKP_radslo[1]))
#	for ii in station_index:   # read file
	for ii in range(event_count):   # read file
		#  find predicted slowness
		arrivals1 = model.get_travel_times(source_depth_in_km=event_dep[ii],distance_in_degree=event_gcdist[ii]-0.5,phase_list=[dphase])
		arrivals2 = model.get_travel_times(source_depth_in_km=event_dep[ii],distance_in_degree=event_gcdist[ii]+0.5,phase_list=[dphase])
		dtime = arrivals2[0].time - arrivals1[0].time
		event_pred_slo[ii]  = dtime/111.  # s/km
		#  find observed slowness
		rad2 = event_PKiKP_radslo[ii]*event_PKiKP_radslo[ii]
		tra2 = event_PKiKP_traslo[ii]*event_PKiKP_traslo[ii]
		event_obs_slo[ii]   = np.sqrt(rad2 + tra2)
		#  find observed back-azimuth
		bazi_rad = np.arctan(event_PKiKP_traslo[ii]/event_PKiKP_radslo[ii])
		event_obs_bazi[ii]  = event_baz[ii] + (bazi_rad * 180 / np.pi)
		#  predicted back-azimuth
		event_pred_bazi[ii] = event_baz[ii]

#		print('t1 t2 dtime: ' + str(arrivals1[0].time) + '  ' + str(arrivals2[0].time) + '  ' + str(dtime))
#		print('rslo tslo in_bazi : ' + str(event_PKiKP_radslo[ii]) + '  ' + str(event_PKiKP_traslo[ii]) + '  ' + str(event_baz[ii]))
#		print('pred: slo bazi  obs: slo bazi : ' + str(event_pred_slo[ii]) + '  ' + str(event_pred_bazi[ii]) + '  ' +
#			 str(event_obs_slo[ii]) + '  ' + str(event_obs_bazi[ii]))

	fig_index = 3
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='polar')
	c = ax.scatter(np.pi * event_baz / 180, event_dist, c=event_ICSflag, s=100, cmap='hsv', alpha=0.75)
	plt.title('LASA events available - map')
	ax.set_theta_zero_location("N")  # theta=0 at the top
	ax.set_theta_direction(-1)  # theta increasing clockwise
	plt.show()

	fig_index = 4
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='polar')

	event_pred_bazi = event_pred_bazi * np.pi / 180.
	event_obs_bazi  = event_obs_bazi * np.pi / 180.

	c = ax.scatter(event_pred_bazi, event_pred_slo, color='blue', s=100, alpha=0.75)
	c = ax.scatter(event_obs_bazi,  event_obs_slo,  color='red', s=100, alpha=0.75)
	for ii in range(event_count):   # read file
		c = ax.plot([event_pred_bazi[ii], event_obs_bazi[ii]], [event_pred_slo[ii], event_obs_slo[ii]], color='black')

	ax.set_rmax(.025)
	ax.grid(True)
	plt.title('Predicted vs observed slowness of PKiKP')
	ax.set_theta_zero_location("N")  # theta=0 at the top
	ax.set_theta_direction(-1)  # theta increasing clockwise
	plt.show()

