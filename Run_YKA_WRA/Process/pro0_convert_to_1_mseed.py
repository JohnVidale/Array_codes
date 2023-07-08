#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 19 14:39:33 2019
Takes individual mseed files and combines them to a single file
@author: vidale
"""
def collect_mseed_files(dir_name):

	from obspy import Stream
	from obspy import read
	import os
	print('Starting' + dir_name)

	#event_yr = '73'
	#event_mo = '12'
	#event_da = '29'
	#event_hrmn = '0019'
	#date_code = 'c' + event_yr + event_mo + event_da + '_' + event_hrmn
	prefix = 'L'
#	date_code = '730425_2134'
	path_name = '/Users/vidale/Documents/Research/IC/Japan/Tanaka_events/' + dir_name
	os.chdir(path_name)

	with open('stations', 'r') as file:
		lines = file.readlines()

	print(str(len(lines)) + ' lines read from file.')

	st   = Stream()
	stfull = Stream()

	station_index = range(len(lines))
	for ii in station_index:
		st = read(lines[ii].rstrip())
		stfull += st

	#nt = len(st71[0].data)
	#dt = st71[0].stats.delta
	#print(str(nt) + ' time pts, time sampling of ' + str(dt) + ' and thus duration of ' + str((nt-1)*dt))

	print('This event: ' + str(len(stfull)) + ' traces.')

	#if do_decimate != 0:
	#	st1.decimate(do_decimate, no_filter=True)

	os.chdir('/Users/vidale/Documents/PyCode/Mseed')

	#fname = 'HD19' + event_yr + '-' + event_mo + '-' + event_da + '.mseed'
	fname = prefix + dir_name + '.mseed'
	stfull.write(fname,format = 'MSEED')