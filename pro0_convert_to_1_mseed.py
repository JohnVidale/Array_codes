#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 19 14:39:33 2019
Takes individual mseed files and combines them to a single file
@author: vidale
"""
from obspy import Stream
from obspy import read
import os

event_no = '83'
event_yr = '72'
event_mo = '01'
event_da = '25'
event_hrmn = '0206'
#date_code = '730701_1333'
date_code = 'c' + event_yr + event_mo + event_da + '_' + event_hrmn
path_name = '/Users/vidale/Documents/PyCode/LASA/Raw_dupes/' + date_code
os.chdir(path_name)

#records_file = '/Users/vidale/Documents/PyCode/LASA/Raw/' + date_code + '/stations'
#records_file = '/Users/vidale/Documents/PyCode/LASA/Raw/stations' + event_no
with open('stations', 'r') as file:
	lines = file.readlines()

print(str(len(lines)) + ' lines read from file, ' + str(len(lines)))

st   = Stream()
stfull = Stream()

station_index = range(len(lines))
for ii in station_index:
#	st = read(date_code + '/' + lines[ii].rstrip())
	st = read(lines[ii].rstrip())
	stfull += st

#nt = len(st71[0].data)
#dt = st71[0].stats.delta
#print(str(nt) + ' time pts, time sampling of ' + str(dt) + ' and thus duration of ' + str((nt-1)*dt))

print('This event: ' + str(len(stfull)) + ' traces.')

#if do_decimate != 0:
#	st1.decimate(do_decimate, no_filter=True)

os.chdir('/Users/vidale/Documents/PyCode/LASA/Mseed')

fname = 'HD19' + event_yr + '-' + event_mo + '-' + event_da + '.mseed'
stfull.write(fname,format = 'MSEED')