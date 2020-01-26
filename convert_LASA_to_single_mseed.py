#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 19 14:39:33 2019

@author: vidale
"""
from obspy import Stream
from obspy import read
import os
os.chdir('/Users/vidale/Documents/PyCode/LASA/Raw_Processed')

#records_file71 = '/Users/vidale/Documents/PyCode/LASA/E71.stations' # get list of files for '71 explosion
records_file = '/Users/vidale/Documents/PyCode/LASA/Raw_Processed/stations12'
with open(records_file, 'r') as file:
	lines = file.readlines()

print(str(len(lines)) + ' lines read from file, ' + str(len(lines)))

st   = Stream()
stfull = Stream()

station_index = range(len(lines))
for ii in station_index:
	st = read('c701217_1605/' + lines[ii].rstrip())
	stfull += st

#nt = len(st71[0].data)
#dt = st71[0].stats.delta
#print(str(nt) + ' time pts, time sampling of ' + str(dt) + ' and thus duration of ' + str((nt-1)*dt))

print('This event: ' + str(len(stfull)) + ' traces.')

#if do_decimate != 0:
#	st1.decimate(do_decimate)

fname = 'event12.mseed'
stfull.write(fname,format = 'MSEED')