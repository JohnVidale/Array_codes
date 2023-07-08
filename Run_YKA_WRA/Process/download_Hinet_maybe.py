#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 26 11:12:51 2022
tested 6/28/21

@author: wwang071
"""

from obspy   import read
from obspy   import Stream
from obspy   import UTCDateTime
from HinetPy import Client
from HinetPy import win32
client = Client('weiwang053', 'Majesty0730')
client = Client('jvidale', 'UmiX7R0f1ZgF')

import os
import numpy as np


    #%% file names
stafile  ='/Users/vidale/Documents/GitHub/Array_codes/Files/sta_hinet.txt';  # list of HiNet stations
dir_data ='/Users/vidale/Documents/Research/IC/YKA_WRA/Mseed'               # output directory
catalog  ='/Users/vidale/Documents/GitHub/Wei/wei_event_file.txt'            # event parameters

    #%% get station names and positions
file = open(stafile,'r')
stalines = file.readlines()
file.close()
stnm = [];
stla = [];
stlo = [];
for line in stalines:
    arr = line.split()
    stnm.append(arr[0]);
    stla.append(float(arr[1]))
    stlo.append(float(arr[2]))

stlo_mean = np.mean(stlo)
stla_mean = np.mean(stla)

    #%% get event catalog
    # example format '2011/04/19      05:50:45.240    -33.942 -72.358 15.30   5.3 155.52'
file = open(catalog,'r')
catalines = file.readlines()
file.close()

num_events = len(catalines)
print('Reading ' + str(num_events) + ' events.')
for ii in range(num_events):
    arr = catalines[ii].split()
    tmpline = arr[0] + 'T' + arr[1]
    datetime_UTC   = UTCDateTime(tmpline)
    datetime_local = datetime_UTC+9*3600
    date_label_UTC = (str(datetime_UTC.year) + '{:0>2}'.format(datetime_UTC.month) + '{:0>2}'.format(datetime_UTC.day) +
                      '{:0>2}'.format(datetime_UTC.hour) + '{:0>2}'.format(datetime_UTC.minute))
    date_label_local = (str(datetime_local.year) + '{:0>2}'.format(datetime_local.month) + '{:0>2}'.format(datetime_local.day)
                        + '{:0>2}'.format(datetime_local.hour) + '{:0>2}'.format(datetime_local.minute))

    data_dir = date_label_UTC

    #if not(os.path.exists(data_dir)):
    #%% get HiNet format files
    data, ctable = client.get_continuous_waveform('0101', date_label_local, 1)     # number is minutes retrieved
    #%% convert to SAC format files
    win32.extract_sac(data, ctable, outdir = data_dir, filter_by_component = 'U')   # U is just vertical, omit for all components

    #%% construct mseed file
    saclist = os.listdir(data_dir)
    st = Stream()
    for sactmp in saclist:
        sacfile = data_dir+'/'+sactmp
        st_tmp = read(sacfile)
        #print(st_tmp[0].stats.station[0:5])
        st_tmp[0].stats.station = st_tmp[0].stats.station[0:5]
        st_tmp[0].stats.starttime -= 9*3600   # adjust JST to UTC
        #st_tmp[0].stats.endtime -= 9*3600
        st += st_tmp

    #%% write mseed file to disk
    seedfile = dir_data + '/' + date_label_UTC +'.mseed'
    st.write(seedfile, format='MSEED')
