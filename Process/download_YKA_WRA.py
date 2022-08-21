#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 26 09:07:16 2022

@author: wwang071
"""

import os
from obspy import read
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
client = Client("IRIS")

import numpy as np
import matplotlib.pyplot as plt

flag1_WRA = True  # get WRA data?
flag2_YKA = True  # get YKA data?

if flag1_WRA:
    # read in catalog
    catalog='/Users/vidale/Documents/Research/IC/YKA_WRA/info/WRA_Aeqlist.txt'
    file=open(catalog,'r')
    catalines=file.readlines()
    file.close()


    network="AU"
    station="WB*,WC*,WR*";
    data_dir="/Users/vidale/Documents/Research/IC/YKA_WRA/Mseed"

    event_beg=1
    event_end=2
    for ii in range(event_beg,event_end):
        line=catalines[ii]
        arr=line.split()
        eq_time=UTCDateTime(arr[1])
        date_label=str(eq_time.year)+'{:0>2}'.format(eq_time.month)+'{:0>2}'.format(eq_time.day)+'_'+'{:0>2}'.format(eq_time.hour)+'{:0>2}'.format(eq_time.minute)

        try:
            st3=client.get_waveforms(network, station , "*", "BHZ", eq_time, eq_time+2600,
                                  attach_response=True)
            st3.remove_response(output="VEL", pre_filt=None,)
        except:
            print(eq_time)
        else:
            st=st3.merge(method=0);

            outfile=data_dir+'/'+date_label+'.mseed'
            st.write(outfile,format='MSEED')
            print(date_label,len(st),len(st3))


if flag2_YKA:
    # read in catalog
    catalog='/Users/vidale/Documents/Research/IC/YKA_WRA/info/YKA_Aeqlist.txt'
    file=open(catalog,'r')
    catalines=file.readlines()
    file.close()

    # Yellowknife
    network="CN"
    station="YK*";
    data_dir="/Users/vidale/Documents/Research/IC/YKA_WRA/Mseed"

    event_beg=30
    event_end=31
    for ii in range(event_beg,event_end):
        line=catalines[ii]
        arr=line.split()
        eq_time=UTCDateTime(arr[1])
        date_label=str(eq_time.year)+'{:0>2}'.format(eq_time.month)+'{:0>2}'.format(eq_time.day)+'_'+'{:0>2}'.format(eq_time.hour)+'{:0>2}'.format(eq_time.minute)

        try:
            st3=client.get_waveforms(network, station , "*", "SHZ", eq_time, eq_time+2600,
                                  attach_response=True)
            st3.remove_response(output="VEL", pre_filt=None,)
        except:
            print(eq_time)
        else:
            st=st3.merge(method=0);

            outfile=data_dir+'/'+date_label+'.mseed'
            st.write(outfile,format='MSEED')
            print(date_label,len(st),len(st3))
