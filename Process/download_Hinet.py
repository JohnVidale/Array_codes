#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 26 11:12:51 2022

@author: wwang071
"""
# import sys
# sys.path.append('/usr/local/bin')

from obspy import read
from obspy import Stream
from obspy import UTCDateTime
from HinetPy import Client
from HinetPy import win32
client = Client('weiwang053', 'Majesty0730')
# client = Client('jvidale', 'UmiX7R0f1ZgF')
client.check_cmd_exists()

import os
import numpy as np
import matplotlib.pyplot as plt

#%%
def main():
    stafile='/Users/vidale/Documents/GitHub/Array_codes/Files/sta_hinet.txt';
    file=open(stafile,'r')
    stalines=file.readlines()
    file.close()
    stnm=[]; stla=[]; stlo=[];
    for line in stalines:
        arr=line.split()
        stnm.append(arr[0]);
        stla.append(float(arr[1]))
        stlo.append(float(arr[2]))
    stlo_mean=np.mean(stlo)
    stla_mean=np.mean(stla)

    catalog='/Users/vidale/Documents/GitHub/Wei/wei_event_file.txt'
    file=open(catalog,'r')
    catalines=file.readlines()
    file.close()
    data_dir2='/Users/vidale/Documents/Research/IC/YKA_WRA/Mseed'


    for ii in range(1):
        arr=catalines[ii].split()
        tmpline=arr[0]+'T'+arr[1]
        datetime_UTC=UTCDateTime(tmpline)
        datetime_local=datetime_UTC+9*3600
        date_label_UTC = str(datetime_UTC.year)+'{:0>2}'.format(datetime_UTC.month)+'{:0>2}'.format(datetime_UTC.day)+'{:0>2}'.format(datetime_UTC.hour)+'{:0>2}'.format(datetime_UTC.minute)
        date_label_local = str(datetime_local.year)+'{:0>2}'.format(datetime_local.month)+'{:0>2}'.format(datetime_local.day)+'{:0>2}'.format(datetime_local.hour)+'{:0>2}'.format(datetime_local.minute)

        data_dir=date_label_UTC

        #if not(os.path.exists(data_dir)):
        data, ctable = client.get_continuous_waveform('0101', date_label_local, 1)
        win32.extract_sac(data, ctable,outdir=data_dir,filter_by_component="U")
    #%%
        plt.close(1); plt.figure(1)
        saclist=os.listdir(data_dir)
        st=Stream()
        indx=0
        for sactmp in saclist:
            sacfile=data_dir+'/'+sactmp
            st_tmp=read(sacfile)
            #print(st_tmp[0].stats.station[0:5])
            st_tmp[0].stats.station=st_tmp[0].stats.station[0:5]
            st_tmp[0].stats.starttime-=9*3600
            #st_tmp[0].stats.endtime-=9*3600
            st+=st_tmp

            ttt=np.arange(st_tmp[0].stats.npts)*st_tmp[0].stats.delta
            arr=st_tmp[0].data/max(abs(st_tmp[0].data))
            if indx<20:
                plt.plot(ttt,arr+indx,'k',linewidth=0.5)
            indx+=1


        seedfile = data_dir2 + '/' + date_label_UTC +'.mseed'
        st.write(seedfile,format='MSEED')

if __name__ == '__main__':
    main()
