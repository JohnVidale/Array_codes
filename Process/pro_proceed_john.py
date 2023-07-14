#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 21:01:52 2020

@author: wwang071
"""

def define_taper_frac(start_buff,end_buff,taper_frac):
    totalt = end_buff - start_buff
    noise_time_skipped = taper_frac * totalt
    if noise_time_skipped >= 0.5 * start_buff:
        taper_frac = 0.5*(-start_buff)/totalt
        if start_buff >= 0:
            taper_frac = 0.05 # pick random minimal window if there is no leader
    return taper_frac

def trim_seismo(st='', eve_file='', sta_file='',
                min_dist=0, max_dist=180, dphase='PKiKP',
                start_buff=-10, end_buff=30,
                lon_ref=0, lat_ref=0, rel_time=0,
                stat_corr=0, corr_threshold=0):

    from obspy import UTCDateTime
    from obspy import Stream
    from obspy.geodetics import gps2dist_azimuth
    from obspy.taup import TauPyModel
    from pro_read_file import read_sta_file, read_sta_corr_file, read_event_file

    model = TauPyModel(model='iasp91')

    # read in event info
    eq_time,ev_lon,ev_lat,ev_dep=read_event_file(eve_file)
    # read in station info
    if stat_corr == 1:  # load static terms, only applies to Hinet and LASA
        sta_num,st_names,st_dist,st_lons,st_lats,st_deps,st_shift,st_corr=read_sta_corr_file(sta_file)
    else: # no static terms, always true for LASA or NORSAR
        sta_num,st_names,st_lons,st_lats,st_deps=read_sta_file(sta_file);


    # compute arrival time at the reference point
    if rel_time==0:
        distance = gps2dist_azimuth(lat_ref,lon_ref,ev_lat,ev_lon)
        dist_ref = distance[0]/(1000*111.19)
        #print(dist_ref)
        arrivals = model.get_travel_times(source_depth_in_km=ev_dep,distance_in_degree=dist_ref,
                                                      phase_list=[dphase])
        atime_ref=arrivals[0].time
        #print("dist_ref=",atime_ref,lat_ref,lon_ref)

    st_tmp=st.copy()
    trace_num_in_range=0
    st_pickalign=Stream()
    for tr in st_tmp:
        if tr.stats.station in st_names:
            if float(eq_time.year) < 1970: # fix the damn 1969 -> 2069 bug in Gibbon's LASA data
                temp_t = str(tr.stats.starttime)
                temp_tt = '19' + temp_t[2:]
                tr.stats.starttime = UTCDateTime(temp_tt)

            indx=st_names.index(tr.stats.station)
            #print(st_lats[indx],st_lons[indx],ev_lat,ev_lon)
            distance = gps2dist_azimuth(st_lats[indx],st_lons[indx],ev_lat,ev_lon)
            tr.stats.distance = distance[0]/1000
            dist = distance[0]/(1000*111.19)
            #print("dist=",dist)
            if min_dist < dist and dist < max_dist: # select distance range from earthquake
                trace_num_in_range += 1
                if rel_time == 1:
                    arrivals = model.get_travel_times(source_depth_in_km=ev_dep,distance_in_degree=dist,
                                                  phase_list=[dphase])
                    atime = arrivals[0].time
                    #print(tr.stats.station,atime)
                    s_t = eq_time+atime+start_buff # start cutting time
                    e_t = eq_time+atime+end_buff # end cutting time
                    #print(s_t,e_t)
                else:
                    s_t = eq_time+atime_ref+start_buff
                    e_t = eq_time+atime_ref+end_buff
                    #print(s_t,e_t,tr.stats.endtime)
                if trace_num_in_range==1:
                    print(s_t,e_t)
                if stat_corr==1:
                    if st_corr[indx] > corr_threshold:
                        #tr.stats.starttime -= float(st_shift[indx])
                        #print(tr.stats.starttime)
                        tr.stats.update({"starttime":tr.stats.starttime-st_shift[indx]})
                        #print(tr.stats.starttime)
                        tr.trim(starttime = s_t,endtime = e_t)
                        tr.starttime=eq_time+start_buff
                        tr.stats.update({"starttime":eq_time+start_buff})
                        st_pickalign += tr
                        #print(st_pickalign[trace_num_in_range-1].stats.starttime)
                elif stat_corr==0:
                    tr.trim(starttime = s_t,endtime = e_t)
                    tr.starttime=eq_time+start_buff
                    tr.stats.update({"starttime":eq_time+start_buff})
                    st_pickalign += tr
                    #print(st_pickalign[trace_num_in_range-1].stats.starttime)
    #print(st_pickalign[0].stats.starttime)
    return trace_num_in_range,st_pickalign


def define_SNR(trace='',signal_twin=10,start_buff=-10,end_buff=30,taper_frac=0.05):
    import numpy as np
    # estimate median noise
    t_noise_start  = int(len(trace.data) * taper_frac)
    t_noise_end    = int(len(trace.data) * (-start_buff)/(end_buff-start_buff))
    noise          = np.median(abs(trace.data[t_noise_start:t_noise_end]))
    # estimate median signal
    t_signal_start = int(len(trace.data) * (-start_buff)/(end_buff-start_buff))
    t_signal_end   = t_signal_start + int(len(trace.data) * signal_twin/(end_buff-start_buff))
    signal         = np.median(abs(trace.data[t_signal_start:t_signal_end]))
    #            test SNR
    SNR = signal/noise;
    return SNR

def derive_cc_trace(st='',eve_file='',sta_file='',tr_ref='',beg_cc_twin=-1,end_cc_twin=20,
                    max_time_shift = 2,dphase='PKIKP',corr_threshold=0):

    from obspy import Stream
    from pro_read_file import read_sta_file, read_event_file
    from obspy.geodetics import gps2dist_azimuth
    from obspy.signal.cross_correlation import xcorr_pick_correction
    from obspy.taup import TauPyModel
    model = TauPyModel(model='iasp91')
    # read in event file
    eq_time,ev_lon,ev_lat,ev_dep=read_event_file(eve_file)
    # read in station file
    sta_num,st_names,st_lons,st_lats,st_deps=read_sta_file(sta_file)

    # initialize lists of statics
    stgood = Stream()
    sta_names   = []
    sta_dists   = []
    sta_lats    = []
    sta_lons    = []
    sta_statics = []
    sta_corrs   = []

    good_corr=0; bad_corr=0;
    for tr in st:
        indx=st_names.index(tr.stats.station)
        distance=gps2dist_azimuth(st_lats[indx],st_lons[indx],ev_lat,ev_lon)
        dist=distance[0]/1000/111.19   # m to deg
        tr.stats.distance=dist
        arrivals=model.get_travel_times(source_depth_in_km=ev_dep,distance_in_degree=dist,phase_list=[dphase])
        dt,coeff = xcorr_pick_correction(eq_time,tr_ref,eq_time,tr,beg_cc_twin,end_cc_twin,max_time_shift,plot=False)
        if coeff <= 1 and coeff > corr_threshold:
            good_corr+=1
            tr.stats.starttime-=dt
            sta_names.append(tr.stats.station)
            sta_dists.append(tr.stats.distance)
            sta_lats.append(st_lats[indx])
            sta_lons.append(st_lons[indx])
            sta_statics.append(dt)
            sta_corrs.append(coeff)
            stgood+=tr
        else:
            bad_corr+=1

    return stgood,sta_names,sta_dists,sta_lons,sta_lats,sta_statics,sta_corrs

def stack_1d(stream='',eq_file='',sta_file='',slow_low=-0.05,slow_high=0.05,slow_delta=0.01,
             ref_lon=0,ref_lat=0,start_buff=-10,end_buff=100,env=0, norm=1):
    from obspy import Stream,Trace
    from pro_read_file import read_sta_file, read_event_file
    from obspy.geodetics import gps2dist_azimuth
    from scipy.signal import hilbert
    import numpy as np

    # read in event file
    eq_time,ev_lon,ev_lat,ev_dep=read_event_file(eq_file)
    # read in station file
    sta_num,st_names,st_lons,st_lats,st_deps=read_sta_file(sta_file)

    dt=stream[0].stats.delta
    nt = len(stream[0].data)
    # Build Stack arrays
    stack = Stream()
    tr = Trace()
    tr.stats.delta = stream[0].stats.delta
    tr.stats.network = 'stack'
    tr.stats.channel = 'BHZ'
    slow_n = int(1 + (slow_high - slow_low)/slow_delta)  # number of slownesses
    stack_nt = int(1 + ((end_buff - start_buff)/dt))  # number of time points
    # In English, stack_slows = range(slow_n) * slow_delta - slowR_lo
    a1 = range(slow_n)
    stack_slows = [(x * slow_delta + slow_low) for x in a1]
    print(str(slow_n) + ' slownesses.')
    tr.stats.starttime = eq_time + start_buff
    tr.data = np.zeros(stack_nt)
    done = 0
    for stack_one in stack_slows:
        tr1 = tr.copy()
        tr1.stats.station = str(int(done))
        stack.extend([tr1])
        done += 1
        #print(done)

    #  Only need to compute ref location to event distance once
    ref_distance = gps2dist_azimuth(ev_lat,ev_lon,ref_lat,ref_lon)

    # Select traces by distance, window and adjust start time to align picked times
    done = 0
    for tr in stream: # traces one by one, find lat-lon by searching entire inventory.  Inefficient but cheap
        if tr.stats.station in st_names:
            indx=st_names.index(tr.stats.station)
            if norm == 1:
                tr.normalize()
            distance = gps2dist_azimuth(st_lats[indx],st_lons[indx],ev_lat,ev_lon) # Get traveltimes again, hard to store
            tr.stats.distance=distance[0] # distance in m
            del_dist = (ref_distance[0] - distance[0])/(1000) # in km
            arr=tr.data;
            if env==1:
                arr=abs(hilbert(arr))
            for slow_i in range(slow_n):  # for this station, loop over slownesses
                time_lag = -del_dist * stack_slows[slow_i]  # time shift due to slowness, flipped to match 2D
                time_correction = ((eq_time-tr.stats.starttime) + (time_lag + start_buff))/dt

                nshift=int(time_correction)
                if time_correction<0:
                    nshift=nshift-1
                if nshift<=0:
                    nbeg1=-nshift; nend1=stack_nt
                    nbeg2=0; nend2=stack_nt+nshift;

                elif nshift>0:
                    nbeg1=0; nend1=stack_nt-nshift
                    nbeg2=nshift; nend2=stack_nt
                stack[slow_i].data[nbeg1:nend1] += arr[nbeg2:nend2]
                #print('nshift=',nshift,nbeg1,nend1,nbeg2,nend2)

#                for it in range(stack_nt):  # check points one at a time
#                    it_in = int(it + time_correction)
#                    #print(it_in,it,int(time_correction),round(time_correction))
#                    if it_in >= 0 and it_in < nt - 1: # does data lie within seismogram?
#                        stack[slow_i].data[it] += arr[it_in]
            done += 1
            if done%50 == 0:
                print('Done stacking ' + str(done) + ' out of ' + str(len(stream)) + ' stations.')
    return stack_slows,stack

def stack_2d(stream='',eq_file='',sta_file='', slowR_low=-0.05, slowR_high=-0.05,
             slowT_low=-0.05, slowT_high=-0.05, slow_delta=0.01, ref_lon=0,ref_lat=0,
             start_buff=-10, end_buff=100,env=0, norm=1, NS=0):
    from obspy import Stream,Trace
    from pro_read_file import read_sta_file, read_event_file
    from obspy.geodetics import gps2dist_azimuth
    from scipy.signal import hilbert
    import numpy as np
    import math


    # read in event file
    eq_time,ev_lon,ev_lat,ev_dep=read_event_file(eq_file)
    # read in station file
    sta_num,st_names,st_lons,st_lats,st_deps=read_sta_file(sta_file)

    dt=stream[0].stats.delta
    nt = len(stream[0].data)

    # initialize parameters
    slowR_n = int(1 + (slowR_high - slowR_low)/slow_delta)    # number of radial slownesses
    slowT_n = int(1 + (slowT_high - slowT_low)/slow_delta)    # number of transverse slownesses
    stack_nt = int(1 + ((end_buff - start_buff)/dt))          # number of time points
    stack_Rslows = [(x * slow_delta + slowR_low) for x in range(slowR_n)]
    stack_Tslows = [(x * slow_delta + slowT_low) for x in range(slowT_n)]

    stack = Stream()
    tr = Trace()
    tr.stats.delta = dt
    tr.stats.starttime = eq_time + start_buff
    tr.stats.npts = stack_nt
    tr.stats.network = 'stack'
    tr.stats.channel = 'BHZ'
    tr.data = np.zeros(stack_nt)
    done=0
    for stackR_one in stack_Rslows:
        for stackT_one in stack_Tslows:
            tr1 = tr.copy()
            tr1.stats.station = str(int(done))
            stack.extend([tr1])
            done+=1

    #  Only need to compute ref location to event distance once
    ref_dist_az = gps2dist_azimuth(ev_lat,ev_lon,ref_lat,ref_lon)
    ref_back_az = ref_dist_az[2]


    done = 0
    for tr in stream: # traces one by one, find lat-lon by searching entire inventory.  Inefficient but cheap
        ii=st_names.index(tr.stats.station)
        stalat=st_lats[ii]
        stalon=st_lons[ii]
        #print(ii,tr.stats.station,stalat,stalon);

        if norm == 1:
            tr.normalize() # trace divided abs(max of trace)
        stalat = float(st_lats[ii])
        stalon = float(st_lons[ii]) # use lat & lon to find distance and back-az
        rel_dist_az = gps2dist_azimuth(stalat,stalon,ref_lat,ref_lon)
        rel_dist    = rel_dist_az[0]/1000  # km
        rel_back_az = rel_dist_az[1]       # radians
        #print(stalat,stalon,rel_back_az,ref_back_az)
        if NS == 0:
            del_distR = rel_dist * math.cos((rel_back_az - ref_back_az)* math.pi/180)
            del_distT = rel_dist * math.sin((rel_back_az - ref_back_az)* math.pi/180)
            # North and east
        else:
            del_distR = rel_dist * math.cos( rel_back_az * math.pi/180)
            del_distT = rel_dist * math.sin( rel_back_az * math.pi/180)

        arr=tr.data;
        if env==1:
            arr=abs(hilbert(arr))
        for slowR_i in range(slowR_n):  # for this station, loop over radial slownesses
            for slowT_i in range(slowT_n):  # loop over transverse slownesses
                time_lag  = del_distR * stack_Rslows[slowR_i]  # time shift due to radial slowness
                time_lag += del_distT * stack_Tslows[slowT_i]  # time shift due to transverse slowness
                time_correction = ((eq_time-tr.stats.starttime) + (time_lag + start_buff))/dt
                indx = int(slowR_i*slowT_n + slowT_i)

                nshift=round(time_correction)
                #if time_correction<0:
                #    nshift=nshift-1
                if nshift<=0:
                    nbeg1=-nshift; nend1=stack_nt
                    nbeg2=0; nend2=stack_nt+nshift;

                elif nshift>0:
                    nbeg1=0; nend1=stack_nt-nshift
                    nbeg2=nshift; nend2=stack_nt
                #print('nshift=',nshift,nbeg1,nend1,nbeg2,nend2)
                stack[indx].data[nbeg1:nend1] += arr[nbeg2:nend2]
                #print('nshift=',nshift,nbeg1,nend1,nbeg2,nend2)

#                for it in range(stack_nt):  # check points one at a time
#                    it_in = int(it + time_correction)
#                    if it_in >= 0 and it_in < nt - 2: # does data lie within seismogram?
#     `               # should be 1, not 2, but 2 prevents the problem "index XX is out of bounds for axis 0 with size XX"
#                        stack[indx].data[it] += arr[it_in]
        done += 1
        if done%20 == 0:
            print('Done stacking ' + str(done) + ' out of ' + str(len(stream)) + ' stations.')

    return stack_Rslows, stack_Tslows,stack

def stack_distance(st='',eve_file='', sta_file='',min_dist=0, max_dist=180, int_dist=1,
                  start_buff=0, end_buff=100):
    from obspy import Stream,Trace
    from obspy.geodetics import gps2dist_azimuth
    from pro_read_file import read_sta_file, read_event_file
    from scipy.signal import hilbert
    import numpy as np
    import matplotlib.pyplot as plt

    # read in event file
    eq_time,ev_lon,ev_lat,ev_dep=read_event_file(eve_file);
    # read in station file
    sta_num,st_names,st_lons,st_lats,st_deps=read_sta_file(sta_file);

    dt=st[0].stats.delta
    nt = len(st[0].data)

    dist_n = int(1 + (max_dist - min_dist)/int_dist);
    stack_nt = nt;
    stack = Stream()
    tr = Trace()
    tr.stats.delta = dt
    tr.stats.starttime = eq_time + start_buff
    tr.stats.npts = stack_nt
    tr.stats.network = 'stack'
    tr.stats.channel = 'BHZ'
    tr.data = np.zeros(stack_nt)
    done=0
    for ii in range(dist_n):
            tr1 = tr.copy()
            tr1.stats.station = str(int(done))
            stack.extend([tr1])
            done+=1


    deg2km=6371*2*np.pi/360;
    dist_indx_arr=[];
    if len(st)>10:
        for tr in st:
            if tr.stats.station in st_names:
                indx=st_names.index(tr.stats.station);
                distance=gps2dist_azimuth(ev_lat,ev_lon,st_lats[indx],st_lons[indx]);
                dist_deg=distance[0]/1000/deg2km;
                dist_indx=int((dist_deg-min_dist)/int_dist)
                dist_indx_arr.append(dist_indx)
                #print(len(tr.data))
                arr=abs(hilbert(tr.data));
                arr=arr/max(arr);
                if len(arr)==stack_nt:
                    stack[dist_indx].data=stack[dist_indx].data+arr;


#    print(min(dist_indx_arr),max(dist_indx_arr))
#    plt.figure(101)
#    plt.hist(dist_indx_arr)
    return stack


def derive_SNR_file(st='',eve_file='', sta_file='',dphase=['P','PKiKP'],start_buff=0,
                    end_buff=60, freq_min=1, freq_max=3, stat_corr=1, corr_threshold='', SNR_file=''):

    from obspy import UTCDateTime
    from obspy import Stream
    from scipy.signal import hilbert
    from obspy.geodetics import gps2dist_azimuth
    import matplotlib.pyplot as plt
    import numpy as np
    from obspy.taup import TauPyModel
    from pro_read_file import read_sta_file, read_sta_corr_file, read_event_file

    model = TauPyModel(model='iasp91')

    # read in event info
    eq_time,ev_lon,ev_lat,ev_dep=read_event_file(eve_file)
    # read in station info
    if stat_corr == 1:  # load static terms, only applies to Hinet and LASA
        sta_num,st_names,st_dist,st_lons,st_lats,st_deps,st_shift,st_corr=read_sta_corr_file(sta_file)
    else: # no static terms, always true for LASA or NORSAR
        sta_num,st_names,st_lons,st_lats,st_deps=read_sta_file(sta_file);

    file=open(SNR_file,'w');
    trace_num_in_range=0

    st_tmp=st;

    st_tmp.detrend(type='simple')
    #st_tmp.taper(0.05)
    st_tmp.filter('bandpass', freqmin=freq_min, freqmax=freq_max, corners=3, zerophase=True)
    #st_tmp.taper(0.05)

    for tr in st_tmp:
        #print(tr.stats.station,st_names)
        if tr.stats.station in st_names:
            if float(eq_time.year) < 1970: # fix the damn 1969 -> 2069 bug in Gibbon's LASA data
                temp_t = str(tr.stats.starttime)
                temp_tt = '19' + temp_t[2:]
                tr.stats.starttime = UTCDateTime(temp_tt)

            indx=st_names.index(tr.stats.station)
            #print(indx,tr.stats.station)
            distance = gps2dist_azimuth(st_lats[indx],st_lons[indx],ev_lat,ev_lon)
            tr.stats.distance = distance[0]/1000
            dist = distance[0]/(1000*111.19)

            arrivals1=model.get_travel_times(source_depth_in_km=ev_dep,distance_in_degree=dist,
                                            phase_list=[dphase[0]])
            tp1=arrivals1[0].time;
            #tp1=tr.stats.sac.a

            arrivals2=model.get_travel_times(source_depth_in_km=ev_dep,distance_in_degree=dist,
                                            phase_list=[dphase[1]])

            tp2=arrivals2[0].time;

            ttp1=eq_time+tp1;
            ttp1_diff=ttp1-tr.stats.starttime

            ttp2=eq_time+tp2;
            ttp2_diff=ttp2-tr.stats.starttime

            nbeg=round((ttp1_diff+start_buff)/tr.stats.delta)
            nend=round((ttp2_diff+end_buff)/tr.stats.delta)
            nn_num=nend-nbeg;
            ttt=np.arange(nn_num)*tr.stats.delta-50;
            #arr=tr.data[nbeg:nend];
            #print(max(abs(arr)))
            arr=(abs((tr.data)));
            #print(ttp1_diff,ttp2_diff)
            # noise twin
            nbeg1=round((ttp1_diff+start_buff-15)/tr.stats.delta)
            nend1=round((ttp1_diff+start_buff-5)/tr.stats.delta)
            noise_lvl=np.nanmedian(arr[nbeg1:nend1])
            # P twin
            nbeg2=round((ttp1_diff)/tr.stats.delta)
            nend2=round((ttp1_diff+50)/tr.stats.delta)
            p_sig_lvl=np.nanmedian(arr[nbeg2:nend2])
            # PKiKP twin
            nbeg3=round((ttp2_diff)/tr.stats.delta)
            nend3=round((ttp2_diff+30)/tr.stats.delta)
            pkikp_sig_lvl=np.nanmedian(arr[nbeg3:nend3])

            file.write('%s %f %f %f\n' % (tr.stats.station,noise_lvl,p_sig_lvl,pkikp_sig_lvl))

            if indx<30:
#                plt.plot(ttt,arr+indx,'k')
#                print(tr.stats.npts,nbeg,nend,ttp1_diff,ttp2_diff)
                #print(indx,noise_lvl,p_sig_lvl,pkikp_sig_lvl)
                indx=indx+1
#            plt.xlim([-70,400])
    file.close();

def separate_LASA_station(sta_file,sta_dir):
    from pro_read_file import read_sta_file
    import numpy as np
    sta_num,st_names,st_lons,st_lats,st_deps=read_sta_file(sta_file);
    lon_median=np.median(st_lons); lat_median=np.median(st_lats);

    inner_file=sta_dir+'/LASA_inner_file.txt';
    outer_file=sta_dir+'/LASA_outer_file.txt';
    file_inner=open(inner_file,'w');
    file_outer=open(outer_file,'w');
    for ii in range(sta_num):
        stlo=st_lons[ii]; stla=st_lats[ii];
        dist_diff=np.sqrt((stlo-lon_median)**2+(stla-lat_median)**2)
        if dist_diff <0.43:
            file_inner.write('%s %8.4f %10.4f %4.1f\n' % (st_names[ii],st_lats[ii],st_lons[ii],0));
        else:
            file_outer.write('%s %8.4f %10.4f %4.1f\n' % (st_names[ii],st_lats[ii],st_lons[ii],0));

        #print(dist_diff)
    file_inner.close()
    file_outer.close()

def find_peak_trough(arr):
    import numpy as np

    peak_arr=[]; peak_indx=[];
    trough_arr=[]; trough_indx=[];

    ii_peak=0; ii_trough=0;
    for ii in range(1,len(arr)-1):
        if arr[ii]<arr[ii-1] and arr[ii]<arr[ii+1]:
            trough_arr.append(arr[ii]);
            trough_indx.append(ii);
            ii_trough+=1;
        if arr[ii]>arr[ii-1] and arr[ii]>arr[ii+1]:
            peak_arr.append(arr[ii]);
            peak_indx.append(ii);
            ii_peak+=1;
#    peak_arr=np.array(peak_arr);
#    peak_indx=np.array(peak_indx);
#    trough_arr=np.array(trough_arr);
#    trough_indx=np.array(trough_indx);
    return peak_arr,peak_indx,trough_arr,trough_indx

def init_grid(lon_ref=0,lat_ref=0,num_x=10,int_x=10,num_y=9,int_y=40,init_option=1):
    from geographiclib.geodesic import Geodesic
    # init_option: 0 for lat & lon grid and 1 for arc and azimuth grid
    if init_option==1:
        lon_arr = [[0 for x in range(num_y)] for y in range(num_x)]
        lat_arr = [[0 for x in range(num_y)] for y in range(num_x)]
        for iarc in range(num_x):
            for iaz in range(num_y):
                arc_tmp=(iarc+1)*int_x;
                az_tmp=iaz*int_y;
                aa=Geodesic.WGS84.ArcDirect(lat_ref, lon_ref, az_tmp, arc_tmp,outmask=1929)
                lon_arr[iarc][iaz]=aa['lon2'];
                lat_arr[iarc][iaz]=aa['lat2'];
    if init_option==0:
        lon_arr = [[0 for x in range(num_y)] for y in range(num_x)]
        lat_arr = [[0 for x in range(num_y)] for y in range(num_x)]
        for iarc in range(num_x):
            for iaz in range(num_y):
                lon_arr[iarc][iaz]=lon_ref+(-int(num_x/2)+iarc)*int_x;
                lat_arr[iarc][iaz]=lat_ref+(-int(num_y/2)+iaz)*int_y;
    return lon_arr, lat_arr

def layerxt_new(p,h,utop,ubot,imth):
    from math import sqrt, log, atan2
    ''' LAYERXT calculates dx and dt for ray in layer with linear velocity gradient

        Inputs:   p     =  horizontal slowness
             h     =  layer thickness
             utop  =  slowness at top of layer
             ubot  =  slowness at bottom of layer
             imth  =  interpolation method
                    imth = 1,  v(z) = 1/sqrt(a - 2*b*z)
                         = 2,  v(z) = a - b*z
                         = 3,  v(z) = a*exp(-b*z)  (better for spherical earth)
        Returns:  dx    =  range offset
             dt    =  travel time
             irtr  =  return code
                   = -1,  zero thickness layer
                   =  0,  ray turned above layer
                   =  1,  ray passed through layer
                   =  2,  ray turned within layer, 1 segment counted in dx,dt'''
# ray turned above layer
    #print(p,h,utop,ubot)
    if p >= utop:
        dx = 0
        dt = 0
        irtr = 0
        return dx,dt,irtr

#Zero thickness layer
    elif h == 0:
        dx = 0
        dt = 0
        irtr = -1
        return dx,dt,irtr

#    if p == 0:
#        dx = 0;
#        dt = 2.*h/(1./utop+1./ubot)
#        irtr = 1
#        return dx,dt,irtr

#Calculate some parameters
    u = utop
    y=u-p

    if (y <= 0.):   #complex vertical slowness
        dx=0.
        dt=0.
        irtr=0
        return dx,dt,irtr
    q=y*(u+p)
    qs=sqrt(q)

    # special function needed for integral at top of layer
    if imth == 2:
        y=u+qs
        if not(p==0):
            y=y/p
        qr=log(y)
    elif imth == 3:
        qr=atan2(qs,p)

    if imth == 1:
        b=-(utop**2-ubot**2)/(2.*h)
    elif imth == 2:
        vtop=1/utop
        vbot=1/ubot
        b=-(vtop-vbot)/h
    else:
        b=-log(ubot/utop)/h

    if b==0:  # constant velocity layer
        b=1./h
        etau=qs
        ex=p/qs
        irtr=1
        dx=ex/b
        dtau=etau/b
        dt=dtau+p*dx
        return dx,dt,irtr

    # intergral at upper limt 1/b factor ommited until end
    if imth == 1:
        etau=-q*qs/3.
        ex=-qs*p
    elif imth == 2:
        ex=qs/u
        etau=qr-ex
        if p!=0:
            ex=ex/p
    else:
        etau=qs-p*qr
        ex=qr

    # Check lower limit to see the turning point
    u=ubot
    if u<=p:
        irtr=2
        dx=ex/b
        dtau=etau/b
        dt=dtau+p*dx
        return dx,dt,irtr

    irtr=1
    q=(u-p)*(u+p)
    qs=sqrt(q)
    if imth == 1:
        etau=etau+q*qs/3.
        ex=ex+qs*p
    elif imth == 2:
        y=u+qs
        z=qs/u
        etau=etau+z
        if p!=0:
            y=y/p
            z=z/p
        qr=log(y)
        etau=etau-qr
        ex=ex-z
    else:
        qr=atan2(qs,p)
        etau=etau-qs+p*qr
        ex=ex-qr
    dx=ex/b
    dtau=etau/b
    dt=dtau+p*dx     #convert tau to t

    return dx,dt,irtr

def flatten_transform(zsph,vsph):
    ''' Calculate the flatten transformation from a spherical Earth.
        zf, vf = flat(zsph,vsph)'''
    import math
# Radius of the Earth
    a = 6371.0

    rsph = a - zsph
    zf = -a*math.log(rsph/float(a))
    vf = (a/float(rsph))*vsph
    return zf, vf

def cc_measure_tshift(tr1='',tr2='',tarr_beg=0,cc_twin=1,cc_len=0.3, cc_delta=0.1, cc_interp1d=10):
    '''
    compute the cc measurement of time shift and cc coefficients

    tr1, tr2:    input traces
    tarr_beg:    trace array time start time
    cc_twin:     time window for cross-correlation (s)
    cc_len:      time window shift to compute CC (fraction of whole time window)
    cc_delta:    time interval for cc (s)
    cc_interp1d: interpolation factor

    return values:
        cc_ttt:    time series of cc series
        cc_coeff:  cc coefficient series
        cc_tshift: cc measured time shift series
    '''
    import numpy as np
    from obspy.signal.cross_correlation import correlate
    from scipy.interpolate import interp1d
    import matplotlib.pyplot as plt

    ttt = np.arange(tr1.stats.npts)*tr1.stats.delta + tarr_beg

    if cc_delta < tr1.stats.delta or cc_delta < tr2.stats.delta:
        print('cc_delta is not greater than or equal to trace dt, quitting')
        return -1
    #%% specify cc time series
    cc_ttt_beg = ttt[0]  + cc_twin/2 + cc_delta;
    cc_ttt_end = ttt[-1] - cc_twin/2 - cc_delta;
    cc_ttt_npts = int((cc_ttt_end - cc_ttt_beg)/cc_delta) + 1;
    cc_ttt = cc_ttt_beg + np.arange(cc_ttt_npts)*cc_delta;

    # transfer cc time series to seismic trace time index
    cc_ttt_indx = np.around((cc_ttt-tarr_beg)/tr1.stats.delta) + 1

    #%% form cc time series & time shift arrays
    cc_coef = np.zeros(len(cc_ttt))
    cc_tshift = np.zeros(len(cc_ttt))

    ntwin = np.around(cc_twin/tr1.stats.delta) # number of cc seismogram samples

    #%% loop through all lags to correlate
    for ii in range(cc_ttt_npts):
        # generate indices of starting and ending points
        nbeg = int(round(cc_ttt_indx[ii] - int(ntwin/2)))
        nend = int(round(cc_ttt_indx[ii] + int(ntwin/2)))
        #print(nbeg,nend)
        ntt = nend - nbeg
        # extract seismogram segment for correlation
        arr1 = tr1.data[nbeg:nend]
        arr2 = tr2.data[nbeg:nend]

        tt = np.arange(ntt)*tr1.stats.delta
        f1 = interp1d(tt, arr1, kind='cubic') # interpolated function
        f2 = interp1d(tt, arr2, kind='cubic')

        cc_ntt = (ntt-1)*cc_interp1d + 1
        cc_tt = np.linspace(tt[0], tt[-1], cc_ntt) # interpolated cc time series
        arr1_new = f1(cc_tt) # interpolate seismogram
        arr2_new = f2(cc_tt) # interpolate seismogram

        nshift = int((cc_twin * cc_len) / tr1.stats.delta) # cc shift number
        nshift_new = nshift * cc_interp1d
        cc = correlate(arr1_new, arr2_new,nshift_new); # cc with interpolation

        cc_coef[ii] = np.amax(  cc)   # peak cc coefficient
        cc_max_indx = np.argmax(cc)   # cc max
        # convert cc shift from index to seconds
        cc_tshift[ii] = -(cc_max_indx-1-int(nshift_new))*tr1.stats.delta / cc_interp1d

    return cc_ttt, cc_coef, cc_tshift
