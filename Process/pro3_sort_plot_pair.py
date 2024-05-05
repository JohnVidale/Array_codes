#!/usr/bin/env python
# input is pair of repeated events from hinet, LASA, or NORSAR
# this program tapers, filters, selects range and SNR
# plots seismograms with traveltime curves, either raw or reduced against traveltimes
# outputs selected traces, "*sel.mseed"
# John Vidale 2/2019

def pro3pair(repeater = '0', stat_corr = 1, simple_taper = False, apply_SNR = False, SNR_thres = 1.7,
            phase1 = 'PKP', phase2 = 'PKiKP', phase3 = 'PKIKP', phase4 = 'pPKP',
            rel_time = 1, start_buff = -200, end_buff = 500, precursor_shift = 0, signal_dur = 0,
            plot_scale_fac = 0.05, corr_threshold = 0, off_center_shift = 0, win_norm = False, wind_buff = 0,
            zoom = False, Zstart_buff = 0, Zend_buff = 0, flip = False, trace_amp = 1,
            freq_min = 1, freq_max = 3, min_dist = 0, max_dist = 180, auto_dist = True, ARRAY = 0,
            ref_loc = False, ref_rad = 0.4, JST = False, fig_index = 100, do_interpolate = False,
            max_taper_length = 5., taper_frac = 0.2):


#%% Import functions
    from obspy import UTCDateTime
    from obspy import Stream
    from obspy import read
    from obspy.geodetics import gps2dist_azimuth
    import numpy as np
    import pandas as pd
    import os
    from obspy.taup import TauPyModel
    import matplotlib.pyplot as plt
    import time
    import math
    model = TauPyModel(model='ak135')

    import sys # don't show any warnings
    import warnings
    from termcolor import colored

    if not sys.warnoptions:
        warnings.simplefilter("ignore")

    print(colored('Running pro3a_sort_plot_pair', 'cyan'))
    start_time_wc = time.time()
    if ARRAY == 7 and stat_corr != 0:
        print(colored('No station corrections for global section', 'red'))
        stat_corr = 0
    if ARRAY == 7 and rel_time != 3:
        print(colored('Rel_time should be 3 for global section', 'red'))
        rel_time = 3

#%%  Set some parameters
    verbose = 0           # more output
    # rel_time   = 1      # timing is relative to a chosen phase, otherwise relative to OT
    # taper_frac = 0.05   # Fraction of window tapered on both ends
    noise_win_max = 10    # maximum length of noise window for SNR estimation, seconds
    # signal_dur = 5.       # signal length used in SNR calculation
    # precursor_shift = 3.  # shift in case SNR is used and signal starts before computed arrival time
    plot_tt = 1           # plot the traveltimes?
    do_decimate = False   # 0 if no decimation desired
    dec_factor = 10
    plot_sta_names = True
    dist_plot = False     # plot seismograms according to distance
    blowup_PKIKP = False
    # if ref_loc ==true,  use ref_rad        to filter station distance
    # if ref_loc ==false, use earthquake loc to filter station distance
    #    ref_rad = 0.4    # ° radius (°) set by input or at top

    if rel_time != 3:  # No reference distance needed if all traces aligned on phase
        if ARRAY == 0:
            ref_lat =    36.00  # °N, HiNet, around middle of Japan
            ref_lon =   139.00  # °E
        elif ARRAY == 1:        #  LASA
            ref_lat =    46.70  # °N keep only inner rings A-D if radius is 0.4°
            ref_lon =  -106.22  # °E
        elif ARRAY == 2:
            ref_lat =    38.00  # °N China
            ref_lon =   104.50  # °E
        elif ARRAY == 4:
            ref_lat =   -19.89  # °N Warramunga
            ref_lon =   134.42  # °E
        elif ARRAY == 5 or ARRAY == 99:
            ref_lat =    62.49  # °N Yellowknife
            ref_lon =  -114.60  # °E
        elif ARRAY == 6:
            ref_lat =    64.77  # °N ILAR
            ref_lon =  -146.89  # °E
        # ARRAY == 7 is global station set

    if rel_time == 0: # SNR requirement not implemented for unaligned traces
        SNR_thres = 0

#%% Get saved event info, also used to name files

    def search_df(df, column, value, partial_match=True):
        df = df.astype({column:'string'})
        if partial_match:
            return df.loc[df[column].str.contains(value, na=False)]
        else:
            return df.loc[df[column] == value]

    # look up pair of earthquakes and time shifts in pairs
    df = pd.read_excel('/Users/vidale/Documents/GitHub/Array_codes/Files/ICevents_full.xlsx', sheet_name='pairs')
    lines0       = search_df(df,'label'      ,repeater,partial_match=True)
    eq_num1      = lines0.index1.iloc[0]
    eq_num2      = lines0.index2.iloc[0]
    tshift       = lines0.tshift.iloc[0]
    shift_both   = lines0.shift_both.iloc[0]

    if ARRAY == 5:  # additional shifts only for YKA
        Y_shift      = lines0.Y_shift.iloc[0]
        shift_bothY  = lines0.shift_bothY.iloc[0]
        tshift     = tshift + Y_shift
        shift_both = shift_both + shift_bothY


    # read origin times for that pair in events
    df = pd.read_excel('/Users/vidale/Documents/GitHub/Array_codes/Files/ICevents_full.xlsx', sheet_name='events')
    lines1 = search_df(df,'INDEX',str(eq_num1),partial_match=True)
    lines2 = search_df(df,'INDEX',str(eq_num2),partial_match=True)

    time1 = lines1.TIME.iloc[0]
    time2 = lines2.TIME.iloc[0]
    t1           = UTCDateTime(time1)
    t2           = UTCDateTime(time2)


    #  new lines to match more specific naming
    date_label1  = time1[0:10]
    date_label2  = time2[0:10]
    year_label1  = time1[0:4]
    year_label2  = time2[0:4]
    year_short_label1  = time1[2:4]

    year_short_label2  = time2[2:4]
    month_label1   = time1[5:7]
    month_label2   = time2[5:7]
    day_label1     = time1[8:10]
    day_label2     = time2[8:10]
    hour_label1    = time1[11:13]
    hour_label2    = time2[11:13]
    minute_label1  = time1[14:16]
    minute_label2  = time2[14:16]

    year1        = time1[0:4]
    year2        = time2[0:4]
    ev_lat1      = float(lines1.LAT)
    ev_lat2      = float(lines2.LAT)
    ev_lon1      = float(lines1.LON)
    ev_lon2      = float(lines2.LON)
    ev_depth1    = float(lines1.DEP)
    ev_depth2    = float(lines2.DEP)
    if abs(ev_depth1 - ev_depth2) > 0.001:  # depths come from 2nd sheet of parameters, should be equal
        print(colored('depth 1 ' + str(ev_depth1) + ' does not equal depth2 ' + str(ev_depth2), 'red'))
    print('1st event: date_label ' + date_label1 + ' time ' + str(t1) + ' lat '
       + str(ev_lat1) + ' lon ' + str( ev_lon1) + ' depth ' + str(ev_depth1))
    print('2nd event: date_label ' + date_label2 + ' time ' + str(t2) + ' lat '
       + str(ev_lat2) + ' lon ' + str( ev_lon2) + ' depth ' + str(ev_depth2))

#%% Get station location file
    if stat_corr == 1 or stat_corr == 2 or stat_corr == 3:  # load static terms

        if ARRAY == 0:
            if stat_corr == 1: # standard set
                sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/Stations/sta_statics_hinet.txt'
            elif stat_corr == 2: # custom set made for SSI to Hi-Net PKiKP
                sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/Stations/SSI_HiNet_statics.txt'
            elif stat_corr == 3: # custom set made for Kawa PKiKP for Hi-NetP
                sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/Stations/Kawa_HiNet_statics.txt'

        elif ARRAY == 1:
            sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/Stations/sta_statics_LASA.txt'
        elif ARRAY == 2:
            sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/Stations/sta_statics_ch.txt'
        elif ARRAY == 5 or ARRAY == 99:
            sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/Stations/sta_statics_YKA.txt'
        elif ARRAY == 6:
            sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/Stations/sta_statics_ILAR.txt'
        else:
            print(colored('ERROR - No station corrections in program for array ' + str(ARRAY), 'cyan'))
            sys.exit(-1)


        with open(sta_file, 'r') as file:
            lines = file.readlines()
        print(str(len(lines)) + ' stations read from ' + sta_file)
        # Load station coords into arrays
        station_index = range(len(lines))
        st_names = []
        st_dist  = []
        st_lats  = []
        st_lons  = []
        st_shift = []
        st_corr  = []
        for ii in station_index:
            line = lines[ii]
            split_line = line.split()
            st_names.append(split_line[0])
            st_dist.append( split_line[1])
            st_lats.append( split_line[2])
            st_lons.append( split_line[3])
            st_shift.append(split_line[4])
            st_corr.append( split_line[5])
    else: # no static terms, always true for NORSAR
        if ARRAY == 0: # Hinet set
            sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/Stations/sta_hinet.txt'
        elif ARRAY == 1: #         LASA set
            sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/Stations/sta_LASA.txt'
        elif ARRAY == 2: #         LASA set
            sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/Stations/sta_ch.txt'
        elif ARRAY == 3: #         NORSAR set
            sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/Stations/sta_NORSAR.txt'
        elif ARRAY == 4: #         Warramunga set
            sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/Stations/sta_AU_WR.txt'
        elif ARRAY == 5 or ARRAY == 99: #         Yellowknife set
            sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/Stations/sta_CN_YK.txt'
        elif ARRAY == 6: #         ILAR set
            sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/Stations/sta_ILAR.txt'
        elif ARRAY == 7: #         Global set
            sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/Stations/sta_global.txt'
        with open(sta_file, 'r') as file:
            lines = file.readlines()

        print(str(len(lines)) + ' stations read from ' + sta_file)
        # Load station coords into arrays
        station_index = range(len(lines))
        st_names = []
        st_lats  = []
        st_lons  = []
        for ii in station_index:
            line = lines[ii]
            split_line = line.split()
            st_names.append(split_line[0])
            st_lats.append( split_line[1])
            st_lons.append( split_line[2])
    if ARRAY == 0:  # shorten and make upper case Hi-net station names to match station list
        for ii in station_index:
            this_name = st_names[ii]
            this_name_truc = this_name[0:5]
            st_names[ii]  = this_name_truc.upper()

#%% Is noise estimation window too short compared to taper?
    totalt = end_buff - start_buff                    # window duration
    noise_time_skipped = taper_frac * totalt          # set ignored SNR interval to be taper
    noise_time_skipped = min(noise_time_skipped,10.0) # cap ignored SNR interval at 10s
    if simple_taper == False:                             # lengthen taper_frac if too short
        if noise_time_skipped >= 0.5 * (-start_buff):
            print('Specified taper of ' + str(taper_frac * totalt) +
               ' is not big enough compared to available noise estimation window ' +
               str(-start_buff - noise_time_skipped) + '. May not work well.')
            old_taper_frac = taper_frac
            taper_frac = -0.5*start_buff/totalt
            print('Taper reset from ' + str(old_taper_frac * totalt) + ' to '
               + str(taper_frac * totalt) + ' seconds.')

#%% Load waveforms and decimate to 10 sps
    st1 = Stream()
    st2 = Stream()

    if ARRAY == 0 or ARRAY == 2:
        mseed_name1 = year_short_label1 + '-' + month_label1 + '-' + day_label1
        mseed_name2 = year_short_label2 + '-' + month_label2 + '-' + day_label2
    elif ARRAY == 1:
        mseed_name1 = year_short_label1 + month_label1 + day_label1 + '_' + hour_label1 + minute_label1
        mseed_name2 = year_short_label2 + month_label2 + day_label2 + '_' + hour_label2 + minute_label2
    else:
        mseed_name1 = year_label1 + month_label1 + day_label1 + '_' + hour_label1 + minute_label1
        mseed_name2 = year_label2 + month_label2 + day_label2 + '_' + hour_label2 + minute_label2

    if ARRAY == 0:
        fname1     = '/Users/vidale/Documents/Research/IC/Mseed/HiNet/HD20' + mseed_name1 + '.mseed'
        fname2     = '/Users/vidale/Documents/Research/IC/Mseed/HiNet/HD20' + mseed_name2 + '.mseed'
    elif ARRAY == 1:
        fname1     = '/Users/vidale/Documents/Research/IC/Mseed/LASA/L' + mseed_name1 + '.mseed'
        fname2     = '/Users/vidale/Documents/Research/IC/Mseed/LASA/L' + mseed_name2 + '.mseed'
    elif ARRAY == 2:
        fname1     = '/Users/vidale/Documents/Research/IC/Mseed/China/HD20' + mseed_name1 + '.mseed'
        fname2     = '/Users/vidale/Documents/Research/IC/Mseed/China/HD20' + mseed_name2 + '.mseed'
    elif ARRAY == 4:
        fname1     = '/Users/vidale/Documents/Research/IC/Mseed/WRA/' + mseed_name1 + '.mseed'
        fname2     = '/Users/vidale/Documents/Research/IC/Mseed/WRA/' + mseed_name2 + '.mseed'
    elif ARRAY == 5:
        fname1     = '/Users/vidale/Documents/Research/IC/Mseed/YKA/' + mseed_name1 + '.mseed'
        fname2     = '/Users/vidale/Documents/Research/IC/Mseed/YKA/' + mseed_name2 + '.mseed'
    elif ARRAY == 6:
        fname1     = '/Users/vidale/Documents/Research/IC/Mseed/ILAR/' + mseed_name1 + '.mseed'
        fname2     = '/Users/vidale/Documents/Research/IC/Mseed/ILAR/' + mseed_name2 + '.mseed'
    elif ARRAY == 7:
        fname1     = '/Users/vidale/Documents/Research/IC/Mseed/Global/' + mseed_name1 + '.mseed'
        fname2     = '/Users/vidale/Documents/Research/IC/Mseed/Global/' + mseed_name2 + '.mseed'
    elif ARRAY == 99: # compare two versions of files
        fname1     = '/Users/vidale/Documents/Research/IC/Mseed/YKA/'  + mseed_name1 + '.mseed'
        fname2     = '/Users/vidale/Documents/Research/IC/Mseed/YKA2/' + mseed_name2 + '.mseed'

    st1=read(fname1)
    st2=read(fname2)
    for tr in st1: # apply a shift to second event, and flip if desired
        tr.stats.starttime = tr.stats.starttime -  shift_both
        if flip:
            tr.data = tr.data * (-1.0)
    for tr in st2: # apply a shift to second event
        if flip:
            tr.data = tr.data * (-1.0)
        tr.stats.starttime = tr.stats.starttime - (shift_both + tshift)
        # if ARRAY == 5 and eq_num2 > 746:  # kludge to minimize effect of YKA time shift in late 2013
        #     tr.stats.starttime = tr.stats.starttime - 0.17

    if do_decimate:
        st1.decimate(dec_factor, no_filter=True)
        st2.decimate(dec_factor, no_filter=True)

    if do_interpolate:
        st1.interpolate(sampling_rate = 40)
        st2.interpolate(sampling_rate = 40)

    print(f'1st trace for event 1 has {len(st1[0].data)} time pts and {st1[0].stats.delta} dt, which is {len(st1[0].data)*st1[0].stats.delta:.1f} s')
    print(f'1st trace for event 2 has {len(st2[0].data)} time pts and {st2[0].stats.delta} dt, which is {len(st2[0].data)*st2[0].stats.delta:.1f} s')
    print('st1 has ' + str(len(st1)) + ' traces')
    print('st2 has ' + str(len(st2)) + ' traces')

    print('1: ' + st1[0].stats.station + ' 2: ' + st2[0].stats.station)
    if JST == False:
        print('1st trace starts at ' + str(st1[0].stats.starttime) + ', event at ' + str(t1))
        print('2nd trace starts at ' + str(st2[0].stats.starttime) + ', event at ' + str(t2))
    if JST == True:  # time shift actually applied in loop below
        temp_time1 = st1[0].stats.starttime - 9*60*60
        temp_time2 = st2[0].stats.starttime - 9*60*60
        print('JST -> UTC 1st trace ' + str(temp_time1) + ', event at ' + str(t1))
        print('JST -> UTC 2nd trace ' + str(temp_time2) + ', event at ' + str(t2))

#%% Select by distance, window and adjust start time to align picked times
    st_pickalign1 = Stream()
    st_pickalign2 = Stream()
    tra1_in_range  = 0
    tra1_sta_found = 0
    nodata1        = 0
    tra2_in_range  = 0
    tra2_sta_found = 0
    nodata2        = 0
    min_dist_auto = 180
    max_dist_auto = 0
    min_time_plot =  1000000
    max_time_plot = -1000000

    if rel_time != 3:  # No reference distance needed if all traces aligned on phase
        ref_distance = gps2dist_azimuth(ref_lat,ref_lon,ev_lat1,ev_lon1)
        ref1_dist  = ref_distance[0]/(1000*111)
        arrivals_ref   = model.get_travel_times(source_depth_in_km=ev_depth1,distance_in_degree=ref1_dist, phase_list=[phase1])

        if(len(arrivals_ref) == 0 and ref1_dist < 10 and phase1 == 'P'):  # in case first arrival might be upgoing P, which is 'p'
            arrivals_ref   = model.get_travel_times(source_depth_in_km=ev_depth1,distance_in_degree=ref1_dist, phase_list='p')
        if(len(arrivals_ref) == 0):
            print('model.get_travel_times failed: dist, phase  ' + str(ref1_dist) + '   ' + phase1)
            sys.exit(-1)
        atime_ref = arrivals_ref[0].time  # phase arrival time at reference distance
        if phase1 == 'PKP':  # use last element; first element might be BC, with limited range and thus step in [0] element
            atime_rayp = arrivals_ref[-1].ray_param  # ray parameter (s/radian) at reference distance
            atime_ref  = arrivals_ref[-1].time       # phase arrival time at reference distance
        else:  # use first element
            atime_rayp = arrivals_ref[0].ray_param
            atime_ref  = arrivals_ref[0].time
        ref_slow  = atime_rayp * 2 * np.pi / 360. # convert to s/°

#%% -- First event
    for tr in st1: # find lat-lon from list, chop, statics, traces one by one
        if float(year1) < 1970: # fix the damn 1969 -> 2069 bug in Gibbon's LASA data
            temp_t = str(tr.stats.starttime)
            temp_tt = '19' + temp_t[2:]
            tr.stats.starttime = UTCDateTime(temp_tt)
        if JST: # if necessary, convert JST -> UTC, time in Greenwich 9 hours earlier than Japan
            tr.stats.starttime = tr.stats.starttime - 9*60*60
        if tr.stats.station in st_names:  # find station in station list
            ii = st_names.index(tr.stats.station)
            tra1_sta_found += 1

            if stat_corr == 0 or float(st_corr[ii]) > corr_threshold: # if using statics, reject low correlations
                stalat = float(st_lats[ii]) # look up lat & lon again to find distance
                stalon = float(st_lons[ii])

                distance = gps2dist_azimuth(stalat,stalon,ev_lat1,ev_lon1) # Get traveltimes again, hard to store
                tr.stats.distance=distance[0] # distance in km
                dist = distance[0]/(1000*111)

                in_range = 0  # flag for whether this trace goes into stack
                rejector = False  # flag in case model.get_travel_times fails

                if ref_loc == False:  # check whether trace is in distance range from earthquake
                    if min_dist < dist and dist < max_dist:
                        in_range = 1
                        tra1_in_range += 1
                elif ref_loc == True:  # alternately, check whether trace is close enough to ref_location
                    ref_distance = gps2dist_azimuth(ref_lat,ref_lon,stalat,stalon)
                    ref_degree_dist = ref_distance[0]/(1000*111)
                    if ref_degree_dist < ref_rad:
                        in_range = 1
                        tra1_in_range += 1

                if in_range == 1:   # trace fulfills the specified criteria for being in range
                    s_t = t1 + start_buff
                    e_t = t1 + end_buff

#%% ---- Apply static
                    if stat_corr == 1 or stat_corr == 2: # apply static station corrections
                        tr.stats.starttime -= float(st_shift[ii])
                    # if rel_time == 0  or rel_time == 3:  #  use absolute time, ignore time of chose phase
                    if rel_time == 0:  #  use absolute time, ignore time of chose phase
                        tr.trim(starttime=s_t,endtime = e_t)
                    else:              # shift relative to a chosen phase
                        arrivals_each = model.get_travel_times(source_depth_in_km=ev_depth1,distance_in_degree=dist,phase_list=[phase1])
                        if(len(arrivals_each) == 0):
                            if(dist < 10 and phase1 == 'P'):  # in case first arrival is upgoing P, try 'p'
                                arrivals_each = model.get_travel_times(source_depth_in_km=ev_depth1,distance_in_degree=dist,phase_list='p')
                        if(len(arrivals_each) == 0):  # did it still fail?
#%% ---- Can't calculate traveltime
                                print('Event 1: model.get_travel_times failed: dist, depth, phase  ' + tr.stats.station + '   ' + str(dist) + '   ' +
                                      '   ' + str(ev_depth1) + '   ' + phase1)
                                tra1_in_range -= 1 # don't count this trace after all
                                rejector = True
                        else:
                            #  reject if there are NANs, except patch last-point NAN
                            count_nan = 0
                            for ii in range(len(tr.data)):  # find where are the NANs, fix if last point, otherwise reject trace
                                if math.isnan(tr.data[ii]):
                                    # print(colored(tr.stats.station + ' has NAN at point ' + str(ii) + ' out of ' + str(len(tr.data)), 'yellow'))
                                    if ii == len(tr.data) - 1:
                                        tr.data[ii] = tr.data[ii-1]
                                        print(colored(tr.stats.station + ' has NAN at last point, fixed', 'yellow'))
                                    else:
                                        count_nan += 1
                            if count_nan != 0:
                                tra1_in_range -= 1 # don't count this trace after all
                                rejector = True
                                print(colored(tr.stats.station + ' has ' + str(count_nan) + ' NANs, trace rejected', 'yellow'))

                            else:
                                # atime_each is phase arrival on this trace
                                # atime_ref  is phase arrival on reference trace
                                if phase1 == 'PKP':
                                    atime_each = arrivals_each[-1].time  # use last instance for AB branch
                                else:
                                    atime_each = arrivals_each[0].time
                                # start time has a shift proportional to ref_dist at phase slowness at ref_dist
                                # preserves slowness of arrivals at ref_dist, sharp even if tt curve has curvature
                                # uniform relative window limits around phase
                                if rel_time == 1:
                                    s_t += atime_each
                                    e_t += atime_each
                                    tr.trim( starttime = s_t, endtime = e_t)
                                    tr.stats.starttime -= atime_each - ((dist-ref1_dist) * ref_slow + off_center_shift)
                                # window time sets 0 at chosen phase only on reference trace
                                # keeps true slowness, but blurs arcuate slowness curves
                                # uniform relative window limits around phase
                                elif rel_time == 2:
                                    s_t += atime_each
                                    e_t += atime_each
                                    tr.trim( starttime = s_t, endtime = e_t)
                                    tr.stats.starttime -= atime_ref
                                # window time sets 0 at chosen phase on all traces
                                # set phase to 0 slowness
                                # uniform relative window limits around phase
                                elif rel_time == 3:
                                    s_t += atime_each
                                    e_t += atime_each
                                    tr.trim( starttime = s_t, endtime = e_t)
                                    tr.stats.starttime -= atime_each
                                    # print(f'data length after {len(tr.data)}')
                                # window time sets 0 at chosen phase only on reference trace
                                # set phase to 0 slowness
                                # use same absolute window around chosen phase for all stations
                                elif rel_time == 4:
                                    s_t += atime_ref
                                    e_t += atime_ref
                                    tr.trim( starttime = s_t, endtime = e_t)
                                    tr.stats.starttime -= atime_ref
                                else:
                                    sys.exit('Invalid rel_time parameter, must be integer 0 to 4')
#%% ---- Reject if not in time window
                    if len(tr.data) > 0 and rejector == False:
                        st_pickalign1 += tr
                    else:
                        nodata1 += 1

#%% ---- Reject if not in station (static) list
        # else:
        #     print(tr.stats.station + ' not found in station list with statics')
        #     print('    ' + str(tra1_in_range) + '  traces in range')
        #     print('    ' + str(len(st_pickalign1)) + '  traces after alignment and correlation selection')
        #     print('    ' + str(nodata1) + '  traces with no data')
    #            sys.exit()

#%% -- Second event
    for tr in st2: # find lat-lon from list, chop, statics, traces one by one
        if float(year2) < 1970: # fix the damn 1969 -> 2069 bug in Gibbon's LASA data
            temp_t = str(tr.stats.starttime)
            temp_tt = '19' + temp_t[2:]
            tr.stats.starttime = UTCDateTime(temp_tt)
        if JST == True and eq_num2 != '180': # if necessary, convert JST -> UTC, time in Greenwich 9 hours earlier than Japan
            tr.stats.starttime = tr.stats.starttime - 9*60*60
        if tr.stats.station in st_names:  # find station in station list
            # print(' station name ' + tr.stats.station)
            ii = st_names.index(tr.stats.station)
            tra2_sta_found += 1

            if stat_corr == 0 or float(st_corr[ii]) > corr_threshold: # if using statics, reject low correlations
                stalat = float(st_lats[ii]) # look up lat & lon again to find distance
                stalon = float(st_lons[ii])

                distance = gps2dist_azimuth(stalat,stalon,ev_lat2,ev_lon2) # Get traveltimes again, hard to store
                tr.stats.distance=distance[0] # distance in km
                dist = distance[0]/(1000*111)

                rejector = False  # flag in case model.get_travel_times fails
                in_range = 0  # flag for whether this trace goes into stack
                if ref_loc == False:  # check whether trace is in distance range from earthquake
                    if min_dist < dist and dist < max_dist:
                        in_range = 1
                        tra2_in_range += 1
                elif ref_loc == True:  # alternately, check whether trace is close enough to ref_location
                    ref_distance = gps2dist_azimuth(ref_lat,ref_lon,stalat,stalon)
                    ref_degree_dist = ref_distance[0]/(1000*111)
                    if ref_degree_dist < ref_rad:
                        in_range = 1
                        tra2_in_range += 1
                if in_range == 1:   # trace fulfills the specified criteria for being in range
                    s_t = t2 + start_buff
                    e_t = t2 + end_buff

#%% ---- Apply static
                    if stat_corr == 1 or stat_corr == 2: # apply static station corrections
                        tr.stats.starttime -= float(st_shift[ii])
                    # if rel_time  == 0 or rel_time == 3:  #  use absolute time, ignore time of chose phase
                    if rel_time  == 0:  #  use absolute time, ignore time of chose phase
                        tr.trim(starttime=s_t,endtime = e_t)
                    else:              # shift relative to a chosen phase
                        arrivals_each = model.get_travel_times(source_depth_in_km=ev_depth2,distance_in_degree=dist,phase_list=[phase1])
                        if(len(arrivals_each) == 0):
                            if(dist < 10 and phase1 == 'P'):  # in case first arrival is upgoing P, try 'p'
                                arrivals_each = model.get_travel_times(source_depth_in_km=ev_depth2,distance_in_degree=dist,phase_list='p')
                        if(len(arrivals_each) == 0):  # did it still fail?
#%% ---- Can't calculate traveltime
                                print('Event 2: model.get_travel_times failed: dist, depth, phase  ' + tr.stats.station + '   ' + str(dist) + '   ' +
                                      '   ' + str(ev_depth2) + '   ' + phase1)
                                tra2_in_range -= 1 # don't count this trace after all
                                rejector = True
                        else:
                            #  reject if there are NANs, except patch last-point NAN
                            count_nan = 0
                            for ii in range(len(tr.data)):  # find where are the NANs, fix if last point, otherwise reject trace
                                if math.isnan(tr.data[ii]):
                                    # print(colored(tr.stats.station + ' has NAN at point ' + str(ii) + ' out of ' + str(len(tr.data)), 'yellow'))
                                    if ii == len(tr.data) - 1:
                                        tr.data[ii] = tr.data[ii-1]
                                        print(colored(tr.stats.station + ' has NAN at last point, fixed', 'yellow'))
                                    else:
                                        count_nan += 1
                            if count_nan != 0:
                                tra1_in_range -= 1 # don't count this trace after all
                                rejector = True
                                print(colored(tr.stats.station + ' has ' + str(count_nan) + ' NANs, trace rejected', 'yellow'))

                            else:
                            # atime_each is phase arrival on this trace
                            # atime_ref  is phase arrival on reference trace
                                if phase1 == 'PKP':
                                    atime_each = arrivals_each[-1].time  # use last instance for AB branch
                                else:
                                    atime_each = arrivals_each[0].time
                            # start time has a shift proportional to ref_dist at phase slowness at ref_dist
                            # preserves slowness of arrivals at ref_dist, sharp even if tt curve has curvature
                            # uniform relative window limits around phase
                            if rel_time == 1:
                                s_t += atime_each
                                e_t += atime_each
                                tr.trim( starttime = s_t, endtime = e_t)
                                tr.stats.starttime -= atime_each - ((dist-ref1_dist) * ref_slow + off_center_shift)
                            # window time sets 0 at chosen phase only on reference trace
                            # keeps true slowness, but blurs arcuate slowness curves
                            # uniform relative window limits around phase
                            elif rel_time == 2:
                                s_t += atime_each
                                e_t += atime_each
                                tr.trim( starttime = s_t, endtime = e_t)
                                tr.stats.starttime -= atime_ref
                            # window time sets 0 at chosen phase on all traces
                            # set phase to 0 slowness
                            # uniform relative window limits around phase
                            elif rel_time == 3:
                                s_t += atime_each
                                e_t += atime_each
                                tr.trim( starttime = s_t, endtime = e_t)
                                tr.stats.starttime -= atime_each
                            # window time sets 0 at chosen phase only on reference trace
                            # set phase to 0 slowness
                            # use same absolute window around chosen phase for all stations
                            elif rel_time == 4:
                                s_t += atime_ref
                                e_t += atime_ref
                                tr.trim( starttime = s_t, endtime = e_t)
                                tr.stats.starttime -= atime_ref
                            else:
                                sys.exit('Invalid rel_time parameter, must be integer 0 to 4')
#%% ---- Reject if not in time window
                    if len(tr.data) > 0 and rejector == False:
                        st_pickalign2 += tr
                    else:
                        nodata2 += 1

#%% ---- Reject if not in station (static) list
        # else:
        #     print(tr.stats.station + ' not found in station list with statics')
        #     print('    ' + str(tra1_in_range) + '  traces in range')
        #     print('    ' + str(len(st_pickalign1)) + '  traces after alignment and correlation selection')
        #     print('    ' + str(nodata1) + '  traces with no data')
#            sys.exit()

    print('After alignment + range and correlation selection')
    print('1st event, Traces found: ' + str(tra1_sta_found) + ', Traces in range: ' + str(tra1_in_range) + ', Traces with no data: ' + str(nodata1))
    print('2nd event, Traces found: ' + str(tra2_sta_found) + ', Traces in range: ' + str(tra2_in_range) + ', Traces with no data: ' + str(nodata2))

    if rel_time != 3:  # No reference distance needed if all traces aligned on phase
        print(f'ref1_distance  {ref1_dist:.3f}  relative start time  {atime_ref:.3f}')
    if ref_loc == True:
        print('ref_loc == True, ref_lat: ' + str(ref_lat) + ' ref_lon: ' + str(ref_lon))
    print(f'last station: distance {dist:.3f}  last station lat: {stalat:.3f}   last station lon: {stalon:.3f}')


    #print(st) # at length
    if verbose:
        print(st1.__str__(extended=True))
        print(st2.__str__(extended=True))
        if rel_time == 1:
            print(st_pickalign1.__str__(extended=True))
            print(st_pickalign2.__str__(extended=True))

#%%  Detrend, taper, filter, optional window normalization
    print('Taper fraction is ' + str(taper_frac) + ' Max taper length is ' + str(max_taper_length) + ' bandpass is ' + str(freq_min) + ' to ' + str(freq_max))
    st_pickalign1.detrend(type='simple')
    st_pickalign2.detrend(type='simple')
    st_pickalign1.taper(taper_frac, max_length = max_taper_length)
    st_pickalign2.taper(taper_frac, max_length = max_taper_length)
    # zerophase filter
    # st_pickalign1.filter('bandpass', freqmin=freq_min, freqmax=freq_max, corners=4, zerophase=True)
    # st_pickalign2.filter('bandpass', freqmin=freq_min, freqmax=freq_max, corners=4, zerophase=True)
    # causal filter
    st_pickalign1.filter('bandpass', freqmin=freq_min, freqmax=freq_max, corners=4, zerophase=False)
    st_pickalign2.filter('bandpass', freqmin=freq_min, freqmax=freq_max, corners=4, zerophase=False)
    # st_pickalign1.taper(taper_frac, max_length = max_taper_length) # double tapering not really necessary
    # st_pickalign2.taper(taper_frac, max_length = max_taper_length)

    if win_norm == True:  # normalize in analysis window
        s_t = t1 + start_buff + wind_buff
        e_t = t1 + end_buff - wind_buff
        for tr in st_pickalign1:
            snippet = tr.copy()
            snippet.trim( starttime = s_t, endtime = e_t)
            smax = abs(snippet.max())
            tr.data = tr.data/smax

        s_t = t2 + start_buff + wind_buff
        e_t = t2 + end_buff - wind_buff
        for tr in st_pickalign2:
            snippet = tr.copy()
            snippet.trim( starttime = s_t, endtime = e_t)
            smax = abs(snippet.max())
            tr.data = tr.data/smax

#%%  Check station is there for both events and impose SNR threshold
    st1good = Stream()
    st2good = Stream()
    if apply_SNR == True:
        print(f'SNR parameters: start buffer {start_buff:.1f} end buffer {end_buff:.1f} pre_shift {precursor_shift:.1f} sig duration {signal_dur:.1f}')
    for tr1 in st_pickalign1:
        # print(tr1.stats.station + ' ' + tr1.stats.station[0:2] + ' ' + tr1.stats.station[3:5])
        if ARRAY == 5 and tr1.stats.station[2] == 'A':  # overcome renaming of YKA stations in late 2013
            tr1.stats.station = tr1.stats.station[0:2] + tr1.stats.station[3:5]
            # tr1.stats.station = tr1.stats.station[0:1] + tr1.stats.station[3:4]
        # print(tr1.stats.station)

        for tr2 in st_pickalign2:
            if ARRAY == 5 and tr2.stats.station[2] == 'A':  # overcome renaming of YKA stations in late 2013
                tr2.stats.station = tr2.stats.station[0:2] + tr2.stats.station[3:5]

            if ((tr1.stats.network  == tr2.stats.network) and
                (tr1.stats.station  == tr2.stats.station)): # if this is a match, then station is on both lists
                if apply_SNR == False:
                    st_name = tr1.stats.station
                    doit = True # weed out bad stations for repeaters
                    if st_name == 'YKB0':
                        doit = False
                    if (eq_num1 == 701 or eq_num2 == 701) and (st_name == 'YKB1'):
                        doit = False
                    if (eq_num1 == 702 or eq_num2 == 702) and (st_name == 'YKR1' or st_name == 'YKR2' or st_name == 'YKR3' or st_name == 'YKR4'):
                        doit = False
                    if (eq_num1 == 705 or eq_num2 == 705) and (st_name == 'YKR2'):
                        doit = False
                    if (eq_num1 == 706 or eq_num2 == 706) and (st_name == 'YKB8'):
                        doit = False
                    if (eq_num1 == 707 or eq_num2 == 707) and (st_name == 'YKB4' or st_name == 'YKR8'):
                        doit = False
                    if (eq_num1 == 709 or eq_num2 == 709) and (st_name == 'YKB1'):
                        doit = False
                    if (eq_num1 == 711 or eq_num2 == 711) and (st_name == 'YKR2'):
                        doit = False
                    if (eq_num1 == 716 or eq_num2 == 716) and (st_name == 'IL13' or st_name == 'IL03' or st_name == 'IL08' or st_name == 'IL18'):
                        doit = False
                    if (eq_num1 == 719 or eq_num2 == 719) and (st_name == 'YKR1'):
                        doit = False
                    if (eq_num1 == 722 or eq_num2 == 722) and (st_name == 'YKR1' or st_name == 'YKR2' or st_name == 'YKR3'):
                        doit = False
                    if (eq_num1 == 723 or eq_num2 == 723) and (st_name == 'IL01' or st_name == 'YKR1' or st_name == 'YKR3'):
                        doit = False
                    if (eq_num1 == 724 or eq_num2 == 724) and (st_name == 'YKB1' or st_name == 'YKR5'):
                        doit = False
                    if (eq_num1 == 734 or eq_num2 == 734) and (st_name == 'YKR6'):
                        doit = False
                    if (eq_num1 == 739 or eq_num2 == 739) and (st_name == 'YKR1' or st_name == 'YKR6' or st_name == 'YKB1' or st_name == 'YKR7'):
                        doit = False
                    if (eq_num1 == 740 or eq_num2 == 740) and (st_name == 'YKB1' or st_name == 'YKB2' or st_name == 'YKR2' or st_name == 'YKR3'
                                                            or st_name == 'IL17' or st_name == 'IL16' or st_name == 'IL14' or st_name == 'IL13'
                                                            or st_name == 'IL02' or st_name == 'YKR1'):
                        doit = False
                    if (eq_num1 == 741 or eq_num2 == 741) and (st_name == 'IL14' or st_name == 'IL15' or st_name == 'YKB1'):
                        doit = False
                    if (eq_num1 == 742 or eq_num2 == 742) and (st_name == 'YKB1'):
                        doit = False
                    if (eq_num1 == 743 or eq_num2 == 743) and (st_name == 'YKB2'):
                        doit = False
                    if (eq_num1 == 745 or eq_num2 == 745) and (st_name == 'YKR2'):
                        doit = False
                    if (eq_num1 == 746 or eq_num2 == 746) and (st_name == 'YKR7'):
                        doit = False
                    if (eq_num1 == 748 or eq_num2 == 748) and (st_name == 'YKB1'):
                        doit = False
                    if (eq_num1 == 750 or eq_num2 == 750) and (st_name == 'YKR7'):
                        doit = False
                    if (eq_num1 == 757 or eq_num2 == 757) and (st_name == 'YKR2'):
                        doit = False
                    if (eq_num1 == 759 or eq_num2 == 759) and (st_name == 'YKB2'):
                        doit = False
                    if (eq_num1 == 800 or eq_num2 == 800) and (st_name == 'YKB9'):
                        doit = False
                    if (eq_num1 == 802 or eq_num2 == 802) and (st_name == 'YKR3' or st_name == 'YKR2' or st_name == 'YKR1'
                                                            or st_name == 'YKB1' or st_name == 'YKB7' or st_name == 'YKB2'):
                        doit = False
                    if (eq_num1 == 808 or eq_num2 == 808) and (st_name == 'YKR9'):
                        doit = False
                    if (eq_num1 == 811 or eq_num2 == 811) and (st_name == 'YKR5' or st_name == 'YKR2'):
                        doit = False
                    if (eq_num1 == 812 or eq_num2 == 812) and (st_name == 'YKB1' or st_name == 'YKB4' or st_name == 'YKR3'):
                        doit = False
                    if (eq_num1 == 813 or eq_num2 == 813) and (st_name == 'IL18'):
                        doit = False
                    if (eq_num1 == 814 or eq_num2 == 814) and (st_name == 'YKR1' or st_name == 'YKR2'):
                        doit = False
                    if (eq_num1 == 816 or eq_num2 == 816) and (st_name == 'YKR3' or st_name == 'YKB1'):
                        doit = False
                    if (eq_num1 == 818 or eq_num2 == 818) and (st_name == 'YKB1'):
                        doit = False
                    if (eq_num1 == 819 or eq_num2 == 819) and (st_name == 'YKB1'):
                        doit = False
                    if (eq_num1 == 823 or eq_num2 == 823) and (st_name == 'YKB1'):
                        doit = False
                    if (eq_num1 == 824 or eq_num2 == 824) and (st_name == 'YKR5' or st_name == 'YKB2'):
                        doit = False
                    if (eq_num1 == 825 or eq_num2 == 825) and (st_name == 'YKB2'):
                        doit = False
                    if (eq_num1 == 826 or eq_num2 == 826) and (st_name == 'YKB2' or st_name == 'YKR1'):
                        doit = False
                    if (eq_num1 == 828 or eq_num2 == 828) and (st_name == 'YKR1' or st_name == 'YKB8'):
                        doit = False
                    if (eq_num1 == 829 or eq_num2 == 829) and (st_name == 'YKR7'):
                        doit = False
                    if (eq_num1 == 830 or eq_num2 == 830) and (st_name == 'YKR1' or st_name == 'YKR2' or st_name == 'YKR3' or st_name == 'YKR5'
                                                            or st_name == 'YKR7'or st_name == 'YKB2'):
                        doit = False
                    if (eq_num1 == 831 or eq_num2 == 831) and (st_name == 'YKB1' or st_name == 'YKR3' or st_name == 'YKR5'):
                        doit = False
                    if (eq_num1 == 837 or eq_num2 == 837) and (st_name == 'YKR9'):
                        doit = False
                    if (eq_num1 == 840 or eq_num2 == 840) and (st_name == 'IL16' or st_name == 'IL17'):
                        doit = False
                    if (eq_num1 == 843 or eq_num2 == 843) and (st_name == 'YKR2' or st_name == 'YKB1' or st_name == 'YKB3'):
                        doit = False
                    if (eq_num1 == 846 or eq_num2 == 846) and (st_name == 'YKB2' or st_name == 'YKR1'):
                        doit = False
                    if (eq_num1 == 847 or eq_num2 == 847) and (st_name == 'YKB8'):
                        doit = False
                    if (eq_num1 == 848 or eq_num2 == 848) and (st_name == 'YKB1'):
                        doit = False
                    if ((eq_num1 >= 852 and eq_num1 <= 883) or (eq_num2 >= 852 and eq_num2 <= 883)) and (st_name == 'YKR4' or st_name == 'IL14'):
                        doit = False
                    if (eq_num1 == 859 or eq_num2 == 859) and (st_name == 'IL08' or st_name == 'IL10' or st_name == 'YKR1'):
                        doit = False
                    if (eq_num1 == 860 or eq_num2 == 860) and (st_name == 'IL08' or st_name == 'YKR1'):
                        doit = False
                    if (eq_num1 == 861 or eq_num2 == 861) and (st_name == 'IL08'):
                        doit = False
                    if (eq_num1 == 862 or eq_num2 == 862) and (st_name == 'IL08' or st_name == 'IL15' or st_name == 'IL16' or st_name == 'IL17'):
                        doit = False
                    if (eq_num1 == 880 or eq_num2 == 880) and (st_name == 'YKR1'):
                        doit = False
                    if (eq_num1 == 884 or eq_num2 == 884) and (st_name == 'IL08'):
                        doit = False
                    if ((eq_num1 >= 884 and eq_num1 <= 888) or (eq_num2 >= 884 and eq_num2 <= 888)) and (st_name == 'IL14'):
                        doit = False
                    if (eq_num1 == 885 or eq_num2 == 885) and (st_name == 'IL16'):
                        doit = False
                    if (eq_num1 == 892 or eq_num2 == 892) and (st_name == 'IL14'):
                        doit = False
                    if (eq_num1 == 894 or eq_num2 == 894) and (st_name == 'YKR2' or st_name == 'YKR1'):
                        doit = False
                    if (eq_num1 == 895 or eq_num2 == 895) and (st_name == 'YKR2' or st_name == 'YKR1' or st_name == 'IL14'):
                        doit = False
                    if (eq_num1 == 896 or eq_num2 == 896) and (st_name == 'IL14'):
                        doit = False
                    if (eq_num1 == 897 or eq_num2 == 897) and (st_name == 'YKR1' or st_name == 'IL14'):
                        doit = False
                    if (eq_num1 == 898 or eq_num2 == 898) and (st_name == 'YKB1' or st_name == 'YKB4' or st_name == 'YKR1' or st_name == 'YKR2' or st_name == 'IL14'):
                        doit = False
                    if (eq_num1 == 899 or eq_num2 == 899) and (st_name == 'YKR1' or st_name == 'YKR2' or st_name == 'IL07' or st_name == 'IL14' or st_name == 'IL08' or st_name == 'IL15'):
                        doit = False
                    if (eq_num1 == 900 or eq_num2 == 900) and (st_name == 'YKR1' or st_name == 'YKR2' or st_name == 'YKB2'):
                        doit = False
                    if (eq_num1 == 901 or eq_num2 == 901) and (st_name == 'IL14'):
                        doit = False
                    if (eq_num1 == 902 or eq_num2 == 902) and (st_name == 'IL14'):
                        doit = False
                    if (eq_num1 == 903 or eq_num2 == 903) and (st_name == 'IL01' or st_name == 'IL13' or st_name == 'IL14' or st_name == 'IL15' or st_name == 'IL16' or st_name == 'IL17'):
                        doit = False
                    if (eq_num1 == 904 or eq_num2 == 904) and (st_name == 'IL08' or st_name == 'IL14'):
                        doit = False
                    if (eq_num1 == 905 or eq_num2 == 905) and (st_name == 'IL08' or st_name == 'IL14'):
                        doit = False
                    if (eq_num1 == 906 or eq_num2 == 906) and (st_name == 'YKR1' or st_name == 'IL08' or st_name == 'IL14'):
                        doit = False
                    if (eq_num1 == 907 or eq_num2 == 907) and (st_name == 'IL08' or st_name == 'IL14'):
                        doit = False
                    if (eq_num1 == 908 or eq_num2 == 908) and (st_name == 'IL08' or st_name == 'IL14'):
                        doit = False
                    if (eq_num1 == 909 or eq_num2 == 909) and (st_name == 'IL08' or st_name == 'IL14'):
                        doit = False
                    if doit == True:
                        st1good += tr1
                        st2good += tr2
                else:
                    # estimate median noise
                    if len(tr1.data) != len(tr2.data):
                        print('length of two traces to be compared is not the same!')

                    buff_time = end_buff - start_buff

                    time_to_signal_start = (precursor_shift - start_buff)
                    if time_to_signal_start - taper_frac * buff_time < noise_win_max: # noise window < max length
                        t_noise_start  = int(tr1.stats.sampling_rate * taper_frac * buff_time) # start just after taper
                        t_noise_end    = int(tr1.stats.sampling_rate * time_to_signal_start)   # end at beam start
                        if t_noise_end - t_noise_start < 20:
                            print('Very short to no window (less than 20 points) to estimate noise')
                    else:  # plenty of leader, set noise window to max length
                        time_to_noise_start = time_to_signal_start - noise_win_max
                        t_noise_start  = int(tr1.stats.sampling_rate * time_to_noise_start)  # noise window = noise_win_max
                        t_noise_end    = int(tr1.stats.sampling_rate * time_to_signal_start) # end at beam start

                    # estimate median noise and signal
                    noise1          = np.median(abs(tr1.data[t_noise_start:t_noise_end]))
                    noise2          = np.median(abs(tr2.data[t_noise_start:t_noise_end]))

                    t_signal_start = int(tr.stats.sampling_rate * time_to_signal_start)
                    t_signal_end   = int(t_signal_start + tr.stats.sampling_rate * signal_dur)
                    signal1         = np.median(abs(tr1.data[t_signal_start:t_signal_end]))
                    signal2         = np.median(abs(tr2.data[t_signal_start:t_signal_end]))
        #            test SNR
                    SNR1 = signal1/noise1;
                    SNR2 = signal2/noise2;
                    if (SNR1 > SNR_thres and SNR2 > SNR_thres):
                        st1good += tr1
                        st2good += tr2
    if apply_SNR == False:
        print('Matches (no SNR test): ' + str(len(st1good)) + ' traces')
    else:
        print('Match and above SNR threshold: ' + str(len(st1good)) + ' traces')

    #%%  get station lat-lon, compute distance for plot, find plot time and distance limits
    min_dist_auto = 180
    max_dist_auto = 0
    min_time_plot =  1000000
    max_time_plot = -1000000

    for tr in st1good:
        if tr.stats.station in st_names:  # find station in station list
            ii = st_names.index(tr.stats.station)
            stalon = float(st_lons[ii]) # look up lat & lon again to find distance
            stalat = float(st_lats[ii])
            distance = gps2dist_azimuth(stalat,stalon,ev_lat1,ev_lon1)
            tr.stats.distance=distance[0]/(1000*111) # distance in km
            if tr.stats.distance < min_dist_auto:
                min_dist_auto = tr.stats.distance
            if tr.stats.distance > max_dist_auto:
                max_dist_auto = tr.stats.distance
            if tr.stats.starttime - t1 < min_time_plot:
                min_time_plot = tr.stats.starttime - t1
            if ((tr.stats.starttime - t1) + ((len(tr.data)-1) * tr.stats.delta)) > max_time_plot:
                max_time_plot =  ((tr.stats.starttime - t1) + ((len(tr.data)-1) * tr.stats.delta))
    for tr in st2good:
        if tr.stats.station in st_names:  # find station in station list
            ii = st_names.index(tr.stats.station)
            stalon = float(st_lons[ii]) # look up lat & lon again to find distance
            stalat = float(st_lats[ii])
            distance = gps2dist_azimuth(stalat,stalon,ev_lat2,ev_lon2)
            tr.stats.distance=distance[0]/(1000*111) # distance in km

    print(f'        Min distance is   {min_dist_auto:.3f}   Max distance is {max_dist_auto:.3f}')
    print(f'        Min time is   {min_time_plot:.2f}   Max time is {max_time_plot:.2f}')
    if min_time_plot > start_buff:
        print(f'Min time {min_time_plot:.2f} > start_buff {start_buff:.2f}')
        print(colored('Warning: write zero-filling into pro3_pair?','red'))
        # sys.exit(-1)
    if max_time_plot < end_buff:
        print(f'Max time {max_time_plot:.2f} < end_buff {end_buff:.2f}')
        print(colored('Warning: write zero-filling into pro3_pair?','red'))
        # sys.exit(-1)

    #%% plot traces evenly spaced
    # fig_index = 3
    plt.figure(figsize=(10,5), num = fig_index)

    plt.xlim(min_time_plot,max_time_plot)

    plt.ylim(-1,len(st1good) + 1)

    tr_cnt = 0
    for tr in st1good:
        ttt = np.arange(len(tr.data)) * tr.stats.delta + (tr.stats.starttime - t1)
        if win_norm == True:
            s_t = t1 + start_buff + wind_buff
            e_t = t1 + end_buff - wind_buff
            snippet = tr.copy()
            snippet.trim( starttime = s_t, endtime = e_t)
            if blowup_PKIKP:
                smax = abs(snippet.max())/20. # blow up PKIKP 20X
            else:
                smax = abs(snippet.max())
        else:
            smax = tr.data.max() - tr.data.min()
        plt.plot(ttt, (tr.data - np.median(tr.data))*trace_amp/smax + tr_cnt, color = 'green')
        # plt.plot(ttt, (tr.data - np.median(tr.data))*trace_amp/(tr.data.max() - tr.data.min()) + tr_cnt, color = 'green')
        tr_cnt = tr_cnt + 1

    tr_cnt = 0
    for tr in st2good:
        ttt = np.arange(len(tr.data)) * tr.stats.delta + (tr.stats.starttime - t2)
        if win_norm == True:
            s_t = t2 + start_buff + wind_buff
            e_t = t2 + end_buff - wind_buff
            snippet = tr.copy()
            snippet.trim( starttime = s_t, endtime = e_t)
            if blowup_PKIKP:
                smax = abs(snippet.max())/20. # blow up PKIKP 20X
            else:
                smax = abs(snippet.max())
        else:
            smax = tr.data.max() - tr.data.min()
        plt.plot(ttt, (tr.data - np.median(tr.data))*trace_amp/smax + tr_cnt, color = 'red')
        # plt.plot(ttt, (tr.data - np.median(tr.data))*trace_amp/(tr.data.max() - tr.data.min()) + tr_cnt, color = 'red')
        if plot_sta_names:
            plt.text(min_time_plot + 0.1 * (max_time_plot - min_time_plot), tr_cnt + 0.03, tr.stats.station, color = 'black')
        tr_cnt = tr_cnt + 1

    if ARRAY == 1:
        plt.title(phase1 + ' for ' + fname1[43:53] + ' vs ' + fname2[43:53] + ' freqs ' + str(freq_min) + '-' + str(freq_max) + ' Hz')
    elif ARRAY == 4 :
        plt.title(phase1 + ' for ' + fname1[46:50] + '-' + fname1[50:52] + '-' + fname1[52:54] + ' vs ' + fname2[46:50] + '-' + fname2[50:52] + '-' + fname2[52:54] + ' for array WRA, ' + str(eq_num1) + ' and ' + str(eq_num2) + ' freqs ' + str(freq_min) + '-' + str(freq_max) + ' Hz')
    elif ARRAY == 5 or ARRAY== 99:
        plt.title(repeater + ' ' + phase1 + ' for ' + fname1[46:50] + '-' + fname1[50:52] + '-' + fname1[52:54] + ' vs ' + fname2[46:50] + '-' + fname2[50:52] + '-' + fname2[52:54] + ' for array YKA, ' + str(eq_num1) + '(g) and ' + str(eq_num2) + '(r), freqs ' + str(freq_min) + '-' + str(freq_max) + ' Hz')
    elif ARRAY == 6:
        plt.title(repeater + ' ' + phase1 + ' for ' + fname1[47:51] + '-' + fname1[51:53] + '-' + fname1[53:55] + ' vs ' + fname2[47:51] + '-' + fname2[51:53] + '-' + fname2[53:55] + ' for array ILAR, ' + str(eq_num1) + '(g) and ' + str(eq_num2) + '(r), freqs ' + str(freq_min) + '-' + str(freq_max) + ' Hz')
    elif ARRAY == 7:
        plt.title(repeater + ' ' + phase1 + ' for ' + fname1[49:53] + '-' + fname1[53:55] + '-' + fname1[55:57] + ' vs ' + fname2[49:53] + '-' + fname2[53:55] + '-' + fname2[55:57] + ' for global array, ' + str(eq_num1) + '(g) and ' + str(eq_num2) + '(r), freqs ' + str(freq_min) + '-' + str(freq_max) + ' Hz', y = 1)

    os.chdir('/Users/vidale/Documents/Research/IC/Plots_hold')
    plt.savefig(repeater + '_Array_' + str(ARRAY) + '_traces')

    #%% plot traces at true distance
    if dist_plot:
        # fig_index = 3
        plt.figure(figsize=(10,5), num = fig_index + 1)

        plt.xlim(min_time_plot,max_time_plot)

        if auto_dist == True:
            if(max_dist_auto == min_dist_auto):
                min_dist_auto -= 1
                max_dist_auto += 1
            dist_diff = max_dist_auto - min_dist_auto # add space at extremes
            plt.ylim(min_dist_auto - 0.1 * dist_diff, max_dist_auto + 0.2 * dist_diff)
        else:
            plt.ylim(min_dist,max_dist)

        for tr in st1good:
            dist_offset = tr.stats.distance # trying for approx degrees
            ttt = np.arange(len(tr.data)) * tr.stats.delta + (tr.stats.starttime - t1)
            plt.plot(ttt, (tr.data - np.median(tr.data))*plot_scale_fac /(tr.data.max()
                - tr.data.min()) + dist_offset, color = 'green')

        for tr in st2good:
            dist_offset = tr.stats.distance # trying for approx degrees
            ttt = np.arange(len(tr.data)) * tr.stats.delta + (tr.stats.starttime - t2)
            plt.plot(ttt, (tr.data - np.median(tr.data))*plot_scale_fac /(tr.data.max()
                - tr.data.min()) + dist_offset, color = 'red')
            if plot_sta_names:
                plt.text(min_time_plot + 0.015 * (max_time_plot - min_time_plot), dist_offset + 0.003 * (max_dist_auto - min_dist_auto), tr.stats.station, color = 'black')

    #%% -- Traveltime curves
        if plot_tt:
            # first traveltime curve
            line_pts = 50
            dist_vec  = np.arange(min_dist_auto, max_dist_auto, (max_dist_auto - min_dist_auto)/line_pts) # distance grid
            time_vec1 = np.arange(min_dist_auto, max_dist_auto, (max_dist_auto - min_dist_auto)/line_pts) # empty time grid of same length (filled with -1000)

            for i in range(0,line_pts):
                arrivals = model.get_travel_times(source_depth_in_km=ev_depth1,distance_in_degree=dist_vec[i],phase_list=[phase1])
                if(len(arrivals) == 0 and dist_vec[i] < 10 and phase1 == 'P'):  # in case first arrival is upgoing P, which is 'p'
                    arrivals = model.get_travel_times(source_depth_in_km=ev_depth1,distance_in_degree=dist_vec[i],phase_list='p')
                    if(len(arrivals) == 0):
                        print('model.get_travel_times failed: dist, phase  ' + str(dist_vec[i]) + '   ' + phase1)
                num_arrivals = len(arrivals)
                found_it = 0
                for j in range(0,num_arrivals):
                    if arrivals[j].name == phase1:
                        time_vec1[i] = arrivals[j].time
                        found_it = 1
                if found_it == 0:
                    time_vec1[i] = np.nan
            # second traveltime curve
            if phase2 != 'no':
                time_vec2 = np.arange(min_dist_auto, max_dist_auto, (max_dist_auto - min_dist_auto)/line_pts) # empty time grid of same length (filled with -1000)
                for i in range(0,line_pts):
                    arrivals = model.get_travel_times(source_depth_in_km=ev_depth1,distance_in_degree=dist_vec[i],phase_list=[phase2])
                    if(len(arrivals) == 0 and dist_vec[i] < 10 and phase2 == 'P'):  # in case first arrival is upgoing P, which is 'p'
                        arrivals = model.get_travel_times(source_depth_in_km=ev_depth1,distance_in_degree=dist_vec[i],phase_list='p')
                        if(len(arrivals) == 0):
                            print('model.get_travel_times failed: dist, phase  ' + str(dist_vec[i]) + '   ' + phase2)
                    num_arrivals = len(arrivals)
                    found_it = 0
                    for j in range(0,num_arrivals):
                        if arrivals[j].name == phase2:
                            time_vec2[i] = arrivals[j].time
                            found_it = 1
                    if found_it == 0:
                        time_vec2[i] = np.nan
                if   rel_time == 1:
                    time_vec2 = time_vec2 - time_vec1 + (dist_vec - ref1_dist) * ref_slow
                if   rel_time == 3:
                    time_vec2 = time_vec2 - time_vec1
                elif rel_time == 2 or rel_time == 4:
                    time_vec2 = time_vec2 - atime_ref
                plt.plot(time_vec2,dist_vec, color = 'orange', label = phase2)
            # third traveltime curve
            if phase3 != 'no':
                time_vec3 = np.arange(min_dist_auto, max_dist_auto, (max_dist_auto - min_dist_auto)/line_pts) # empty time grid of same length (filled with -1000)
                for i in range(0,line_pts):
                    arrivals = model.get_travel_times(source_depth_in_km=ev_depth1,distance_in_degree=dist_vec[i],phase_list=[phase3])
                    if(len(arrivals) == 0 and dist_vec[i] < 10 and phase3 == 'P'):  # in case first arrival is upgoing P, which is 'p'
                        arrivals = model.get_travel_times(source_depth_in_km=ev_depth1,distance_in_degree=dist_vec[i],phase_list='p')
                        if(len(arrivals) == 0):
                            print('model.get_travel_times failed: dist, phase  ' + str(dist_vec[i]) + '   ' + phase3)
                    num_arrivals = len(arrivals)
                    found_it = 0
                    for j in range(0,num_arrivals):
                        if arrivals[j].name == phase3:
                            time_vec3[i] = arrivals[j].time
                            found_it = 1
                    if found_it == 0:
                        time_vec3[i] = np.nan
                if   rel_time == 1:
                    time_vec3 = time_vec3 - time_vec1 + (dist_vec - ref1_dist) * ref_slow
                if   rel_time == 3:
                    time_vec3 = time_vec3 - time_vec1
                elif rel_time == 2 or rel_time == 4:
                    time_vec3 = time_vec3 - atime_ref
                plt.plot(time_vec3,dist_vec, color = 'yellow', label = phase3)
            # fourth traveltime curve
            if phase4 != 'no':
                time_vec4 = np.arange(min_dist_auto, max_dist_auto, (max_dist_auto - min_dist_auto)/line_pts) # empty time grid of same length (filled with -1000)
                for i in range(0,line_pts):
                    arrivals = model.get_travel_times(source_depth_in_km=ev_depth1,distance_in_degree=dist_vec[i],phase_list=[phase4])
                    if(len(arrivals) == 0 and dist_vec[i] < 10 and phase4 == 'P'):  # in case first arrival is upgoing P, which is 'p'
                        arrivals = model.get_travel_times(source_depth_in_km=ev_depth1,distance_in_degree=dist_vec[i],phase_list='p')
                        if(len(arrivals) == 0):
                            print('model.get_travel_times failed: dist, phase  ' + str(dist_vec[i]) + '   ' + phase4)
                    num_arrivals = len(arrivals)
                    found_it = 0
                    for j in range(0,num_arrivals):
                        if arrivals[j].name == phase4:
                            time_vec4[i] = arrivals[j].time
                            found_it = 1
                    if found_it == 0:
                        time_vec4[i] = np.nan
                if   rel_time == 1:
                    time_vec4 = time_vec4 - time_vec1 + (dist_vec - ref1_dist) * ref_slow
                if   rel_time == 3:
                    time_vec4 = time_vec4 - time_vec1
                elif rel_time == 2 or rel_time == 4:
                    time_vec4 = time_vec4 - atime_ref
                plt.plot(time_vec4,dist_vec, color = 'purple', label = phase4)
            if   rel_time == 1:
                time_vec1 = (dist_vec - ref1_dist) * ref_slow
            if   rel_time == 3:
                time_vec1 = time_vec1 - time_vec1
            elif rel_time == 2 or rel_time == 4:
                time_vec1 = time_vec1 - atime_ref
            plt.plot(time_vec1,dist_vec, color = 'blue', label = phase1)

        if zoom: # draw gray lines to show zoom area
            plt.plot((Zstart_buff, Zstart_buff), (min_dist,max_dist), color = 'gray', label = 'beam start')
            plt.plot((Zend_buff, Zend_buff), (min_dist,max_dist), color = 'lightgray', label = 'beam end')

        plt.xlabel('Time (s)')
        plt.ylabel('Epicentral distance from event (°)')
        plt.legend(loc="upper left")
        if ARRAY == 1:
            plt.title(phase1 + ' for ' + fname1[43:53] + ' vs ' + fname2[43:53] + ' freqs ' + str(freq_min) + '-' + str(freq_max) + ' Hz')
        elif ARRAY == 4 :
            plt.title(phase1 + ' for ' + fname1[46:50] + '-' + fname1[50:52] + '-' + fname1[52:54] + ' vs ' + fname2[46:50] + '-' + fname2[50:52] + '-' + fname2[52:54] + ' for array WRA, ' + str(eq_num1) + '(g) and ' + str(eq_num2) + '(r) freqs ' + str(freq_min) + '-' + str(freq_max) + ' Hz')
        elif ARRAY == 5 or ARRAY== 99:
            plt.title(repeater + ' ' + phase1 + ' for ' + fname1[46:50] + '-' + fname1[50:52] + '-' + fname1[52:54] + ' vs ' + fname2[46:50] + '-' + fname2[50:52] + '-' + fname2[52:54] + ' for array YKA, ' + str(eq_num1) + '(g) and ' + str(eq_num2) + '(r) freqs ' + str(freq_min) + '-' + str(freq_max) + ' Hz')
        elif ARRAY == 6:
            plt.title(repeater + ' ' + phase1 + ' for ' + fname1[47:51] + '-' + fname1[51:53] + '-' + fname1[53:55] + ' vs ' + fname2[47:51] + '-' + fname2[51:53] + '-' + fname2[53:55] + ' for array ILAR, ' + str(eq_num1) + '(g) and ' + str(eq_num2) + '(r) freqs ' + str(freq_min) + '-' + str(freq_max) + ' Hz')
        elif ARRAY == 7:
            plt.title(repeater + ' ' + phase1 + ' for ' + fname1[49:53] + '-' + fname1[53:55] + '-' + fname1[55:57] + ' vs ' + fname2[49:53] + '-' + fname2[53:55] + '-' + fname2[55:57] + ' for global array, ' + str(eq_num1) + '(g) and ' + str(eq_num2) + '(r), freqs ' + str(freq_min) + '-' + str(freq_max) + ' Hz', y = 1)
        else:
            plt.title(phase1 + ' for ' + fname1[42:52] + ' vs ' + fname2[42:52] + ' freqs ' + str(freq_min) + '-' + str(freq_max) + ' Hz')
        ax = plt.gca()
        ax.get_yaxis().get_major_formatter().set_useOffset(False)
        # plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)  # turn off conversion of y-axis to offset numbering
        # plt.rcParams['axes.formatter.useoffset'] = False #  these two commented commands worked for some but not all cases
        os.chdir('/Users/vidale/Documents/Research/IC/Plots_hold')
        plt.savefig(repeater + '_Array_' + str(ARRAY) + '_section')

#%%  Save processed files
    cnt1 = len(st1good)
    cnt2 = len(st2good)
    print('Matches:  st1good' + str(cnt1) + ' st2good ' + str(cnt2))
    if cnt1 == cnt2 and cnt1 != 0:
        fname1 = '/Users/vidale/Documents/Research/IC/Pro_Files/HD' + date_label1 + 'sel.mseed'
        fname2 = '/Users/vidale/Documents/Research/IC/Pro_Files/HD' + date_label2 + 'sel.mseed'
        st1good.write(fname1,format = 'MSEED')
        st2good.write(fname2,format = 'MSEED')
    else:
        print('No pairs found\n\n')

    elapsed_time_wc = time.time() - start_time_wc
    print(f'This job took   {elapsed_time_wc:.1f}   seconds')
    os.system('say "three done"')
