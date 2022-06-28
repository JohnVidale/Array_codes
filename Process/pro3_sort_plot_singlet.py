#!/usr/bin/env python
# reads in raw mseed traces from a single event
# this program tapers, filters, selects range and SNR
# plots against traveltime curves, either raw or reduced against traveltimes
# writes out "*sel.mseed" file
# plot lines are blue, orange, yellow, purple for phases 1 through 4
# John Vidale 2/2019

def pro3singlet(eq_num, stat_corr = 1, corr_threshold = 0, rel_time = 1, shift_tt = 0,
            max_taper_length = 3., simple_taper = True, apply_SNR = False, SNR_thres = 0,
            dphase = 'P', dphase2 = '', dphase3 = '', dphase4 = '',
            start_buff      =   -10, end_buff   =    10,
            start_beam      =     0, end_beam   =     0,
            Zstart_buff     =   -10, Zend_buff  =    10, zoom = 0,
            precursor_shift = -1000, signal_dur = -1000,
            freq_min = 0.25, freq_max = 1, do_filt = 1, zerophase = True, plot_scale_fac = 1,
            min_dist = 0, max_dist = 180, plot_auto_dist = True, do_decimate = False, decimate_fac = 1,
            alt_statics = 0, statics_file = 'nothing', ARRAY = 0,
            ref_rad = 180, ref_loc = False, ref_lat = 0, ref_lon = 0,
            verbose = False, fig_index = 1, JST = False):
# 0 is Hinet, 1 is LASA, 2 is NORSAR

#%% Import
#%% -- Functions
    from obspy import UTCDateTime
    from obspy import Stream
    from obspy import read
    from obspy.geodetics import gps2dist_azimuth
    import numpy as np
    import os
    import sys
    from obspy.taup import TauPyModel
    import matplotlib.pyplot as plt
    import time
    from termcolor import colored
    model = TauPyModel(model='iasp91')

    if precursor_shift == -1000: # if precursor_ shift is not specified, use beam time for SNR decisions
        precursor_shift = start_beam
        signal_dur = end_beam - start_beam

#    import sys # don't show any warnings
#    import warnings
#
#    if not sys.warnoptions:
#        warnings.simplefilter("ignore")

    print(colored('Running pro3_sort_plot_singlet', 'cyan'))
    start_time_wc = time.time()
    plot_sta_names = True

#%% -- Event info
    #  input event data with 1-line file of format
    #  event 2016-05-28T09:47:00.000 -56.241 -26.935 78
    fname = '/Users/vidale/Documents/Research/IC/EvLocs/event' + str(eq_num) + '.txt'
    print('Opening event file ' + fname)
    file = open(fname, 'r')
    lines=file.readlines()
    split_line = lines[0].split()
#            ids.append(split_line[0])  ignore label for now
    t           = UTCDateTime(split_line[1])
    date_label  = split_line[1][0:10]
    year_label  = split_line[1][0:4]
    year_short_label  = split_line[1][2:4]
    month_label   = split_line[1][5:7]
    day_label     = split_line[1][8:10]
    hour_label    = split_line[1][11:13]
    minute_label  = split_line[1][14:16]
    ev_lat      = float(      split_line[2])
    ev_lon      = float(      split_line[3])
    ev_depth    = float(      split_line[4])
    print('    event parameters: ' + str(t) + ' lat ' + str(ev_lat) + ' lon ' + str( ev_lon) + ' depth ' + str(ev_depth))

#%% -- Station locations and statics
    if stat_corr != 0:  # load static terms, only applies to Hinet and LASA
        if ARRAY == 0:
            if   stat_corr == 1: # standard set
                sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/sta_statics_hinet.txt'
            elif stat_corr == 2: # custom set made for SSI to Hi-Net PKiKP
                sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/SSI_HiNet_statics.txt'
            elif stat_corr == 3: # custom set made for Kawa PKiKP for Hi-Net
                sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/Kawa_HiNet_statics.txt'
        elif ARRAY == 1:
            sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/sta_statics_LASA.txt'
        elif ARRAY == 2:
            sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/sta_statics_ch.txt'
        elif ARRAY == 4:
            sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/sta_statics_WR.txt'
        elif ARRAY == 5:
            sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/sta_statics_YK.txt'
        with open(sta_file, 'r') as file:
            lines = file.readlines()
        print('    ' + str(len(lines)) + ' station statics read from ' + sta_file)
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
            if ARRAY == 0 or ARRAY == 1:
                st_dist.append( split_line[1])
                st_lats.append( split_line[2])
                st_lons.append( split_line[3])
                st_shift.append(split_line[4])
                st_corr.append( split_line[5])
            elif ARRAY == 2:
                st_lats.append( split_line[1])
                st_lons.append( split_line[2])
                st_shift.append(split_line[3])
                st_corr.append( split_line[4]) # but really std dev
            elif ARRAY == 4:
                st_lats.append( split_line[2])
                st_lons.append( split_line[3])
                st_shift.append(split_line[4])
                st_corr.append( split_line[5]) # but really std dev
            else:
                exit('Failure: Asking for non-existent statics.')

    else: # no static terms, always true for NORSAR
        if   ARRAY == 0: # Hinet set
            sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/sta_hinet.txt'
        elif ARRAY == 1: #         LASA set
            sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/sta_LASA.txt'
        elif ARRAY == 2: #         China set
            sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/sta_ch.txt'
        elif ARRAY == 3: #         NORSAR set
            sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/sta_NORSAR.txt'
        elif ARRAY == 4: #         Warramunga set
            sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/sta_AU_WR.txt'
        elif ARRAY == 5: #         Yellowknife set
            sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/sta_CN_YK.txt'
        with open(sta_file, 'r') as file:
            lines = file.readlines()
        print('    ' + str(len(lines)) + ' stations read from ' + sta_file)
        # Load station coords into arrays
        station_index = range(len(lines))
        st_names = []
        st_lats  = []
        st_lons  = []
        for ii in station_index:
            line = lines[ii]
            split_line = line.split()
            # print('input station ' + split_line[0])
            st_names.append(split_line[0])
            st_lats.append( split_line[1])
            st_lons.append( split_line[2])
    if ARRAY == 0:  # shorten and make upper case Hi-net station names to match station list
        for ii in station_index:
            this_name = st_names[ii]
            this_name_truc = this_name[0:5]
            st_names[ii]  = this_name_truc.upper()

#%%  -- Rest of parameters
#    stat_corr = 1, 2, or 3 # apply station static corrections
#    rel_time = 1           # timing is relative to a chosen phase, otherwise relative to OT
#    dphase  = 'PKIKP'      # phase to be aligned
#    dphase2 = 'PKiKP'      # another phase to have traveltime plotted
#    dphase3 = 'PKP'        # another phase to have traveltime plotted
#    dphase4 = 'pP'        # another phase to have traveltime plotted
    taper_frac = .05      # Fraction of window tapered on both ends
    noise_win_max = 20    # maximum length of noise window for SNR estimation, seconds
#    plot_scale_fac = 0.5    #  Bigger numbers make each trace amplitude bigger on plot
#    SNR_thres      = 2 # minimum SNR
#    corr_threshold = 0.7  # minimum correlation in measuring shift to use station
    plot_tt = 1           # plot the traveltimes?
    plot_auto_dist = 1    # plot just the traces, not potentially bigger dist_min to dist_max
    # if ref_loc == True,  use distance to (ref_lat,ref_lon) to filter stations
    # if ref_loc == False, use min_dist and max_dist to earthquake loc to filter stations
    #    ref_rad = 0.4    # ° radius (°) set by input or at top
    if ref_loc == False:  # use default array centers unless specified
        if ARRAY == 0:
            ref_lat = 36  # °N, around middle of Japan
            ref_lon = 139 # °E
        elif ARRAY == 1:
            ref_lat = 46.7      # °N keep only inner rings A-D if radius is 0.4°
            ref_lon = -106.22   # °E
        elif ARRAY == 2:
            ref_lat = 38      # °N
            ref_lon = 104.5   # °E
        elif ARRAY == 4:
            ref_lat = -19.89  # °N Warramunga
            ref_lon = 134.42  # °E
        elif ARRAY == 5:
            ref_lat =  62.49  # °N Yellowknife
            ref_lon = -114.6  # °E

#%% -- Test taper, needs adjustment?
#   Is taper too long compared to noise estimation window?
    totalt = end_buff - start_buff
    noise_time_skipped = taper_frac * totalt
    noise_time_skipped = min(noise_time_skipped,10.0) # set max of 10s to taper length
    print(f'        noise time at start skipped due to taper length: {noise_time_skipped:.2f} seconds')
    if simple_taper == False:
        if noise_time_skipped >= -0.5 * start_buff:
            print('        ' + 'Specified taper of ' + str(taper_frac * totalt) +
               ' is not big enough compared to available noise estimation window ' +
               str(-start_buff - noise_time_skipped) + '. May not work well.')
            old_taper_frac = taper_frac
            taper_frac = -0.5*start_buff/totalt
            if start_buff > 0:
                    taper_frac = 0.05 # pick random minimal window if there is no leader
            print('        ' + 'Taper reset from ' + str(old_taper_frac * totalt) + ' to '
               + str(taper_frac * totalt) + ' seconds.')

    if rel_time == 0: # SNR requirement not implemented for unaligned traces
        SNR_thres = 0

#%% Load waveforms
    st = Stream()
    # fname     = '/Users/vidale/Documents/GitHub/LASA_data/HD' + date_label + '.mseed'
    if ARRAY == 1:
        mseed_name = year_short_label + month_label + day_label + '_' + hour_label + minute_label
        fname     = '/Users/vidale/Documents/Research/IC/Mseed/L' + mseed_name + '.mseed'
    elif ARRAY == 4 or ARRAY == 5:
        mseed_name = year_label + month_label + day_label + '_' + hour_label + minute_label
        fname     = '/Users/vidale/Documents/Research/IC/Mseed/' + mseed_name + '.mseed'
    else:
        mseed_name = year_short_label + '-' +month_label + '-' + day_label
        fname     = '/Users/vidale/Documents/Research/IC/Mseed/HD20' + mseed_name + '.mseed'

    print('Opening data file: ' + fname)
    st=read(fname)
    print('data len ' + str(len(st[0].data)) + ' # traces ' + str(len(st)) + ' delta ' + str(st[0].stats.delta))
#%% -- Ensure 10 sps
    if do_decimate == True:
        st.decimate(decimate_fac, no_filter=True)

    print('    ' + str(len(st)) + '  traces read in')
    if(len(st) == 0):
        exit('Failure: No traces read.')
    print('        1st trace parameters')
    print('        Start time: ' + str(st[0].stats.starttime) + '  event time : ' + str(t))
    if JST: # if necessary, convert JST -> UTC, time in Greenwich 9 hours earlier than Japan
        print('        Start lag: ' + str((st[0].stats.starttime - 9*60*60) - t))
    else:
        print('        Start lag: ' + str(st[0].stats.starttime - t))
    # print('    ' + str(len(st)) + '  traces after decimation')
    nt = len(st[0].data)
    dt = st[0].stats.delta
    print('        Time pts: ' + str(nt) + ' time sampling: ' + str(dt) + ' duration: ' + str((nt-1)*dt))

#%% Select and process traces
    st_pickalign = Stream()

    tra_in_range  = 0
    tra_sta_found = 0
    nodata        = 0
    no_tt_calc    = 0
    min_dist_auto = 180
    max_dist_auto = 0
    min_time_plot =  1000000
    max_time_plot = -1000000

    # find reference distance and slowness at array center
    ref_distance = gps2dist_azimuth(ref_lat,ref_lon,ev_lat,ev_lon)
    ref1_dist  = ref_distance[0]/(1000*111)
    arrivals_ref   = model.get_travel_times(source_depth_in_km=ev_depth,distance_in_degree=ref1_dist, phase_list=[dphase])
    if(len(arrivals_ref) == 0 and ref1_dist < 10 and dphase == 'P'):  # in case first arrival might be upgoing P, which is 'p'
        arrivals_ref   = model.get_travel_times(source_depth_in_km=ev_depth,distance_in_degree=ref1_dist, phase_list='p')
    if(len(arrivals_ref) == 0):
        print('model.get_travel_times failed: dist, phase  ' + str(ref1_dist) + '   ' + dphase)
    if dphase == 'PKP':  # use last element; first element might be BC, with limited range and thus step in [0] element
        atime_rayp = arrivals_ref[-1].ray_param  # ray parameter (s/radian) at reference distance
        atime_ref  = arrivals_ref[-1].time       # phase arrival time at reference distance
    else:  # use first element
        atime_rayp = arrivals_ref[0].ray_param
        atime_ref  = arrivals_ref[0].time
    ref_slow  = atime_rayp * 2 * np.pi / 360. # convert to s/°
    # print(f'ref_slow {ref_slow:.2f}')

    for tr in st: # traces one by one, find lat-lon
        print('station ' + tr.stats.station + ' trace start time relative to event time ' + str(tr.stats.starttime  - t) + ' duration ' + str(len(tr.data)*tr.stats.delta))
        if float(year_label) < 1970: # fix the damn 1969 -> 2069 bug in Gibbon's LASA data
            temp_t = str(tr.stats.starttime)
            temp_tt = '19' + temp_t[2:]
            tr.stats.starttime = UTCDateTime(temp_tt)
        if JST: # if necessary, convert JST -> UTC, time in Greenwich 9 hours earlier than Japan
            tr.stats.starttime = tr.stats.starttime - 9*60*60
        if tr.stats.station in st_names:  # find station in station list
            ii = st_names.index(tr.stats.station)
            tra_sta_found += 1

            corr = 1
            if stat_corr != 0:
                corr = float(st_corr[ii])

#%% -- Reject low correlation
            if corr > corr_threshold: # if using statics, reject low correlations
                stalat = float(st_lats[ii]) # look up lat & lon again to find distance
                stalon = float(st_lons[ii])

#%% -- Reject out of range
                distance = gps2dist_azimuth(stalat,stalon,ev_lat,ev_lon) # Get traveltimes again, hard to store
                tr.stats.distance=distance[0] # distance in km
                dist = distance[0]/(1000*111)

                in_range = 0  # flag for whether this trace goes into stack
                rejector = False  # flag in case model.get_travel_times fails

                if ref_loc == False:  # check whether trace is in distance range from earthquake
                    # if min_dist < dist and dist < max_dist and tr.stats.station[2] == 'B':
                    if min_dist < dist and dist < max_dist:
                        in_range = 1
                        tra_in_range += 1
                elif ref_loc == True:  # alternately, check whether trace is close enough to ref_location
                    ref_distance = gps2dist_azimuth(ref_lat,ref_lon,stalat,stalon)
                    ref_degree_dist = ref_distance[0]/(1000*111)
                    if ref_degree_dist < ref_rad:
                        in_range = 1
                        tra_in_range += 1

                if in_range == 1:   # trace fulfills the specified criteria for being in range
                    s_t = t  + start_buff
                    e_t = t  + end_buff

#%% -- Apply static
                    if stat_corr != 0: # apply static station corrections
                        tr.stats.starttime -= float(st_shift[ii])
                        # print('static correction is' + str(st_shift[ii]) + ' for station ' + tr.stats.station)
                    if rel_time == 0:  #  use absolute time, ignore time of chose phase
                        tr.trim(starttime=s_t,endtime = e_t)
                    else:              # shift relative to a chosen phase
                        arrivals_each = model.get_travel_times(source_depth_in_km=ev_depth, distance_in_degree=dist,phase_list=[dphase])
                        if(len(arrivals_each) == 0):
                            if(dist < 10 and dphase == 'P'):  # in case first arrival is upgoing P, try 'p'
                                arrivals_each = model.get_travel_times(source_depth_in_km=ev_depth, distance_in_degree=dist,phase_list='p')
                        if(len(arrivals_each) == 0):  # did it still fail?
#%% -- Can't calculate traveltime
                                # print(f'travel_time failed: station {tr.stats.station} dist {dist:.2f} phase {dphase}')
                                tra_in_range -= 1 # don't count this trace after all
                                rejector = True
                                no_tt_calc += 1
                        else:
                            # atime_each is phase arrival on this trace
                            # atime_ref  is phase arrival on reference trace

                            if dphase == 'PKP':
                                atime_each = arrivals_each[-1].time  # use last instance for AB branch
                            else:
                                atime_each = arrivals_each[0].time
                            # start time has a shift proportional to ref_dist at phase slowness at ref_dist
                            # preserves slowness of arrivals at ref_dist, sharp even if tt curve has curvature
                            # uniform relative window limits around phase
                            if rel_time == 1:
                                s_t += atime_each
                                e_t += atime_each
                                tr.trim( starttime=s_t, endtime = e_t)
                                tr.stats.starttime -= atime_each - (dist-ref1_dist) * ref_slow
                            # window time sets 0 at chosen phase only on reference trace
                            # keeps true slowness, but blurs arcuate slowness curves
                            # uniform relative window limits around phase
                            elif rel_time == 2:
                                s_t += atime_each
                                e_t += atime_each
                                tr.trim( starttime=s_t, endtime = e_t)
                                tr.stats.starttime -= atime_ref
                            # window time sets 0 at chosen phase on all traces
                            # set phase to 0 slowness
                            # uniform relative window limits around phase
                            elif rel_time == 3:
                                s_t += atime_each
                                e_t += atime_each
                                tr.trim( starttime=s_t, endtime = e_t)
                                tr.stats.starttime -= atime_each
                            # window time sets 0 at chosen phase only on reference trace
                            # set phase to 0 slowness
                            # use same absolute window around chosen phase for all stations
                            elif rel_time == 4:
                                s_t += atime_ref
                                e_t += atime_ref
                                tr.trim( starttime=s_t, endtime = e_t)
                                tr.stats.starttime -= atime_ref
                            else:
                                sys.exit('Invalid rel_time parameter, must be integer 0 to 4')
#%% -- Reject if not in time window
                    if len(tr.data) > 0 and rejector == False:
                        st_pickalign += tr
                    else:
                        nodata += 1
#%% -- Reject if not in station (static) list
        # else:
        #     print(tr.stats.station + ' not found in station list with statics')
    print('    ' + str(tra_in_range) + '  traces in range')
    print('        ' + str(nodata) + '  traces with no data,  ' + str(no_tt_calc) + '  traces tt calc failed')
    print('    ' + str(len(st_pickalign)) + '  traces after alignment, check for data, tt calc, correlation selection')

    #print(st) # at length
    if verbose == True:
        print(st.__str__(extended=True))
        if rel_time == 1:
            print(st_pickalign.__str__(extended=True))

#%%  -- Detrend, taper, filter, taper
    st_pickalign.detrend(type='simple')
    st_pickalign.taper(taper_frac, max_length = max_taper_length)
    if do_filt == 1:
        print('Filter' + ' min ' + str(freq_min) + ' max ' + str(freq_max))
        print('data len ' + str(len(st_pickalign[0].data)) + ' # traces ' + str(len(st_pickalign)) + ' delta ' + str(st_pickalign[0].stats.delta))
        st_pickalign.filter('bandpass', freqmin=freq_min, freqmax=freq_max, corners=4, zerophase=zerophase)
    st_pickalign.taper(taper_frac, max_length = max_taper_length)

#%%  -- Cull by SNR threshold
    if apply_SNR == False:
        stgood = st_pickalign.copy()
    else:
        stgood = Stream()
        buff_time = end_buff - start_buff
        time_to_signal_start = (precursor_shift - start_buff)
        print(f'SNR parameters: start buffer {start_buff:.1f} end buffer {end_buff:.1f} pre_shift {precursor_shift:.1f} sig duration {signal_dur:.1f} SNR onset offset {precursor_shift:.1f}')
        for tr in st_pickalign:
            # estimate median noise
            if time_to_signal_start - taper_frac * buff_time < noise_win_max: # noise window < max length
                t_noise_start  = int(tr.stats.sampling_rate * taper_frac * buff_time) # start just after taper
                t_noise_end    = int(tr.stats.sampling_rate * time_to_signal_start)   # end at beam start
                if t_noise_end - t_noise_start < 20:
                    print('Very short to no window (less than 20 points) to estimate noise')
            else:  # plenty of leader, set noise window to max length
                time_to_noise_start = time_to_signal_start - noise_win_max
                t_noise_start  = int(tr.stats.sampling_rate * time_to_noise_start)  # noise window = noise_win_max
                t_noise_end    = int(tr.stats.sampling_rate * time_to_signal_start) # end at beam start
            noise = np.median(abs(tr.data[t_noise_start:t_noise_end]))

            # estimate median signal
            t_signal_start = int(tr.stats.sampling_rate * time_to_signal_start)
            t_signal_end   = int(t_signal_start + tr.stats.sampling_rate * signal_dur)
            signal = np.median(abs(tr.data[t_signal_start:t_signal_end]))

            # test SNR
            SNR = signal/noise;
            if (SNR > SNR_thres):
                stgood += tr
            # else:
            #     print('Failed SNR test: SNR ' + str(SNR))

    print('    ' + str(len(stgood)) + '  traces above SNR threshold')
    if verbose == True:
        for tr in stgood:
            print('        Distance is ' + str(tr.stats.distance/(1000*111)) + ' for station ' + tr.stats.station)

#%% Plot traces
#%% -- Get station lat-lon
    for tr in stgood:
        if tr.stats.station in st_names:  # find station in station list
            ii = st_names.index(tr.stats.station)
            stalon = float(st_lons[ii]) # look up lat & lon again to find distance
            stalat = float(st_lats[ii])
#%% -- Get distance and arrival time
            distance = gps2dist_azimuth(stalat,stalon,ev_lat,ev_lon)
            # atime   = model.get_travel_times(source_depth_in_km=ev_depth,distance_in_degree=distance[0], phase_list=[dphase])
            # tr.stats.atime = atime
            tr.stats.distance=distance[0]/(1000*111) # distance in degrees
            # print(tr.stats.station + '  ' + str(tr.stats.distance))
            if tr.stats.distance < min_dist_auto:
                min_dist_auto = tr.stats.distance
            if tr.stats.distance > max_dist_auto:
                max_dist_auto = tr.stats.distance
            if tr.stats.starttime - t < min_time_plot:
                min_time_plot = tr.stats.starttime - t
            if ((tr.stats.starttime - t) + ((len(tr.data)-1) * tr.stats.delta)) > max_time_plot:
                max_time_plot =  ((tr.stats.starttime - t) + ((len(tr.data)-1) * tr.stats.delta))
    print(f'        Min distance is   {min_dist_auto:.3f}   Max distance is {max_dist_auto:.3f}')
    print(f'        Min time is   {min_time_plot:.2f}   Max time is {max_time_plot:.2f}')
    if min_time_plot > start_buff:
        print(f'Min time {min_time_plot:.2f} > start_buff {start_buff:.2f}')
        print(colored('Warning: write zero-filling into pro3?','red'))
    #     sys.exit(-1)
    if max_time_plot < end_buff:
        print(f'Max time {max_time_plot:.2f} < end_buff {end_buff:.2f}')
        print(colored('Warning: write zero-filling into pro3?','red'))
    #     sys.exit(-1)

#%% -- Set x, y limits
    # plt.close(fig_index)

    # fig_tit = str(eq_num) + '_traces'
    # plt.figure(fig_tit,figsize=(10,10), num = fig_index)
    plt.figure(figsize=(10,10), num = fig_index)
    plt.xlim(min_time_plot,max_time_plot)

#%%  -- Set distance limits
    if plot_auto_dist == True:
        dist_diff = max_dist_auto - min_dist_auto # add space at extremes
        plt.ylim(min_dist_auto - 0.1 * dist_diff, max_dist_auto + 0.1 * dist_diff)
        max_dist = max_dist_auto
        min_dist = min_dist_auto
    else:
        plt.ylim(min_dist,max_dist)

    for tr in stgood:
        # list stations with distance to ID best traces for correlations
        # print('        Distance is ' + str(tr.stats.distance) + ' for station ' + tr.stats.station)
        dist_offset = tr.stats.distance
        ttt = np.arange(len(tr.data)) * tr.stats.delta + (tr.stats.starttime - t)

#%%  -- Set time limits
        if len(tr.data) > 0:  # plot if data exists
            if tr.data.max() - tr.data.min() > 0:
                if zoom:  # normalize plot amp to zoom window if zooming
                    get_window_max = tr.copy()
                    get_window_max.trim(starttime = Zstart_buff,endtime = Zend_buff)
                    if len(get_window_max.data) > 0:  # only norm if zoom window exists
                        if get_window_max.data.max() - get_window_max.data.min() > 0:
                            win_max = max(abs(get_window_max.data))
                            tr.data = tr.data/(1.5*win_max)
                plt.plot(ttt, (tr.data - np.median(tr.data))*plot_scale_fac /(tr.data.max() - tr.data.min()) + dist_offset, color = 'black')
            else:
                print('Max ' + str(tr.data.max()) + ' equals min ' + str(tr.data.min()) + ', no trace, skip plotting')
            if plot_sta_names:
                plt.text(min_time_plot + 0.015 * (max_time_plot - min_time_plot), dist_offset + 0.003 * (max_dist - min_dist), tr.stats.station)
        else:
            nodata += 1
            print('Trace ' + tr.stats.station + ' has : ' + str(len(tr.data)) + ' time pts, skip plotting')
#%% -- Traveltime curves
    if plot_tt:
        # first traveltime curve
        line_pts = 50
        dist_vec  = np.arange(min_dist_auto, max_dist_auto, (max_dist_auto - min_dist_auto)/line_pts) # distance grid
        time_vec1 = np.arange(min_dist_auto, max_dist_auto, (max_dist_auto - min_dist_auto)/line_pts) # empty time grid of same length (filled with -1000)

        for i in range(0,line_pts):
            arrivals = model.get_travel_times(source_depth_in_km=ev_depth,distance_in_degree=dist_vec[i],phase_list=[dphase])
            if(len(arrivals) == 0 and dist_vec[i] < 10 and dphase == 'P'):  # in case first arrival is upgoing P, which is 'p'
                arrivals = model.get_travel_times(source_depth_in_km=ev_depth,distance_in_degree=dist_vec[i],phase_list='p')
                if(len(arrivals) == 0):
                    print(f'travel_time failed: blue ttime1 dist {dist_vec[i]:.2f} phase {dphase}')
            num_arrivals = len(arrivals)
            found_it = 0
            for j in range(0,num_arrivals):
                if arrivals[j].name == dphase:
                    time_vec1[i] = arrivals[j].time
                    found_it = 1
            if found_it == 0:
                time_vec1[i] = np.nan
        # second traveltime curve
        if dphase2 != 'no':
            time_vec2 = np.arange(min_dist_auto, max_dist_auto, (max_dist_auto - min_dist_auto)/line_pts) # empty time grid of same length (filled with -1000)
            for i in range(0,line_pts):
                arrivals = model.get_travel_times(source_depth_in_km=ev_depth,distance_in_degree=dist_vec[i],phase_list=[dphase2])
                if(len(arrivals) == 0 and dist_vec[i] < 10 and dphase2 == 'P'):  # in case first arrival is upgoing P, which is 'p'
                    arrivals = model.get_travel_times(source_depth_in_km=ev_depth,distance_in_degree=dist_vec[i],phase_list='p')
                    if(len(arrivals) == 0):
                        print(f'travel_time failed: orange ttime2 dist {dist_vec[i]:.2f} phase {dphase2}')
                num_arrivals = len(arrivals)
                found_it = 0
                for j in range(0,num_arrivals):
                    if arrivals[j].name == dphase2:
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
            plt.plot(time_vec2,dist_vec, color = 'orange', label = dphase2)
        # third traveltime curve
        if dphase3 != 'no':
            time_vec3 = np.arange(min_dist_auto, max_dist_auto, (max_dist_auto - min_dist_auto)/line_pts) # empty time grid of same length (filled with -1000)
            for i in range(0,line_pts):
                arrivals = model.get_travel_times(source_depth_in_km=ev_depth,distance_in_degree=dist_vec[i],phase_list=[dphase3])
                if(len(arrivals) == 0 and dist_vec[i] < 10 and dphase3 == 'P'):  # in case first arrival is upgoing P, which is 'p'
                    arrivals = model.get_travel_times(source_depth_in_km=ev_depth,distance_in_degree=dist_vec[i],phase_list='p')
                    if(len(arrivals) == 0):
                        print(f'travel_time failed: yellow ttime3 dist {dist_vec[i]:.2f} phase {dphase3}')
                num_arrivals = len(arrivals)
                found_it = 0
                for j in range(0,num_arrivals):
                    if arrivals[j].name == dphase3:
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
            plt.plot(time_vec3,dist_vec, color = 'yellow', label = dphase3)
        # fourth traveltime curve
        if dphase4 != 'no':
            time_vec4 = np.arange(min_dist_auto, max_dist_auto, (max_dist_auto - min_dist_auto)/line_pts) # empty time grid of same length (filled with -1000)
            for i in range(0,line_pts):
                arrivals = model.get_travel_times(source_depth_in_km=ev_depth,distance_in_degree=dist_vec[i],phase_list=[dphase4])
                if(len(arrivals) == 0 and dist_vec[i] < 10 and dphase4 == 'P'):  # in case first arrival is upgoing P, which is 'p'
                    arrivals = model.get_travel_times(source_depth_in_km=ev_depth,distance_in_degree=dist_vec[i],phase_list='p')
                    if(len(arrivals) == 0):
                        print(f'travel_time failed: purple ttime4 dist {dist_vec[i]:.2f} phase {dphase4}')
                num_arrivals = len(arrivals)
                found_it = 0
                for j in range(0,num_arrivals):
                    if arrivals[j].name == dphase4:
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
            plt.plot(time_vec4,dist_vec, color = 'purple', label = dphase4)
        if   rel_time == 1:
            time_vec1 = (dist_vec - ref1_dist) * ref_slow
        if   rel_time == 3:
            time_vec1 = time_vec1 - time_vec1
        elif rel_time == 2 or rel_time == 4:
            time_vec1 = time_vec1 - atime_ref
        plt.plot(time_vec1 + shift_tt, dist_vec, color = 'blue', label = dphase)
        plt.plot((Zstart_buff, Zstart_buff), (min_dist,max_dist), color = 'red', label = 'beam start')
        plt.plot((Zend_buff, Zend_buff), (min_dist,max_dist), color = 'red', label = 'beam end')

    plt.xlabel('Time (s)')
    plt.ylabel('Epicentral distance from event (°)')
    plt.title(date_label + ' event #' + str(eq_num) + ' ' + dphase + ' freq ' + str(freq_min) + ' to ' + str(freq_max) + ' Hz')
    plt.legend(loc="upper left")
    os.chdir('/Users/vidale/Documents/Research/IC/Plots_hold')
    plt.savefig('section_' + date_label + '_' + str(eq_num) + '_' + str(int(freq_min)) + '-' + str(int(freq_max)))
    plt.show()


#%%  Save processed seismograms
    fname3 = '/Users/vidale/Documents/Research/IC/Pro_Files/HD' + date_label + 'sel.mseed'

    stgood.write(fname3,format = 'MSEED')

    elapsed_time_wc = time.time() - start_time_wc
    print(f'    This job took   {elapsed_time_wc:.1f}   seconds')
    os.system('say "three"')
