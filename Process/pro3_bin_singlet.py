#!/usr/bin/env python
# reads in raw mseed traces from a single event
# bins Hinet seismograms onto lat-lon squares
# this program tapers, filters, selects range and SNR
# plots against traveltime curves, either raw or reduced against traveltimes
# writes out "*sel_bin.mseed" file
# plot lines are blue, orange, yellow, purple for phases 1 through 4
# John Vidale 10/2021

def pro3bin(eq_num, stat_corr = 1, corr_threshold = 0, rel_time = 1, norm = 1,
            max_taper_length = 5., simple_taper = 1, skip_SNR = True, SNR_thres = 0,
            dphase = 'P', dphase2 = '', dphase3 = '', dphase4 = '',
            start_buff = -10, end_buff = 10,
            freq_min = 0.25, freq_max = 1, do_filt = 1, plot_scale_fac = 0.2,
            min_dist = 0, max_dist = 180, plot_auto_dist = True, do_decimate = 0,
            alt_statics = 0, statics_file = 'nothing', ARRAY = 0, JST = False, ref_loc = False, ref_rad = 0.4,
            verbose = 0, fig_index = 1):
# 0 is Hinet, 1 is LASA, 2 is NORSAR

#%% Import
#%% -- Functions
    from obspy import UTCDateTime
    from obspy import Stream, Trace
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

    print(colored('Running pro3b_sort_plot_singlet', 'cyan'))
    start_time_wc = time.time()
    if ARRAY != 0:
        print(colored('Only runs on Hinet - ARRAY must be 0!', 'red'))
        exit()

#%% -- Event info
    #  input event data with 1-line file of format
    #  event 2016-05-28T09:47:00.000 -56.241 -26.935 78
    fname = '/Users/vidale/Documents/Research/IC/EvLocs/event' + str(eq_num) + '.txt'
    print('Opening ' + fname)
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
    print(date_label + ' year_label ' + year_label + ' hour_label ' + hour_label + ' min_label ' + minute_label)
    ev_lat      = float(      split_line[2])
    ev_lon      = float(      split_line[3])
    ev_depth    = float(      split_line[4])
    print('        date_label ' + date_label + ' time ' + str(t) + ' lat ' + str(ev_lat) + ' lon ' + str( ev_lon) + ' depth ' + str(ev_depth))

#%% -- Station locations and statics
    if stat_corr == 1:  # load static terms, only applies to Hinet, LASA, and China
        if ARRAY == 0:
            if alt_statics == 0: # standard set
                sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/sta_statics_hinet.txt'
            else: # custom set made by this event for this event
                sta_file = ('/Users/vidale/Documents/GitHub/Array_codes/Files/Tuned_Statics/' + 'HA' +
                   date_label[:10] + 'pro4_' + dphase + '.statics')
        elif ARRAY == 1:
            sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/sta_statics_LASA.txt'
        elif ARRAY == 2:
            sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/sta_statics_ch.txt'
        with open(sta_file, 'r') as file:
            lines = file.readlines()
        print('    ' + str(len(lines)) + ' coarse station statics read from ' + sta_file)
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
                st_dist.append(split_line[1])
                st_lats.append( split_line[2])
                st_lons.append( split_line[3])
                st_shift.append(split_line[4])
                st_corr.append(split_line[5])
            elif ARRAY == 2:
                st_lats.append( split_line[1])
                st_lons.append( split_line[2])
                st_shift.append(split_line[3])
                st_corr.append(split_line[4]) # but really std dev

    else: # no static terms, always true for NORSAR
        if ARRAY == 0: # Hinet set
            sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/sta_hinet.txt'
        elif ARRAY == 1: #         LASA set
            sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/sta_LASA.txt'
        elif ARRAY == 2: #         China set
            sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/sta_ch.txt'
        else: #         NORSAR set
            sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/sta_NORSAR.txt'
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
            st_names.append(split_line[0])
            st_lats.append( split_line[1])
            st_lons.append( split_line[2])
    if ARRAY == 0:  # shorten and make upper case Hi-net station names to match station list
        for ii in station_index:
            this_name = st_names[ii]
            this_name_truc = this_name[0:5]
            st_names[ii]  = this_name_truc.upper()

#%%  -- Rest of parameters
#    stat_corr = 1 # apply station static corrections
#    rel_time = 1          # timing is relative to a chosen phase, otherwise relative to OT
#    dphase  = 'PKIKP'       # phase to be aligned
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
    plot_map           = False
    plot_input_traces  = True
    plot_binned_traces = True

    # if ref_loc ==true,  use distance to (ref_lat,ref_lon) to filter stations
    # if ref_loc ==false, use distance to earthquake loc to filter stations
    #    ref_rad = 0.4    # ° radius (°) set by input or at top
    if ARRAY == 0:
        ref_lat = 36  # °N, around middle of Japan
        ref_lon = 139 # °E
    if ARRAY == 1:
        ref_lat = 46.7      # °N keep only inner rings A-D if radius is 0.4°
        ref_lon = -106.22   # °E
    if ARRAY == 2:
        ref_lat = 38      # °N
        ref_lon = 104.5   # °E
#        ref_rad = 0.4    # ° radius (°) set by input or at top

#%% -- Test taper, needs adjustment?
#   Is taper too long compared to noise estimation window?
    totalt = end_buff - start_buff
    noise_time_skipped = taper_frac * totalt
    noise_time_skipped = min(noise_time_skipped,10.0) # set max of 10s to taper length
    if simple_taper == 0:
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
        skip_SNR = True
        SNR_thres = 0

#%% Load waveforms
    st = Stream()
    # fname     = '/Users/vidale/Documents/GitHub/LASA_data/HD' + date_label + '.mseed'
    if ARRAY == 1:
        mseed_name = year_short_label + month_label + day_label + '_' + hour_label + minute_label
        fname     = '/Users/vidale/Documents/Research/IC/Mseed/L' + mseed_name + '.mseed'
    else:
        mseed_name = year_short_label + '-' +month_label + '-' + day_label
        fname     = '/Users/vidale/Documents/Research/IC/Mseed/HD20' + mseed_name + '.mseed'

    print('Reading seismogram file: ' + fname)
    st=read(fname)
#%% -- Ensure 10 sps
    if do_decimate != 0:
        st.decimate(do_decimate, no_filter=True)

    print('        ' + fname)
    print('    ' + str(len(st)) + '  traces read in')
    if(len(st) == 0):
        exit('No traces means read failed')
    print('        First trace has : ' + str(len(st[0].data)) + ' time pts ')
    print('        Start time : ' + str(st[0].stats.starttime) + '  event time : ' + str(t))
    # print('    ' + str(len(st)) + '  traces after decimation')
    nt = len(st[0].data)
    dt = st[0].stats.delta
    print('        First trace has : ' + str(nt) + ' time pts, time sampling of ' + str(dt) + ' and thus duration of ' + str((nt-1)*dt))

#%%  -- Detrend, taper, filter, taper
    st.detrend(type='simple')
    st.taper(taper_frac, max_length = max_taper_length)
    if do_filt == 1:
        st.filter('bandpass', freqmin=freq_min, freqmax=freq_max, corners=4, zerophase=True)
    st.taper(taper_frac, max_length = max_taper_length)

#%% Make bins for stacking
#%% -- Define grid of bins
    lat_lo =  30  # rounds off lat_hi and lon_hi to a round increment of grid_delta
    lat_hi =  46
    lon_lo = 129
    lon_hi = 146
    grid_delta = 1 # mesh spacing, in °
    lat_n = int(round(1 + (lat_hi - lat_lo)/grid_delta))  # number of latitude nodes
    lon_n = int(round(1 + (lon_hi - lon_lo)/grid_delta))  # number of longitude nodes
    stack_nt = int(round(1 + ((end_buff - start_buff)/dt)))  # number of time points

    # In English, stack_slows = range(slow_n) * slow_delta - slow_lo
    a1lat = range(lat_n)
    a1lon = range(lon_n)
    stack_lats = [(x * grid_delta + lat_lo) for x in a1lat]
    stack_lons = [(x * grid_delta + lon_lo) for x in a1lon]

    # testing slownesses in indexing
    print(str(lat_n) + ' lats, ' + str(lon_n) + ' lons')
    print('Lats 0' + ' ' + str(stack_lats[0]) + '   ' 'end' + ' ' + str(stack_lats[-1]))
    print('Lons 1' + ' ' + str(stack_lons[0]) + '   ' 'end' + ' ' + str(stack_lons[-1]))

#%% -- Template empty-bin seismogram
    stack = Stream()
    tr = Trace()
    tr.stats.delta = dt
    tr.stats.starttime = t + start_buff
    tr.stats.npts = stack_nt
    tr.stats.network = 'stack'
    tr.stats.station = ''
    tr.stats.channel = 'BHZ'
    tr.data = np.zeros(stack_nt)
    done1 = 0
    done2 = 0
    bin_count = np.zeros(lat_n * lon_n)
    bin_lat   = np.zeros(lat_n * lon_n)
    bin_lon   = np.zeros(lat_n * lon_n)

#%% -- Seed bins with empty seismograms, initialize sgram count, compute (lat,lon) of each bin center
    for lat_one in stack_lats: # make a trace per bin, and find vector index for each (lat_index,lon_index) pair
        for lon_one in stack_lons:
            tr1 = tr.copy()
            tr1.stats.station = 'A' + str(int(round(lat_one))) +  str(int(round(lon_one-100)))
            stack.extend([tr1])
            rel_lat = (lat_one - lat_lo)/grid_delta  # grid point offset from lat at origin
            rel_lon = (lon_one - lon_lo)/grid_delta  # grid point offset from lon at origin
            index_lat = int(round(rel_lat))  # lat index
            index_lon = int(round(rel_lon))  # lat index
            index = int(round(index_lat*lon_n + index_lon))
            # print(tr1.stats.station + ' tr1.stats.station ' + str(rel_lat) + ' rel_lat ' + str(rel_lon) + ' rel_lon ' + str(lat_one) + ' lat_one ' + str(lon_one) + ' lon_one ' + str(index_lat) + ' index_lat ' + str(index_lon) + ' index_lon ' + str(index) + ' index')
            bin_lat[index] = lat_one
            bin_lon[index] = lon_one
            done1 += 1
    print(str(len(stack_lats)) + ' lats in grid ' + str(len(stack_lons)) + ' lons in grid')
    print(str(done1) + ' bins zeroed ')

#%% Loop over input traces
    binned_file = Stream()
    st_trimmed  = Stream()

    tra_in_range  = 0
    tra_sta_found = 0
    nodata        = 0
    min_dist_auto = 180
    max_dist_auto = 0
    min_time_plot =  1000000
    max_time_plot = -1000000

#%% -- Normalize, align chosen phase, add trace to bin sgram, keep count
    for tr in st: # stack traces one by one
        if JST: # if necessary, convert JST -> UTC, time in Greenwich 9 hours earlier than Japan
            tr.stats.starttime = tr.stats.starttime - 9*60*60
        if tr.stats.station in st_names:  # find station in station list
            ii = st_names.index(tr.stats.station)
            if stat_corr == 1: # or use value from correlation file
                corr = float(st_corr[ii])
            else:
                corr = 1 # 1 says correlation is high
            if corr > corr_threshold: # if using statics, reject low correlations

                if norm == 1:
                    tr.normalize() # trace divided abs(max of trace)

    #%% -- -- Find bin index and its coordinates
                stalat = float(st_lats[ii])
                stalon = float(st_lons[ii])
                dist_lat = (stalat - lat_lo)/float(grid_delta)
                dist_lon = (stalon - lon_lo)/float(grid_delta)
                index_lat = int(round(dist_lat))  # lat index
                index_lon = int(round(dist_lon))  # lon index
                index = int(round(index_lat*lon_n + index_lon)) # compute index of bin
                local_bin_lat = bin_lat[index]
                local_bin_lon = bin_lon[index]

    #%% -- -- Find bin distance, arrival time, slowness, and increment count
                ref_distance = gps2dist_azimuth(local_bin_lat,local_bin_lon,ev_lat,ev_lon)
                sta_distance = gps2dist_azimuth(stalat       ,stalon       ,ev_lat,ev_lon)
                ref_dist  = ref_distance[0]/(1000*111)
                dist      = sta_distance[0]/(1000*111)
                dist_minus = ref_dist - 0.5
                dist_plus  = ref_dist + 0.5

                arrivals_ref   = model.get_travel_times(source_depth_in_km=ev_depth,distance_in_degree=ref_dist  ,phase_list=[dphase])
                arrivals_minus = model.get_travel_times(source_depth_in_km=ev_depth,distance_in_degree=dist_minus,phase_list=[dphase])
                arrivals_plus  = model.get_travel_times(source_depth_in_km=ev_depth,distance_in_degree=dist_plus ,phase_list=[dphase])
                atime_ref = arrivals_ref[0].time  # phase arrival time at reference distance
                ref_slow = arrivals_plus[0].time - arrivals_minus[0].time  # dt over 1 degree at ref distance

    # distance = gps2dist_azimuth(stalat,stalon,ev_lat,ev_lon) # Get traveltimes again, hard to store
    # tr.stats.distance=sta_distance[0] # distance in km
    # dist = sta_distance[0]/(1000*111)

                if(len(arrivals_ref) == 0 or len(arrivals_minus) == 0 or len(arrivals_plus) == 0):
                    print('model.get_travel_times failed: dist, phase  ' + str(ref_dist) + '   ' + dphase)
                in_range = 0  # flag for whether this trace goes into stack
                rejector = False  # flag in case model.get_travel_times fails

    #%% -- -- Check whether in distance range or within distance circle
                if ref_loc == False:  # check whether trace is in distance range from earthquake
                    if min_dist < ref_dist and ref_dist < max_dist:
                        in_range = 1
                elif ref_loc == True:  # alternately, check whether trace is close enough to the ref_location
                    ref_distance = gps2dist_azimuth(ref_lat,ref_lon,stalat,stalon)
                    ref2_dist = ref_distance[0]/(1000*111)
                    if ref2_dist < ref_rad:
                        in_range = 1

                if in_range == 1:   # trace fulfills the specified criteria for being in range
                    tra_in_range += 1
                    s_t = t + start_buff
                    e_t = t + end_buff

#%% -- Apply static
                    if stat_corr == 1: # apply static station corrections
                        tr.stats.starttime -= float(st_shift[ii])

#%% -- Four rel_time options
                    if rel_time == 0:  #  use absolute time, ignore time of chosen phase
                        tr.trim(starttime=s_t,endtime = e_t)
                    else:              # shift relative to a chosen phase
                        # print(f'        dist_m is   {dist:.2f}   ev_depth is {ev_depth:.2f}   dphase is {dphase}')
                        arrivals_each = model.get_travel_times(source_depth_in_km=ev_depth,distance_in_degree=dist,phase_list=[dphase])
                        # print(f'        cal_time is   {arrivals_each[0].time:.2f}')
                        if(len(arrivals_each) == 0):
                            if(dist < 10 and dphase == 'P'):  # in case first arrival is upgoing P, try 'p'
                                arrivals_each = model.get_travel_times(source_depth_in_km=ev_depth,distance_in_degree=dist,phase_list='p')
                        if(len(arrivals_each) == 0):  # did it still fail?
                            print('model.get_travel_times failed: dist, depth, phase  ' + tr.stats.station + '   ' + str(ref_dist) + '   ' +
                                  '   ' + str(ev_depth) + '   ' + dphase)
                            tra_in_range -= 1 # don't count this trace after all
                            rejector = True
                        else:
                            # atime_each is phase arrival on this trace
                            # atime_ref  is phase arrival on reference trace

                            atime_each = arrivals_each[0].time
                            # start time has a shift proportional to ref_dist at phase slowness at ref_dist
                            # preserves slowness of arrivals at ref_dist
                            # uniform relative window limits around phase
                            if rel_time == 1:
                                s_t += atime_each
                                e_t += atime_each
                                tr.trim( starttime=s_t, endtime = e_t)
                                tr.stats.starttime -= atime_each - (dist-ref_dist) * ref_slow
                                # print('station ' + tr.stats.station + ' Start ' + str(s_t)  + ' End ' + str(e_t))
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

    #%% -- Window, shift, and add to stack
                    #%% -- -- Stack
                    if len(tr.data) >= stack_nt and rejector == False:
                        for it in range(stack_nt):  # check points one at a time
                            # adjust time to stack chosen phase to its time at bin center
                            # start with shifts of zero
                            # time_lag  = radial distance from center * chosen phase slowness  # time shift due to radial slowness

                            # time_lag  = del_distR * stack_Rslows[slowR_i]  # time shift due to radial slowness
                            # time_lag = 0
                            # time_correction = ((t-tr.stats.starttime) + (time_lag + start_buff))/dt
                            # it_in = int(round(it + time_correction))

                            # if it_in >= 0 and it_in < nt - 1: # does data lie within seismogram?
                            #     stack[indx].data[it] += tr[it_in]

                            time_correction = ((t + start_buff) - tr.stats.starttime)/dt
                            it_in = int(round(it + time_correction))

                            if it_in >= 0 and it_in < stack_nt - 1: # does data lie within seismogram?
                                stack[index].data[it] += tr[it_in]
                            # stack[index].data[it] += tr[it]
                        st_trimmed += tr
                    #%% -- -- Reject if not in time window
                    else:
                        nodata += 1 # one could complicate this code to allow partial traces into stack

                    bin_count[index] += 1
                    done2 += 1
                    if done2 % 100 == 0:
                        print('Done with ' + str(done2))
                    # add trace to stack
                    # print('stalat,lon ' + str(stalat) + ' ' + str(stalon) + ' index_lat,lon vector index ' + str(index_lat) + ' ' + str(index_lon) + ' ' + str(index))
        else:
            print(tr.stats.station + ' not found in station list')
#%% -- Reject if not in station (static) list
        # else:
        #     print(tr.stats.station + ' not found in station list with statics')
    # print('    ' + str(tra_in_range) + '  traces in range')
    # print('    ' + str(len(st_pickalign)) + '  traces after alignment and correlation selection')
    # print('    ' + str(nodata) + '  traces with no data')

    non_zero = 0
    for i in range(int(lat_n * lon_n)):
        if bin_count[i] > 0:
            non_zero += 1
            if bin_count[i] > 0:
                stack[i].data /= float(bin_count[i])
            binned_file += stack[i]

    print(str(done2) + ' input traces counted')
    print(str(non_zero) + ' grids have traces, and ' + str(len(binned_file)) + ' traces in collection')

#%%  Save binned seismograms
    fname3 = '/Users/vidale/Documents/Research/IC/Pro_Files/HD' + date_label + 'binned.mseed'

    binned_file.write(fname3,format = 'MSEED')

    elapsed_time_wc = time.time() - start_time_wc
    print(f'    To start of plotting took   {elapsed_time_wc:.1f}   seconds')
    os.system('say "binning done"')

#%% Plot map of station density
    if plot_map:
        Z = np.reshape(bin_count,(lat_n,lon_n))

        plt.close('all')
        plt.figure(1,figsize=(10,10))
        plt.imshow(Z, extent=[129,146,30,46], cmap='jet',
               vmin = 0, vmax = 10, origin='lowest', aspect='auto')
        plt.colorbar()
        plt.show()

#%% Plot input traces
    if plot_input_traces == True:
#%% -- Re-find station lat-lon
        for tr in st_trimmed:
            if tr.stats.station in st_names:  # find station in station list
                ii = st_names.index(tr.stats.station)
                stalon = float(st_lons[ii]) # look up lat & lon again to find distance
                stalat = float(st_lats[ii])
    #%% -- Re-calculate distance and arrival time
                distance = gps2dist_azimuth(stalat,stalon,ev_lat,ev_lon)
                # atime   = model.get_travel_times(source_depth_in_km=ev_depth,distance_in_degree=distance[0], phase_list=[dphase])
                # tr.stats.atime = atime
                tr.stats.distance=distance[0]/(1000*111) # distance in degrees
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
            print(colored('Write zero-filling into pro3 for this code to work','red'))
            sys.exit(-1)
        if max_time_plot < end_buff:
            print(f'Max time {max_time_plot:.2f} < end_buff {end_buff:.2f}')
            print(colored('Write zero-filling into pro3 for this code to work','red'))
            sys.exit(-1)

    #%% -- Set x, y limits
        # plt.close(fig_index)
        fig_tit = str(eq_num) + '_input'
        plt.figure(fig_tit,figsize=(10,10))
        plt.xlim(min_time_plot,max_time_plot)

    #%%  -- Set distance limits
        if plot_auto_dist == True:
            dist_diff = max_dist_auto - min_dist_auto # add space at extremes
            min_dist = min_dist_auto - 0.1 * dist_diff
            max_dist = max_dist_auto + 0.1 * dist_diff
        plt.ylim(min_dist,max_dist)

        print('            Plotting - First trace has : ' + str(len(st_trimmed[0].data)) + ' time pts ')
        print('            Start time : ' + str(st_trimmed[0].stats.starttime) + '  event time : ' + str(t))
        nt = len(st_trimmed[0].data)
        dt = st_trimmed[0].stats.delta
        print('            First trace has : ' + str(nt) + ' time pts, time sampling of ' + str(dt) + ' and thus duration of ' + str((nt-1)*dt))
        print('            ' + str(len(st_trimmed)) + '  traces after decimation')

        for tr in st_trimmed:
            dist_offset = tr.stats.distance

            atime   = model.get_travel_times(source_depth_in_km=ev_depth,distance_in_degree=tr.stats.distance, phase_list=[dphase])
            ttt = np.arange(len(tr.data)) * tr.stats.delta + (tr.stats.starttime - t)
            # print('first ttt = ' + str(ttt[0]) + '  last ttt = ' + str(ttt[-1]) + '  ' + tr.stats.station)

    #%%  -- Set time limits
            if len(tr.data) > 0:
                if tr.data.max() - tr.data.min() > 0:
                    plt.plot(ttt, (tr.data - np.median(tr.data))*plot_scale_fac /(tr.data.max() - tr.data.min()) + dist_offset, color = 'black')
                else:
                    print('Max ' + str(tr.data.max()) + ' equals min ' + str(tr.data.min()) + ', skip plotting')
            else:
                nodata += 1
                print('Trace ' + tr.stats.station + ' has : ' + str(len(tr.data)) + ' time pts, skip plotting')
        #%% -- Traveltime curves
        if plot_tt:
            # first traveltime curve
            line_pts = 50
            dist_vec  = np.arange(min_dist_auto, max_dist_auto, (max_dist_auto - min_dist_auto)/line_pts) # distance grid
            time_vec1 = np.arange(min_dist_auto, max_dist_auto, (max_dist_auto - min_dist_auto)/line_pts) # empty time grid of same length (filled with -1000)

            arrivals_ref   = model.get_travel_times(source_depth_in_km=ev_depth,distance_in_degree=ref_dist, phase_list=[dphase])
            if(len(arrivals_ref) == 0 and ref_dist < 10 and dphase == 'P'):  # in case first arrival is upgoing P, which is 'p'
                arrivals_ref   = model.get_travel_times(source_depth_in_km=ev_depth,distance_in_degree=ref_dist, phase_list='p')
            atime_ref = arrivals_ref[0].time  # phase arrival time at reference distance


            for i in range(0,line_pts):
                arrivals = model.get_travel_times(source_depth_in_km=ev_depth,distance_in_degree=dist_vec[i],phase_list=[dphase])
                if(len(arrivals) == 0 and dist_vec[i] < 10 and dphase == 'P'):  # in case first arrival is upgoing P, which is 'p'
                    arrivals = model.get_travel_times(source_depth_in_km=ev_depth,distance_in_degree=dist_vec[i],phase_list='p')
                    if(len(arrivals) == 0):
                        print('model.get_travel_times failed: dist, phase  ' + str(dist_vec[i]) + '   ' + dphase)
                num_arrivals = len(arrivals)
                found_it = 0
                for j in range(0,num_arrivals):
                    if arrivals[j].name == dphase:
                        time_vec1[i] = arrivals[j].time
                        found_it = 1
                if found_it == 0:
                    time_vec1[i] = np.nan
            if   rel_time == 3 or rel_time == 4 or rel_time == 1:
                time_vec1 = time_vec1 - time_vec1
            elif rel_time == 2:
                time_vec1 = time_vec1 - atime_ref
            plt.plot(time_vec1,dist_vec, color = 'blue')
            # second traveltime curve
            if dphase2 != 'no':
                time_vec2 = np.arange(min_dist_auto, max_dist_auto, (max_dist_auto - min_dist_auto)/line_pts) # empty time grid of same length (filled with -1000)
                for i in range(0,line_pts):
                    arrivals = model.get_travel_times(source_depth_in_km=ev_depth,distance_in_degree=dist_vec[i],phase_list=[dphase2])
                    if(len(arrivals) == 0 and dist_vec[i] < 10 and dphase2 == 'P'):  # in case first arrival is upgoing P, which is 'p'
                        arrivals = model.get_travel_times(source_depth_in_km=ev_depth,distance_in_degree=dist_vec[i],phase_list='p')
                        if(len(arrivals) == 0):
                            print('model.get_travel_times failed: dist, phase  ' + str(dist_vec[i]) + '   ' + dphase2)
                    num_arrivals = len(arrivals)
                    found_it = 0
                    for j in range(0,num_arrivals):
                        if arrivals[j].name == dphase2:
                            time_vec2[i] = arrivals[j].time
                            found_it = 1
                    if found_it == 0:
                        time_vec2[i] = np.nan
                if   rel_time == 3 or rel_time == 4 or rel_time == 1:
                    time_vec2 = time_vec2 - time_vec1
                elif rel_time == 2:
                    time_vec2 = time_vec2 - atime_ref
                plt.plot(time_vec2,dist_vec, color = 'orange')
            # third traveltime curve
            if dphase3 != 'no':
                time_vec3 = np.arange(min_dist_auto, max_dist_auto, (max_dist_auto - min_dist_auto)/line_pts) # empty time grid of same length (filled with -1000)
                for i in range(0,line_pts):
                    arrivals = model.get_travel_times(source_depth_in_km=ev_depth,distance_in_degree=dist_vec[i],phase_list=[dphase3])
                    if(len(arrivals) == 0 and dist_vec[i] < 10 and dphase3 == 'P'):  # in case first arrival is upgoing P, which is 'p'
                        arrivals = model.get_travel_times(source_depth_in_km=ev_depth,distance_in_degree=dist_vec[i],phase_list='p')
                        if(len(arrivals) == 0):
                            print('model.get_travel_times failed: dist, phase  ' + str(dist_vec[i]) + '   ' + dphase3)
                    num_arrivals = len(arrivals)
                    found_it = 0
                    for j in range(0,num_arrivals):
                        if arrivals[j].name == dphase3:
                            time_vec3[i] = arrivals[j].time
                            found_it = 1
                    if found_it == 0:
                        time_vec3[i] = np.nan
                if   rel_time == 3 or rel_time == 4 or rel_time == 1:
                    time_vec3 = time_vec3 - time_vec1
                elif rel_time == 2:
                    time_vec3 = time_vec3 - atime_ref
                plt.plot(time_vec3,dist_vec, color = 'yellow')
            # fourth traveltime curve
            if dphase4 != 'no':
                time_vec4 = np.arange(min_dist_auto, max_dist_auto, (max_dist_auto - min_dist_auto)/line_pts) # empty time grid of same length (filled with -1000)
                for i in range(0,line_pts):
                    arrivals = model.get_travel_times(source_depth_in_km=ev_depth,distance_in_degree=dist_vec[i],phase_list=[dphase4])
                    if(len(arrivals) == 0 and dist_vec[i] < 10 and dphase4 == 'P'):  # in case first arrival is upgoing P, which is 'p'
                        arrivals = model.get_travel_times(source_depth_in_km=ev_depth,distance_in_degree=dist_vec[i],phase_list='p')
                        if(len(arrivals) == 0):
                            print('model.get_travel_times failed: dist, phase  ' + str(dist_vec[i]) + '   ' + dphase4)
                    num_arrivals = len(arrivals)
                    found_it = 0
                    for j in range(0,num_arrivals):
                        if arrivals[j].name == dphase4:
                            time_vec4[i] = arrivals[j].time
                            found_it = 1
                    if found_it == 0:
                        time_vec4[i] = np.nan
                if   rel_time == 3 or rel_time == 4 or rel_time == 1:
                    time_vec4 = time_vec4 - time_vec1
                elif rel_time == 2:
                    time_vec4 = time_vec4 - atime_ref
                plt.plot(time_vec4,dist_vec, color = 'purple')

            plt.xlabel('Time (s)')
            plt.ylabel('Epicentral distance from event (°)')
            plt.title('Input: ' + date_label + ' event #' + str(eq_num))
        #    os.chdir('/Users/vidale/Documents/PyCode/Plots')
        #    plt.savefig(date_label + '_' + str(event_no) + '_raw.png')
            plt.show()

    nodata        = 0
    min_dist_auto = 180
    max_dist_auto = 0
    min_time_plot =  1000000
    max_time_plot = -1000000

#%% Plot binned traces
    if plot_binned_traces:
#%% -- Re-find station lat-lon
        for tr in binned_file:
            # if tr.stats.station in st_names:  # find station in station list
            #     ii = st_names.index(tr.stats.station)
            #     stalon = float(st_lons[ii]) # look up lat & lon again to find distance
            #     stalat = float(st_lats[ii])
            stalat = float(tr.stats.station[1:3]) + grid_delta/2.
            stalon = float(tr.stats.station[3:5]) + grid_delta/2. + 100
            # print(tr.stats.station + ' lat ' + str(stalat) + ' lon ' + str(stalon))

    #%% -- Re-calculate distance and arrival time
            distance = gps2dist_azimuth(stalat,stalon,ev_lat,ev_lon)
            # atime   = model.get_travel_times(source_depth_in_km=ev_depth,distance_in_degree=distance[0], phase_list=[dphase])
            # tr.stats.atime = atime
            tr.stats.distance=distance[0]/(1000*111) # distance in degrees
            # print(f'tr.stats.station - t {t:.2f} < tr.stats.starttime {tr.stats.starttime:.2f}')
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
            print(colored('Write zero-filling into pro3 for this code to work','red'))
            sys.exit(-1)
        if max_time_plot < end_buff:
            print(f'Max time {max_time_plot:.2f} < end_buff {end_buff:.2f}')
            print(colored('Write zero-filling into pro3 for this code to work','red'))
            sys.exit(-1)

    #%% -- Set x, y limits
        # plt.close(fig_index)
        fig_tit = str(eq_num) + '_binned'
        plt.figure(fig_tit,figsize=(10,10))
        plt.xlim(min_time_plot,max_time_plot)

    #%%  -- Set distance limits
        print('min_dist is ' + str(min_dist) + ' max_dist is ' + str(max_dist))
        if plot_auto_dist == True:
            dist_diff = max_dist_auto - min_dist_auto # add space at extremes
            min_dist = min_dist_auto - 0.1 * dist_diff
            max_dist = max_dist_auto + 0.1 * dist_diff
        plt.ylim(min_dist,max_dist)
        print('After : min_dist is ' + str(min_dist) + ' max_dist is ' + str(max_dist))

        print('            Plotting - First trace has : ' + str(len(binned_file[0].data)) + ' time pts ')
        print('            Start time : ' + str(binned_file[0].stats.starttime) + '  event time : ' + str(t))
        nt = len(binned_file[0].data)
        dt = binned_file[0].stats.delta
        print('            First trace has : ' + str(nt) + ' time pts, time sampling of ' + str(dt) + ' and thus duration of ' + str((nt-1)*dt))
        print('            ' + str(len(binned_file)) + '  traces after decimation')

        for tr in binned_file:
            dist_offset = tr.stats.distance

            atime   = model.get_travel_times(source_depth_in_km=ev_depth,distance_in_degree=tr.stats.distance, phase_list=[dphase])
            ttt = np.arange(len(tr.data)) * tr.stats.delta + (tr.stats.starttime - t)
            # print('first ttt = ' + str(ttt[0]) + '  last ttt = ' + str(ttt[-1]) + '  ' + tr.stats.station)

    #%%  -- Set time limits
            if len(tr.data) > 0:
                if tr.data.max() - tr.data.min() > 0:
                    plt.plot(ttt, (tr.data - np.median(tr.data))*plot_scale_fac /(tr.data.max() - tr.data.min()) + dist_offset, color = 'black')
                else:
                    print('Max ' + str(tr.data.max()) + ' equals min ' + str(tr.data.min()) + ', skip plotting')
            else:
                nodata += 1
                print('Trace ' + tr.stats.station + ' has : ' + str(len(tr.data)) + ' time pts, skip plotting')
        #%% -- Traveltime curves
        if plot_tt:
            # first traveltime curve
            line_pts = 50
            dist_vec  = np.arange(min_dist_auto, max_dist_auto, (max_dist_auto - min_dist_auto)/line_pts) # distance grid
            time_vec1 = np.arange(min_dist_auto, max_dist_auto, (max_dist_auto - min_dist_auto)/line_pts) # empty time grid of same length (filled with -1000)

            arrivals_ref   = model.get_travel_times(source_depth_in_km=ev_depth,distance_in_degree=ref_dist, phase_list=[dphase])
            if(len(arrivals_ref) == 0 and ref_dist < 10 and dphase == 'P'):  # in case first arrival is upgoing P, which is 'p'
                arrivals_ref   = model.get_travel_times(source_depth_in_km=ev_depth,distance_in_degree=ref_dist, phase_list='p')
            atime_ref = arrivals_ref[0].time  # phase arrival time at reference distance


            for i in range(0,line_pts):
                arrivals = model.get_travel_times(source_depth_in_km=ev_depth,distance_in_degree=dist_vec[i],phase_list=[dphase])
                if(len(arrivals) == 0 and dist_vec[i] < 10 and dphase == 'P'):  # in case first arrival is upgoing P, which is 'p'
                    arrivals = model.get_travel_times(source_depth_in_km=ev_depth,distance_in_degree=dist_vec[i],phase_list='p')
                    if(len(arrivals) == 0):
                        print('model.get_travel_times failed: dist, phase  ' + str(dist_vec[i]) + '   ' + dphase)
                num_arrivals = len(arrivals)
                found_it = 0
                for j in range(0,num_arrivals):
                    if arrivals[j].name == dphase:
                        time_vec1[i] = arrivals[j].time
                        found_it = 1
                if found_it == 0:
                    time_vec1[i] = np.nan
            if   rel_time == 1:
                time_vec1 = time_vec1 - time_vec1 + (dist_vec - ref_dist) * ref_slow
            if   rel_time == 3:
                time_vec1 = time_vec1 - time_vec1
            elif rel_time == 2 or rel_time == 4:
                time_vec1 = time_vec1 - atime_ref
            plt.plot(time_vec1,dist_vec, color = 'blue')
            # second traveltime curve
            if dphase2 != 'no':
                time_vec2 = np.arange(min_dist_auto, max_dist_auto, (max_dist_auto - min_dist_auto)/line_pts) # empty time grid of same length (filled with -1000)
                for i in range(0,line_pts):
                    arrivals = model.get_travel_times(source_depth_in_km=ev_depth,distance_in_degree=dist_vec[i],phase_list=[dphase2])
                    if(len(arrivals) == 0 and dist_vec[i] < 10 and dphase2 == 'P'):  # in case first arrival is upgoing P, which is 'p'
                        arrivals = model.get_travel_times(source_depth_in_km=ev_depth,distance_in_degree=dist_vec[i],phase_list='p')
                        if(len(arrivals) == 0):
                            print('model.get_travel_times failed: dist, phase  ' + str(dist_vec[i]) + '   ' + dphase2)
                    num_arrivals = len(arrivals)
                    found_it = 0
                    for j in range(0,num_arrivals):
                        if arrivals[j].name == dphase2:
                            time_vec2[i] = arrivals[j].time
                            found_it = 1
                    if found_it == 0:
                        time_vec2[i] = np.nan
                if   rel_time == 1:
                    time_vec2 = time_vec2 - time_vec1 + (dist_vec - ref_dist) * ref_slow
                if   rel_time == 3:
                    time_vec2 = time_vec2 - time_vec1
                elif rel_time == 2 or rel_time == 4:
                    time_vec2 = time_vec2 - atime_ref
                plt.plot(time_vec2,dist_vec, color = 'orange')
            # third traveltime curve
            if dphase3 != 'no':
                time_vec3 = np.arange(min_dist_auto, max_dist_auto, (max_dist_auto - min_dist_auto)/line_pts) # empty time grid of same length (filled with -1000)
                for i in range(0,line_pts):
                    arrivals = model.get_travel_times(source_depth_in_km=ev_depth,distance_in_degree=dist_vec[i],phase_list=[dphase3])
                    if(len(arrivals) == 0 and dist_vec[i] < 10 and dphase3 == 'P'):  # in case first arrival is upgoing P, which is 'p'
                        arrivals = model.get_travel_times(source_depth_in_km=ev_depth,distance_in_degree=dist_vec[i],phase_list='p')
                        if(len(arrivals) == 0):
                            print('model.get_travel_times failed: dist, phase  ' + str(dist_vec[i]) + '   ' + dphase3)
                    num_arrivals = len(arrivals)
                    found_it = 0
                    for j in range(0,num_arrivals):
                        if arrivals[j].name == dphase3:
                            time_vec3[i] = arrivals[j].time
                            found_it = 1
                    if found_it == 0:
                        time_vec3[i] = np.nan
                if   rel_time == 1:
                    time_vec3 = time_vec3 - time_vec1 + (dist_vec - ref_dist) * ref_slow
                if   rel_time == 3:
                    time_vec3 = time_vec3 - time_vec1
                elif rel_time == 2 or rel_time == 4:
                    time_vec3 = time_vec3 - atime_ref
                plt.plot(time_vec3,dist_vec, color = 'yellow')
            # fourth traveltime curve
            if dphase4 != 'no':
                time_vec4 = np.arange(min_dist_auto, max_dist_auto, (max_dist_auto - min_dist_auto)/line_pts) # empty time grid of same length (filled with -1000)
                for i in range(0,line_pts):
                    arrivals = model.get_travel_times(source_depth_in_km=ev_depth,distance_in_degree=dist_vec[i],phase_list=[dphase4])
                    if(len(arrivals) == 0 and dist_vec[i] < 10 and dphase4 == 'P'):  # in case first arrival is upgoing P, which is 'p'
                        arrivals = model.get_travel_times(source_depth_in_km=ev_depth,distance_in_degree=dist_vec[i],phase_list='p')
                        if(len(arrivals) == 0):
                            print('model.get_travel_times failed: dist, phase  ' + str(dist_vec[i]) + '   ' + dphase4)
                    num_arrivals = len(arrivals)
                    found_it = 0
                    for j in range(0,num_arrivals):
                        if arrivals[j].name == dphase4:
                            time_vec4[i] = arrivals[j].time
                            found_it = 1
                    if found_it == 0:
                        time_vec4[i] = np.nan
                if   rel_time == 1:
                    time_vec4 = time_vec4 - time_vec1 + (dist_vec - ref_dist) * ref_slow
                if   rel_time == 3:
                    time_vec4 = time_vec4 - time_vec1
                elif rel_time == 2 or rel_time == 4:
                    time_vec4 = time_vec4 - atime_ref
                plt.plot(time_vec4,dist_vec, color = 'purple')

            plt.xlabel('Time (s)')
            plt.ylabel('Epicentral distance from event (°)')
            plt.title(dphase + ' binned ' + date_label + ' event #' + str(eq_num))
        #    os.chdir('/Users/vidale/Documents/PyCode/Plots')
        #    plt.savefig(date_label + '_' + str(event_no) + '_raw.png')
            plt.show()

    elapsed_time_wc = time.time() - start_time_wc
    print(f'    Whole binning code took   {elapsed_time_wc:.1f}   seconds')
    os.system('say "Whole binning code is done"')
