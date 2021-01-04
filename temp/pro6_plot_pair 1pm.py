#!/usr/bin/env python
# Read in 2D stacks for two events
# Compute tdiff, ave_amp, amp_ratio
# Plot radial and transverse cuts through stack, plus beam sum
# Write out tdiff, ave_amp, amp_ratio results
# John Vidale 3/2019

def pro6stacked_pair(eq_file1, eq_file2, plot_scale_fac = 0.03, slow_delta = 0.0005,
              slowR_lo = -0.1, slowR_hi = 0.1, slowT_lo = -0.1, slowT_hi = 0.1,
              start_buff = -50, end_buff = 50, norm = 0, freq_corr = 1.0,
              plot_dyn_range = 1000, fig_index = 401, get_stf = 0, start_beam = 0, end_beam = 0,
              ARRAY = 0, max_rat = 1.8, min_amp = 0.2, turn_off_black = 0,
              R_slow_plot = 0, T_slow_plot = 0, tdiff_clip = 0.5, decimate_fac = 0,
              ref_loc = 0, ref_lat = 36.3, ref_lon = 138.5, dphase = 'PKiKP'):

    import obspy
    import obspy.signal
    from obspy import UTCDateTime
    from obspy import Stream, Trace
    from obspy import read
    from obspy.geodetics import gps2dist_azimuth
    import numpy as np
    import os
    from obspy.taup import TauPyModel
    import obspy.signal as sign
    import matplotlib.pyplot as plt
    model = TauPyModel(model='iasp91')
    from scipy.signal import hilbert
    import math
    import time
    import statistics

#%% Get info
    #%% get locations
    print('Running pro6_plot_stacked_seis')
    start_time_wc = time.time()

    sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/events_good.txt'
    with open(sta_file, 'r') as file:
        lines = file.readlines()
    event_count = len(lines)

    print(str(event_count) + ' lines read from ' + sta_file)
    # Load station coords into arrays
    # Only lat, lon, depth, baz, radslo's are used
    station_index = range(event_count)
    event_names        = []

    event_index = np.zeros(event_count)
    event_year  = np.zeros(event_count)
    event_mo    = np.zeros(event_count)
    event_day   = np.zeros(event_count)
    event_hr    = np.zeros(event_count)
    event_min   = np.zeros(event_count)
    event_sec   = np.zeros(event_count)
    event_lat   = np.zeros(event_count)
    event_lon   = np.zeros(event_count)
    event_dep   = np.zeros(event_count)
    event_mb    = np.zeros(event_count)
    event_ms    = np.zeros(event_count)
    event_tstart       = np.zeros(event_count)
    event_tend         = np.zeros(event_count)
    event_gcdist       = np.zeros(event_count)
    event_dist         = np.zeros(event_count)
    event_baz          = np.zeros(event_count)
    event_SNR          = np.zeros(event_count)
    event_Sflag        = np.zeros(event_count)
    event_PKiKPflag    = np.zeros(event_count)
    event_ICSflag      = np.zeros(event_count)
    event_PKiKP_radslo = np.zeros(event_count)
    event_PKiKP_traslo = np.zeros(event_count)
    event_PKiKP_qual   = np.zeros(event_count)
    event_ICS_qual     = np.zeros(event_count)

    iii = 0
    for ii in station_index:   # read file
        line = lines[ii]
        split_line = line.split()

        event_index[ii]  = float(split_line[0])
        event_names.append(split_line[1])
        event_year[ii]   = float(split_line[2])
        event_mo[ii]     = float(split_line[3])
        event_day[ii]    = float(split_line[4])
        event_hr[ii]     = float(split_line[5])
        event_min[ii]    = float(split_line[6])
        event_sec[ii]    = float(split_line[7])
        event_lat[ii]    = float(split_line[8])
        event_lon[ii]    = float(split_line[9])
        event_dep[ii]    = float(split_line[10])
        event_mb[ii]     = float(split_line[11])
        event_ms[ii]     = float(split_line[12])
        event_tstart[ii] = float(split_line[13])
        event_tend[ii]   = float(split_line[14])
        event_gcdist[ii] = float(split_line[15])
        event_dist[ii]   = float(split_line[16])
        event_baz[ii]    = float(split_line[17])
        event_SNR[ii]    = float(split_line[18])
        event_Sflag[ii]  = float(split_line[19])
        event_PKiKPflag[ii]     = float(split_line[20])
        event_ICSflag[ii]       = float(split_line[21])
        event_PKiKP_radslo[ii]  = float(split_line[22])
        event_PKiKP_traslo[ii]  = float(split_line[23])
        event_PKiKP_qual[ii]    = float(split_line[24])
        event_ICS_qual[ii]      = float(split_line[25])
#        print('Event ' + str(ii) + ' is ' + str(event_index[ii]))

    #  find predicted slowness
    if ref_loc == 0:
        if ARRAY == 0:
            ref_lat = 36.3  # °N, around middle of Japan
            ref_lon = 138.5 # °E
        elif ARRAY == 1:
            ref_lat = 46.7  # °N keep only inner rings A-D
            ref_lon = -106.22   # °E
        elif ARRAY == 2: # China set and center
            ref_lat = 38      # °N
            ref_lon = 104.5   # °E
    ref_dist_az = gps2dist_azimuth(event_lat[iii],event_lon[iii],ref_lat,ref_lon)
    ref_back_az = ref_dist_az[2]
    ref_dist    = ref_dist_az[0]/(1000*111)  # distance in °

    arrivals1 = model.get_travel_times(source_depth_in_km=event_dep[iii],distance_in_degree=ref_dist-0.5,phase_list=[dphase])
    arrivals2 = model.get_travel_times(source_depth_in_km=event_dep[iii],distance_in_degree=ref_dist+0.5,phase_list=[dphase])
    dtime = arrivals2[0].time - arrivals1[0].time
    event_pred_slo  = dtime/111.  # s/km

    # convert to pred rslo and tslo
    sin_baz = np.sin(ref_back_az * np.pi /180)
    cos_baz = np.cos(ref_back_az * np.pi /180)
    pred_Nslo = event_pred_slo * cos_baz
    pred_Eslo = event_pred_slo * sin_baz

    #  rotate observed slowness to N and E
    obs_Nslo = (event_PKiKP_radslo[iii] * cos_baz) - (event_PKiKP_traslo[iii] * sin_baz)
    obs_Eslo = (event_PKiKP_radslo[iii] * sin_baz) + (event_PKiKP_traslo[iii] * cos_baz)

    print('PR '+ str(pred_Nslo) + ' PT ' + str(pred_Eslo) + ' OR ' + str(obs_Nslo) + ' OT ' + str(obs_Eslo))
    #  find observed back-azimuth
    #    bazi_rad = np.arctan(event_PKiKP_traslo[ii]/event_PKiKP_radslo[ii])
    #    event_obs_bazi  = event_baz[ii] + (bazi_rad * 180 / np.pi)

    goto = '/Users/vidale/Documents/Research/IC/EvLocs'
    os.chdir(goto)

    file = open(eq_file1, 'r')
    lines=file.readlines()
    split_line = lines[0].split()
    t1           = UTCDateTime(split_line[1])
    date_label1  = split_line[1][0:10]

    file = open(eq_file2, 'r')
    lines=file.readlines()
    split_line = lines[0].split()
    t2           = UTCDateTime(split_line[1])
    date_label2  = split_line[1][0:10]

    #%% read files
    # Get saved event info, also used to name files
    # date_label = '2018-04-02' # date for filename
    goto = '/Users/vidale/Documents/Research/IC/Pro_files'
    os.chdir(goto)
    fname1 = 'HD' + date_label1 + '_2dstack.mseed'
    fname2 = 'HD' + date_label2 + '_2dstack.mseed'
    st1 = Stream()
    st2 = Stream()
    st1 = read(fname1)
    st2 = read(fname2)

    tshift    = st1.copy()  # make array for time shift
    amp_ratio = st1.copy()  # make array for relative amplitude
    amp_ave   = st1.copy()  # make array for average amplitude

    print('Read in: event 1 ' + str(len(st1)) + ' event 2 ' + str(len(st2)) + ' traces')
    nt1 = len(st1[0].data)
    nt2 = len(st2[0].data)
    dt1 = st1[0].stats.delta
    dt2 = st2[0].stats.delta
    print('Event 1 - First trace has ' + str(nt1) + ' time pts, time sampling of '
          + str(dt1) + ' and thus duration of ' + str((nt1-1)*dt1))
    print('Event 2 - First trace has ' + str(nt2) + ' time pts, time sampling of '
          + str(dt2) + ' and thus duration of ' + str((nt2-1)*dt2))
    if nt1 != nt2 or dt1 != dt2:
        print('nt or dt not does not match')
        exit(-1)

    #%% Make grid of slownesses
    slowR_n = int(1 + (slowR_hi - slowR_lo)/slow_delta)  # number of slownesses
    slowT_n = int(1 + (slowT_hi - slowT_lo)/slow_delta)  # number of slownesses
    print(str(slowT_n) + ' trans slownesses, hi and lo are ' + str(slowT_hi) + '  ' + str(slowT_lo))
    # In English, stack_slows = range(slow_n) * slow_delta - slow_lo
    a1R = range(slowR_n)
    a1T = range(slowT_n)
    stack_Rslows = [(x * slow_delta + slowR_lo) for x in a1R]
    stack_Tslows = [(x * slow_delta + slowT_lo) for x in a1T]
    print(str(slowR_n) + ' radial slownesses, ' + str(slowT_n) + ' trans slownesses, ')

#%%  Loop over slowness
    total_slows = slowR_n * slowT_n
    global_max = 0
    for slow_i in range(total_slows): # find envelope, phase, tshift, and global max
        if slow_i % 200 == 0:
            print('At line 101, ' +str(slow_i) + ' slowness out of ' + str(total_slows))
        if len(st1[slow_i].data) == 0: # test for zero-length traces
                print('%d data has zero length ' % (slow_i))

        seismogram1 = hilbert(st1[slow_i].data)  # make analytic seismograms
        seismogram2 = hilbert(st2[slow_i].data)

        env1 = np.abs(seismogram1) # amplitude
        env2 = np.abs(seismogram2)
        amp_ave[slow_i].data    = 0.5 * (env1 + env2)
        amp_ratio[slow_i].data  = env1/env2

        angle1 = np.angle(seismogram1) # time shift
        angle2 = np.angle(seismogram2)
        phase1 = np.unwrap(angle1)
        phase2 = np.unwrap(angle2)
        d_phase = (angle1 - angle2)
#        d_phase = phase1 - phase2
        for it in range(nt1):
            if d_phase[it] > math.pi:
                d_phase[it] -= 2 * math.pi
            elif d_phase[it] < -1 * math.pi:
                d_phase[it] += 2 * math.pi
            if d_phase[it] > math.pi or d_phase[it] < -math.pi:
                print(f'Bad d_phase value {d_phase[it]:.2f}  {it:4d}')
        # freq1 = np.diff(phase1) # observed freq in radians/sec
        # freq2 = np.diff(phase2)
        # ave_freq = 0.5*(freq1 + freq2)
        # ave_freq_plus = np.append(ave_freq,[1]) # ave_freq one element too short
        # tshift[slow_i].data     = d_phase / ave_freq_plus # 2*pi top and bottom cancels
        tshift[slow_i].data     = d_phase/(2*math.pi*freq_corr)

        local_max = max(abs(amp_ave[slow_i].data))
        if local_max > global_max:
            global_max = local_max
#%% Extract slices
    tshift_full = tshift.copy()  # make array for time shift
    for slow_i in range(total_slows): # ignore less robust points
        if slow_i % 200 == 0:
            print('At line 140, ' +str(slow_i) + ' slowness out of ' + str(total_slows))
        for it in range(nt1):
            if ((amp_ratio[slow_i].data[it] < (1/max_rat)) or (amp_ratio[slow_i].data[it] > max_rat) or (amp_ave[slow_i].data[it] < (min_amp * global_max))):
                tshift[slow_i].data[it] = np.nan
    #%% If desired, find transverse slowness nearest T_slow_plot
    lowest_Tslow = 1000000
    for slow_i in range(slowT_n):
        if abs(stack_Tslows[slow_i] - T_slow_plot) < lowest_Tslow:
            lowest_Tindex = slow_i
            lowest_Tslow = abs(stack_Tslows[slow_i] - T_slow_plot)

    print(str(slowT_n) + ' T slownesses, index ' + str(lowest_Tindex) + ' is closest to input parameter ' + str(T_slow_plot) + ', slowness diff there is ' + str(lowest_Tslow) + ' and slowness is ' + str(stack_Tslows[lowest_Tindex]))
    # Select only stacks with that slowness for radial plot
    centralR_st1 = Stream()
    centralR_st2 = Stream()
    centralR_amp   = Stream()
    centralR_ampr  = Stream()
    centralR_tdiff = Stream()
    for slowR_i in range(slowR_n):
        ii = slowR_i*slowT_n + lowest_Tindex
        centralR_st1 += st1[ii]
        centralR_st2 += st2[ii]
        centralR_amp   += amp_ave[ii]
        centralR_ampr  += amp_ratio[ii]
        centralR_tdiff += tshift[ii]

    #%% If desired, find radial slowness nearest R_slow_plot
    lowest_Rslow = 1000000
    for slow_i in range(slowR_n):
        if abs(stack_Rslows[slow_i] - R_slow_plot) < lowest_Rslow:
            lowest_Rindex = slow_i
            lowest_Rslow = abs(stack_Rslows[slow_i] - R_slow_plot)

    print(str(slowR_n) + ' R slownesses, index ' + str(lowest_Rindex) + ' is closest to input parameter ' + str(R_slow_plot) + ', slowness diff there is ' + str(lowest_Rslow) + ' and slowness is ' + str(stack_Rslows[lowest_Rindex]))

    # Select only stacks with that slowness for transverse plot
    centralT_st1 = Stream()
    centralT_st2 = Stream()
    centralT_amp   = Stream()
    centralT_ampr  = Stream()
    centralT_tdiff = Stream()

    #%% to extract stacked time functions
    event1_sample = Stream()
    event2_sample = Stream()

    for slowT_i in range(slowT_n):
        ii = lowest_Rindex*slowT_n + slowT_i
        centralT_st1 += st1[ii]
        centralT_st2 += st2[ii]
        centralT_amp   += amp_ave[ii]
        centralT_ampr  += amp_ratio[ii]
        centralT_tdiff += tshift[ii]

    #%% compute timing time series
    ttt = (np.arange(len(st1[0].data)) * st1[0].stats.delta + start_buff) # in units of seconds

#%% Plot radial amp and tdiff vs time plots
    fig_index = 6
#    plt.close(fig_index)
    plt.figure(fig_index,figsize=(30,10))
    plt.xlim(start_buff,end_buff)
    plt.ylim(stack_Rslows[0], stack_Rslows[-1])
    for slowR_i in range(slowR_n):  # loop over radial slownesses
        dist_offset = stack_Rslows[slowR_i] # trying for approx degrees
        ttt = (np.arange(len(centralR_st1[slowR_i].data)) * centralR_st1[slowR_i].stats.delta
         + (centralR_st1[slowR_i].stats.starttime - t1))
        plt.plot(ttt, (centralR_st1[slowR_i].data - np.median(centralR_st1[slowR_i].data))*plot_scale_fac /global_max + dist_offset, color = 'green')
        plt.plot(ttt, (centralR_st2[slowR_i].data - np.median(centralR_st2[slowR_i].data))*plot_scale_fac /global_max + dist_offset, color = 'red')
        # extract stacked time functions
        if get_stf != 0:
            if np.abs(stack_Rslows[slowR_i]- 0.005) < 0.000001: # kludge, not exactly zero when desired
                event1_sample = centralR_st1[slowR_i].copy()
                event2_sample = centralR_st2[slowR_i].copy()
#        plt.plot(ttt, (centralR_amp[slowR_i].data)  *plot_scale_fac/global_max + dist_offset, color = 'purple')
        if turn_off_black == 0:
            plt.plot(ttt, (centralR_tdiff[slowR_i].data)*plot_scale_fac/1 + dist_offset, color = 'black')
            plt.plot(ttt, (centralR_amp[slowR_i].data)*0.0 + dist_offset, color = 'lightgray') # reference lines
    plt.xlabel('Time (s)')
    plt.ylabel('R Slowness (s/km)')
    plt.title(dphase + ' seismograms and tdiff at ' + str(T_slow_plot) + ' T slowness, green is event1, red is event2')
    # Plot transverse amp and tdiff vs time plots
    fig_index = 7
#    plt.close(fig_index)
    plt.figure(fig_index,figsize=(30,10))
    plt.xlim(start_buff,end_buff)
    plt.ylim(stack_Tslows[0], stack_Tslows[-1])

    for slowT_i in range(slowT_n):  # loop over transverse slownesses
        dist_offset = stack_Tslows[slowT_i] # trying for approx degrees
        ttt = (np.arange(len(centralT_st1[slowT_i].data)) * centralT_st1[slowT_i].stats.delta
         + (centralT_st1[slowT_i].stats.starttime - t1))
        plt.plot(ttt, (centralT_st1[slowT_i].data - np.median(centralT_st1[slowT_i].data))*plot_scale_fac /global_max + dist_offset, color = 'green')
        plt.plot(ttt, (centralT_st2[slowT_i].data - np.median(centralT_st2[slowT_i].data))*plot_scale_fac /global_max + dist_offset, color = 'red')
#        plt.plot(ttt, (centralT_amp[slowT_i].data)  *plot_scale_fac/global_max + dist_offset, color = 'purple')
        if turn_off_black == 0:
            plt.plot(ttt, (centralT_tdiff[slowT_i].data)*plot_scale_fac/1 + dist_offset, color = 'black')
            plt.plot(ttt, (centralT_amp[slowT_i].data)*0.0 + dist_offset, color = 'lightgray') # reference lines
    plt.xlabel('Time (s)')
    plt.ylabel('T Slowness (s/km)')
    plt.title(date_label1 + '  ' + dphase + ' seismograms and tdiff ' + str(R_slow_plot) + ' R slowness, green is event1, red is event2')
    os.chdir('/Users/vidale/Documents/Research/IC/Plots')
#    plt.savefig(date_label1 + '_' + str(start_buff) + '_' + str(end_buff) + '_stack.png')

#%% R-T tshift averaged over time window
    fig_index = 8
    stack_slice = np.zeros((slowR_n,slowT_n))

    if start_beam == 0 and end_beam == 0:
        full_beam = 1
    else:  # beam just part of stack volume
        full_beam = 0
        start_index = int((start_beam - start_buff) / dt1)
        end_index   = int((end_beam   - start_buff) / dt1)
        print('beam is ' + str(start_beam) + ' to ' + str(end_beam) + 's, out of ' + str(start_buff)
            + ' to ' + str(end_buff) + 's, dt is ' + str(dt1)  + 's, and indices are '+ str(start_index) + ' ' + str(end_index))

    for slowR_i in range(slowR_n):  # loop over radial slownesses
        for slowT_i in range(slowT_n):  # loop over transverse slownesses
            index = slowR_i*slowT_n + slowT_i
            if full_beam == 1:
                num_val = np.nanmedian(tshift[index].data)
            else:
                num_val = np.nanmedian(tshift[index].data[start_index:end_index])
#            num_val = statistics.median(tshift_full[index].data)
            stack_slice[slowR_i, slowT_i] = num_val # adjust for dominant frequency of 1.2 Hz, not 1 Hz
#    stack_slice[0,0] = -0.25
#    stack_slice[0,1] =  0.25
#    tdiff_clip = 0.4/1.2
    tdiff_clip_max =  tdiff_clip  # DO NOT LEAVE COMMENTED OUT!!
    tdiff_clip_min = -tdiff_clip

    y1, x1 = np.mgrid[slice(stack_Rslows[0], stack_Rslows[-1] + slow_delta, slow_delta),
                 slice(stack_Tslows[0], stack_Tslows[-1] + slow_delta, slow_delta)]

    fig, ax = plt.subplots(1, figsize=(7,0.8*7*(slowR_n/slowT_n)))  # try to make correct aspect ratio plot
#        fig, ax = plt.subplots(1, figsize=(9,2))
#        fig.subplots_adjust(bottom=0.3)
#    c = ax.pcolormesh(x1, y1, stack_slice, cmap=plt.cm.bwr,      vmin = tdiff_clip_min, vmax = tdiff_clip_max)
    c = ax.pcolormesh(x1, y1, stack_slice, cmap=plt.cm.coolwarm, vmin = tdiff_clip_min, vmax = tdiff_clip_max)
    ax.axis([x1.min(), x1.max(), y1.min(), y1.max()])
    circle1 = plt.Circle((0, 0), 0.019, color='black', fill=False)
    ax.add_artist(circle1)
    circle2 = plt.Circle((0, 0), 0.040, color='black', fill=False)
    ax.add_artist(circle2)  #outer core limit
    fig.colorbar(c, ax=ax)
    plt.ylabel('R Slowness (s/km)')
    plt.xlabel('Transverse Slowness (s/km)')
    plt.title(dphase + ' time shift ' + date_label1 + ' ' + date_label2)
    os.chdir('/Users/vidale/Documents/Research/IC/Plots')
    plt.savefig(date_label1 + '_' + date_label2 + '_' + str(start_buff) + '_' + str(end_buff) + '_tshift.png')
    plt.show()

#%% R-T amplitude averaged over time window
    fig_index = 9
    stack_slice = np.zeros((slowR_n,slowT_n))
    smax = 0
    for slowR_i in range(slowR_n):  # loop over radial slownesses
        for slowT_i in range(slowT_n):  # loop over transverse slownesses
            index = slowR_i*slowT_n + slowT_i
            if full_beam == 1:
                num_val = np.nanmedian(amp_ave[index].data)
            else:
                num_val = np.nanmedian(amp_ave[index].data[start_index:end_index])
            stack_slice[slowR_i, slowT_i] = num_val
            if num_val > smax:
                smax = num_val
#    stack_slice[0,0] = 0

    y1, x1 = np.mgrid[slice(stack_Rslows[0], stack_Rslows[-1] + slow_delta, slow_delta),
                 slice(stack_Tslows[0], stack_Tslows[-1] + slow_delta, slow_delta)]

#    fig, ax = plt.subplots(1)
    fig, ax = plt.subplots(1, figsize=(7,0.8*7*(slowR_n/slowT_n)))
#    c = ax.pcolormesh(x1, y1, stack_slice/smax, cmap=plt.cm.gist_yarg, vmin = 0.5)
    c = ax.pcolormesh(x1, y1, stack_slice/smax, cmap=plt.cm.gist_rainbow_r, vmin = 0)
#    c = ax.pcolormesh(x1, y1, stack_slice, cmap=plt.cm.gist_rainbow_r, vmin = 0)
    fig.colorbar(c, ax=ax, label='linear amplitude')
    ax.axis([x1.min(), x1.max(), y1.min(), y1.max()])
    circle1 = plt.Circle((0, 0), 0.019, color='black', fill=False)
    ax.add_artist(circle1)  #inner core limit
    circle2 = plt.Circle((0, 0), 0.040, color='black', fill=False)
    ax.add_artist(circle2)  #outer core limit

    c = ax.scatter(pred_Eslo, pred_Nslo, color='black'  , s=50, alpha=0.75)
#    c = ax.scatter( obs_Eslo,  obs_Nslo, color='purple', s=100, alpha=0.75)
    c = ax.scatter(        0,         0, color='black' , s=50,  alpha=0.75)

    plt.xlabel('Transverse Slowness (s/km)')
    plt.ylabel('Radial Slowness (s/km)')
    plt.title(date_label1 + ' ' + date_label2 + '  ' + dphase + ' beam amplitude')
    os.chdir('/Users/vidale/Documents/Research/IC/Plots')
    plt.savefig(date_label1 + '_' + date_label2 + '_' + str(start_buff) + '_' + str(end_buff) + '_beam.png')
    plt.show()

#%%  Save processed files
    goto = '/Users/vidale/Documents/Research/IC/Pro_Files'
    os.chdir(goto)

    if decimate_fac != 0:
        for i in range(len(tshift_full)):  # loop over traces
            tshift_full[i].decimate(decimate_fac, no_filter=True)
            amp_ave[i].decimate(decimate_fac, no_filter=True)
            amp_ratio[i].decimate(decimate_fac, no_filter=True)

    fname = 'HD' + date_label1 + '_' + date_label2 + '_tshift.mseed'
    tshift_full.write(fname,format = 'MSEED')

    fname = 'HD' + date_label1 + '_' + date_label2 + '_amp_ave.mseed'
    amp_ave.write(fname,format = 'MSEED')

    fname = 'HD' + date_label1 + '_' + date_label2 + '_amp_ratio.mseed'
    amp_ratio.write(fname,format = 'MSEED')

#%% Option to write out stf
    if get_stf != 0:
        event1_sample.taper(0.1)
        event2_sample.taper(0.1)
        fname = 'HD' + date_label1 + '_stf.mseed'
        event1_sample.write(fname,format = 'MSEED')
        fname = 'HD' + date_label2 + '_stf.mseed'
        event2_sample.write(fname,format = 'MSEED')

    elapsed_time_wc = time.time() - start_time_wc
    print(f'This job took   {elapsed_time_wc:.1f}   seconds')
    os.system('say "Done"')
