#!/usr/bin/env python
# looks at envelope stack differential parameters from pro6
# Reads in tdiff, ave_amp, cc computed from a pair of events
# window by signal quality
# eleven_slice == T produces sections at 4 radial slownesses (0.0, 0.05, 0.01, 0.015)
# plus plots sections at 7 transv slownesses (-0.015 to 0.015)
# eleven_slice == F plots one radial and one transverse plot, plus snaps
# John Vidale 2/2019, overhauled 1/2021

def pro7_pair_scan(eq_file1, eq_file2, slow_delta = 0.0005, turn_off_black = 0,
              slowR_lo = -0.1, slowR_hi = 0.1, slowT_lo = -0.1, slowT_hi = 0.1,
              start_buff = 50, end_buff = 50,fig_index = 401, do_T = False, do_R = False,
              tdiff_clip = 1, cc_thres = 0.8, min_amp = 0.2, plot_scale_fac = 1,
              R_slow_plot = 0.06, T_slow_plot = 0.0, snaptime = 8, snaps = 10,
              nR_plots  = 3, nT_plots = 3, slow_incr = 0.01, NS = False, dphase = 'PKiKP',
              log_plot = False,log_plot_range = 2):

    from obspy import read
    from obspy.taup import TauPyModel
    import numpy as np
    import os
    import matplotlib.pyplot as plt
    import time
    import sys
    import math
    from obspy import Stream
    from termcolor import colored
    from obspy.geodetics import gps2dist_azimuth
    model = TauPyModel(model='iasp91')

    print(colored('Running pro7b_plot_stack', 'cyan'))
    start_time_wc = time.time()

    if NS == True:
        print('NS == True, these are not radial and transverse coordinates')
        sys.exit()

    #%% Input parameters and computed files
    folder_name = '/Users/vidale/Documents/Research/IC/'
    file1 = open(folder_name + 'EvLocs/' + eq_file1, 'r')
    file2 = open(folder_name + 'EvLocs/' + eq_file2, 'r')
    lines1=file1.readlines()
    lines2=file2.readlines()
    split_line1 = lines1[0].split()
    split_line2 = lines2[0].split()
    # t2 = UTCDateTime(split_line2[1])
    date_label1  = split_line1[1][0:10]
    date_label2  = split_line2[1][0:10]
    ev_lat      = float(      split_line1[2])
    ev_lon      = float(      split_line1[3])
    ev_depth    = float(      split_line1[4])
    # date_label = '2018-04-02' # dates in filename

    ref_lat = 46.7  # °N keep only inner rings A-D
    ref_lon = -106.22   # °E
    ref_distance = gps2dist_azimuth(ref_lat,ref_lon,ev_lat,ev_lon)
    ref_back_az = ref_distance[2]

    ref1_dist  = ref_distance[0]/(1000*111)

    arrivals1 = model.get_travel_times(source_depth_in_km=ev_depth,distance_in_degree=ref1_dist-0.5,phase_list=[dphase])
    arrivals2 = model.get_travel_times(source_depth_in_km=ev_depth,distance_in_degree=ref1_dist+0.5,phase_list=[dphase])
    dtime = arrivals2[0].time - arrivals1[0].time
    event_pred_slo  = dtime/111.  # s/km
    # convert to pred rslo and tslo
    if NS == True:
        sin_baz = np.sin(ref_back_az * np.pi /180)
        cos_baz = np.cos(ref_back_az * np.pi /180)
    #  rotate predicted slowness to N and E
        pred_Nslo = event_pred_slo * cos_baz
        pred_Eslo = event_pred_slo * sin_baz
    else:
        pred_Nslo = event_pred_slo
        pred_Eslo = 0

    name_str = folder_name + 'Pro_files/HD' + date_label1 + '_' + date_label2 + '_'
    fname1  = name_str + 'tshift.mseed'
    fname2  = name_str + 'amp_ave.mseed'
    fname3  = name_str + 'cc.mseed'
    tdiff   = Stream()
    amp_ave = Stream()
    cc      = Stream()
    tdiff   = read(fname1)
    amp_ave = read(fname2)
    cc      = read(fname3)

    nt = len(tdiff[0].data)

    #%% Make grid of slownesses
    slowR_n = int(round(1 + (slowR_hi - slowR_lo)/slow_delta))  # number of slownesses
    slowT_n = int(round(1 + (slowT_hi - slowT_lo)/slow_delta))  # number of slownesses
    print(str(slowT_n) + ' trans slownesses, hi and lo are ' + str(slowT_hi) + '  ' + str(slowT_lo))
    print(str(slowR_n) + ' radial slownesses, hi and lo are ' + str(slowR_hi) + '  ' + str(slowR_lo))
    # In English, stack_slows = range(slow_n) * slow_delta - slow_lo
    a1R = range(slowR_n)
    a1T = range(slowT_n)
    stack_Rslows = [(x * slow_delta + slowR_lo) for x in a1R]
    stack_Tslows = [(x * slow_delta + slowT_lo) for x in a1T]

    #%% Select subset if Zoomed

    #%% Mask out weak and/or less correlated points

    global_max = 0
    for slow_i in range(len(amp_ave)): # find global max of ave_amp
        for data_i in range(len(amp_ave[slow_i].data)): # find global max of ave_amp
            if log_plot == True:
                amp_ave[slow_i].data[data_i] = math.log10(amp_ave[slow_i].data[data_i])
        local_max = max(amp_ave[slow_i].data)
        if local_max > global_max:
            global_max = local_max

    amp_ave_thres = amp_ave.copy()  # copy amp array so that thresholding can normalize properly
    for slow_i in range(len(tdiff)): # ignore less robust points
        for it in range(nt):
            if (cc[slow_i].data[it] < cc_thres) or (amp_ave[slow_i].data[it] < (min_amp * global_max)):
                tdiff[        slow_i].data[it] = np.nan
                amp_ave_thres[slow_i].data[it] = np.nan

    #%% Beam sum plots
    #%% -- R-T tdiff amp-normed
    stack_slice = np.zeros((slowR_n,slowT_n))

    for slowR_i in range(slowR_n):  # loop over radial slownesses
        for slowT_i in range(slowT_n):  # loop over transverse slownesses
            index = slowR_i*slowT_n + slowT_i
            num_val = np.nansum(np.multiply(tdiff[index].data, amp_ave_thres[index].data))/np.nansum(amp_ave_thres[index].data)
            stack_slice[slowR_i, slowT_i] = num_val

    y1, x1 = np.mgrid[slice(stack_Rslows[0], stack_Rslows[-1] + slow_delta, slow_delta),
                 slice(stack_Tslows[0], stack_Tslows[-1] + slow_delta, slow_delta)]

    fig, ax = plt.subplots(1, figsize=(7,0.8*7*(slowR_n/slowT_n)))  # try to make correct aspect ratio plot
    c = ax.pcolormesh(x1, y1, stack_slice, cmap=plt.cm.coolwarm, vmin = -tdiff_clip, vmax = tdiff_clip)
    fig.colorbar(c, ax=ax, label='time lag (s)')
    ax.axis([x1.min(), x1.max(), y1.min(), y1.max()])
    circle1 = plt.Circle((0, 0), 0.019, color='black', fill=False)
    ax.add_artist(circle1)
    circle2 = plt.Circle((0, 0), 0.040, color='black', fill=False)
    ax.add_artist(circle2)  #outer core limit

    c = ax.scatter(pred_Eslo, pred_Nslo, color='black'  , s=50, alpha=0.75)
    c = ax.scatter(        0,         0, color='black' , s=50,  alpha=0.75)

    plt.ylabel('R Slowness (s/km)')
    plt.xlabel('Transverse Slowness (s/km)')
    plt.title(dphase + ' time shift ' + date_label1 + ' ' + date_label2 + ' amp weighted')
    os.chdir('/Users/vidale/Documents/Research/IC/Plots')
    plt.savefig(date_label1 + '_' + date_label2 + '_' + str(start_buff) + '_' + str(end_buff) + '_tdiff.png')
    plt.show()

#%% -- R-T amplitude averaged over time window
    stack_slice = np.zeros((slowR_n,slowT_n))
    smax = 0
    for slowR_i in range(slowR_n):  # loop over radial slownesses
        for slowT_i in range(slowT_n):  # loop over transverse slownesses
            index = slowR_i*slowT_n + slowT_i
            num_val = np.nanmean(amp_ave[index].data)
            stack_slice[slowR_i, slowT_i] = num_val
            if num_val > smax:
                smax = num_val

    y1, x1 = np.mgrid[slice(stack_Rslows[0], stack_Rslows[-1] + slow_delta, slow_delta),
                 slice(stack_Tslows[0], stack_Tslows[-1] + slow_delta, slow_delta)]

    fig, ax = plt.subplots(1, figsize=(7,0.8*7*(slowR_n/slowT_n)))
    c = ax.pcolormesh(x1, y1, stack_slice/smax, cmap=plt.cm.gist_rainbow_r, vmin = 0)
    if log_plot == True:
        fig.colorbar(c, ax=ax, label='log amplitude')
    else:
        fig.colorbar(c, ax=ax, label='linear amplitude')
    ax.axis([x1.min(), x1.max(), y1.min(), y1.max()])
    circle1 = plt.Circle((0, 0), 0.019, color='black', fill=False)
    ax.add_artist(circle1)  #inner core limit
    circle2 = plt.Circle((0, 0), 0.040, color='black', fill=False)
    ax.add_artist(circle2)  #outer core limit

    c = ax.scatter(pred_Eslo, pred_Nslo, color='black' , s=50, alpha=0.75)
    c = ax.scatter(        0,         0, color='black' , s=50, alpha=0.75)

    plt.xlabel('Transverse Slowness (s/km)')
    plt.ylabel('Radial Slowness (s/km)')
    plt.title(date_label1 + ' ' + date_label2 + '  ' + dphase + ' beam amplitude')
    plt.show()

    os.system('say "Done"')