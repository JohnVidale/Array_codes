#!/usr/bin/env python
# looks at envelope stack differential parameters from pro6
# Reads in tdiff, ave_amp, cc computed from a pair of events
# window by signal quality
# eleven_slice == T produces sections at 4 radial slownesses (0.0, 0.05, 0.01, 0.015)
# plus plots sections at 7 transv slownesses (-0.015 to 0.015)
# eleven_slice == F plots one radial and one transverse plot, plus snaps
# John Vidale 2/2019, overhauled 1/2021

def pro7_pair_scan(eq_num1, eq_num2, slow_delta = 0.0005, turn_off_black = 1,
              slowR_lo = -0.1, slowR_hi = 0.1, slowT_lo = -0.1, slowT_hi = 0.1,
              start_buff = 50, end_buff = 50,fig_index = 401, do_T = False, do_R = False,
              ZslowR_lo = -0.1, ZslowR_hi = 0.1, ZslowT_lo = -0.1, ZslowT_hi = 0.1,
              Zstart_buff = 50, Zend_buff = 50, zoom = False, tdiff_clip = 1,
              ref_phase = 'blank', cc_thres = 0.8, min_amp = 0.2,
              R_slow_plot = 0.06, T_slow_plot = 0.0, no_plots = False,
              snaptime = 8, snaps = 10, snap_depth = 0,
              nR_plots  = 3, nT_plots = 3, slow_incr = 0.01, NS = False,
              ARRAY = 0, auto_slice = True, two_slice_plots = False, beam_sums = True,
              wiggly_plots = 0, start_beam = 0, end_beam = 0, log_plot = False,
              log_plot_range = 2, no_tdiff_plot = False,
              wig_scale_fac = 1, tdiff_scale_fac = 1):

    from obspy import read
    from obspy.taup import TauPyModel
    import numpy as np
    import os
    import matplotlib.pyplot as plt
    import time
    import sys
    import math
    from obspy import UTCDateTime
    from obspy import Stream
    from termcolor import colored
    from obspy.geodetics import gps2dist_azimuth
    model = TauPyModel(model='iasp91')

    print(colored('Running pro7b_pair_scan', 'cyan'))
    start_time_wc = time.time()

    if zoom == True:
        if Zstart_buff  < start_buff:
            print(f'Zstart_buff of {Zstart_buff:.1f} cannot be < start_buff of {start_buff:.1f}')
            Zstart_buff = start_buff
        if Zend_buff    > end_buff:
            print(f'Zend_buff of {Zend_buff:.1f} cannot be > end_buff of {end_buff:.1f}')
            Zend_buff   = end_buff

    #%% Input parameters and computed files
    fname1 = '/Users/vidale/Documents/Research/IC/EvLocs/event' + str(eq_num1) + '.txt'
    fname2 = '/Users/vidale/Documents/Research/IC/EvLocs/event' + str(eq_num2) + '.txt'
    file1 = open(fname1, 'r')
    file2 = open(fname2, 'r')
    lines1=file1.readlines()
    lines2=file2.readlines()
    split_line1 = lines1[0].split()
    split_line2 = lines2[0].split()
    t1          = UTCDateTime(split_line1[1])
    # t2 = UTCDateTime(split_line2[1])
    date_label1  = split_line1[1][0:10]
    date_label2  = split_line2[1][0:10]
    save_name = '/Users/vidale/Documents/Research/IC/Plots/' + date_label1 + '_' + date_label2 + '_'
    ev_lat      = float(      split_line1[2])
    ev_lon      = float(      split_line1[3])
    ev_depth    = float(      split_line1[4])
    # date_label = '2018-04-02' # dates in filename

    if ARRAY == 0:
        ref_lat = 36.3  # °N, around middle of Japan
        ref_lon = 138.5 # °E
    elif ARRAY == 1:
        ref_lat = 46.7  # °N keep only inner rings A-D
        ref_lon = -106.22   # °E
    elif ARRAY == 2: # China set and center
        ref_lat = 38      # °N
        ref_lon = 104.5   # °E
    ref_distance = gps2dist_azimuth(ref_lat,ref_lon,ev_lat,ev_lon)
    ref_dist     = ref_distance[0]/(1000*111)
    ref_az       = ref_distance[1]
    ref_back_az  = ref_distance[2]

    arrivals_ref   = model.get_travel_times(source_depth_in_km=ev_depth,distance_in_degree=ref_dist, phase_list=[ref_phase])
    arrival_time = arrivals_ref[0].time

    arrivals1 = model.get_travel_times(source_depth_in_km=ev_depth,distance_in_degree=ref_dist-0.5,phase_list=[ref_phase])
    arrivals2 = model.get_travel_times(source_depth_in_km=ev_depth,distance_in_degree=ref_dist+0.5,phase_list=[ref_phase])
    dtime = arrivals2[0].time - arrivals1[0].time
    event_pred_slo  = dtime/111.  # s/km
    # convert to pred rslo and tslo
    if NS == True:    #  rotate predicted slowness to N and E
        print(f'Array  lat {ref_lat:.0f}, lon  {ref_lon:.0f}, Event lat {ev_lat:.0f}, lon {ev_lon:.0f}, az {ref_az:.0f}, baz {ref_back_az:.0f}')
        sin_baz = np.sin(ref_az * np.pi /180)
        cos_baz = np.cos(ref_az * np.pi /180)
        pred_Nslo = event_pred_slo * cos_baz
        pred_Eslo = event_pred_slo * sin_baz
    else:
        pred_Nslo = event_pred_slo
        pred_Eslo = 0

    name_str = '/Users/vidale/Documents/Research/IC/Pro_files/HD' + date_label1 + '_' + date_label2 + '_'
    fname1  = name_str + 'tshift.mseed'
    fname2  = name_str + 'amp_ave.mseed'
    fname3  = name_str + 'cc.mseed'
    tdiff   = Stream()
    amp_ave = Stream()
    cc      = Stream()
    tdiff   = read(fname1)
    amp_ave = read(fname2)
    cc      = read(fname3)

    dt     = tdiff[  0].stats.delta
    dt_cc  = cc[     0].stats.delta
    dt_amp = amp_ave[0].stats.delta
    nt     = len(tdiff[  0].data)
    nt_cc  = len(cc[     0].data)
    nt_amp = len(amp_ave[0].data)
    print(f'cc      data length is {nt_cc } time pts, dt is {dt_cc :.2f}, so record length is {dt_cc  * nt_cc:.0f} seconds')
    print(f'tdiff   data length is {nt    } time pts, dt is {dt    :.2f}, so record length is {dt     * nt   :.0f} seconds')
    print(f'amp_ave data length is {nt_amp} time pts, dt is {dt_amp:.2f}, so record length is {dt_amp * nt   :.0f} seconds')
    print(f'slowness grid has {len(tdiff)} elements')

    #%% Make grid of slownesses
    slowR_n = int(round(1 + (slowR_hi - slowR_lo)/slow_delta))  # number of slownesses
    slowT_n = int(round(1 + (slowT_hi - slowT_lo)/slow_delta))  # number of slownesses
    stack_nt = int(round(1 + (end_buff - start_buff)/dt))  # number of time points
    print(f'{slowT_n} trans slownesses, hi and lo are {slowT_hi} and {slowT_lo}')
    print(f'{slowR_n} trans slownesses, hi and lo are {slowR_hi} and {slowR_lo}')
    # In English, stack_slows = range(slow_n) * slow_delta - slow_lo
    a1R = range(slowR_n)
    a1T = range(slowT_n)
    stack_Rslows = [(x * slow_delta + slowR_lo) for x in a1R]
    stack_Tslows = [(x * slow_delta + slowT_lo) for x in a1T]

    #%% Select subset if Zoomed
    if zoom == True:
        Ztdiff   = Stream()
        Zamp_ave = Stream()
        Zcc      = Stream()
        print(f'before calculation, tdiff[0] has length {len(tdiff[0])})')
        for slowR_i in range(slowR_n):  # loop over radial slownesses
            for slowT_i in range(slowT_n):  # loop over transverse slownesses, kludge to evade rounding error
                if ((stack_Rslows[slowR_i] >= ZslowR_lo - 0.000001) and (stack_Rslows[slowR_i] <= ZslowR_hi + 0.000001) and
                    (stack_Tslows[slowT_i] >= ZslowT_lo - 0.000001) and (stack_Tslows[slowT_i] <= ZslowT_hi + 0.000001)):
                    index = slowR_i*slowT_n + slowT_i
                    s_t = t1 + Zstart_buff
                    e_t = t1 + Zend_buff
                    Ztdiff   += tdiff[  index].trim(starttime=s_t, endtime=e_t)
                    Zamp_ave += amp_ave[index].trim(starttime=s_t, endtime=e_t)
                    Zcc      += cc[     index].trim(starttime=s_t, endtime=e_t)
                            #tr.trim(starttime=s_t,endtime = e_t)
        tdiff   = Ztdiff
        amp_ave = Zamp_ave
        cc      = Zcc
        nt = len(tdiff[0].data)
        start_buff = Zstart_buff
        # make time series
        print(f'after calculation, Ztdiff[0] has length {len(Ztdiff[0])}')
        print(f'after calculation, tdiff[0] has length {len(tdiff[0])}')
        print(f'slowR_lo  is {slowR_lo}  and slowR_hi  is {slowR_hi}  and slowT_lo  is {slowT_lo}  and slowT_hi is {slowT_hi}')
        print(f'ZslowR_lo is {ZslowR_lo} and ZslowR_hi is {ZslowR_hi} and ZslowT_lo is {ZslowT_lo} and ZslowT_hi is {ZslowT_hi}')

        #%% -- Re-make finer grid of slownesses
        slowR_lo   = ZslowR_lo
        slowR_hi   = ZslowR_hi
        slowT_lo   = ZslowT_lo
        slowT_hi   = ZslowT_hi
        end_buff   = Zend_buff
        start_buff = Zstart_buff
        end_buff   = Zend_buff
        slowR_n = int(round(1 + (slowR_hi - slowR_lo)/slow_delta))  # number of slownesses
        slowT_n = int(round(1 + (slowT_hi - slowT_lo)/slow_delta))  # number of slownesses
        stack_nt = int(round(1 + ((end_buff - start_buff)/dt)))  # number of time points
        print('After zoom ' + str(slowT_n) + ' trans slownesses, hi and lo are ' + str(slowT_hi) + '  ' + str(slowT_lo) + ' stack_nt is ' + str(stack_nt))
        # In English, stack_slows = range(slow_n) * slow_delta - slow_lo
        a1R = range(slowR_n)
        a1T = range(slowT_n)
        stack_Rslows = [(x * slow_delta + slowR_lo) for x in a1R]
        stack_Tslows = [(x * slow_delta + slowT_lo) for x in a1T]
        print('After zoom ' + str(slowR_n) + ' radial slownesses, ' + str(slowT_n) + ' trans slownesses, ')
        print('Output trace starttime ' + str(Ztdiff[0].stats.starttime))

    ttt = (np.arange(len(tdiff[0].data)) * tdiff[0].stats.delta + start_buff) # in units of seconds

    global_max = 0  # find global_max, largest amplitude in amp_ave beam array envelopes
    for slow_i in range(len(amp_ave)): # find global max of ave_amp
        local_max = max(amp_ave[slow_i].data)
        if local_max > global_max:
            global_max = local_max

    #%% Mask out weak and/or less correlated points
    amp_ave_thres = amp_ave.copy()  # copy amp envelope array, set amps and tdiff below thresholds to NaN using global_max
    nt = len(tdiff[0].data)
    for slow_i in range(len(tdiff)): # don't plot less robust points
        for it in range(nt):
            if (cc[slow_i].data[it] < cc_thres) or (amp_ave[slow_i].data[it] < (min_amp * global_max)):
                tdiff[        slow_i].data[it] = np.nan
                amp_ave_thres[slow_i].data[it] = np.nan

    for slow_i in range(len(tdiff)): # set NaN to avoid plotting phase wrapping when tdiff jumps too much
        for it in range(nt-1):
            if (abs(tdiff[slow_i].data[it+1] - tdiff[slow_i].data[it]) > tdiff_clip):
                tdiff[slow_i].data[it] = np.nan

    if log_plot == True:  # convert amp envelope array to log amp and record global_max of logs
        global_max = -100  # different global max if converting to plotting log amp
                           # remember logs can be negative
        for slow_i in range(len(amp_ave)): # find global max of ave_amp
            for data_i in range(len(amp_ave[slow_i].data)): # find global max of ave_amp
                amp_ave[slow_i].data[data_i] = math.log10(amp_ave[slow_i].data[data_i])
            local_max = max(amp_ave[slow_i].data)
            if local_max > global_max:
                global_max = local_max

#%% Auto slice option
    if auto_slice == True:

#%% -- compute timing time series
        #%% -- R slices
        ttt = (np.arange(stack_nt) * dt + start_buff)
        if do_R == True:  # remember plots scanning R are those at constant T
            for T_cnt in range(-nR_plots, nR_plots + 1):
                if nR_plots * slow_incr > slowT_hi:
                    print('nR_plots * slow_incr > slowT_hi, out of range')
                    sys.exit()
                if -(nR_plots * slow_incr) < slowT_lo:
                    print('-nR_plots * slow_incr < slowT_lo, out of range')
                    sys.exit()
                #%% -- -- gather R data
                lowest_Tslow = 1000000  # find index of row with T_cnt slowness, awkward coding
                target_slow = (T_cnt * slow_incr)
                for slow_i in range(slowT_n):
                    if abs(stack_Tslows[slow_i] - target_slow) < lowest_Tslow:
                        lowest_Tindex = slow_i
                        lowest_Tslow = abs(stack_Tslows[slow_i] - target_slow)
                print(f'For R plot {T_cnt:2d}, {lowest_Tindex:3d} is T slow nearest {target_slow:.3f}, difference is {lowest_Tslow:.3f}')

                # Collect data with that slowness for R (T=const) plot
                Rcentral_st = Stream()
                Rcentral_am = Stream()
                for slowR_i in range(slowR_n):
                    Rcentral_st += tdiff[  slowR_i*slowT_n + lowest_Tindex]
                    Rcentral_am += amp_ave[slowR_i*slowT_n + lowest_Tindex]


                #%% -- -- plot R tdiff
                if no_tdiff_plot == False:
                    stack_arrayR_Tdf = np.zeros((slowR_n,stack_nt))
                    for it in range(stack_nt):  # check points one at a time
                        for slowR_i in range(slowR_n):  # for this station, loop over slownesses
                            num_val = Rcentral_st[slowR_i].data[it]
                            stack_arrayR_Tdf[slowR_i, it] = num_val

                    y, x = np.mgrid[slice(stack_Rslows[0], stack_Rslows[-1] + slow_delta, slow_delta),
                                 slice(ttt[0], ttt[-1] + dt, dt)]

                    fig, ax = plt.subplots(1, figsize=(10,3))
                    c = ax.pcolormesh(x, y, stack_arrayR_Tdf, cmap=plt.cm.coolwarm, vmin= -tdiff_clip, vmax=tdiff_clip)
                    fig.subplots_adjust(bottom=0.2)
                    ax.axis([x.min(), x.max(), y.min(), y.max()])
                    fig.colorbar(c, ax=ax, label='time shift (s)')
                    c = ax.scatter(arrival_time, event_pred_slo, color='black'  , s=50, alpha=0.75)
                    plt.xlabel('Time (s)')
                    plt.ylabel('Radial slowness (s/km)')
                    plt.title(f'{ref_phase} Tdiff at {target_slow:.3f} s/km T slowness, {fname1[48:58]}  {fname1[59:69]}  min amp {min_amp:.1f}  cc_thres {cc_thres:.2f}')
                    if no_plots == False:
                        plt.show()
                    fig_index += 1

                #%% -- -- plot R amp
                stack_arrayR_Amp = np.zeros((slowR_n,stack_nt))
                for it in range(stack_nt):  # check points one at a time
                    for slowR_i in range(slowR_n):  # for this station, loop over slownesses
                        num_val = Rcentral_am[slowR_i].data[it]
                        stack_arrayR_Amp[slowR_i, it] = num_val

                y, x = np.mgrid[slice(stack_Rslows[0], stack_Rslows[-1] + slow_delta, slow_delta),
                             slice(ttt[0], ttt[-1] + dt, dt)]

                fig, ax = plt.subplots(1, figsize=(10,3))
                if log_plot == True:
                    c = ax.pcolormesh(x, y, stack_arrayR_Amp - global_max, cmap=plt.cm.gist_rainbow_r, vmin= - log_plot_range, vmax=0)
                else:
                    c = ax.pcolormesh(x, y, stack_arrayR_Amp, cmap=plt.cm.gist_rainbow_r, vmin= 0, vmax=global_max)
                fig.subplots_adjust(bottom=0.2)
                ax.axis([x.min(), x.max(), y.min(), y.max()])
                if log_plot == True:
                    fig.colorbar(c, ax=ax, label='log amplitude')
                else:
                    fig.colorbar(c, ax=ax, label='linear amplitude')
                c = ax.scatter(arrival_time, event_pred_slo, color='black'  , s=50, alpha=0.75)
                plt.xlabel('Time (s)')
                if NS == True:
                    plt.ylabel('North Slowness (s/km)')
                    plt.title(f'Amp at {target_slow:.3f} s/km E slowness, {fname1[48:58]}  {fname1[59:69]}')
                else:
                    plt.ylabel('Radial Slowness (s/km)')
                    plt.title(f'Amp at {target_slow:.3f} s/km T slowness, {fname1[48:58]}  {fname1[59:69]}')
                if no_plots == False:
                    plt.show()
                fig_index += 1

        #%% -- T slices
        if do_T == True:  # remember plots scanning T are those at constant R
            for R_cnt in range(-nT_plots, nT_plots + 1):
                if nT_plots * slow_incr > slowR_hi:
                    print('nT_plots * slow_incr > slowR_hi, out of range')
                    sys.exit()
                if -(nT_plots * slow_incr) < slowR_lo:
                    print('-nT_plots * slow_incr < slowR_lo, out of range')
                    sys.exit()
                #%% -- -- gather T data
                lowest_Rslow = 1000000  # find index of row with R_cnt slowness, awkward coding
                target_slow = (R_cnt * slow_incr)
                for slow_i in range(slowR_n):
                    if abs(stack_Rslows[slow_i] - target_slow) < lowest_Rslow:
                        lowest_Rindex = slow_i
                        lowest_Rslow = abs(stack_Rslows[slow_i] - target_slow)

                print(f'For T plot {R_cnt:2d}, {lowest_Rindex:3d} is R slow nearest {target_slow:.3f}, difference is {lowest_Rslow:.3f}')

                # Collect data with that slowness for T (R=const) plot
                Tcentral_st = Stream()
                Tcentral_am = Stream()
                for slowT_i in range(slowT_n):
                    Tcentral_st += tdiff[  lowest_Rindex*slowT_n + slowT_i]
                    Tcentral_am += amp_ave[lowest_Rindex*slowT_n + slowT_i]

                #%% -- -- plot T tdiff
                if no_tdiff_plot == False:
                    stack_arrayT_Tdf = np.zeros((slowT_n,stack_nt))
                    for it in range(stack_nt):  # check points one at a time
                        for slowT_i in range(slowT_n):  # for this station, loop over slownesses
                            num_val = Tcentral_st[slowT_i].data[it]
                            stack_arrayT_Tdf[slowT_i, it] = num_val

                    y, x = np.mgrid[slice(stack_Tslows[0], stack_Tslows[-1] + slow_delta, slow_delta),
                                 slice(ttt[0], ttt[-1] + dt, dt)]

                    fig, ax = plt.subplots(1, figsize=(10,3))
                    c = ax.pcolormesh(x, y, stack_arrayT_Tdf, cmap=plt.cm.coolwarm, vmin= -tdiff_clip, vmax=tdiff_clip)
                    fig.subplots_adjust(bottom=0.2)
                    ax.axis([x.min(), x.max(), y.min(), y.max()])
                    fig.colorbar(c, ax=ax, label='time shift (s)')
                    c = ax.scatter(arrival_time, 0, color='black'  , s=50, alpha=0.75)
                    plt.xlabel('Time (s)')
                    plt.ylabel('Transverse slowness (s/km)')
                    plt.title(f'{ref_phase} Tdiff at {target_slow:.3f} s/km R slowness, {fname1[48:58]}  {fname1[59:69]}  min amp {min_amp:.1f}  cc_thres {cc_thres:.2f}')
                    if no_plots == False:
                        plt.show()
                    fig_index += 1

                #%% -- -- plot T amp
                stack_arrayT_Amp = np.zeros((slowT_n,stack_nt))
                for it in range(stack_nt):  # check points one at a time
                    for slowT_i in range(slowT_n):  # for this station, loop over slownesses
                        num_val = Tcentral_am[slowT_i].data[it]
                        stack_arrayT_Amp[slowT_i, it] = num_val

                y, x = np.mgrid[slice(stack_Tslows[0], stack_Tslows[-1] + slow_delta, slow_delta),
                             slice(ttt[0], ttt[-1] + dt, dt)]

                fig, ax = plt.subplots(1, figsize=(10,3))
                if log_plot == True:
                    c = ax.pcolormesh(x, y, stack_arrayT_Amp - global_max, cmap=plt.cm.gist_rainbow_r, vmin= - log_plot_range, vmax=0)
                else:
                    c = ax.pcolormesh(x, y, stack_arrayT_Amp, cmap=plt.cm.gist_rainbow_r, vmin= 0, vmax=global_max)
                fig.subplots_adjust(bottom=0.2)
                ax.axis([x.min(), x.max(), y.min(), y.max()])
                if log_plot == True:
                    fig.colorbar(c, ax=ax, label='log amplitude')
                else:
                    fig.colorbar(c, ax=ax, label='linear amplitude')
                c = ax.scatter(arrival_time, 0, color='black'  , s=50, alpha=0.75)
                plt.xlabel('Time (s)')
                if NS == True:
                    plt.ylabel('East Slowness (s/km)')
                    plt.title(f'Amp at {target_slow:.3f} s/km N slowness, {fname1[48:58]}  {fname1[59:69]}')
                else:
                    plt.ylabel('Transverse Slowness (s/km)')
                    plt.title(f'Amp at {target_slow:.3f} s/km R slowness, {fname1[48:58]}  {fname1[59:69]}')
                if no_plots == False:
                    plt.show()
                fig_index += 1

    #%% 2-slices-plus-snaps option
    if two_slice_plots == True:
        #%% -- Find slowness arrays for R and T slices
        lowest_Tslow = 1000000
        for slow_i in range(slowT_n):
            if abs(stack_Tslows[slow_i] - T_slow_plot) < lowest_Tslow:
                lowest_Tindex = slow_i
                lowest_Tslow = abs(stack_Tslows[slow_i] - T_slow_plot)

        print(f'{lowest_Tindex:4d} is T slow nearest {T_slow_plot:.3f}, difference is {lowest_Tslow:.3f}')

        lowest_Rslow = 1000000
        for slow_i in range(slowR_n):
            if abs(stack_Rslows[slow_i] - R_slow_plot) < lowest_Rslow:
                lowest_Rindex = slow_i
                lowest_Rslow = abs(stack_Rslows[slow_i] - R_slow_plot)

        print(f'{lowest_Rindex:4d} is R slow nearest 0.005, difference is {lowest_Rslow:.3f}')

        #%% -- Extract and sum tdiff and amp slices from beam matrix
        # Select only stacks with that slowness for Transverse plot
        centralR_Dst = Stream()
        for slowR_i in range(slowR_n):
            centralR_Dst += tdiff[slowR_i*slowT_n + lowest_Tindex]

        # Select only stacks with that slowness for Radial plot
        centralT_Dst = Stream()
        for slowT_i in range(slowT_n):
            centralT_Dst += tdiff[lowest_Rindex*slowT_n + slowT_i]

        # Select only stacks with that slowness for Transverse plot
        centralR_Ast = Stream()
        for slowR_i in range(slowR_n):
            centralR_Ast += amp_ave[slowR_i*slowT_n + lowest_Tindex]

        # Select only stacks with that slowness for Radial plot
        centralT_Ast = Stream()
        for slowT_i in range(slowT_n):
            centralT_Ast += amp_ave[lowest_Rindex*slowT_n + slowT_i]

        #%% -- Stack plots
        stack_array = np.zeros((slowR_n,stack_nt))

        for it in range(stack_nt):  # check points one at a time
            for slowR_i in range(slowR_n):  # for this station, loop over slownesses
                num_val = centralR_Dst[slowR_i].data[it]
                stack_array[slowR_i, it] = num_val

        y, x = np.mgrid[slice(stack_Rslows[0], stack_Rslows[-1] + slow_delta, slow_delta),
                     slice(ttt[0], ttt[-1] + dt, dt)]

        fig, ax = plt.subplots(1, figsize=(10,3))
        print(f'len(x) is {len(x)} and len(y) is {len(y)}')
        print(f'len(stack_Rslows) is {len(stack_Rslows)} and len(ttt) is {len(ttt)}')
        print(f'slowR_n is {slowR_n} and stack_nt is {stack_nt}')
        c = ax.pcolormesh(x, y, stack_array, cmap=plt.cm.coolwarm, vmin=-tdiff_clip, vmax=tdiff_clip)
        fig.colorbar(c, ax=ax, label='time lag (s)')
        fig.subplots_adjust(bottom=0.2)
        ax.axis([x.min(), x.max(), y.min(), y.max()])
        c = ax.scatter(arrival_time, event_pred_slo, color='black'  , s=50, alpha=0.75)
        plt.xlabel('Time (s)')
        if NS == True:
            plt.ylabel('N slowness (s/km)')
            plt.title('Time lag at ' + str(T_slow_plot) + ' s/km E slowness, ' + date_label1 + ' ' + date_label2)
            plt.savefig(save_name + str(start_buff) + '_' + str(end_buff) + '_NtdiffSection.png')
        else:
            plt.ylabel('R slowness (s/km)')
            plt.title('Time lag at ' + str(T_slow_plot) + ' s/km T slowness, ' + date_label1 + ' ' + date_label2)
            plt.savefig(save_name + str(start_buff) + '_' + str(end_buff) + '_RtdiffSection.png')
        if no_plots == False:
            plt.show()

        fig_index += 1

        stack_array = np.zeros((slowT_n,stack_nt))

        for it in range(stack_nt):  # check points one at a time
            for slowT_i in range(slowT_n):  # for this station, loop over slownesses
                num_val = centralT_Dst[slowT_i].data[it]
                stack_array[slowT_i, it] = num_val

        y, x = np.mgrid[slice(stack_Tslows[0], stack_Tslows[-1] + slow_delta, slow_delta),
                     slice(ttt[0], ttt[-1] + dt, dt)]

        fig, ax = plt.subplots(1, figsize=(10,3))
        fig.subplots_adjust(bottom=0.2)
        c = ax.pcolormesh(x, y, stack_array, cmap=plt.cm.coolwarm, vmin=-tdiff_clip, vmax=tdiff_clip)
        fig.colorbar(c, ax=ax, label='time lag (s)')
        ax.axis([x.min(), x.max(), y.min(), y.max()])
        c = ax.scatter(arrival_time, 0, color='black'  , s=50, alpha=0.75)
        plt.xlabel('Time (s)')
        if NS == True:
            plt.ylabel('E slowness (s/km)')
            plt.title('Time lag at ' + str(R_slow_plot) + ' s/km N slowness, ' + date_label1 + ' ' + date_label2)
            plt.savefig(save_name + str(start_buff) + '_' + str(end_buff) + '_EtdiffSection.png')
        else:
            plt.ylabel('T slowness (s/km)')
            plt.title('Time lag at ' + str(R_slow_plot) + ' s/km R slowness, ' + date_label1 + ' ' + date_label2)
            plt.savefig(save_name + str(start_buff) + '_' + str(end_buff) + '_TtdiffSection.png')
        if no_plots == False:
            plt.show()

        fig_index += 1

        stack_array = np.zeros((slowR_n,stack_nt))

        for it in range(stack_nt):  # check points one at a time
            for slowR_i in range(slowR_n):  # for this station, loop over slownesses
                num_val = centralR_Ast[slowR_i].data[it]
                stack_array[slowR_i, it] = num_val

        y, x = np.mgrid[slice(stack_Rslows[0], stack_Rslows[-1] + slow_delta, slow_delta),
                     slice(ttt[0], ttt[-1] + dt, dt)]

        fig, ax = plt.subplots(1, figsize=(10,3))
        print(f'len(x) is {len(x)} and len(y) is {len(y)}')
        print(f'len(stack_Rslows) is {len(stack_Rslows)} and len(ttt) is {len(ttt)}')
        print(f'slowR_n is {slowR_n} and stack_nt is {stack_nt}')
        c = ax.pcolormesh(x, y, stack_array, cmap=plt.cm.gist_rainbow_r, vmin=0)
        fig.colorbar(c, ax=ax, label='linear amplitude')
        fig.subplots_adjust(bottom=0.2)
        ax.axis([x.min(), x.max(), y.min(), y.max()])
        c = ax.scatter(arrival_time, event_pred_slo, color='black'  , s=50, alpha=0.75)
        plt.xlabel('Time (s)')
        if NS == True:
            plt.ylabel('N slowness (s/km)')
            plt.title('Amp at ' + str(T_slow_plot) + ' s/km E slowness, ' + date_label1 + ' ' + date_label2)
            plt.savefig(save_name + str(start_buff) + '_' + str(end_buff) + '_NampSection.png')
        else:
            plt.ylabel('R slowness (s/km)')
            plt.title('Amp at ' + str(T_slow_plot) + ' s/km T slowness, ' + date_label1 + ' ' + date_label2)
            plt.savefig(save_name + str(start_buff) + '_' + str(end_buff) + '_RampSection.png')
        if no_plots == False:
            plt.show()

        fig_index += 1

        stack_array = np.zeros((slowT_n,stack_nt))

        for it in range(stack_nt):  # check points one at a time
            for slowT_i in range(slowT_n):  # for this station, loop over slownesses
                num_val = centralT_Ast[slowT_i].data[it]
                stack_array[slowT_i, it] = num_val

        y, x = np.mgrid[slice(stack_Tslows[0], stack_Tslows[-1] + slow_delta, slow_delta),
                     slice(ttt[0], ttt[-1] + dt, dt)]

        fig, ax = plt.subplots(1, figsize=(10,3))
        fig.subplots_adjust(bottom=0.2)
        c = ax.pcolormesh(x, y, stack_array, cmap=plt.cm.gist_rainbow_r, vmin=0)
        fig.colorbar(c, ax=ax, label='linear amplitude')
        ax.axis([x.min(), x.max(), y.min(), y.max()])
        c = ax.scatter(arrival_time, 0, color='black'  , s=50, alpha=0.75)
        plt.xlabel('Time (s)')
        if NS == True:
            plt.ylabel('E slowness (s/km)')
            plt.title('Amp at ' + str(R_slow_plot) + ' s/km N slowness, ' + date_label1 + ' ' + date_label2)
            plt.savefig(save_name + str(start_buff) + '_' + str(end_buff) + '_EampSection.png')
        else:
            plt.ylabel('T slowness (s/km)')
            plt.title('Amp at ' + str(R_slow_plot) + ' s/km R slowness, ' + date_label1 + ' ' + date_label2)
            plt.savefig(save_name + str(start_buff) + '_' + str(end_buff) + '_TampSection.png')
        if no_plots == False:
            plt.show()

        fig_index += 1

    #%% -- Snap plots
        stack_slice = np.zeros((slowR_n,slowT_n))
        if snaps > 0:
            # check for impossible parameters
            if (start_buff + snaptime) < start_buff:
                print(f'snaptime {start_buff + snaptime:.0f} is earlier than start_buff of {start_buff:.0f}')
                sys.exit(-1)
            last_snap = snaptime + snaps*snap_depth
            if (start_buff + last_snap) > end_buff:
                print(f'last snap {last_snap:.0f} is later than end_buff of {end_buff:.0f}')
                sys.exit(-1)

            for snap_num in range(snaps):
                snap_start = start_buff + snaptime + (snap_num  ) * snap_depth
                snap_end   = start_buff + snaptime + (snap_num+1) * snap_depth
                fig_index += 1
                it_start = int(round((snap_start - start_buff)/dt))
                it_end   = int(round((snap_end   - start_buff)/dt))
                # plt.savefig(save_name + str(start_buff + snap_start) + '_' + str(start_buff + snap_end) + '_Abeam.png')
                for slowR_i in range(slowR_n):  # loop over radial slownesses
                    for slowT_i in range(slowT_n):  # loop over transverse slownesses
                        index = slowR_i*slowT_n + slowT_i
                        num_val = np.nanmean(tdiff[index].data[it_start:it_end])
                        stack_slice[slowR_i, slowT_i] = num_val

                y1, x1 = np.mgrid[slice(stack_Rslows[0], stack_Rslows[-1] + slow_delta, slow_delta),
                             slice(stack_Tslows[0], stack_Tslows[-1] + slow_delta, slow_delta)]

                fig, ax = plt.subplots(1, figsize=(7,0.8*7))
                c = ax.pcolormesh(x1, y1, stack_slice, cmap=plt.cm.coolwarm, vmin=-tdiff_clip, vmax=tdiff_clip)
                ax.axis([x1.min(), x1.max(), y1.min(), y1.max()])
                fig.colorbar(c, ax=ax, label = 'time shift (s)')
                c = ax.scatter(pred_Eslo, pred_Nslo, color='black'  , s=50, alpha=0.75)
                c = ax.scatter(        0,         0, color='black' , s=50,  alpha=0.75)
                circle1 = plt.Circle((0, 0), 0.019, color='black', fill=False)
                ax.add_artist(circle1)  # inner core limit
                circle2 = plt.Circle((0, 0), 0.040, color='black', fill=False)
                ax.add_artist(circle2)  # outer core limit
                plt.title(f'Tdiff {snap_start:.0f} to {snap_end:.0f}s  {date_label1} {date_label2}  events {eq_num1}&{eq_num2}')
                if NS == True:
                    plt.xlabel('E Slowness (s/km)')
                    plt.ylabel('N Slowness (s/km)')
                else:
                    plt.xlabel('T Slowness (s/km)')
                    plt.ylabel('R Slowness (s/km)')
                plt.savefig(save_name + str(snap_start) + '_' + str(snap_end) + '_Tbeam.png')
                if no_plots == False:
                    plt.show()

            for snap_num in range(snaps):
                snap_start = start_buff + snaptime + (snap_num  ) * snap_depth
                snap_end   = start_buff + snaptime + (snap_num+1) * snap_depth
                fig_index += 1
                it_start = int(round((snap_start - start_buff)/dt))
                it_end   = int(round((snap_end   - start_buff)/dt))
                for slowR_i in range(slowR_n):  # loop over radial slownesses
                    for slowT_i in range(slowT_n):  # loop over transverse slownesses
                        index = slowR_i*slowT_n + slowT_i
                        num_val = np.nanmean(amp_ave[index].data[it_start:it_end])
                        stack_slice[slowR_i, slowT_i] = num_val

                y1, x1 = np.mgrid[slice(stack_Rslows[0], stack_Rslows[-1] + slow_delta, slow_delta),
                             slice(stack_Tslows[0], stack_Tslows[-1] + slow_delta, slow_delta)]

                fig, ax = plt.subplots(1, figsize=(7,0.8*7))
                c = ax.pcolormesh(x1, y1, stack_slice, cmap=plt.cm.gist_rainbow_r)
                ax.axis([x1.min(), x1.max(), y1.min(), y1.max()])
                if log_plot == True:
                    fig.colorbar(c, ax=ax, label = 'log amp')
                else:
                    fig.colorbar(c, ax=ax, label = 'linear amp')
                c = ax.scatter(pred_Eslo, pred_Nslo, color='black'  , s=50, alpha=0.75)
                c = ax.scatter(        0,         0, color='black' , s=50,  alpha=0.75)
                circle1 = plt.Circle((0, 0), 0.019, color='black', fill=False)
                ax.add_artist(circle1)  # inner core limit
                circle2 = plt.Circle((0, 0), 0.040, color='black', fill=False)
                ax.add_artist(circle2)  # outer core limit
                plt.title(f'Amp {snap_start:.0f} to {snap_end:.0f}s  {date_label1} {date_label2}  events {eq_num1}&{eq_num2}')
                if NS == True:
                    plt.xlabel('E Slowness (s/km)')
                    plt.ylabel('N Slowness (s/km)')
                else:
                    plt.xlabel('T Slowness (s/km)')
                    plt.ylabel('R Slowness (s/km)')
                plt.savefig(save_name + str(snap_start) + '_' + str(snap_end) + '_Abeam.png')
                if no_plots == False:
                    plt.show()

    #%% Wiggly plots
    if wiggly_plots == True:

        #%% -- read wiggle beams
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

        #%% -- Extract slices to wiggle plot, cumbersome, every slowness along slice is now selected
        #%% -- -- Collect T slowness nearest T_slow
        lowest_Tslow = 1000000
        for slow_i in range(slowT_n):
            if abs(stack_Tslows[slow_i] - T_slow_plot) < lowest_Tslow:
                lowest_Tindex = slow_i
                lowest_Tslow = abs(stack_Tslows[slow_i] - T_slow_plot)

        print(f'{slowT_n} T slownesses, index {lowest_Tindex} is closest to requested plot T slowness {T_slow_plot:.4f}, slowness diff there is {lowest_Tslow:.4f} and slowness is {stack_Tslows[lowest_Tindex]:.4f}')
        # Select only stacks with that slowness for radial plot
        centralR_st1 = Stream()
        centralR_st2 = Stream()
        centralR_amp   = Stream()
        centralR_tdiff = Stream()
        for slowR_i in range(slowR_n):
            ii = slowR_i*slowT_n + lowest_Tindex
            centralR_st1 += st1[ii]
            centralR_st2 += st2[ii]
            centralR_amp   += amp_ave[ii]
            centralR_tdiff += tdiff[ii]

        #%% -- -- Collect R slowness nearest R_slow
        lowest_Rslow = 1000000
        for slow_i in range(slowR_n):
            if abs(stack_Rslows[slow_i] - R_slow_plot) < lowest_Rslow:
                lowest_Rindex = slow_i
                lowest_Rslow = abs(stack_Rslows[slow_i] - R_slow_plot)

        print(f'{slowR_n} R slownesses, index {lowest_Rindex} is closest to requested plot R slowness {R_slow_plot:.4f}, slowness diff there is {lowest_Rslow:.4f} and slowness is {stack_Rslows[lowest_Rindex]:.4f}')

        # Select only stacks with that slowness for transverse plot
        centralT_st1 = Stream()
        centralT_st2 = Stream()
        centralT_amp   = Stream()
        centralT_tdiff = Stream()

        for slowT_i in range(slowT_n):
            ii = lowest_Rindex*slowT_n + slowT_i
            centralT_st1 += st1[ii]
            centralT_st2 += st2[ii]
            centralT_amp   += amp_ave[ii]
            centralT_tdiff += tdiff[ii]

    #%% -- Plot wiggles
        #%% -- -- Compute timing time series
        ttt_dec = (np.arange(len(tdiff[0].data)) * tdiff[0].stats.delta + start_buff) # in units of seconds

        #%% -- -- R amp and tdiff vs time plots with black line for time shift
        scale_plot_wig = wig_scale_fac / (200 * global_max)
        scale_plot_tdiff = tdiff_scale_fac / 500.
        if log_plot == True:
            scale_plot_wig /= 30  # not quite sure why this renormalization works
            # scale_plot_tdiff = plot_scale_fac / 500.
        fig_index = 116
        plt.figure(fig_index,figsize=(30,10))
        plt.xlim(start_buff,end_buff)
        plt.ylim(stack_Rslows[0], stack_Rslows[-1])
        for slowR_i in range(slowR_n):  # loop over radial slownesses
            dist_offset = stack_Rslows[slowR_i] # trying for approx degrees
            ttt1 = (np.arange(len(centralR_st1[slowR_i].data)) * centralR_st1[slowR_i].stats.delta
              + (centralR_st1[slowR_i].stats.starttime - t1))
            plt.plot(ttt1, ((centralR_st1[slowR_i].data - np.median(centralR_st1[slowR_i].data)) * scale_plot_wig) + dist_offset, color = 'green')
            plt.plot(ttt1, ((centralR_st2[slowR_i].data - np.median(centralR_st2[slowR_i].data)) * scale_plot_wig) + dist_offset, color = 'red')
            if turn_off_black == 0:
                plt.plot(ttt1,     (centralT_st1[slowT_i].data)*0.0 + dist_offset, color = 'gray') # reference lines
                plt.plot(ttt_dec, (centralR_tdiff[slowR_i].data) * scale_plot_tdiff + dist_offset, color = 'black')

        plt.xlabel('Time (s)')
        if NS == True:
            plt.ylabel('N Slowness (s/km)')
            plt.title(date_label1 + '  ' + date_label2 + '  ' + ' seismograms ' + str(T_slow_plot) + ' E slowness, green is event1, red is event2')
            plt.savefig(save_name + str(start_buff) + '_' + str(end_buff) + '_N_pro_wig.png')
        else:
            plt.ylabel('R Slowness (s/km)')
            plt.title(date_label1 + '  ' + date_label2 + '  ' + ' seismograms ' + str(T_slow_plot) + ' T slowness, green is event1, red is event2')
            plt.savefig(save_name + str(start_buff) + '_' + str(end_buff) + '_R_pro_wig.png')
        #%% -- -- T amp and tdiff vs time plots with black line for time shift
        fig_index = 117
        plt.figure(fig_index,figsize=(30,10))
        plt.xlim(start_buff,end_buff)
        plt.ylim(stack_Tslows[0], stack_Tslows[-1])

        for slowT_i in range(slowT_n):  # loop over transverse slownesses
            dist_offset = stack_Tslows[slowT_i] # trying for approx degrees
            ttt2 = (np.arange(len(centralT_st1[slowT_i].data)) * centralT_st1[slowT_i].stats.delta
              + (centralT_st1[slowT_i].stats.starttime - t1))
            plt.plot(ttt2, ((centralT_st1[slowT_i].data - np.median(centralT_st1[slowT_i].data)) * scale_plot_wig) + dist_offset, color = 'green')
            plt.plot(ttt2, ((centralT_st2[slowT_i].data - np.median(centralT_st2[slowT_i].data)) * scale_plot_wig) + dist_offset, color = 'red')
            if turn_off_black == 0:
                plt.plot(ttt2,     (centralT_st1[slowT_i].data)*0.0 + dist_offset, color = 'gray') # reference lines
                plt.plot(ttt_dec, (centralT_tdiff[slowT_i].data) * scale_plot_tdiff + dist_offset, color = 'black')
        plt.xlabel('Time (s)')
        if NS == True:
            plt.ylabel('E Slowness (s/km)')
            plt.title(date_label1 + '  ' + date_label2 + '  ' + ' seismograms ' + str(R_slow_plot) + ' N slowness, green is event1, red is event2')
            plt.savefig(save_name + str(start_buff) + '_' + str(end_buff) + '_E_pro_wig.png')
        else:
            plt.ylabel('T Slowness (s/km)')
            plt.title(date_label1 + '  ' + date_label2 + '  ' + ' seismograms ' + str(R_slow_plot) + ' R slowness, green is event1, red is event2')
            plt.savefig(save_name + str(start_buff) + '_' + str(end_buff) + '_T_pro_wig.png')

    #%% Beam sum plots
    if beam_sums == True:
    #%% -- R-T tdiff amp-normed
        stack_slice = np.zeros((slowR_n,slowT_n))

        if start_beam == 0 and end_beam == 0:
            full_beam = 1
            print('Full beam is specified.')
        else:  # beam just part of stack volume
            full_beam = 0
            start_index = int((start_beam - start_buff) / dt)
            end_index   = int((end_beam   - start_buff) / dt)
            print(f'Beam is {start_beam:.4f} to {end_beam:.4f}s, out of {start_buff:.4f} to {end_buff:.4f}s, dt is {dt:.4f}s, and indices are {start_index} {end_index}')

        for slowR_i in range(slowR_n):  # loop over radial slownesses
            for slowT_i in range(slowT_n):  # loop over transverse slownesses
                index = slowR_i*slowT_n + slowT_i
                if full_beam == 1: # using elementwise multiplication, amplitude weighted
                    num_val = np.nansum(np.multiply(tdiff[index].data, amp_ave_thres[index].data))/np.nansum(amp_ave_thres[index].data)
                else:
                    num_val = np.nansum(np.multiply(tdiff[start_index:end_index].data, amp_ave_thres[start_index:end_index].data
                                                     ))/np.nansum(amp_ave_thres[start_index:end_index].data)
                stack_slice[slowR_i, slowT_i] = num_val

        y1, x1 = np.mgrid[slice(stack_Rslows[0], stack_Rslows[-1] + slow_delta, slow_delta),
                     slice(stack_Tslows[0], stack_Tslows[-1] + slow_delta, slow_delta)]

        fig, ax = plt.subplots(1, figsize=(7,0.8*7*(slowR_n/slowT_n)))  # try to make correct aspect ratio plot
        c = ax.pcolormesh(x1, y1, stack_slice, cmap=plt.cm.coolwarm, vmin = -tdiff_clip, vmax = tdiff_clip)
        fig.colorbar(c, ax=ax, label='time shift (s)')
        ax.axis([x1.min(), x1.max(), y1.min(), y1.max()])
        circle1 = plt.Circle((0, 0), 0.019, color='black', fill=False)
        ax.add_artist(circle1)
        circle2 = plt.Circle((0, 0), 0.040, color='black', fill=False)
        ax.add_artist(circle2)  #outer core limit

        c = ax.scatter(pred_Eslo, pred_Nslo, color='black'  , s=50, alpha=0.75)
        c = ax.scatter(        0,         0, color='black' , s=50,  alpha=0.75)

        if NS == True:
            plt.ylabel('N Slowness (s/km)')
            plt.xlabel('E Slowness (s/km)')
        else:
            plt.ylabel('R Slowness (s/km)')
            plt.xlabel('T Slowness (s/km)')
        plt.title(f'{date_label1} {date_label2} {start_buff:.0f} to {end_buff:.0f} time shift')
        plt.savefig(save_name + str(start_buff) + '_' + str(end_buff) + '_Tbeam.png')
        if no_plots == False:
            plt.show()

    #%% -- R-T amplitude averaged over time window
        stack_slice = np.zeros((slowR_n,slowT_n))
        for slowR_i in range(slowR_n):  # loop over radial slownesses
            for slowT_i in range(slowT_n):  # loop over transverse slownesses
                index = slowR_i*slowT_n + slowT_i
                if full_beam == 1:
                    num_val = np.nanmean(amp_ave[index].data)
                else:
                    num_val = np.nanmean(amp_ave[index].data[start_index:end_index])
                stack_slice[slowR_i, slowT_i] = num_val

        y1, x1 = np.mgrid[slice(stack_Rslows[0], stack_Rslows[-1] + slow_delta, slow_delta),
                     slice(stack_Tslows[0], stack_Tslows[-1] + slow_delta, slow_delta)]

        fig, ax = plt.subplots(1, figsize=(7,0.8*7*(slowR_n/slowT_n)))
        smax = np.max(stack_slice)
        smin = np.min(stack_slice)
        if log_plot == True:
            if (smax - smin) < log_plot_range:  # use full color scale even if range is less than specified
                log_plot_range = smax - smin
            c = ax.pcolormesh(x1, y1, stack_slice - smax, cmap=plt.cm.gist_rainbow_r, vmin= - log_plot_range, vmax=0)
        else:
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

        c = ax.scatter(pred_Eslo, pred_Nslo, color='black'  , s=50, alpha=0.75)
        c = ax.scatter(        0,         0, color='black' , s=50,  alpha=0.75)

        if NS == True:
            plt.xlabel('E Slowness (s/km)')
            plt.ylabel('N Slowness (s/km)')
        else:
            plt.xlabel('T Slowness (s/km)')
            plt.ylabel('R Slowness (s/km)')
        plt.title(f'{date_label1} {date_label2} {start_buff:.0f} to {end_buff:.0f} beam amp')
        plt.savefig(save_name + str(start_buff) + '_' + str(end_buff) + '_Abeam.png')
        if no_plots == False:
            plt.show()

    #  Save processed files
#    fname = 'HD' + date_label + '_slice.mseed'
#    stack.write(fname,format = 'MSEED')

    elapsed_time_wc = time.time() - start_time_wc
    print(f'This job took {elapsed_time_wc:.1f} seconds')
    os.system('say "Done"')