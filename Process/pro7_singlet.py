#!/usr/bin/env python
# looks at envelope stacks from pro6
# window by signal quality
# John Vidale 2/2019, overhauled 1/2021, revisited 9/2021

def pro7_singlet(eq_num, slow_delta = 0.0005, turn_off_black = 0,
              slowR_lo = -0.1, slowR_hi = 0.1, slowT_lo = -0.1, slowT_hi = 0.1,
              start_buff = 50, end_buff = 50, fig_index = 401, do_T = False, do_R = False,
              ZslowR_lo = -0.1, ZslowR_hi = 0.1, ZslowT_lo = -0.1, ZslowT_hi = 0.1,
              Zstart_buff = 50, Zend_buff = 50, zoom = False, ref_phase = 'blank', min_amp = 0.2,
              R_slow_plot = 0.004, T_slow_plot = -0.00,
              snaptime = 0, snaps = 0, snap_depth = 1,
              nR_plots  = 1, nT_plots = 1, slow_incr = 0.01, NS = False,
              ARRAY = 0, auto_slice = True, two_slice_plots = True, beam_sums = 1,
              wiggly_plots = False, max_wiggly_plot = True, start_beam = 0, end_beam = 0, log_plot = False,
              log_plot_range = 2, wig_scale_fac = 1, ref_loc = True, ref_lat = 0, ref_lon = 0):

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

    print(colored('Running pro7b_singlet', 'cyan'))
    start_time_wc = time.time()

    if zoom == True:
        orig_start_buff = start_buff
        orig_end_buff   = end_buff
        if Zstart_buff  < start_buff:
            print(f'Zstart_buff of {Zstart_buff:.1f} cannot be < start_buff of {start_buff:.1f}')
            Zstart_buff = start_buff
        if Zend_buff    > end_buff:
            print(f'Zend_buff of {Zend_buff:.1f} cannot be > end_buff of {end_buff:.1f}')
            Zend_buff   = end_buff

    #%% Input parameters and computed files
    fname = '/Users/vidale/Documents/Research/IC/EvLocs/event' + str(eq_num) + '.txt'
    folder_name = '/Users/vidale/Documents/Research/IC/'
    file = open(fname, 'r')
    lines=file.readlines()
    split_line = lines[0].split()
    t          = UTCDateTime(split_line[1])
    date_label  = split_line[1][0:10]
    ev_lat      = float(      split_line[2])
    ev_lon      = float(      split_line[3])
    ev_depth    = float(      split_line[4])
    # date_label = '2018-04-02' # dates in filename

    if ref_loc == False:
        if ARRAY == 0: # Hinet set
            ref_lat = 36.3
            ref_lon = 138.5
        elif ARRAY == 1:         # LASA set
            ref_lat = 46.69
            ref_lon = -106.22
        elif ARRAY == 2: # China set and center
            ref_lat = 38      # °N
            ref_lon = 104.5   # °E
    ref_distance = gps2dist_azimuth(ref_lat, ref_lon, ev_lat, ev_lon)
    ref_dist     = ref_distance[0]/(1000*111)
    ref_az       = ref_distance[1]
    ref_back_az  = ref_distance[2]

    print(ref_phase + ' is ref phase, depth is ' + str(ev_depth) + ' at distance ' + str(ref_dist))

    arrivals_ref   = model.get_travel_times(source_depth_in_km=ev_depth,distance_in_degree=ref_dist, phase_list=[ref_phase])
    arrival_time = arrivals_ref[0].time
    print(ref_phase + ' arrival time is ' + str(arrival_time) + ' at distance ' + str(ref_dist))

    atime_rayp = arrivals_ref[0].ray_param
    # atime_ref  = arrivals_ref[0].time
    event_pred_slo  = atime_rayp * 2 * np.pi / (360. * 111) # convert to s/°

    # convert to pred rslo and tslo
    if NS == True:     #  rotate predicted slowness to N and E
        print(f'Array  lat {ref_lat:.0f}, lon  {ref_lon:.0f}, Event lat {ev_lat:.0f}, lon {ev_lon:.0f}, az {ref_az:.0f}, baz {ref_back_az:.0f}')
        sin_baz = np.sin(ref_az * np.pi /180)
        cos_baz = np.cos(ref_az * np.pi /180)
        pred_Nslo = event_pred_slo * cos_baz
        pred_Eslo = event_pred_slo * sin_baz
        print(f'pred slowness is {pred_Nslo:.3f} N and {pred_Eslo:.3f} E')
    else:
        pred_Nslo = event_pred_slo
        pred_Eslo = 0
        print(f'pred slowness is {pred_Nslo:.3f} R and {pred_Eslo:.3f} T')

    name_str = folder_name + 'Pro_files/HD' + date_label + '_'
    fname  = name_str + 'amp_ave.mseed'
    amp_ave = Stream()
    amp_ave = read(fname)

    dt     =     amp_ave[0].stats.delta
    nt     = len(amp_ave[0].data)
    nt_amp = len(amp_ave[0].data)
    print(f'amp_ave data length is {nt_amp} time pts, dt is {dt:.2f}, so record length is {dt * nt   :.0f} seconds')
    print(f'slowness grid has {len(amp_ave)} elements')

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
        Zamp_ave = Stream()
        print(f'before trim, amp_ave[0] has length {len(amp_ave[0])}')
        for slowR_i in range(slowR_n):  # loop over radial slownesses
            for slowT_i in range(slowT_n):  # loop over transverse slownesses, kludge to evade rounding error
                if ((stack_Rslows[slowR_i] >= ZslowR_lo - 0.000001) and (stack_Rslows[slowR_i] <= ZslowR_hi + 0.000001) and
                    (stack_Tslows[slowT_i] >= ZslowT_lo - 0.000001) and (stack_Tslows[slowT_i] <= ZslowT_hi + 0.000001)):
                    index = slowR_i*slowT_n + slowT_i
                    s_t = t + Zstart_buff
                    e_t = t + Zend_buff
                    Zamp_ave += amp_ave[index].trim(starttime=s_t, endtime=e_t)
                            #tr.trim(starttime=s_t,endtime = e_t)
        amp_ave = Zamp_ave
        nt = len(amp_ave[0].data)
        start_buff = Zstart_buff
        # make time series
        print(f'after trim, amp_ave[0] has length {len(amp_ave[0])}')
        print(f'input slowR_lo  is {slowR_lo}  and slowR_hi  is {slowR_hi}  and slowT_lo  is {slowT_lo}  and slowT_hi is {slowT_hi}')
        print(f'trimmed ZslowR_lo is {ZslowR_lo} and ZslowR_hi is {ZslowR_hi} and ZslowT_lo is {ZslowT_lo} and ZslowT_hi is {ZslowT_hi}')

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
        print('Output trace starttime ' + str(amp_ave[0].stats.starttime))

    ttt = (np.arange(len(amp_ave[0].data)) * amp_ave[0].stats.delta + start_buff) # in units of seconds

    global_max = 0  # find global_max, largest amplitude in amp_ave beam array envelopes
    for slow_i in range(len(amp_ave)): # find global max of ave_amp
        local_max = max(amp_ave[slow_i].data)
        if local_max > global_max:
            global_max = local_max

    #%% Mask out weak and/or less correlated slowness points
    amp_ave_thres = amp_ave.copy()  # copy amp envelope array, set amps and tdiff below thresholds to NaN using global_max
    nt = len(amp_ave[0].data)
    for slow_i in range(len(amp_ave)): # don't plot less robust points
        for it in range(nt):
            if (amp_ave[slow_i].data[it] < (min_amp * global_max)):
                amp_ave_thres[slow_i].data[it] = np.nan

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
                Rcentral_am = Stream()
                for slowR_i in range(slowR_n):
                    Rcentral_am += amp_ave[slowR_i*slowT_n + lowest_Tindex]

                #%% -- -- plot R amp
                stack_arrayR_Amp = np.zeros((slowR_n,stack_nt))
                for it in range(stack_nt):  # check points one at a time
                    for slowR_i in range(slowR_n):  # for this station, loop over slownesses
                        num_val = Rcentral_am[slowR_i].data[it]
                        stack_arrayR_Amp[slowR_i, it] = num_val

                y, x = np.mgrid[slice(stack_Rslows[0], stack_Rslows[-1] + slow_delta, slow_delta),
                             slice(ttt[0], ttt[-1] + dt, dt)]

                fig_index += 1
                fig, ax = plt.subplots(1, figsize=(10,3), num = fig_index)
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
                    plt.title(f'Amp at {target_slow:.3f} s/km E slowness, {fname[48:58]} #{eq_num}')
                else:
                    plt.ylabel('Radial Slowness (s/km)')
                    plt.title(f'Amp at {target_slow:.3f} s/km T slowness, {fname[48:58]} #{eq_num}')
                plt.show()

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
                Tcentral_am = Stream()
                for slowT_i in range(slowT_n):
                    Tcentral_am += amp_ave[lowest_Rindex*slowT_n + slowT_i]

                #%% -- -- plot T amp
                stack_arrayT_Amp = np.zeros((slowT_n,stack_nt))
                for it in range(stack_nt):  # check points one at a time
                    for slowT_i in range(slowT_n):  # for this station, loop over slownesses
                        num_val = Tcentral_am[slowT_i].data[it]
                        stack_arrayT_Amp[slowT_i, it] = num_val

                y, x = np.mgrid[slice(stack_Tslows[0], stack_Tslows[-1] + slow_delta, slow_delta),
                             slice(ttt[0], ttt[-1] + dt, dt)]

                fig_index += 1
                fig, ax = plt.subplots(1, figsize=(10,3), num = fig_index)
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
                    plt.title(f'Amp at {target_slow:.3f} s/km N slowness, {fname[48:58]} #{eq_num}')
                else:
                    plt.ylabel('Transverse Slowness (s/km)')
                    plt.title(f'Amp at {target_slow:.3f} s/km R slowness, {fname[48:58]} #{eq_num}')
                plt.show()

    #%% 2-slices-plus-snaps option
    if two_slice_plots == True:
    #%% -- R-T amplitude snap plots
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

        #%% -- Extract and sum amp slices from beam matrix

        # Select only stacks with that slowness for Transverse plot
        centralR_Ast = Stream()
        for slowR_i in range(slowR_n):
            centralR_Ast += amp_ave[slowR_i*slowT_n + lowest_Tindex]

        # Select only stacks with that slowness for Radial plot
        centralT_Ast = Stream()
        for slowT_i in range(slowT_n):
            centralT_Ast += amp_ave[lowest_Rindex*slowT_n + slowT_i]

        #%% -- Amp slice plots
        #%% -- -- Radial slice plot
        stack_array = np.zeros((slowR_n,stack_nt))

        for it in range(stack_nt):  # check points one at a time
            for slowR_i in range(slowR_n):  # for this station, loop over slownesses
                num_val = centralR_Ast[slowR_i].data[it]
                if log_plot == True:
                    stack_array[slowR_i, it] = num_val - global_max
                else:
                    stack_array[slowR_i, it] = num_val

        y, x = np.meshgrid(stack_Rslows,ttt)
        # if x.shape[1] != len(ttt):  # does array come out too long in time steps?
        #     y, x = np.mgrid[slice(stack_Rslows[0], stack_Rslows[-1] + slow_delta, slow_delta),
        #                  slice(ttt[0], ttt[-1] + (dt/1.5), dt)]

        fig_index += 1
        fig, ax = plt.subplots(1, figsize=(10,3), num = fig_index)
        print(f'len(x) is {len(x)} and len(y) is {len(y)}')
        print(f'x1 is {x.shape[0]} and x2 is {x.shape[1]}')
        print(f'len(stack_Rslows) is {len(stack_Rslows)} and len(ttt) is {len(ttt)}')
        print(f'slowR_n is {slowR_n} and stack_nt is {stack_nt}')
        if log_plot == True:
            c = ax.pcolormesh(x, y, np.transpose(stack_array), cmap=plt.cm.gist_rainbow_r, vmin = -log_plot_range, vmax=0)
            fig.colorbar(c, ax=ax, label='log amplitude')
        else:
            c = ax.pcolormesh(x, y, np.transpose(stack_array), cmap=plt.cm.gist_rainbow_r, vmin=0)
            fig.colorbar(c, ax=ax, label='Linear amplitude')
        fig.subplots_adjust(bottom=0.2)
        ax.axis([x.min(), x.max(), y.min(), y.max()])
        c = ax.scatter(arrival_time, event_pred_slo, color='black'  , s=50, alpha=0.75)
        plt.xlabel('Time (s)')
        if NS == True:
            plt.ylabel('North Slowness (s/km)')
            plt.title('Amp at ' + str(T_slow_plot) + ' s/km East slowness, ' + date_label + ' #' + str(eq_num))
        else:
            plt.ylabel('Radial slowness (s/km)')
            plt.title('Amp at ' + str(T_slow_plot) + ' s/km transverse slowness, ' + date_label + ' #' + str(eq_num))
        plt.show()

        #%% -- -- Transverse slice plot

        stack_array = np.zeros((slowT_n,stack_nt))

        for it in range(stack_nt):  # check points one at a time
            for slowT_i in range(slowT_n):  # for this station, loop over slownesses
                num_val = centralT_Ast[slowT_i].data[it]
                if log_plot == True:
                    stack_array[slowT_i, it] = num_val - global_max
                else:
                    stack_array[slowT_i, it] = num_val

        y, x = np.meshgrid(stack_Tslows,ttt)

        fig_index += 1
        fig, ax = plt.subplots(1, figsize=(10,3), num = fig_index)
        if log_plot == True:
            c = ax.pcolormesh(x, y, np.transpose(stack_array), cmap=plt.cm.gist_rainbow_r, vmin = -log_plot_range, vmax=0)
            fig.colorbar(c, ax=ax, label='log amplitude')
        else:
            c = ax.pcolormesh(x, y, np.transpose(stack_array), cmap=plt.cm.gist_rainbow_r, vmin=0)
            fig.colorbar(c, ax=ax, label='Linear amplitude')
        fig.subplots_adjust(bottom=0.2)
        ax.axis([x.min(), x.max(), y.min(), y.max()])
        c = ax.scatter(arrival_time, 0, color='black'  , s=50, alpha=0.75)
        plt.xlabel('Time (s)')
        if NS == True:
            plt.ylabel('East Slowness (s/km)')
            plt.title('Amp at ' + str(R_slow_plot) + ' s/km North slowness, ' + date_label + ' #' + str(eq_num))
        else:
            plt.ylabel('Transverse slowness (s/km)')
            plt.title('Amp at ' + str(R_slow_plot) + ' s/km radial slowness, ' + date_label + ' #' + str(eq_num))
        plt.show()

        #%% -- Snap plots
        stack_slice = np.zeros((slowR_n,slowT_n))
        if snaps > 0 :
            # check for impossible parameters
            if snaptime < start_buff:
                print(f'snaptime {snaptime:.0f} is earlier than start_buff of {start_buff:.0f}')
                sys.exit(-1)
            last_snap = snaptime + snaps*snap_depth
            if last_snap > end_buff:
                print(f'last snap {last_snap:.0f} is later than end_buff of {end_buff:.0f}')
                sys.exit(-1)

            for snap_num in range(snaps):
                snap_start = snaptime + (snap_num  ) * snap_depth
                snap_end   = snaptime + (snap_num+1) * snap_depth
                it_start = int(round((snap_start - start_buff)/dt))
                it_end   = int(round((snap_end   - start_buff)/dt))
                print(f'snap start {snap_start:.1f} snap end {snap_end:.1f} it start {it_start:.1f} it end {it_end:.1f} snap num {snap_num:.1f} start buff {start_buff:.1f}')
                for slowR_i in range(slowR_n):  # loop over radial slownesses
                    for slowT_i in range(slowT_n):  # loop over transverse slownesses
                        index = slowR_i*slowT_n + slowT_i
                        num_val = np.nanmean(amp_ave[index].data[it_start:it_end])
                        stack_slice[slowR_i, slowT_i] = num_val

                y1, x1 = np.mgrid[slice(stack_Rslows[0], stack_Rslows[-1] + slow_delta, slow_delta),
                             slice(stack_Tslows[0], stack_Tslows[-1] + slow_delta, slow_delta)]

                fig_index += 1
                fig, ax = plt.subplots(1, figsize=(7,0.8*7), num = fig_index)
                c = ax.pcolormesh(x1, y1, stack_slice, cmap=plt.cm.gist_rainbow_r)
                ax.axis([x1.min(), x1.max(), y1.min(), y1.max()])
                if log_plot == True:
                    fig.colorbar(c, ax=ax, label = 'log amp')
                else:
                    fig.colorbar(c, ax=ax, label = 'linear amp')
                c = ax.scatter(pred_Eslo, pred_Nslo, color='purple'  , s=50, alpha=0.75)
                c = ax.scatter(        0,         0, color='black' , s=50,  alpha=0.75)
                circle1 = plt.Circle((0, 0), 0.019, color='black', fill=False)
                ax.add_artist(circle1)  # inner core limit
                circle2 = plt.Circle((0, 0), 0.040, color='black', fill=False)
                ax.add_artist(circle2)  # outer core limit
                plt.title(f'Amplitude from time {snap_start:.1f} to {snap_end:.1f}s  {date_label}  event #{eq_num}')
                if NS == True:
                    plt.xlabel('E Slowness (s/km)')
                    plt.ylabel('N Slowness (s/km)')
                else:
                    plt.xlabel('T Slowness (s/km)')
                    plt.ylabel('R Slowness (s/km)')
                plt.show()

    #%% Beam sum plots
    if beam_sums == True or max_wiggly_plot == True:
    #%% -- R-T amplitude averaged over time window
        stack_slice = np.zeros((slowR_n,slowT_n))
        if start_beam == 0 and end_beam == 0: # entire volume is summed
            full_beam = 1
        else:  # beam just part of stack volume
            full_beam = 0
            start_index = int((start_beam - start_buff) / dt)
            end_index   = int((end_beam   - start_buff) / dt)
            print('beam is ' + str(start_beam) + ' to ' + str(end_beam) + 's, out of ' + str(start_buff)
                + ' to ' + str(end_buff) + 's, dt is ' + str(dt)  + 's, and indices are '+ str(start_index) + ' ' + str(end_index))
            print(f'Beam is {start_beam:.4f} to {end_beam:.4f}s, out of {start_buff:.4f} to {end_buff:.4f}s, dt is {dt:.4f}s, and indices are {start_index} {end_index}')

        for slowR_i in range(slowR_n):  # loop over radial slownesses
            for slowT_i in range(slowT_n):  # loop over transverse slownesses
                index = slowR_i*slowT_n + slowT_i
                if full_beam == 1:
                    num_val = np.nanmean(amp_ave[index].data)
                else:
                    num_val = np.nanmean(amp_ave[index].data[start_index:end_index])
                stack_slice[slowR_i, slowT_i] = num_val

        y1, x1 = np.meshgrid(stack_Rslows,stack_Tslows)

        # y1, x1 = np.mgrid[slice(stack_Rslows[0], stack_Rslows[-1] + slow_delta, slow_delta),
        #              slice(stack_Tslows[0], stack_Tslows[-1] + slow_delta, slow_delta)]

        fig_index += 1
        fig, ax = plt.subplots(1, figsize=(7,0.8*7*(slowR_n/slowT_n)), num = fig_index)
        smax = np.max(stack_slice)
        smin = np.min(stack_slice)

        max_xy = np.where(stack_slice == stack_slice.max() ) # find indices of slowness with max amplitude, needed for max_wiggle_plot
        print('Beam max - (x,y):' + str(max_xy[0]) + '  ' + str(max_xy[1]) + ' out of nR, nT slownesses ' + str(slowR_n) + '  ' + str(slowR_n))
        Rslow_max = int(max_xy[0]) * slow_delta + slowR_lo
        Tslow_max = int(max_xy[1]) * slow_delta + slowT_lo

        if log_plot == True:
            if (smax - smin) < log_plot_range:  # use full color scale even if range is less than specified
                log_plot_range = smax - smin
            c = ax.pcolormesh(x1, y1, np.transpose(stack_slice - smax), cmap=plt.cm.gist_rainbow_r, vmin= - log_plot_range, vmax=0)
        else:
            c = ax.pcolormesh(x1, y1, np.transpose(stack_slice/smax), cmap=plt.cm.gist_rainbow_r, vmin = 0)
        if log_plot == True:
            fig.colorbar(c, ax=ax, label='log amplitude')
        else:
            fig.colorbar(c, ax=ax, label='linear amplitude')
        ax.axis([x1.min(), x1.max(), y1.min(), y1.max()])
        circle1 = plt.Circle((0, 0), 0.019, color='black', fill=False)
        ax.add_artist(circle1)  #inner core limit
        circle2 = plt.Circle((0, 0), 0.040, color='black', fill=False)
        ax.add_artist(circle2)  #outer core limit

        c = ax.scatter(pred_Eslo, pred_Nslo, color='purple' , s=50, alpha=1) # predicted slowness
        c = ax.scatter(Tslow_max, Rslow_max, color='white' , s=30, alpha=1)  # max amp observed
        c = ax.scatter(        0,         0, color='black' , s=50, alpha=1)  # (0,0)

        if NS == True:
            plt.xlabel('East Slowness (s/km)')
            plt.ylabel('North Slowness (s/km)')
        else:
            plt.xlabel('Transverse Slowness (s/km)')
            plt.ylabel('Radial Slowness (s/km)')
        plt.title(f'{date_label} {ref_phase} {start_buff:.1f} to {end_buff:.1f} beam amp #{eq_num}')
        os.chdir('/Users/vidale/Documents/Research/IC/Plots_hold')
        # plt.savefig(date_label + '_' + str(start_buff) + '_' + str(end_buff) + '_beam.png')
        plt.savefig(f'{eq_num:02}_{date_label}_{start_buff:.0f}-{end_buff:.0f}_beam')
        plt.show()

    #%% Wiggly plots
    if wiggly_plots == True or max_wiggly_plot == True:

        #%% -- read wiggle beams
        # Get saved event info, also used to name files
        # date_label = '2018-04-02' # date for filename
        goto = '/Users/vidale/Documents/Research/IC/Pro_files'
        os.chdir(goto)
        fname = 'HD' + date_label + '_2dstack.mseed'
        st = Stream()
        st = read(fname)

        #%% -- Extract slices to wiggle plot, cumbersome, every slowness along slice is now selected
        #%% -- -- Collect T slowness nearest T_slow
        lowest_Tslow = 1000000
        for slow_i in range(slowT_n):
            if abs(stack_Tslows[slow_i] - T_slow_plot) < lowest_Tslow:
                lowest_Tindex = slow_i
                lowest_Tslow = abs(stack_Tslows[slow_i] - T_slow_plot)

        print(f'{slowT_n} T slownesses, index {lowest_Tindex} is closest to requested plot T slowness {T_slow_plot:.4f}, slowness diff there is {lowest_Tslow:.4f} and slowness is {stack_Tslows[lowest_Tindex]:.4f}')
        # Select only stacks with that slowness for radial plot
        centralR_st    = Stream()
        centralR_amp   = Stream()
        for slowR_i in range(slowR_n):
            ii = slowR_i*slowT_n + lowest_Tindex
            centralR_st  += st[ii]
            centralR_amp += amp_ave[ii]

        #%% -- -- Collect R slowness nearest R_slow
        lowest_Rslow = 1000000
        for slow_i in range(slowR_n):
            if abs(stack_Rslows[slow_i] - R_slow_plot) < lowest_Rslow:
                lowest_Rindex = slow_i
                lowest_Rslow = abs(stack_Rslows[slow_i] - R_slow_plot)

        print(f'{slowR_n} R slownesses, index {lowest_Rindex} is closest to requested plot R slowness {R_slow_plot:.4f}, slowness diff there is {lowest_Rslow:.4f} and slowness is {stack_Rslows[lowest_Rindex]:.4f}')

        # Select only stacks with that slowness for transverse plot
        centralT_st = Stream()
        centralT_amp   = Stream()

        for slowT_i in range(slowT_n):
            ii = lowest_Rindex*slowT_n + slowT_i
            centralT_st += st[ii]
            centralT_amp   += amp_ave[ii]

    #%% -- Plot wiggles
    if wiggly_plots == True:
        #%% -- -- R amp and tdiff vs time plots with black line for time shift
        scale_plot_wig = wig_scale_fac / (200 * global_max)
        if log_plot == True:
            scale_plot_wig /= 30  # not quite sure why this renormalization works
            # scale_plot_tdiff = plot_scale_fac / 500.
        fig_tit = str(eq_num) + '_Ramp'
        plt.figure(fig_tit,figsize=(30,10))
        plt.xlim(start_buff,end_buff)
        plt.ylim(stack_Rslows[0], stack_Rslows[-1])
        for slowR_i in range(slowR_n):  # loop over radial slownesses
            dist_offset = stack_Rslows[slowR_i] # trying for approx degrees
            ttt1 = (np.arange(len(centralR_st[slowR_i].data)) * centralR_st[slowR_i].stats.delta
              + (centralR_st[slowR_i].stats.starttime - t))
            plt.plot(ttt1, ((centralR_st[slowR_i].data - np.median(centralR_st[slowR_i].data)) * scale_plot_wig) + dist_offset, color = 'black')
            if turn_off_black == 0:
                plt.plot(ttt1,     (centralT_st[slowT_i].data)*0.0 + dist_offset, color = 'gray') # reference lines

        plt.xlabel('Time (s)')
        if NS == True:
            plt.ylabel('N Slowness (s/km)')
            plt.title(ref_phase + '  ' + date_label + '  ' + ' seismograms ' + str(T_slow_plot) + ' E slowness' + ' #' + str(eq_num))
        else:
            plt.ylabel('R Slowness (s/km)')
            plt.title(ref_phase + '  ' + date_label + '  ' + ' seismograms ' + str(T_slow_plot) + ' T slowness' + ' #' + str(eq_num))
        #%% -- -- T amp and tdiff vs time plots with black line for time shift
        fig_tit = str(eq_num) + '_Tamp'
        plt.figure(fig_tit,figsize=(30,10))
        plt.xlim(start_buff,end_buff)
        plt.ylim(stack_Tslows[0], stack_Tslows[-1])

        for slowT_i in range(slowT_n):  # loop over transverse slownesses
            dist_offset = stack_Tslows[slowT_i] # trying for approx degrees
            ttt2 = (np.arange(len(centralT_st[slowT_i].data)) * centralT_st[slowT_i].stats.delta
              + (centralT_st[slowT_i].stats.starttime - t))
            plt.plot(ttt2, ((centralT_st[slowT_i].data - np.median(centralT_st[slowT_i].data)) * scale_plot_wig) + dist_offset, color = 'black')
            if turn_off_black == 0:
                plt.plot(ttt2,     (centralT_st[slowT_i].data)*0.0 + dist_offset, color = 'gray') # reference lines
        plt.xlabel('Time (s)')
        if NS == True:
            plt.ylabel('E Slowness (s/km)')
            plt.title(ref_phase + '  ' + date_label + '  ' + ' seismograms ' + str(R_slow_plot) + ' N slowness' + ' #' + str(eq_num))
        else:
            plt.ylabel('T Slowness (s/km)')
            plt.title(ref_phase + '  ' + date_label + '  ' + ' seismograms ' + str(R_slow_plot) + ' R slowness' + ' #' + str(eq_num))
        os.chdir('/Users/vidale/Documents/Research/IC/Plots_hold')
    #    plt.savefig(date_label1 + '_' + str(start_buff) + '_' + str(end_buff) + '_stack.png')

    #%% Wiggly plot at max amplitude
    if max_wiggly_plot == True:
        # have already read wiggle beams
        # if beam_sums == False or wiggly_plots == False:
        #     print(colored('Need beam_sums and wiggly_plots True to make max_wiggly_plot', 'magenta'))
        # else:
        #%% -- Select stack at max amp
        max_trace = Stream()
        # iii = 84*slowT_n + 81
        iii = int(max_xy[0])*slowT_n + int(max_xy[1])
        max_trace = st[iii]

        max_trace.data /= max(abs(max_trace.data))
        if log_plot == True:
            scale_plot_wig /= 30  # not quite sure why this renormalization works
        # fig_tit = str(eq_num) + '_maxamp'
        fig_index += 1
        plt.figure(figsize=(10,5), num = fig_index)
        plt.xlim(orig_start_buff,orig_end_buff)
        plt.ylim(-1, 1)
        ttt = (np.arange(len(max_trace)) * centralR_st[0].stats.delta
          + (centralR_st[0].stats.starttime - t))
        print('Length of ttt and max_trace:  ' + str(len(ttt)) + '  ' +  str(len(max_trace)))
        plt.plot(ttt, max_trace, color = 'black')

        plt.xlabel('Time (s)')
        plt.ylabel('Amplitude')

        plt.title(f'{ref_phase} {date_label}  max amp stack in event {eq_num}  Rslow  {Rslow_max:.4f}  Tslow {Tslow_max:.4f}')
        # plt.title(date_label + '  ' + ' max amp stack in event ' + str(eq_num) + ' Rslow ' + str(Rslow) + ' Tslow ' + str(Tslow))
        os.chdir('/Users/vidale/Documents/Research/IC/Plots_hold')
    #    plt.savefig(date_label1 + '_' + str(start_buff) + '_' + str(end_buff) + '_stack.png')

    #  Save processed files
#    fname = 'HD' + date_label + '_slice.mseed'
#    stack.write(fname,format = 'MSEED')

    elapsed_time_wc = time.time() - start_time_wc
    print(f'This job took {elapsed_time_wc:.1f} seconds')
    os.system('say "Seven"')