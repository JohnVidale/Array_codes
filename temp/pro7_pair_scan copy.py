#!/usr/bin/env python
# looks at envelope stack differential parameters from pro6
# Reads in tdiff, ave_amp, cc computed from a pair of events
# window by signal quality
# eleven_slice == T produces sections at 4 radial slownesses (0.0, 0.05, 0.01, 0.015)
# plus plots sections at 7 transv slownesses (-0.015 to 0.015)
# eleven_slice == F plots one radial and one transverse plot, plus snaps
# John Vidale 2/2019, overhauled 1/2021

def pro7_pair_scan(eq_file1, eq_file2, slow_delta = 0.0005,
              slowR_lo = -0.1, slowR_hi = 0.1, slowT_lo = -0.1, slowT_hi = 0.1,
              start_buff = 50, end_buff = 50,fig_index = 401, skip_T = 0, skip_R = 0,
              ZslowR_lo = -0.1, ZslowR_hi = 0.1, ZslowT_lo = -0.1, ZslowT_hi = 0.1,
              Zstart_buff = 50, Zend_buff = 50, zoom = 0, tdiff_clip = 1,
              ref_phase = 'blank', cc_thres = 0.8, min_amp = 0.2, auto_slice = True,
              R_slow_plot = 0.06, T_slow_plot = 0.0, snaptime = 8, snaps = 10, skip_snaps = 0,
              nR_plots  = 3, nT_plots = 3, slow_incr = 0.01, NS = False, dphase = 'PKiKP',
              ARRAY = 0):

    from obspy import read
    from obspy.taup import TauPyModel
    import numpy as np
    import os
    import matplotlib.pyplot as plt
    import time
    import sys
    from obspy import UTCDateTime
    from obspy import Stream
    from termcolor import colored
    from obspy.geodetics import gps2dist_azimuth
    model = TauPyModel(model='iasp91')

    print(colored('Running pro7b_plot_stack', 'cyan'))
    start_time_wc = time.time()

    if zoom == 1:
        if Zstart_buff  < start_buff:
            print(f'Zstart_buff of {Zstart_buff:.1f} cannot be < start_buff of {start_buff:.1f}')
            Zstart_buff = start_buff
        if Zend_buff    > end_buff:
            print(f'Zend_buff of {Zend_buff:.1f} cannot be > end_buff of {end_buff:.1f}')
            Zend_buff   = end_buff

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
    t1          = UTCDateTime(split_line1[1])
    # t2 = UTCDateTime(split_line2[1])
    date_label1  = split_line1[1][0:10]
    date_label2  = split_line2[1][0:10]
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
    ref1_dist  = ref_distance[0]/(1000*111)
    arrivals_ref   = model.get_travel_times(source_depth_in_km=ev_depth,distance_in_degree=ref1_dist, phase_list=[dphase])
    arrival_time = arrivals_ref[0].time

    arrivals1 = model.get_travel_times(source_depth_in_km=ev_depth,distance_in_degree=ref1_dist-0.5,phase_list=[dphase])
    arrivals2 = model.get_travel_times(source_depth_in_km=ev_depth,distance_in_degree=ref1_dist+0.5,phase_list=[dphase])
    dtime = arrivals2[0].time - arrivals1[0].time
    event_pred_slo  = dtime/111.  # s/km

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
    print('Read in: ')
    print('Input trace starttime ' + str(tdiff[0].stats.starttime))
    print(str(len(tdiff)) + '  ' + str(len(amp_ave)) + '  ' + str(len(cc)) + ' traces for tdiff, amp_ave, cc')
    print(str(len(tdiff[0].data)) + '  ' + str(len(amp_ave[0].data)) + '  ' + str(len(cc[0].data)) + ' time pts for tdiff, amp_ave, cc')
    print(str(tdiff[0].stats.delta) + '  ' + str(amp_ave[0].stats.delta) + '  ' + str(cc[0].stats.delta) + ' dt for tdiff, amp_ave, cc')

    dt = tdiff[0].stats.delta
    nt = len(tdiff[0].data)

    #%% Make grid of slownesses
    slowR_n = int(round(1 + (slowR_hi - slowR_lo)/slow_delta))  # number of slownesses
    slowT_n = int(round(1 + (slowT_hi - slowT_lo)/slow_delta))  # number of slownesses
    stack_nt = int(round(1 + (end_buff - start_buff)/dt))  # number of time points
    print(str(slowT_n) + ' trans slownesses, hi and lo are ' + str(slowT_hi) + '  ' + str(slowT_lo))
    print(str(slowR_n) + ' radial slownesses, hi and lo are ' + str(slowR_hi) + '  ' + str(slowR_lo))
    # In English, stack_slows = range(slow_n) * slow_delta - slow_lo
    a1R = range(slowR_n)
    a1T = range(slowT_n)
    stack_Rslows = [(x * slow_delta + slowR_lo) for x in a1R]
    stack_Tslows = [(x * slow_delta + slowT_lo) for x in a1T]

    #%% Select subset if Zoomed
    if zoom == 1:
        Ztdiff   = Stream()
        Zamp_ave = Stream()
        Zcc      = Stream()
        for slowR_i in range(slowR_n):  # loop over radial slownesses
            for slowT_i in range(slowT_n):  # loop over transverse slownesses, kludge to evade rounding error
                if ((stack_Rslows[slowR_i] >= ZslowR_lo - 0.000001) and (stack_Rslows[slowR_i] <= ZslowR_hi + 0.000001) and
                    (stack_Tslows[slowT_i] >= ZslowT_lo - 0.000001) and (stack_Tslows[slowT_i] <= ZslowT_hi + 0.000001)):
                    index = slowR_i*slowT_n + slowT_i
                    s_t = t1 - Zstart_buff
                    e_t = t1 + Zend_buff
                    Ztdiff   += tdiff[  index].trim(starttime=s_t, endtime=e_t)
                    Zamp_ave += amp_ave[index].trim(starttime=s_t, endtime=e_t)
                    Zcc      += cc[     index].trim(starttime=s_t, endtime=e_t)
                            #tr.trim(starttime=s_t,endtime = e_t)
        tdiff   = Ztdiff
        amp_ave = Zamp_ave
        cc      = Zcc
        nt = len(tdiff[0].data)

        #%% -- Re-make finer grid of slownesses
        slowR_lo   = ZslowR_lo
        slowR_hi   = ZslowR_hi
        slowT_lo   = ZslowT_lo
        slowT_hi   = ZslowT_hi
        end_buff   = Zend_buff
        start_buff = Zstart_buff
        slowR_n = int(1 + (slowR_hi - slowR_lo)/slow_delta)  # number of slownesses
        slowT_n = int(1 + (slowT_hi - slowT_lo)/slow_delta)  # number of slownesses
        stack_nt = int(1 + ((end_buff - start_buff)/dt))  # number of time points
        print('After zoom ' + str(slowT_n) + ' trans slownesses, hi and lo are ' + str(slowT_hi) + '  ' + str(slowT_lo) + ' stack_nt is ' + str(stack_nt))
        # In English, stack_slows = range(slow_n) * slow_delta - slow_lo
        a1R = range(slowR_n)
        a1T = range(slowT_n)
        stack_Rslows = [(x * slow_delta + slowR_lo) for x in a1R]
        stack_Tslows = [(x * slow_delta + slowT_lo) for x in a1T]
        print('After zoom ' + str(slowR_n) + ' radial slownesses, ' + str(slowT_n) + ' trans slownesses, ')
        print('Output trace starttime ' + str(Ztdiff[0].stats.starttime))

    #%% Mask out weak and/or less correlated points
    global_max = 0
    for slow_i in range(len(amp_ave)): # find global max of ave_amp
        local_max = max(amp_ave[slow_i].data)
        if local_max > global_max:
            global_max = local_max

    for slow_i in range(len(tdiff)): # ignore less robust points
        for it in range(nt):
            if (cc[slow_i].data[it] < cc_thres) or (amp_ave[slow_i].data[it] < (min_amp * global_max)):
                tdiff[slow_i].data[it] = np.nan
#%% auto slice option
    if auto_slice == True:

#%% -- compute timing time series
        ttt = (np.arange(len(tdiff[0].data)) * tdiff[0].stats.delta + start_buff) # in units of seconds

        #%% -- R slices
        if skip_R != 1:  # remember plots scanning R are those at constant T
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

                stack_array = np.zeros((slowR_n,stack_nt))

                #%% -- -- plot R tdiff
                for it in range(stack_nt):  # check points one at a time
                    for slowR_i in range(slowR_n):  # for this station, loop over slownesses
                        num_val = Rcentral_st[slowR_i].data[it]
                        stack_array[slowR_i, it] = num_val

                y, x = np.mgrid[slice(stack_Rslows[0], stack_Rslows[-1] + slow_delta, slow_delta),
                             slice(ttt[0], ttt[-1] + dt, dt)]

                fig, ax = plt.subplots(1, figsize=(10,3))
                c = ax.pcolormesh(x, y, stack_array, cmap=plt.cm.coolwarm, vmin= -tdiff_clip, vmax=tdiff_clip)
                fig.subplots_adjust(bottom=0.2)
                ax.axis([x.min(), x.max(), y.min(), y.max()])
                fig.colorbar(c, ax=ax)
                c = ax.scatter(arrival_time, event_pred_slo, color='black'  , s=50, alpha=0.75)
                plt.xlabel('Time (s)')
                plt.ylabel('Radial slowness (s/km)')
                plt.title(f'{ref_phase} Time lag at {target_slow:.3f} s/km T slowness, {fname1[48:58]}  {fname1[59:69]}  min amp {min_amp:.3f}  cc_thres {cc_thres:.3f}')
                plt.show()
                fig_index += 1

                #%% -- -- plot R amp
                for it in range(stack_nt):  # check points one at a time
                    for slowR_i in range(slowR_n):  # for this station, loop over slownesses
                        num_val = Rcentral_am[slowR_i].data[it]
                        stack_array[slowR_i, it] = num_val

                y, x = np.mgrid[slice(stack_Rslows[0], stack_Rslows[-1] + slow_delta, slow_delta),
                             slice(ttt[0], ttt[-1] + dt, dt)]

                fig, ax = plt.subplots(1, figsize=(10,3))
                c = ax.pcolormesh(x, y, stack_array, cmap=plt.cm.gist_rainbow_r, vmin= 0, vmax=global_max)
                fig.subplots_adjust(bottom=0.2)
                ax.axis([x.min(), x.max(), y.min(), y.max()])
                fig.colorbar(c, ax=ax)
                c = ax.scatter(arrival_time, event_pred_slo, color='black'  , s=50, alpha=0.75)
                plt.xlabel('Time (s)')
                plt.ylabel('Radial slowness (s/km)')
                plt.title(f'{ref_phase} Amp at {target_slow:.3f} s/km T slowness, {fname1[48:58]}  {fname1[59:69]}')
                plt.show()
                fig_index += 1

        #%% -- T slices
        if skip_T != 1:  # remember plots scanning T are those at constant R
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

                stack_array = np.zeros((slowT_n,stack_nt))

                #%% -- -- plot T tdiff
                for it in range(stack_nt):  # check points one at a time
                    for slowT_i in range(slowT_n):  # for this station, loop over slownesses
                        num_val = Tcentral_st[slowT_i].data[it]
                        stack_array[slowT_i, it] = num_val

                y, x = np.mgrid[slice(stack_Tslows[0], stack_Tslows[-1] + slow_delta, slow_delta),
                             slice(ttt[0], ttt[-1] + dt, dt)]

                fig, ax = plt.subplots(1, figsize=(10,3))
                c = ax.pcolormesh(x, y, stack_array, cmap=plt.cm.coolwarm, vmin= -tdiff_clip, vmax=tdiff_clip)
                fig.subplots_adjust(bottom=0.2)
                ax.axis([x.min(), x.max(), y.min(), y.max()])
                fig.colorbar(c, ax=ax)
                c = ax.scatter(arrival_time, 0, color='black'  , s=50, alpha=0.75)
                plt.xlabel('Time (s)')
                plt.ylabel('Transverse slowness (s/km)')
                plt.title(f'{ref_phase} Time lag at {target_slow:.3f} s/km R slowness, {fname1[48:58]}  {fname1[59:69]}  min amp {min_amp:.3f}  cc_thres {cc_thres:.3f}')
                plt.show()
                fig_index += 1

                #%% -- -- plot T amp
                for it in range(stack_nt):  # check points one at a time
                    for slowT_i in range(slowT_n):  # for this station, loop over slownesses
                        num_val = Tcentral_am[slowT_i].data[it]
                        stack_array[slowT_i, it] = num_val

                y, x = np.mgrid[slice(stack_Tslows[0], stack_Tslows[-1] + slow_delta, slow_delta),
                             slice(ttt[0], ttt[-1] + dt, dt)]

                fig, ax = plt.subplots(1, figsize=(10,3))
                c = ax.pcolormesh(x, y, stack_array, cmap=plt.cm.gist_rainbow_r, vmin= 0, vmax=global_max)
                fig.subplots_adjust(bottom=0.2)
                ax.axis([x.min(), x.max(), y.min(), y.max()])
                fig.colorbar(c, ax=ax)
                c = ax.scatter(arrival_time, 0, color='black'  , s=50, alpha=0.75)
                plt.xlabel('Time (s)')
                plt.ylabel('Transverse slowness (s/km)')
                plt.title(f'{ref_phase} Amp at {target_slow:.3f} s/km R slowness, {fname1[48:58]}  {fname1[59:69]}')
                plt.show()
                fig_index += 1

    #%% 2-slices-plus-snaps option
    else:
        #%% -- R Slice near T slowness T_slow
        if skip_R != 1:
            lowest_Tslow = 1000000
            for slow_i in range(slowT_n):
                if abs(stack_Tslows[slow_i] - T_slow_plot) < lowest_Tslow:
                    lowest_Tindex = slow_i
                    lowest_Tslow = abs(stack_Tslows[slow_i] - T_slow_plot)

            print(f'{lowest_Tindex:4d} is T slow nearest {T_slow_plot:.3f}, difference is {lowest_Tslow:.3f}')

            # Select only stacks with that slowness for Transverse plot
            centralR_st = Stream()
            for slowR_i in range(slowR_n):
                centralR_st += tdiff[slowR_i*slowT_n + lowest_Tindex]

        #%% -- T Slice near R slowness R_slow
        if skip_T != 1:
            lowest_Rslow = 1000000
            for slow_i in range(slowR_n):
                if abs(stack_Rslows[slow_i] - R_slow_plot) < lowest_Rslow:
                    lowest_Rindex = slow_i
                    lowest_Rslow = abs(stack_Rslows[slow_i] - R_slow_plot)

            print(f'{lowest_Rindex:4d} is R slow nearest 0.005, difference is {lowest_Rslow:.3f}')

            # Select only stacks with that slowness for Radial plot
            centralT_st = Stream()
            for slowT_i in range(slowT_n):
                centralT_st += tdiff[lowest_Rindex*slowT_n + slowT_i]

        ttt = (np.arange(len(tdiff[0].data)) * tdiff[0].stats.delta + start_buff) # in units of seconds

        #%% -- R-time stack plot
        if skip_R != 1:
            stack_array = np.zeros((slowR_n,stack_nt))

            for it in range(stack_nt):  # check points one at a time
                for slowR_i in range(slowR_n):  # for this station, loop over slownesses
                    num_val = centralR_st[slowR_i].data[it]
                    stack_array[slowR_i, it] = num_val

            y, x = np.mgrid[slice(stack_Rslows[0], stack_Rslows[-1] + slow_delta, slow_delta),
                         slice(ttt[0], ttt[-1] + dt, dt)]

            fig, ax = plt.subplots(1, figsize=(10,3))
            c = ax.pcolormesh(x, y, stack_array, cmap=plt.cm.coolwarm, vmin=-tdiff_clip, vmax=tdiff_clip)
            fig.subplots_adjust(bottom=0.2)
            ax.axis([x.min(), x.max(), y.min(), y.max()])
            fig.colorbar(c, ax=ax)
            plt.xlabel('Time (s)')
            plt.ylabel('Radial slowness (s/km)')
            plt.title(ref_phase + ' Time lag at ' + str(T_slow_plot) + ' s/km transverse slowness, ' + fname1[12:22] + ' ' + fname1[23:33])
            plt.show()

            fig_index += 1
        #%% -- T-time stack plot
        if skip_T != 1:
            stack_array = np.zeros((slowT_n,stack_nt))

            for it in range(stack_nt):  # check points one at a time
                for slowT_i in range(slowT_n):  # for this station, loop over slownesses
                    num_val = centralT_st[slowT_i].data[it]
                    stack_array[slowT_i, it] = num_val

            y, x = np.mgrid[slice(stack_Tslows[0], stack_Tslows[-1] + slow_delta, slow_delta),
                         slice(ttt[0], ttt[-1] + dt, dt)]

            fig, ax = plt.subplots(1, figsize=(10,3))
            fig.subplots_adjust(bottom=0.2)
            c = ax.pcolormesh(x, y, stack_array, cmap=plt.cm.coolwarm, vmin=-tdiff_clip, vmax=tdiff_clip)
            ax.axis([x.min(), x.max(), y.min(), y.max()])
            fig.colorbar(c, ax=ax)
            plt.xlabel('Time (s)')
            plt.ylabel('Transverse slowness (s/km)')
            plt.title(ref_phase + ' Time lag at ' + str(R_slow_plot) + ' s/km radial slowness, ' + fname1[12:22] + ' ' + fname1[23:33])
            plt.show()

            fig_index += 1

    #%% -- R-T stack time difference snap plots
        if skip_snaps == 0:
            stack_slice = np.zeros((slowR_n,slowT_n))
            for snap_num in range(snaps):
                fig_index += 1
                it = int((snaptime - start_buff)/dt) + snap_num
                for slowR_i in range(slowR_n):  # loop over radial slownesses
                    for slowT_i in range(slowT_n):  # loop over transverse slownesses
                        index = slowR_i*slowT_n + slowT_i
                        num_val = tdiff[index].data[it]
                        stack_slice[slowR_i, slowT_i] = num_val

                y1, x1 = np.mgrid[slice(stack_Rslows[0], stack_Rslows[-1] + slow_delta, slow_delta),
                             slice(stack_Tslows[0], stack_Tslows[-1] + slow_delta, slow_delta)]

                fig, ax = plt.subplots(1)
                c = ax.pcolormesh(x1, y1, stack_slice, cmap=plt.cm.gist_rainbow_r, vmin=-tdiff_clip, vmax=tdiff_clip)
                ax.axis([x1.min(), x1.max(), y1.min(), y1.max()])
                fig.colorbar(c, ax=ax)
                circle1 = plt.Circle((0, 0), 0.019, color='black', fill=False)
                ax.add_artist(circle1)
                plt.xlabel('T Slowness (s/km)')
                plt.ylabel('R Slowness (s/km)')
                plt.title(ref_phase + ' T-R plot of time lag at rel time ' + str(snaptime + snap_num*dt) + '  ' + fname1[12:22] + ' ' + fname1[23:33])
                plt.show()

    #%% -- R-T stack amplitude snap plots
        if skip_snaps == 0:
            stack_slice = np.zeros((slowR_n,slowT_n))
            for snap_num in range(snaps):
                fig_index += 1
                it = int((snaptime - start_buff)/dt) + snap_num
                for slowR_i in range(slowR_n):  # loop over radial slownesses
                    for slowT_i in range(slowT_n):  # loop over transverse slownesses
                        index = slowR_i*slowT_n + slowT_i
                        num_val = amp_ave[index].data[it]
                        stack_slice[slowR_i, slowT_i] = num_val

                y1, x1 = np.mgrid[slice(stack_Rslows[0], stack_Rslows[-1] + slow_delta, slow_delta),
                             slice(stack_Tslows[0], stack_Tslows[-1] + slow_delta, slow_delta)]

                fig, ax = plt.subplots(1)
                c = ax.pcolormesh(x1, y1, stack_slice, cmap=plt.cm.gist_rainbow_r)
                ax.axis([x1.min(), x1.max(), y1.min(), y1.max()])
                fig.colorbar(c, ax=ax)
                circle1 = plt.Circle((0, 0), 0.019, color='black', fill=False)
                ax.add_artist(circle1)
                plt.xlabel('T Slowness (s/km)')
                plt.ylabel('R Slowness (s/km)')
                plt.title(ref_phase + ' T-R plot of amplitude at rel time ' + str(snaptime + snap_num*dt) + '  ' + fname1[12:22] + ' ' + fname1[23:33])
                plt.show()

    #  Save processed files
#    fname = 'HD' + date_label + '_slice.mseed'
#    stack.write(fname,format = 'MSEED')

    elapsed_time_wc = time.time() - start_time_wc
    print('This job took ' + str(elapsed_time_wc) + ' seconds')
    os.system('say "Done"')