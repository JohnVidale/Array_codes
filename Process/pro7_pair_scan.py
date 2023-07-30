#!/usr/bin/env python
# inputs phase stack comparison for a doublet from pro6
# Reads in tdiff, ave_amp, cc computed from a pair of events
# window by signal quality
# John Vidale 2/2019, overhauled 1/2021, reviewed 6/2021

def pro7_pair_scan(eq_num1, eq_num2, repeater = '0', slow_delta = 0.0005, turn_off_black = 1,
              slowR_lo = -0.1, slowR_hi = 0.1, slowT_lo = -0.1, slowT_hi = 0.1,
              start_buff = 50, end_buff = 50,fig_index = 401, do_T = False, do_R = False,
              ZslowR_lo = -0.1, ZslowR_hi = 0.1, ZslowT_lo = -0.1, ZslowT_hi = 0.1,
              Zstart_buff = 50, Zend_buff = 50, zoom = False, tdiff_clip = 1,
              phase1 = 'PKiKP', phase2 = 'no', phase3 = 'no', phase4 = 'no',
              cc_thres = 0.8, min_amp = 0.2, cc_twin = 1, cc_len = 1,
              R_slow_plot = 0.06, T_slow_plot = 0.0,
              snaptime = 8, snaps = 10, snap_depth = 0, freq_min = 0, freq_max = 0,
              nR_plots  = 3, nT_plots = 3, slow_incr = 0.01, NS = False,
              ARRAY = 0, auto_slice = True, two_slice_plots = False, beam_sums = True,
              wiggly_plots = True, start_beam = 0, end_beam = 0, log_plot = False,
              log_plot_range = 2, tdiff_plots_too = False, pred_wiggles = True,
              wig_scale_fac = 1, tdiff_scale_fac = 1, do_trans = False,
              pair_name = '', plot_peak = 1):

    from obspy import read
    from obspy.taup import TauPyModel
    import numpy as np
    import os
    import matplotlib.pyplot as plt
    import time
    import sys
    import math
    import copy
    from obspy import UTCDateTime
    from obspy import Stream
    from termcolor import colored
    from obspy.geodetics import gps2dist_azimuth
    model = TauPyModel(model='ak135')

    print(colored('Running pro7_pair_scan', 'cyan'))
    tdiff_clip = cc_twin * cc_len
    phasePKP_single = False
    phasePKP_double = False
    # if phase1 == 'PKP':
        # print(colored('code not configured to handle double PKP entries as phase1', 'yellow'))
        # sys.exit(-1)
    if phase2 == 'PKP':
        phasePKP_double = True
        phase2 = 'no'
    if phase3 == 'PKP':
        phasePKP_double = True
        phase3 = 'no'
    if phase4 == 'PKP':
        phasePKP_double = True
        phase4 = 'no'

    if ARRAY == 0:
        arrayname = 'HiNet '
    elif ARRAY == 1:
        arrayname = 'LASA '
    elif ARRAY == 2:
        arrayname = 'China '
    elif ARRAY == 3:
        arrayname = 'NORSAR '
    elif ARRAY == 4:
        arrayname = 'WRA '
    elif ARRAY == 5:
        arrayname = 'YKA '
    elif ARRAY == 6:
        arrayname = 'ILAR '

    start_time_wc = time.time()
    beam_env_plot   = True
    max_wiggly_plot = False

    IC_beam = False
    beam_stack_rad = 0.01
    folder_name = '/Users/vidale/Documents/Research/IC/'

    if zoom == True:
        if Zstart_buff  < start_buff:
            print(colored(f'Zstart_buff of {Zstart_buff:.1f} cannot be < start_buff of {start_buff:.1f}', 'red'))
            Zstart_buff = start_buff
            exit()
        if Zend_buff    > end_buff:
            print(colored(f'Zend_buff of {Zend_buff:.1f} cannot be < end_buff of {end_buff:.1f}', 'red'))
            Zend_buff   = end_buff
            exit()

    #%% Input parameters and computed files
    fname1 = '/Users/vidale/Documents/Research/IC/EvLocs/event' + str(eq_num1) + '.txt'
    fname2 = '/Users/vidale/Documents/Research/IC/EvLocs/event' + str(eq_num2) + '.txt'
    file1 = open(fname1, 'r')
    file2 = open(fname2, 'r')
    lines1=file1.readlines()
    lines2=file2.readlines()
    split_line1  = lines1[0].split()
    split_line2  = lines2[0].split()
    t1           = UTCDateTime(split_line1[1])
    # t2 = UTCDateTime(split_line2[1])  # not needed
    # date_label = '2018-04-02' # dates in filename
    date_label1  = split_line1[1][0:10]
    date_label2  = split_line2[1][0:10]
    save_name = '/Users/vidale/Documents/Research/IC/Plots_hold/' + repeater + '_Array_' + str(ARRAY)
    ev_lat       = float(split_line1[2])
    ev_lon       = float(split_line1[3])
    ev_depth     = float(split_line1[4])

    if ARRAY == 0:
        ref_lat =   36.30  # °N, around middle of Japan
        ref_lon =  138.50 # °E
    elif ARRAY == 1:
        ref_lat =   46.70  # °N keep only inner rings A-D
        ref_lon = -106.22   # °E
    elif ARRAY == 2:
        ref_lat =   38.00  # °N China
        ref_lon =  104.50  # °E
    elif ARRAY == 4:
        ref_lat =  -19.89  # °N Warramunga
        ref_lon =  134.42  # °E
    elif ARRAY == 5:
        ref_lat =   62.49  # °N Yellowknife
        ref_lon = -114.60  # °E
    elif ARRAY == 6:
        ref_lat =   64.77  # °N ILAR
        ref_lon = -146.89  # °E

    ref_distance = gps2dist_azimuth(ref_lat,ref_lon,ev_lat,ev_lon)
    ref_dist     = ref_distance[0]/(1000*111)
    ref_az       = ref_distance[1]
    ref_back_az  = ref_distance[2]

    # Estimate slowness of reference phases
    arrivals_ref   = model.get_travel_times(source_depth_in_km=ev_depth,distance_in_degree=ref_dist, phase_list=[phase1])
    print(str(len(arrivals_ref)))
    if len(arrivals_ref) == 0:
        print('reference phase ' + phase1 + ' does not exist at distance ' + str(ref_distance))
        exit(-1)
    if len(arrivals_ref) == 2:
        print('reference phase ' + phase1 + ' has two arrivals at distance ' + str(ref_distance) + ', may not work')
        exit(-1)
    arrival_time = arrivals_ref[0].time
    arrival_time1 = 0 # for
    atime_rayp = arrivals_ref[0].ray_param
    event_pred_slo  = atime_rayp * 2 * np.pi / (111. * 360.) # convert to s/km

    if phase2 != 'no':
        arrivals_ref2   = model.get_travel_times(source_depth_in_km=ev_depth,distance_in_degree=ref_dist, phase_list=[phase2])
        if len(arrivals_ref2) == 0:
            print('phase2 ' + phase2 + ' does not exist at distance ' + str(ref_distance[0]/(111. * 1000)))
            phase2 = 'no'
        else:
            arrival_time2 = arrivals_ref2[0].time - arrival_time
            atime_rayp2 = arrivals_ref2[0].ray_param
            event_pred_slo2  = atime_rayp2 * 2 * np.pi / (111. * 360.) # convert to s/km

    if phase3 != 'no':
        arrivals_ref3   = model.get_travel_times(source_depth_in_km=ev_depth,distance_in_degree=ref_dist, phase_list=[phase3])
        if len(arrivals_ref3) == 0:
            print('phase3 ' + phase3 + ' does not exist at distance ' + str(ref_distance[0]/(111. * 1000)))
            phase3 = 'no'
        else:
            arrival_time3 = arrivals_ref3[0].time - arrival_time
            atime_rayp3 = arrivals_ref3[0].ray_param
            event_pred_slo3  = atime_rayp3 * 2 * np.pi / (111. * 360.) # convert to s/km

    if phase4 != 'no':
        arrivals_ref4   = model.get_travel_times(source_depth_in_km=ev_depth,distance_in_degree=ref_dist, phase_list=[phase4])
        if len(arrivals_ref4) == 0:
            print('phase4 ' + phase4 + ' does not exist at distance ' + str(ref_distance[0]/(111. * 1000)))
            phase4 = 'no'
        else:
            arrival_time4 = arrivals_ref4[0].time - arrival_time
            atime_rayp4 = arrivals_ref4[0].ray_param
            event_pred_slo4  = atime_rayp4 * 2 * np.pi / (111. * 360.) # convert to s/km

    if phasePKP_double:  # 0, 1, or 2 PKP arrivals from TauP
        arrivals_refPKP   = model.get_travel_times(source_depth_in_km=ev_depth,distance_in_degree=ref_dist, phase_list=['PKP'])
        PKP_count = len(arrivals_refPKP)
        if PKP_count == 0:
            print('phase PKP does not compute at distance ' + str(ref_distance[0]/(111. * 1000)))
            phasePKP_double = False
        elif PKP_count == 1:
            print('phase PKP produced just 1 arrival at distance ' + str(ref_distance[0]/(111. * 1000)))
            phasePKP_double = False
            phasePKP_single = True
            arrival_timePKP1 = arrivals_refPKP[0].time - arrival_time
            atime_raypPKP1 = arrivals_refPKP[0].ray_param
            event_pred_sloPKP1  = atime_raypPKP1 * 2 * np.pi / (111. * 360.) # convert to s/km
        elif PKP_count == 2:
            arrival_timePKP1 = arrivals_refPKP[0].time - arrival_time
            atime_raypPKP1 = arrivals_refPKP[0].ray_param
            event_pred_sloPKP1  = atime_raypPKP1 * 2 * np.pi / (111. * 360.) # convert to s/km
            arrival_timePKP2 = arrivals_refPKP[1].time - arrival_time
            atime_raypPKP2 = arrivals_refPKP[1].ray_param
            event_pred_sloPKP2  = atime_raypPKP2 * 2 * np.pi / (111. * 360.) # convert to s/km
        else:
            print(colored('PKP_count is ' + str(PKP_count) + ', too high!?', 'yellow'))
            sys.exit(-1)

    # convert to pred rslo and tslo
    if NS:    #  rotate predicted slowness to N and E
        print(f'Array  lat {ref_lat:.0f}, lon  {ref_lon:.0f}, Event lat {ev_lat:.0f}, lon {ev_lon:.0f}, az {ref_az:.0f}, baz {ref_back_az:.0f}')
        sin_baz = np.sin(ref_az * np.pi /180)
        cos_baz = np.cos(ref_az * np.pi /180)
        pred_Nslo = event_pred_slo * cos_baz
        pred_Eslo = event_pred_slo * sin_baz
        if phase2 != 'no':
            pred_Nslo2 = event_pred_slo2 * cos_baz
            pred_Eslo2 = event_pred_slo2 * sin_baz
        if phase3 != 'no':
            pred_Nslo3 = event_pred_slo3 * cos_baz
            pred_Eslo3 = event_pred_slo3 * sin_baz
        if phase4 != 'no':
            pred_Nslo4 = event_pred_slo4 * cos_baz
            pred_Eslo4 = event_pred_slo4 * sin_baz
        if phasePKP_single or phasePKP_double:
            pred_NsloPKP1 = event_pred_sloPKP1 * cos_baz
            pred_EsloPKP1 = event_pred_sloPKP1 * sin_baz
        if phasePKP_double:
            pred_NsloPKP2 = event_pred_sloPKP2 * cos_baz
            pred_EsloPKP2 = event_pred_sloPKP2 * sin_baz
    else:
        pred_Nslo  = event_pred_slo
        pred_Eslo  = 0
        if phase2 != 'no':
            pred_Nslo2 = event_pred_slo2
            pred_Eslo2 = 0
        if phase3 != 'no':
            pred_Nslo3 = event_pred_slo3
            pred_Eslo3 = 0
        if phase4 != 'no':
            pred_Nslo4 = event_pred_slo4
            pred_Eslo4 = 0
        if phasePKP_single or phasePKP_double:
            pred_NsloPKP1 = event_pred_sloPKP1
            pred_EsloPKP1 = 0
        if phasePKP_double:
            pred_NsloPKP2 = event_pred_sloPKP2
            pred_EsloPKP2 = 0

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
    print(f'input grids: tdiff {len(tdiff)} elements, cc {len(cc)} elements, amp_ave {len(amp_ave)} elements')

    #%% Make grid of slownesses
    slowR_n = int(round(1 + (slowR_hi - slowR_lo)/slow_delta))  # number of slownesses
    slowT_n = int(round(1 + (slowT_hi - slowT_lo)/slow_delta))  # number of slownesses
    stack_nt = int(round(1 + (end_buff - start_buff)/dt))  # number of time points
    print(f'{slowR_n} radial slownesses, low is {slowR_lo}, high is {slowR_hi}')
    print(f'{slowT_n} transv slownesses, low is {slowT_lo}, high is {slowT_hi}')
    print(f'{stack_nt} time points, low is {start_buff} s, high is {end_buff} s, dt is {dt}')
    # In English, stack_slows = range(slow_n) * slow_delta - slow_lo
    a1R = range(slowR_n)
    a1T = range(slowT_n)
    stack_Rslows = [(x * slow_delta + slowR_lo) for x in a1R]
    stack_Tslows = [(x * slow_delta + slowT_lo) for x in a1T]

    #%% plot beam envelope time series

    if beam_env_plot == True:

        beam_env_trace = Stream()
        beam_env_trace = amp_ave[0].copy() #just chosen to get right dimension and dt
        beam_env_trace.data *= 0

        for slowR_i in range(slowR_n):  # loop over radial slownesses
            for slowT_i in range(slowT_n):  # loop over transverse slownesses
                index = slowR_i*slowT_n + slowT_i
                slowR_actual = stack_Rslows[slowR_i]
                slowT_actual = stack_Tslows[slowT_i]
                if IC_beam:
                    slow_amp = np.sqrt((slowR_actual * slowR_actual) + (slowT_actual * slowT_actual))
                    if slow_amp < 0.019:
                        beam_env_trace.data += amp_ave[index].data
                        # print(f'included slowness Rindex {slowR_i} Tindex {slowT_i}  Rslow  {slowR_actual:.4f}  Tslow {slowT_actual:.4f} slow_amp {slow_amp:.4f}')
                else:
                    slow_anomaly = np.sqrt(((slowR_actual - pred_Nslo) * (slowR_actual - pred_Nslo)) +
                                           ((slowT_actual - pred_Eslo) * (slowT_actual - pred_Eslo)))
                    # print(f'included Rslow  {slowR_actual:.4f}  Tslow {slowT_actual:.4f} pred_Nslo  {pred_Nslo:.4f}  pred_Eslo {pred_Eslo:.4f} slow_anomaly {slow_anomaly:.4f}')
                    # if slow_anomaly < 0.01:
                    if slow_anomaly < beam_stack_rad:
                        beam_env_trace.data += amp_ave[index].data
                        # print(f'included Rslow  {slowR_actual:.4f}  Tslow {slowT_actual:.4f} pred_Nslo  {pred_Nslo:.4f}  pred_Eslo {pred_Eslo:.4f} slow_anomaly {slow_anomaly:.4f}')

        beam_env_trace.data /= max(abs(beam_env_trace.data))
        fig_index += 1
        plt.figure(figsize=(10,5), num = fig_index)
        plt.xlim(start_buff,end_buff)
        plt.ylim(0, 1.2)
        ttt = (np.arange(stack_nt)*dt + start_buff)

        # normalize in zoom window
        questor = (ttt >= Zstart_buff) & (ttt < Zend_buff) # identify zoom window
        ts_sel = beam_env_trace[questor]  #extract zoom window
        max_env = max(ts_sel)
        beam_env_trace.data = beam_env_trace.data/(max_env)

        print('diff ' + str(start_buff) + ' stack_nt ' + str(stack_nt))

        print('Length of ttt and beam envelope:  ' + str(len(ttt)) + '  ' +  str(len(beam_env_trace)))
        plt.plot(ttt, beam_env_trace, color = 'black')
        plt.plot((Zstart_buff, Zstart_buff), (0, 1), color = 'red')
        plt.text(Zstart_buff, 0, 'start', color = 'black')
        plt.plot((Zend_buff,     Zend_buff), (0, 1), color = 'red')
        plt.text(Zend_buff  , 0,   'end', color = 'black')
        plt.plot((0,0), (0, 1), color = 'black')
        plt.text(arrival_time1, 0.9, phase1, color = 'black')
        if phase2 != 'no':
            plt.plot((arrival_time2, arrival_time2), (0, 1), color = 'gray')
            plt.text(arrival_time2, 1, phase2, color = 'black')
        if phase3 != 'no':
            plt.text(arrival_time3, 1, phase3, color = 'black')
            plt.plot((arrival_time3, arrival_time3), (0, 1), color = 'gray')
        if phase4 != 'no':
            plt.text(arrival_time4, 1, phase4, color = 'black')
            plt.plot((arrival_time4, arrival_time4), (0, 1), color = 'gray')
        if phasePKP_single or phasePKP_double:
            plt.text(arrival_timePKP1, 0.95, 'PKP1', color = 'black')
            plt.plot((arrival_timePKP1, arrival_timePKP1), (0, 1), color = 'gray')
        if phasePKP_double:
            plt.text(arrival_timePKP2, 0.95, 'PKP2', color = 'black')
            plt.plot((arrival_timePKP2, arrival_timePKP2), (0, 1), color = 'gray')

        plt.xlabel('Time (s)')
        plt.ylabel('Amplitude')
        print('beam_stack_rad is ' + str(beam_stack_rad))

        if IC_beam:
            plt.title(f'{date_label1} and {date_label2} IC stack in events {eq_num1} and {eq_num2} sum inside {beam_stack_rad} s/° , {freq_min}-{freq_max} Hz')
        else:
            plt.title(f'{pair_name} {arrayname} {eq_num1} and {eq_num2} beam env stack, sum within {beam_stack_rad} s/° of prediction, {freq_min}-{freq_max} Hz')
        os.chdir('/Users/vidale/Documents/Research/IC/Plots_hold')
        plt.savefig(save_name + '_envelope.png')

    #%% Select subset if Zoomed
    orig_start_buff = start_buff
    orig_end_buff   = end_buff
    if zoom == True:
        Ztdiff   = Stream()
        Zamp_ave = Stream()
        Zcc      = Stream()
        print(f'before calculation, tdiff[0] has length {len(tdiff[0])}')
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
        tdiff   = Ztdiff  # tdiff might be one element shorter than expected
        amp_ave = Zamp_ave
        cc      = Zcc
        nt = len(tdiff[0].data)
        start_buff = Zstart_buff
        # make time series
        print(f'after calculation, Ztdiff[0] has length {len(Ztdiff[0])}')
        print(f'after calculation, tdiff[0] has length {len(tdiff[0])}')
        print(f'slowR_lo  is {slowR_lo:.4f}  and slowR_hi  is {slowR_hi:.4f}  and slowT_lo  is {slowT_lo:.4f}  and slowT_hi is {slowT_hi:.4f}')
        print(f'ZslowR_lo is {ZslowR_lo:.4f} and ZslowR_hi is {ZslowR_hi:.4f} and ZslowT_lo is {ZslowT_lo:.4f} and ZslowT_hi is {ZslowT_hi:.4f}')

        #%% -- Re-make subset with more limited grid of slownesses and times
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
        if stack_nt != len(Ztdiff[0]):
            print(f'Array length rounding clash, tdiff[0] has length {len(tdiff[0])}, stack_nt is {stack_nt}')
            stack_nt = len(tdiff[0])
        print('After zoom ' + str(slowT_n) + ' trans slownesses, hi and lo are ' + str(slowT_hi) + '  ' + str(slowT_lo) + ' stack_nt is ' + str(stack_nt))
        # In English, stack_slows = range(slow_n) * slow_delta - slow_lo
        a1R = range(slowR_n)
        a1T = range(slowT_n)
        stack_Rslows = [(x * slow_delta + slowR_lo) for x in a1R]
        stack_Tslows = [(x * slow_delta + slowT_lo) for x in a1T]
        print('After zoom ' + str(slowR_n) + ' radial slownesses, ' + str(slowT_n) + ' trans slownesses.')
        print('Output trace starttime ' + str(Ztdiff[0].stats.starttime))

    ttt = (np.arange(len(tdiff[0].data)) * tdiff[0].stats.delta + start_buff) # in units of seconds

    global_max = 0  # find global_max, largest amplitude in amp_ave beam array envelopes
    for slow_i in range(len(amp_ave)): # find global max of ave_amp
        local_max = max(amp_ave[slow_i].data)
        if local_max > global_max:
            global_max = local_max
    if global_max == 0:
        print(colored('Global max is ' + str(global_max), 'yellow'))
        sys.exit(-1)

    #%% Mask out weak and/or less correlated points
    amp_ave_thres    = amp_ave.copy()  # copy amp envelope array, set amps and tdiff below thresholds to NaN using global_max
    nt = len(tdiff[0].data)
    for slow_i in range(len(tdiff)): # don't plot less robust points, change them to NANs
        for it in range(nt):
            if (cc[slow_i].data[it] < cc_thres) or (amp_ave[slow_i].data[it] < (min_amp * global_max)):
                tdiff[        slow_i].data[it] = np.nan
                amp_ave_thres[slow_i].data[it] = np.nan

    for slow_i in range(len(tdiff)): # set NaNs to avoid including (errant?) large time shifts
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
                if tdiff_plots_too == True:
                    stack_arrayR_Tdf = np.zeros((slowR_n,stack_nt))
                    for it in range(stack_nt):  # check points one at a time
                        for slowR_i in range(slowR_n):  # loop over slownesses
                            num_val = Rcentral_st[slowR_i].data[it]
                            stack_arrayR_Tdf[slowR_i, it] = num_val

                    y, x = np.meshgrid(stack_Rslows,ttt)
                    fig, ax = plt.subplots(1, figsize=(10,3))
                    c = ax.pcolormesh(x, y, np.transpose(stack_arrayR_Tdf), cmap=plt.cm.coolwarm, vmin= -tdiff_clip, vmax=tdiff_clip)
                    fig.subplots_adjust(bottom=0.2)
                    ax.axis([x.min(), x.max(), y.min(), y.max()])
                    fig.colorbar(c, ax=ax, label='time shift (s)')
                    print(str(len(arrival_time1)))
                    c = ax.scatter(arrival_time1, event_pred_slo, color='black'  , s=50, alpha=0.75)
                    if phase2 != 'no':
                        c = ax.scatter(arrival_time2, event_pred_slo2, color='gray'  , s=50, alpha=0.75)
                    if phase3 != 'no':
                        c = ax.scatter(arrival_time3, event_pred_slo3, color='gray'  , s=50, alpha=0.75)
                    if phase4 != 'no':
                        c = ax.scatter(arrival_time4, event_pred_slo4, color='gray'  , s=50, alpha=0.75)
                    plt.xlabel('Time (s)')
                    plt.ylabel('Radial slowness (s/km)')
                    plt.title(f'Tdiff at {target_slow:.3f} s/km T slowness, {fname1[48:58]}  {fname1[59:69]}  min amp {min_amp:.1f}  cc_thres {cc_thres:.2f}')
                    fig_index += 1

                #%% -- -- plot R amp
                stack_arrayR_Amp = np.zeros((slowR_n,stack_nt))
                for it in range(stack_nt):  # check points one at a time
                    for slowR_i in range(slowR_n):  # for this station, loop over slownesses
                        num_val = Rcentral_am[slowR_i].data[it]
                        stack_arrayR_Amp[slowR_i, it] = num_val

                y, x = np.meshgrid(stack_Rslows,ttt)
                fig, ax = plt.subplots(1, figsize=(10,3))
                if log_plot == True:
                    c = ax.pcolormesh(x, y, np.transpose(stack_arrayR_Amp - global_max), cmap=plt.cm.gist_rainbow_r, vmin= - log_plot_range, vmax=0)
                else:
                    c = ax.pcolormesh(x, y, np.transpose(stack_arrayR_Amp), cmap=plt.cm.gist_rainbow_r, vmin= 0, vmax=global_max)
                fig.subplots_adjust(bottom=0.2)
                ax.axis([x.min(), x.max(), y.min(), y.max()])
                if log_plot == True:
                    fig.colorbar(c, ax=ax, label='log amplitude')
                else:
                    fig.colorbar(c, ax=ax, label='linear amplitude')
                c = ax.scatter(arrival_time, event_pred_slo, color='black'  , s=50, alpha=0.75)
                plt.xlabel('Time (s)')
                if NS:
                    plt.ylabel('North Slowness (s/km)')
                    plt.title(f'{pair_name}  {arrayname}  Amp at {target_slow:.3f} s/km E slowness, {fname1[48:58]}  {fname1[59:69]}')
                else:
                    plt.ylabel('Radial Slowness (s/km)')
                    plt.title(f'{pair_name}  {arrayname}  Amp at {target_slow:.3f} s/km T slowness, {fname1[48:58]}  {fname1[59:69]}')
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
                lowest_Rslow = 1000000  # find index of row closest to R_cnt slowness
                target_slow = (R_cnt * slow_incr) # radial slowness of this slice in s/°
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
                if tdiff_plots_too == True:
                    stack_arrayT_Tdf = np.zeros((slowT_n,stack_nt))
                    for it in range(stack_nt):  # check points one at a time
                        for slowT_i in range(slowT_n):  # for this station, loop over slownesses
                            num_val = Tcentral_st[slowT_i].data[it]
                            stack_arrayT_Tdf[slowT_i, it] = num_val

                    y, x = np.meshgrid(stack_Tslows,ttt)
                    fig, ax = plt.subplots(1, figsize=(10,3))
                    c = ax.pcolormesh(x, y, np.transpose(stack_arrayT_Tdf), cmap=plt.cm.coolwarm, vmin= -tdiff_clip, vmax=tdiff_clip)
                    fig.subplots_adjust(bottom=0.2)
                    ax.axis([x.min(), x.max(), y.min(), y.max()])
                    fig.colorbar(c, ax=ax, label='time shift (s)')
                    c = ax.scatter(arrival_time, 0, color='black'  , s=50, alpha=0.75)
                    plt.xlabel('Time (s)')
                    plt.ylabel('Transverse slowness (s/km)')
                    plt.title(f'{phase1} Tdiff at {target_slow:.3f} s/km R slowness, {fname1[48:58]}  {fname1[59:69]}  min amp {min_amp:.1f}  cc_thres {cc_thres:.2f}')
                    fig_index += 1

                #%% -- -- plot T amp
                stack_arrayT_Amp = np.zeros((slowT_n,stack_nt))
                for it in range(stack_nt):  # check points one at a time
                    for slowT_i in range(slowT_n):  # for this station, loop over slownesses
                        num_val = Tcentral_am[slowT_i].data[it]
                        stack_arrayT_Amp[slowT_i, it] = num_val

                y, x = np.meshgrid(stack_Tslows,ttt)
                fig, ax = plt.subplots(1, figsize=(10,3))
                if log_plot == True:
                    c = ax.pcolormesh(x, y, np.transpose(stack_arrayT_Amp - global_max), cmap=plt.cm.gist_rainbow_r, vmin= - log_plot_range, vmax=0)
                else:
                    c = ax.pcolormesh(x, y, np.transpose(stack_arrayT_Amp), cmap=plt.cm.gist_rainbow_r, vmin= 0, vmax=global_max)
                fig.subplots_adjust(bottom=0.2)
                ax.axis([x.min(), x.max(), y.min(), y.max()])
                if log_plot == True:
                    fig.colorbar(c, ax=ax, label='log amplitude')
                else:
                    fig.colorbar(c, ax=ax, label='linear amplitude')
                c = ax.scatter(arrival_time, 0, color='black'  , s=50, alpha=0.75)
                plt.xlabel('Time (s)')
                if NS:
                    plt.ylabel('East Slowness (s/km)')
                    plt.title(f'{pair_name}  {arrayname}  Amp at {target_slow:.3f} s/km N slowness, {fname1[48:58]}  {fname1[59:69]}')
                else:
                    plt.ylabel('Transverse Slowness (s/km)')
                    plt.title(f'{pair_name}  {arrayname}  Amp at {target_slow:.3f} s/km R slowness, {fname1[48:58]}  {fname1[59:69]}')
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
        centralR_Ast = Stream()
        centralR_cc  = Stream()
        for slowR_i in range(slowR_n):
            centralR_Dst +=   tdiff[slowR_i*slowT_n + lowest_Tindex]
            centralR_Ast += amp_ave[slowR_i*slowT_n + lowest_Tindex]
            centralR_cc  +=      cc[slowR_i*slowT_n + lowest_Tindex]

        # Select only stacks with that slowness for Radial plot
        centralT_Dst = Stream()
        centralT_Ast = Stream()
        centralT_cc  = Stream()
        for slowT_i in range(slowT_n):
            centralT_Dst +=   tdiff[lowest_Rindex*slowT_n + slowT_i]
            centralT_Ast += amp_ave[lowest_Rindex*slowT_n + slowT_i]
            centralT_cc  +=      cc[lowest_Rindex*slowT_n + slowT_i]

        #%% -- Primary stack plots
        #%% -- -- Time lag plots
        #%% -- -- -- R/N plot
        stack_array = np.zeros((slowR_n,stack_nt))

        for it in range(stack_nt):  # check points one at a time
            for slowR_i in range(slowR_n):  # for this station, loop over slownesses
                num_val = centralR_Dst[slowR_i].data[it]
                stack_array[slowR_i, it] = num_val

        y, x = np.meshgrid(stack_Rslows,ttt)
        fig, ax = plt.subplots(1, figsize=(10,3))
        print(f'len(x) is {len(x)} and len(y) is {len(y)}')
        print(f'len(stack_Rslows) is {len(stack_Rslows)} and len(ttt) is {len(ttt)}')
        print(f'slowR_n is {slowR_n} and stack_nt is {stack_nt}')
        c = ax.pcolormesh(x, y, np.transpose(stack_array), cmap=plt.cm.coolwarm, vmin=-tdiff_clip, vmax=tdiff_clip)
        fig.colorbar(c, ax=ax, label='time lag (s)')
        fig.subplots_adjust(bottom=0.2)
        ax.axis([x.min(), x.max(), y.min(), y.max()])

        c = ax.scatter(arrival_time1, pred_Nslo, color='black'  , s=50, alpha=0.75)
        plt.text(arrival_time1, pred_Nslo, phase1, color = 'black')
        if phase2 != 'no':
            c = ax.scatter(arrival_time2, pred_Nslo2, color='gray'  , s=50, alpha=0.75)
            plt.text(arrival_time2, pred_Nslo2, phase2, color = 'black')
        if phase3 != 'no':
            c = ax.scatter(arrival_time3, pred_Nslo3, color='gray'  , s=50, alpha=0.75)
            plt.text(arrival_time3, pred_Nslo3, phase3, color = 'black')
        if phase4 != 'no':
            c = ax.scatter(arrival_time4, pred_Nslo4, color='gray'  , s=50, alpha=0.75)
            plt.text(arrival_time4, pred_Nslo4, phase4, color = 'black')
        if phasePKP_single or phasePKP_double:
            c = ax.scatter(arrival_timePKP1, pred_NsloPKP1, color='purple'  , s=50, alpha=0.75)
            plt.text(arrival_timePKP1, pred_NsloPKP1, 'PKP1', color = 'black')
        if phasePKP_double:
            c = ax.scatter(arrival_timePKP2, pred_NsloPKP2, color='purple'  , s=50, alpha=0.75)
            plt.text(arrival_timePKP2, pred_NsloPKP2, 'PKP2', color = 'black')

        plt.xlabel('Time (s)')
        if NS:
            plt.ylabel('N slowness (s/km)')
            plt.title(pair_name + ' ' + arrayname + ' ' + ' Time lag at ' + str(T_slow_plot) + ' s/km E slowness, ' + date_label1 + ' ' + date_label2)
            plt.savefig(save_name + '_Ntdiff_hist.png')
        else:
            plt.ylabel('R slowness (s/km)')
            plt.title(pair_name + ' ' + arrayname + ' ' + ' Time lag at ' + str(T_slow_plot) + ' s/km T slowness, ' + date_label1 + ' ' + date_label2)
            plt.savefig(save_name + '_Rtdiff_hist.png')
        fig_index += 1

        #%% -- -- -- T/E plot
        if do_trans:
            stack_array = np.zeros((slowT_n,stack_nt))

            for it in range(stack_nt):  # check points one at a time
                for slowT_i in range(slowT_n):  # for this station, loop over slownesses
                    num_val = centralT_Dst[slowT_i].data[it]
                    stack_array[slowT_i, it] = num_val

            y, x = np.meshgrid(stack_Tslows,ttt)
            fig, ax = plt.subplots(1, figsize=(10,3))
            fig.subplots_adjust(bottom=0.2)
            # c = ax.pcolormesh(x, y, np.transpose(stack_array), cmap=plt.cm.coolwarm, vmin=-tdiff_clip, vmax=tdiff_clip)
            c = ax.pcolormesh(x, y, np.transpose(stack_array), cmap=plt.cm.coolwarm)
            fig.colorbar(c, ax=ax, label='time lag (s)')
            ax.axis([x.min(), x.max(), y.min(), y.max()])

            c = ax.scatter(arrival_time1, pred_Eslo, color='black'  , s=50, alpha=0.75)
            plt.text(arrival_time1, pred_Eslo, phase1, color = 'black')
            if phase2 != 'no':
                c = ax.scatter(arrival_time2, pred_Eslo2, color='gray'  , s=50, alpha=0.75)
                plt.text(arrival_time2, pred_Eslo2, phase2, color = 'black')
            if phase3 != 'no':
                c = ax.scatter(arrival_time3, pred_Eslo3, color='gray'  , s=50, alpha=0.75)
                plt.text(arrival_time3, pred_Eslo3, phase3, color = 'black')
            if phase4 != 'no':
                c = ax.scatter(arrival_time4, pred_Eslo4, color='gray'  , s=50, alpha=0.75)
                plt.text(arrival_time4, pred_Eslo4, phase4, color = 'black')
            if phasePKP_single or phasePKP_double:
                c = ax.scatter(arrival_timePKP1, pred_EsloPKP1, color='purple'  , s=50, alpha=0.75)
                plt.text(arrival_timePKP1, pred_EsloPKP1, 'PKP1', color = 'black')
            if phasePKP_double:
                c = ax.scatter(arrival_timePKP2, pred_EsloPKP2, color='purple'  , s=50, alpha=0.75)
                plt.text(arrival_timePKP2, pred_EsloPKP2, 'PKP2', color = 'black')

            plt.xlabel('Time (s)')
            if NS:
                plt.ylabel('E slowness (s/km)')
                plt.title(pair_name + ' ' + arrayname + ' ' + ' Time lag at ' + str(R_slow_plot) + ' s/km N slowness, ' + date_label1 + ' ' + date_label2)
                plt.savefig(save_name + '_Etdiff_hist.png')
            else:
                plt.ylabel('T slowness (s/km)')
                plt.title(pair_name + ' ' + arrayname + ' ' + ' Time lag at ' + str(R_slow_plot) + ' s/km R slowness, ' + date_label1 + ' ' + date_label2)
                plt.savefig(save_name + '_Ttdiff_hist.png')
            fig_index += 1

        #%% -- -- Correlation plots
        #%% -- -- -- R/N plot
        stack_array = np.zeros((slowR_n,stack_nt))

        for it in range(stack_nt):  # check points one at a time
            for slowR_i in range(slowR_n):  # for this station, loop over slownesses
                num_val = centralR_cc[slowR_i].data[it]
                stack_array[slowR_i, it] = num_val

        y, x = np.meshgrid(stack_Rslows,ttt)
        fig, ax = plt.subplots(1, figsize=(10,3))
        print(f'len(x) is {len(x)} and len(y) is {len(y)}')
        print(f'len(stack_Rslows) is {len(stack_Rslows)} and len(ttt) is {len(ttt)}')
        print(f'slowR_n is {slowR_n} and stack_nt is {stack_nt}')
        c = ax.pcolormesh(x, y, np.transpose(stack_array), cmap=plt.cm.coolwarm, vmin=0, vmax=1)
        fig.colorbar(c, ax=ax, label='Correlation')
        fig.subplots_adjust(bottom=0.2)
        ax.axis([x.min(), x.max(), y.min(), y.max()])

        c = ax.scatter(arrival_time1, pred_Nslo, color='black'  , s=50, alpha=0.75)
        plt.text(arrival_time1, pred_Nslo, phase1, color = 'black')
        if phase2 != 'no':
            c = ax.scatter(arrival_time2, pred_Nslo2, color='gray'  , s=50, alpha=0.75)
            plt.text(arrival_time2, pred_Nslo2, phase2, color = 'black')
        if phase3 != 'no':
            c = ax.scatter(arrival_time3, pred_Nslo3, color='gray'  , s=50, alpha=0.75)
            plt.text(arrival_time3, pred_Nslo3, phase3, color = 'black')
        if phase4 != 'no':
            c = ax.scatter(arrival_time4, pred_Nslo4, color='gray'  , s=50, alpha=0.75)
            plt.text(arrival_time4, pred_Nslo4, phase4, color = 'black')
        if phasePKP_single or phasePKP_double:
            c = ax.scatter(arrival_timePKP1, pred_NsloPKP1, color='purple'  , s=50, alpha=0.75)
            plt.text(arrival_timePKP1, pred_NsloPKP1, 'PKP1', color = 'black')
        if phasePKP_double:
            c = ax.scatter(arrival_timePKP2, pred_NsloPKP2, color='purple'  , s=50, alpha=0.75)
            plt.text(arrival_timePKP2, pred_NsloPKP2, 'PKP2', color = 'black')

        plt.xlabel('Time (s)')
        if NS:
            plt.ylabel('N slowness (s/km)')
            plt.title(pair_name + ' ' + arrayname + ' ' + ' Correlation at ' + str(T_slow_plot) + ' E slowness, ' + date_label1 + ' ' + date_label2)
            plt.savefig(save_name + '_Ncc_hist.png')
        else:
            plt.ylabel('R slowness (s/km)')
            plt.title(pair_name + ' ' + arrayname + ' ' + ' Correlation at ' + str(T_slow_plot) + ' T slowness, ' + date_label1 + ' ' + date_label2)
            plt.savefig(save_name + '_Rcc_hist.png')
        fig_index += 1

        #%% -- -- -- T/E plot
        if do_trans:
            stack_array = np.zeros((slowT_n,stack_nt))

            for it in range(stack_nt):  # check points one at a time
                for slowT_i in range(slowT_n):  # for this station, loop over slownesses
                    num_val = centralT_cc[slowT_i].data[it]
                    stack_array[slowT_i, it] = num_val

            y, x = np.meshgrid(stack_Tslows,ttt)
            fig, ax = plt.subplots(1, figsize=(10,3))
            fig.subplots_adjust(bottom=0.2)
            c = ax.pcolormesh(x, y, np.transpose(stack_array), cmap=plt.cm.coolwarm, vmin=0, vmax=1)
            fig.colorbar(c, ax=ax, label='Correlation')
            ax.axis([x.min(), x.max(), y.min(), y.max()])

            c = ax.scatter(arrival_time1, pred_Eslo, color='black'  , s=50, alpha=0.75)
            plt.text(arrival_time1, pred_Eslo, phase1, color = 'black')
            if phase2 != 'no':
                c = ax.scatter(arrival_time2, pred_Eslo2, color='gray'  , s=50, alpha=0.75)
                plt.text(arrival_time2, pred_Eslo2, phase2, color = 'black')
            if phase3 != 'no':
                c = ax.scatter(arrival_time3, pred_Eslo3, color='gray'  , s=50, alpha=0.75)
                plt.text(arrival_time3, pred_Eslo3, phase3, color = 'black')
            if phase4 != 'no':
                c = ax.scatter(arrival_time4, pred_Eslo4, color='gray'  , s=50, alpha=0.75)
                plt.text(arrival_time4, pred_Eslo4, phase4, color = 'black')
            if phasePKP_single or phasePKP_double:
                c = ax.scatter(arrival_timePKP1, pred_EsloPKP1, color='purple'  , s=50, alpha=0.75)
                plt.text(arrival_timePKP1, pred_EsloPKP1, 'PKP1', color = 'black')
            if phasePKP_double:
                c = ax.scatter(arrival_timePKP2, pred_EsloPKP2, color='purple'  , s=50, alpha=0.75)
                plt.text(arrival_timePKP2, pred_EsloPKP2, 'PKP2', color = 'black')

            plt.xlabel('Time (s)')
            if NS:
                plt.ylabel('E slowness (s/km)')
                plt.title(pair_name + ' ' + arrayname + ' ' + ' Correlation at ' + str(R_slow_plot) + ' N slowness, ' + date_label1 + ' ' + date_label2)
                plt.savefig(save_name + '_EtdiffSection.png')
            else:
                plt.ylabel('T slowness (s/km)')
                plt.title(pair_name + ' ' + arrayname + ' ' + ' Correlation at ' + str(R_slow_plot) + ' R slowness, ' + date_label1 + ' ' + date_label2)
                plt.savefig(save_name + '_TtdiffSection.png')
            fig_index += 1

        #%% -- -- Amp plots
        #%% -- -- -- R/N plot
        stack_array = np.zeros((slowR_n,stack_nt))

        for it in range(stack_nt):  # check points one at a time
            for slowR_i in range(slowR_n):  # for this station, loop over slownesses
                num_val = centralR_Ast[slowR_i].data[it]
                stack_array[slowR_i, it] = num_val

        y, x = np.meshgrid(stack_Rslows,ttt)
        fig, ax = plt.subplots(1, figsize=(10,3))
        print(f'len(x) is {len(x)} and len(y) is {len(y)}')
        print(f'len(stack_Rslows) is {len(stack_Rslows)} and len(ttt) is {len(ttt)}')
        print(f'slowR_n is {slowR_n} and stack_nt is {stack_nt}')
        c = ax.pcolormesh(x, y, np.transpose(stack_array), cmap=plt.cm.gist_rainbow_r, vmin=0)
        fig.colorbar(c, ax=ax, label='linear amplitude')
        fig.subplots_adjust(bottom=0.2)
        ax.axis([x.min(), x.max(), y.min(), y.max()])

        c = ax.scatter(arrival_time1, pred_Nslo, color='black'  , s=50, alpha=0.75)
        plt.text(arrival_time1, pred_Nslo, phase1, color = 'black')
        if phase2 != 'no':
            c = ax.scatter(arrival_time2, pred_Nslo2, color='gray'  , s=50, alpha=0.75)
            plt.text(arrival_time2, pred_Nslo2, phase2, color = 'black')
        if phase3 != 'no':
            c = ax.scatter(arrival_time3, pred_Nslo3, color='gray'  , s=50, alpha=0.75)
            plt.text(arrival_time3, pred_Nslo3, phase3, color = 'black')
        if phase4 != 'no':
            c = ax.scatter(arrival_time4, pred_Nslo4, color='gray'  , s=50, alpha=0.75)
            plt.text(arrival_time4, pred_Nslo4, phase4, color = 'black')
        if phasePKP_single or phasePKP_double:
            c = ax.scatter(arrival_timePKP1, pred_NsloPKP1, color='purple'  , s=50, alpha=0.75)
            plt.text(arrival_timePKP1, pred_NsloPKP1, 'PKP1', color = 'black')
        if phasePKP_double:
            c = ax.scatter(arrival_timePKP2, pred_NsloPKP2, color='purple'  , s=50, alpha=0.75)
            plt.text(arrival_timePKP2, pred_NsloPKP2, 'PKP2', color = 'black')

        plt.xlabel('Time (s)')
        if NS:
            plt.ylabel('N slowness (s/km)')
            plt.title(pair_name + ' ' + arrayname + ' ' + ' Amp at ' + str(T_slow_plot) + ' s/km E slowness, ' + date_label1 + ' ' + date_label2)
            plt.savefig(save_name + '_Namp_hist.png')
        else:
            plt.ylabel('R slowness (s/km)')
            plt.title(pair_name + ' ' + arrayname + ' ' + ' Amp at ' + str(T_slow_plot) + ' s/km T slowness, ' + date_label1 + ' ' + date_label2)
            plt.savefig(save_name + '_Ramp_hist.png')

        fig_index += 1

        #%% -- -- -- T/E plot
        if do_trans:
            stack_array = np.zeros((slowT_n,stack_nt))

            for it in range(stack_nt):  # check points one at a time
                for slowT_i in range(slowT_n):  # for this station, loop over slownesses
                    num_val = centralT_Ast[slowT_i].data[it]
                    stack_array[slowT_i, it] = num_val

            y, x = np.meshgrid(stack_Tslows,ttt)
            fig, ax = plt.subplots(1, figsize=(10,3))
            fig.subplots_adjust(bottom=0.2)
            c = ax.pcolormesh(x, y, np.transpose(stack_array), cmap=plt.cm.gist_rainbow_r, vmin=0)
            fig.colorbar(c, ax=ax, label='linear amplitude')
            ax.axis([x.min(), x.max(), y.min(), y.max()])

            plt.text(arrival_time1, pred_Eslo, phase1, color = 'black')
            c = ax.scatter(arrival_time1, pred_Eslo, color='black'  , s=50, alpha=0.75)
            if phase2 != 'no':
                plt.text(arrival_time2, pred_Eslo2, phase2, color = 'black')
                c = ax.scatter(arrival_time2, pred_Eslo2, color='gray'  , s=50, alpha=0.75)
            if phase3 != 'no':
                plt.text(arrival_time3, pred_Eslo3, phase3, color = 'black')
                c = ax.scatter(arrival_time3, pred_Eslo3, color='gray'  , s=50, alpha=0.75)
            if phase4 != 'no':
                plt.text(arrival_time4, pred_Eslo4, phase4, color = 'black')
                c = ax.scatter(arrival_time4, pred_Eslo4, color='gray'  , s=50, alpha=0.75)
            if phasePKP_single or phasePKP_double:
                c = ax.scatter(arrival_timePKP1, pred_EsloPKP1, color='purple'  , s=50, alpha=0.75)
                plt.text(arrival_timePKP1, pred_EsloPKP1, 'PKP1', color = 'black')
            if phasePKP_double:
                c = ax.scatter(arrival_timePKP2, pred_EsloPKP2, color='purple'  , s=50, alpha=0.75)
                plt.text(arrival_timePKP2, pred_EsloPKP2, 'PKP2', color = 'black')

            plt.xlabel('Time (s)')
            if NS:
                plt.ylabel('E slowness (s/km)')
                plt.title(pair_name + ' ' + arrayname + ' ' + ' Amp at ' + str(R_slow_plot) + ' s/km N slowness, ' + date_label1 + ' ' + date_label2)
                plt.savefig(save_name + '_Eamp_hist.png')
            else:
                plt.ylabel('T slowness (s/km)')
                plt.title(pair_name + ' ' + arrayname + ' ' + ' Amp at ' + str(R_slow_plot) + ' s/km R slowness, ' + date_label1 + ' ' + date_label2)
                plt.savefig(save_name + '_Tamp_hist.png')

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
                for slowR_i in range(slowR_n):  # loop over radial slownesses
                    for slowT_i in range(slowT_n):  # loop over transverse slownesses
                        index = slowR_i*slowT_n + slowT_i
                        num_val = np.nanmean(tdiff[index].data[it_start:it_end])
                        stack_slice[slowR_i, slowT_i] = num_val

                y1, x1 = np.meshgrid(stack_Rslows,stack_Tslows)
                fig, ax = plt.subplots(1, figsize=(7,0.8*7))
                c = ax.pcolormesh(x1, y1, np.transpose(stack_slice), cmap=plt.cm.coolwarm, vmin=-tdiff_clip, vmax=tdiff_clip)
                ax.axis([x1.min(), x1.max(), y1.min(), y1.max()])
                fig.colorbar(c, ax=ax, label = 'time shift (s)')
                c = ax.scatter(        0,         0, color='black' , s=50,  alpha=0.75)
                c = ax.scatter(pred_Eslo , pred_Nslo , color='black' , s=50, alpha=0.75)
                if phase2 != 'no':
                    c = ax.scatter(pred_Eslo2, pred_Nslo2, color='black'  , s=50, alpha=0.75)
                if phase3 != 'no':
                    c = ax.scatter(pred_Eslo3, pred_Nslo3, color='black'  , s=50, alpha=0.75)
                if phase4 != 'no':
                    c = ax.scatter(pred_Eslo4, pred_Nslo4, color='black'  , s=50, alpha=0.75)
                circle1 = plt.Circle((0, 0), 0.019, color='black', fill=False)
                ax.add_artist(circle1)  # inner core limit
                circle2 = plt.Circle((0, 0), 0.040, color='black', fill=False)
                ax.add_artist(circle2)  # outer core limit
                plt.title(f'Tdiff {snap_start:.0f} to {snap_end:.0f}s  {date_label1} {date_label2}  events {eq_num1}&{eq_num2}')
                if NS:
                    plt.xlabel('E Slowness (s/km)')
                    plt.ylabel('N Slowness (s/km)')
                else:
                    plt.xlabel('T Slowness (s/km)')
                    plt.ylabel('R Slowness (s/km)')
                plt.savefig(save_name + str(snap_start) + '_' + str(snap_end) + '_Tbeam.png')

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

                y1, x1 = np.meshgrid(stack_Rslows,stack_Tslows)
                fig, ax = plt.subplots(1, figsize=(7,0.8*7))
                c = ax.pcolormesh(x1, y1, np.transpose(stack_slice), cmap=plt.cm.gist_rainbow_r)
                ax.axis([x1.min(), x1.max(), y1.min(), y1.max()])
                if log_plot == True:
                    fig.colorbar(c, ax=ax, label = 'log amp')
                else:
                    fig.colorbar(c, ax=ax, label = 'linear amp')
                c = ax.scatter(        0,         0, color='black' , s=50,  alpha=0.75)
                c = ax.scatter(pred_Eslo , pred_Nslo , color='black' , s=50, alpha=0.75)
                if phase2 != 'no':
                    c = ax.scatter(pred_Eslo2, pred_Nslo2, color='gray'  , s=50, alpha=0.75)
                if phase3 != 'no':
                    c = ax.scatter(pred_Eslo3, pred_Nslo3, color='gray'  , s=50, alpha=0.75)
                if phase4 != 'no':
                    c = ax.scatter(pred_Eslo4, pred_Nslo4, color='gray'  , s=50, alpha=0.75)
                circle1 = plt.Circle((0, 0), 0.019, color='black', fill=False)
                ax.add_artist(circle1)  # inner core limit
                circle2 = plt.Circle((0, 0), 0.040, color='black', fill=False)
                ax.add_artist(circle2)  # outer core limit
                plt.title(f'Amp {snap_start:.0f} to {snap_end:.0f}s  {date_label1} {date_label2}  events {eq_num1}&{eq_num2}')
                if NS:
                    plt.xlabel('E Slowness (s/km)')
                    plt.ylabel('N Slowness (s/km)')
                else:
                    plt.xlabel('T Slowness (s/km)')
                    plt.ylabel('R Slowness (s/km)')
                plt.savefig(save_name + str(snap_start) + '_' + str(snap_end) + '_Abeam.png')

    #%% Wiggly plots
    if wiggly_plots or pred_wiggles or max_wiggly_plot:

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

    if wiggly_plots:
        #%% -- -- R amp and tdiff vs time plots with black line for time shift
        scale_plot_wig = wig_scale_fac / (200 * global_max)
        scale_plot_tdiff = tdiff_scale_fac / 500.
        if log_plot == True:
            scale_plot_wig /= 30  # not quite sure why this renormalization works
            # scale_plot_tdiff = plot_scale_fac / 500.
        fig_index += 1
        plt.figure(fig_index,figsize=(15,6))
        plt.xlim(start_buff,end_buff)
        del_y = stack_Rslows[-1] - stack_Rslows[0]
        plt.ylim(stack_Rslows[0] - (del_y * 0.05), stack_Rslows[-1] + (del_y * 0.05))
        for slowR_i in range(slowR_n):  # loop over radial slownesses
            dist_offset = stack_Rslows[slowR_i] # trying for approx degrees
            ttt1 = (np.arange(len(centralR_st1[slowR_i].data)) * centralR_st1[slowR_i].stats.delta
              + (centralR_st1[slowR_i].stats.starttime - t1))
            plt.plot(ttt1, ((centralR_st1[slowR_i].data - np.median(centralR_st1[slowR_i].data)) * scale_plot_wig) + dist_offset, color = 'green')
            plt.plot(ttt1, ((centralR_st2[slowR_i].data - np.median(centralR_st2[slowR_i].data)) * scale_plot_wig) + dist_offset, color = 'red')
            if turn_off_black == 0:
                plt.plot(ttt1,     (centralT_st1[slowT_i].data)*0.0 + dist_offset, color = 'gray') # reference lines
                plt.plot(ttt_dec, (centralR_tdiff[slowR_i].data) * scale_plot_tdiff + dist_offset, color = 'black')

        plt.plot((Zstart_buff, Zstart_buff), (-1, 1), color = 'gray')
        plt.plot((Zend_buff,     Zend_buff), (-1, 1), color = 'lightgray')
        plt.plot((arrival_time1, arrival_time1), (-1, 1), color = 'black')
        plt.text(arrival_time1, 0, phase1, color = 'black')
        if phase2 != 'no':
            plt.plot((arrival_time2, arrival_time2), (-1, 1), color = 'gray')
            plt.text(arrival_time2, 0, phase2, color = 'black')
        if phase3 != 'no':
            plt.plot((arrival_time3, arrival_time3), (-1, 1), color = 'gray')
            plt.text(arrival_time3, 0, phase3, color = 'black')
        if phase4 != 'no':
            plt.plot((arrival_time4, arrival_time4), (-1, 1), color = 'gray')
            plt.text(arrival_time4, 0, phase4, color = 'black')

        plt.xlabel('Time (s)')

        if NS:
            plt.ylabel('N Slowness (s/km)')
            plt.title(date_label1 + '  ' + date_label2 + '  ' + ' seismograms ' + str(T_slow_plot) + ' E slowness, green is event1, red is event2')
            plt.savefig(save_name + str(start_buff) + '_' + str(end_buff) + '_N_pro_wig.png')
        else:
            plt.ylabel('R Slowness (s/km)')
            plt.title(date_label1 + '  ' + date_label2 + '  ' + ' seismograms ' + str(T_slow_plot) + ' T slowness, green is event1, red is event2')
            plt.savefig(save_name + str(start_buff) + '_' + str(end_buff) + '_R_pro_wig.png')
        #%% -- -- T amp and tdiff vs time plots with black line for time shift
        if do_trans:
            fig_index += 1
            plt.figure(fig_index,figsize=(15,6))
            plt.xlim(start_buff,end_buff)
            del_y = stack_Tslows[-1] - stack_Tslows[0]
            plt.ylim(stack_Tslows[0] - (del_y * 0.05), stack_Tslows[-1] + (del_y * 0.05))

            for slowT_i in range(slowT_n):  # loop over transverse slownesses
                dist_offset = stack_Tslows[slowT_i] # trying for approx degrees
                ttt2 = (np.arange(len(centralT_st1[slowT_i].data)) * centralT_st1[slowT_i].stats.delta
                  + (centralT_st1[slowT_i].stats.starttime - t1))
                plt.plot(ttt2, ((centralT_st1[slowT_i].data - np.median(centralT_st1[slowT_i].data)) * scale_plot_wig) + dist_offset, color = 'green')
                plt.plot(ttt2, ((centralT_st2[slowT_i].data - np.median(centralT_st2[slowT_i].data)) * scale_plot_wig) + dist_offset, color = 'red')
                if turn_off_black == 0:
                    plt.plot(ttt2,     (centralT_st1[slowT_i].data)*0.0 + dist_offset, color = 'gray') # reference lines
                    plt.plot(ttt_dec, (centralT_tdiff[slowT_i].data) * scale_plot_tdiff + dist_offset, color = 'black')

            plt.plot((Zstart_buff, Zstart_buff), (-1, 1), color = 'gray')
            plt.plot((Zend_buff,     Zend_buff), (-1, 1), color = 'lightgray')
            plt.plot((arrival_time1, arrival_time1), (-1, 1), color = 'black')
            plt.text(arrival_time1, 0, phase1, color = 'black')
            if phase2 != 'no':
                plt.plot((arrival_time2, arrival_time2), (-1, 1), color = 'gray')
                plt.text(arrival_time2, 0, phase2, color = 'black')
            if phase3 != 'no':
                plt.plot((arrival_time3, arrival_time3), (-1, 1), color = 'gray')
                plt.text(arrival_time3, 0, phase3, color = 'black')
            if phase4 != 'no':
                plt.plot((arrival_time4, arrival_time4), (-1, 1), color = 'gray')
                plt.text(arrival_time4, 0, phase4, color = 'black')
            plt.xlabel('Time (s)')

            if NS:
                plt.ylabel('E Slowness (s/km)')
                plt.title(repeater + '_' + str(ARRAY) + '_' + date_label1 + '  ' + date_label2 + '  ' + ' seismograms ' + str(R_slow_plot) + ' N slowness, green is event1, red is event2')
                plt.savefig(save_name + str(start_buff) + '_' + str(end_buff) + '_E_pro_wig.png')
            else:
                plt.ylabel('T Slowness (s/km)')
                plt.title(repeater + '_'  + str(ARRAY) + '_' ++ date_label1 + '  ' + date_label2 + '  ' + ' seismograms ' + str(R_slow_plot) + ' R slowness, green is event1, red is event2')
                plt.savefig(save_name + str(start_buff) + '_' + str(end_buff) + '_T_pro_wig.png')

    #%% Beam sum plots
    if beam_sums or max_wiggly_plot:
    #%% -- R-T tdiff amp-normed
        fig_index += 1
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

        y1, x1 = np.meshgrid(stack_Rslows,stack_Tslows)
        fig, ax = plt.subplots(1, figsize=(7,0.8*7*(slowR_n/slowT_n)))  # try to make correct aspect ratio plot
        c = ax.pcolormesh(x1, y1, np.transpose(stack_slice), cmap=plt.cm.coolwarm, vmin = -tdiff_clip, vmax = tdiff_clip)
        fig.colorbar(c, ax=ax, label='time shift (s)')
        ax.axis([x1.min(), x1.max(), y1.min(), y1.max()])
        circle1 = plt.Circle((0, 0), 0.019, color='black', fill=False)
        ax.add_artist(circle1)
        circle2 = plt.Circle((0, 0), 0.040, color='black', fill=False)
        ax.add_artist(circle2)  #outer core limit

        c = ax.scatter(pred_Eslo, pred_Nslo, color='black'  , s=50, alpha=0.75)
        plt.text(pred_Eslo, pred_Nslo, phase1, color = 'black')
        c = ax.scatter(        0,         0, color='black' , s=50,  alpha=0.75)
        plt.text(0, 0, '0', color = 'black')
        if phase2 != 'no':
            plt.text(pred_Eslo2, pred_Nslo2, phase2, color = 'black')
            c = ax.scatter(pred_Eslo2, pred_Nslo2, color='gray'  , s=50, alpha=0.75)
        if phase3 != 'no':
            plt.text(pred_Eslo3, pred_Nslo3, phase3, color = 'black')
            c = ax.scatter(pred_Eslo3, pred_Nslo3, color='gray'  , s=50, alpha=0.75)
        if phase4 != 'no':
            plt.text(pred_Eslo4, pred_Nslo4, phase4, color = 'black')
            c = ax.scatter(pred_Eslo4, pred_Nslo4, color='gray'  , s=50, alpha=0.75)
        if phasePKP_single or phasePKP_double:
            plt.text(pred_EsloPKP1, pred_NsloPKP1, 'PKP1', color = 'black')
            c = ax.scatter(pred_EsloPKP1, pred_NsloPKP1, color='gray'  , s=50, alpha=0.75)
        if phasePKP_double:
            plt.text(pred_EsloPKP2, pred_NsloPKP2, 'PKP2', color = 'black')
            c = ax.scatter(pred_EsloPKP2, pred_NsloPKP2, color='gray'  , s=50, alpha=0.75)

        if NS:
            plt.ylabel('N Slowness (s/km)')
            plt.xlabel('E Slowness (s/km)')
        else:
            plt.ylabel('R Slowness (s/km)')
            plt.xlabel('T Slowness (s/km)')
        plt.title(f'{pair_name} {arrayname} {date_label1} {date_label2} {start_buff:.0f} to {end_buff:.0f} time shift')
        plt.savefig(save_name + '_Tbeam.png')

    #%% -- R-T amplitude averaged over time window
        fig_index += 1
        stack_slice = np.zeros((slowR_n,slowT_n))
        for slowR_i in range(slowR_n):  # loop over radial slownesses
            for slowT_i in range(slowT_n):  # loop over transverse slownesses
                index = slowR_i*slowT_n + slowT_i
                if full_beam == 1:
                    num_val = np.nanmean(amp_ave[index].data)
                else:
                    num_val = np.nanmean(amp_ave[index].data[start_index:end_index])
                stack_slice[slowR_i, slowT_i] = num_val
                # print('slowR_n is ' + str(slowR_n) + ' slowT_n is ' + str(slowT_n) + ' index is ' + str(index) + ' num_val is ' + str(num_val))

        # print('stack_slice[8,8] is ' + str(stack_slice[8,8]))
        y1, x1 = np.meshgrid(stack_Rslows,stack_Tslows)
        fig, ax = plt.subplots(1, figsize=(7,0.8*7*(slowR_n/slowT_n)))
        smax = np.max(stack_slice)
        smin = np.min(stack_slice)

        max_xy = np.where(stack_slice == stack_slice.max() ) # find indices of slowness with max amplitude, needed for max_wiggle_plot
        print('len of max_xy[0] is ' + str(len(max_xy[0])) + ' max is ' + str(stack_slice.max()) + ' min is ' + str(stack_slice.min()))
        if len(max_xy[0]) == 2:  # rare case of two identical maxima, just use first one
               max_xy = max_xy[0]
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

        c = ax.scatter(pred_Eslo, pred_Nslo, color='black'  , s=50, alpha=0.75)
        plt.text(pred_Eslo, pred_Nslo, phase1, color = 'black')
        c = ax.scatter(        0,         0, color='black' , s=50,  alpha=0.75)
        text_offset = (x1.max() - x1.min())/50
        plt.text(text_offset, 0, '0', color = 'black')
        if phase2 != 'no':
            plt.text(pred_Eslo2 + text_offset, pred_Nslo2, phase2, color = 'black')
            c = ax.scatter(pred_Eslo2, pred_Nslo2, color='gray'  , s=50, alpha=0.75)
        if phase3 != 'no':
            plt.text(pred_Eslo3 + text_offset, pred_Nslo3, phase3, color = 'black')
            c = ax.scatter(pred_Eslo3, pred_Nslo3, color='gray'  , s=50, alpha=0.75)
        if phase4 != 'no':
            plt.text(pred_Eslo4 + text_offset, pred_Nslo4, phase4, color = 'black')
            c = ax.scatter(pred_Eslo4, pred_Nslo4, color='gray'  , s=50, alpha=0.75)
        if phasePKP_single or phasePKP_double:
            plt.text(pred_EsloPKP1 + text_offset, pred_NsloPKP1, 'PKP1', color = 'black')
            c = ax.scatter(pred_EsloPKP1, pred_NsloPKP1, color='gray'  , s=50, alpha=0.75)
        if phasePKP_double:
            plt.text(pred_EsloPKP2 + text_offset, pred_NsloPKP2, 'PKP2', color = 'black')
            c = ax.scatter(pred_EsloPKP2, pred_NsloPKP2, color='gray'  , s=50, alpha=0.75)

        if NS:
            plt.xlabel('E Slowness (s/km)')
            plt.ylabel('N Slowness (s/km)')
        else:
            plt.xlabel('T Slowness (s/km)')
            plt.ylabel('R Slowness (s/km)')
        plt.title(f'{pair_name} {arrayname} {date_label1} {date_label2} {start_buff:.0f} to {end_buff:.0f} beam amp')
        plt.savefig(save_name + '_Abeam.png')

    #%% Wiggly plot at max amplitude
    if max_wiggly_plot and zoom:
        #%% -- Select stack at max amp
        max_trace1 = Stream()
        max_trace2 = Stream()
        iii = int(max_xy[0])*slowT_n + int(max_xy[1])
        max_trace1 = st1[iii]
        max_trace2 = st2[iii]

        max_trace1.data /= max(abs(max_trace1.data))
        max_trace2.data /= max(abs(max_trace2.data))
        if log_plot:
            scale_plot_wig /= 30  # not quite sure why this renormalization works
        # fig_tit = str(eq_num) + '_maxamp'
        fig_index += 1
        plt.figure(figsize=(10,5), num = fig_index)
        plt.xlim(orig_start_buff,orig_end_buff)
        plt.ylim(-1.1, 1.1)
        ttt = (np.arange(len(max_trace1)) * centralR_st1[0].stats.delta
          + (centralR_st1[0].stats.starttime - t1))

        # normalize in zoom window
        questor = (ttt >= Zstart_buff) & (ttt < Zend_buff) # identify zoom window
        ts_sel1 = max_trace1[questor]  #extract zoom wondow
        ts_sel2 = max_trace2[questor]  #extract zoom wondow
        max_env1 = max(abs(ts_sel1))
        max_env2 = max(abs(ts_sel2))
        max_trace1.data = max_trace1.data/(max_env1)
        max_trace2.data = max_trace2.data/(max_env2)

        print('Length of ttt and max_trace1 and 2:  ' + str(len(ttt)) + ' ' +  str(len(max_trace1)) + ' ' +  str(len(max_trace2)))
        plt.plot(ttt, max_trace1, color = 'green')
        plt.plot(ttt, max_trace2, color = 'red')
        plt.plot((Zstart_buff, Zstart_buff), (-1, 1), color = 'gray')
        plt.text(Zstart_buff, -1, 'start', color = 'black')
        plt.plot((Zend_buff,     Zend_buff), (-1, 1), color = 'lightgray')
        plt.text(Zend_buff  , -1,   'end', color = 'black')
        plt.plot((arrival_time1,arrival_time1), (-1, 1), color = 'black')
        plt.text(arrival_time1, 0.8, phase1, color = 'black')
        if phase2 != 'no':
            plt.plot((arrival_time2, arrival_time2), (-1, 1), color = 'gray')
            plt.text(arrival_time2, 1, phase2, color = 'black')
        if phase3 != 'no':
            plt.text(arrival_time3, 1, phase3, color = 'black')
            plt.plot((arrival_time3, arrival_time3), (-1, 1), color = 'gray')
        if phase4 != 'no':
            plt.text(arrival_time4, 1, phase4, color = 'black')
            plt.plot((arrival_time4, arrival_time4), (-1, 1), color = 'gray')
        if phasePKP_single or phasePKP_double:
            plt.text(arrival_timePKP1, 0.9, 'PKP1', color = 'black')
            plt.plot((arrival_timePKP1, arrival_timePKP1), (-1, 1), color = 'gray')
        if phasePKP_double:
            plt.text(arrival_timePKP2, 0.9, 'PKP2', color = 'black')
            plt.plot((arrival_timePKP2, arrival_timePKP2), (-1, 1), color = 'gray')

        plt.xlabel('Time (s)')
        plt.ylabel('Amplitude')

        plt.title(f'{arrayname}  {phase1} {date_label1} {date_label2}  max amp stack \n in events {eq_num1} and {eq_num2}  Rslow  {Rslow_max:.4f}  Tslow {Tslow_max:.4f}')
        # os.chdir('/Users/vidale/Documents/Research/IC/Plots_hold')
        # if zoom:
        #     plt.savefig('wig_' + date_label + '_' + str(eq_num) + '_' + str(start_buff) + '_' + str(end_buff) + '_' + str(fig_index) + '.png')

    #%% Wiggly plot at predicted slowness
    if pred_wiggles and zoom:

        misfit = 1000000  # find the grid point closest to the predicted slowness of phase1
        for slowR_i in range(slowR_n):  # loop over radial slownesses
            for slowT_i in range(slowT_n):  # loop over transverse slownesses
                index = slowR_i*slowT_n + slowT_i
                slowR_actual = stack_Rslows[slowR_i]
                slowT_actual = stack_Tslows[slowT_i]
                slow_anomaly = np.sqrt(((slowR_actual - pred_Nslo) * (slowR_actual - pred_Nslo)) +
                                           ((slowT_actual - pred_Eslo) * (slowT_actual - pred_Eslo)))
                if slow_anomaly < misfit:
                    misfit = slow_anomaly
                    closest_slowR_i = slowR_i
                    closest_slowT_i = slowT_i

        iii = closest_slowR_i * slowT_n + closest_slowT_i

        Rslow_pred = closest_slowR_i * slow_delta + slowR_lo
        Tslow_pred = closest_slowT_i * slow_delta + slowT_lo

        pred_trace1 = Stream()
        pred_trace2 = Stream()

        pred_trace1 = st1[iii]
        pred_trace2 = st2[iii]

        pred_trace1.data /= max(abs(pred_trace1.data))
        pred_trace2.data /= max(abs(pred_trace2.data))
        if log_plot:
            scale_plot_wig /= 30  # not quite sure why this renormalization works
        fig_index += 1
        plt.figure(figsize=(10,5), num = fig_index)
        # plt.xlim(orig_start_buff,orig_end_buff)
        plt.xlim(start_buff,end_buff)
        edge_height = plot_peak # clipping amplitude of plot
        plt.ylim(-edge_height * 1.1, edge_height * 1.1)
        ttt = (np.arange(len(pred_trace1)) * centralR_st1[0].stats.delta
          + (centralR_st1[0].stats.starttime - t1))

        print('Length of ttt and pred_trace1 and 2:  ' + str(len(ttt)) + ' ' +  str(len(pred_trace1)) + ' ' +  str(len(pred_trace2)))
        plt.plot(ttt, pred_trace1, color = 'green')
        plt.plot(ttt, pred_trace2, color = 'red')
        # plt.plot((Zstart_buff, Zstart_buff), (-edge_height, edge_height), color = 'gray')
        # plt.text(Zstart_buff, -edge_height, 'start', color = 'black')
        # plt.plot((Zend_buff,     Zend_buff), (-edge_height, edge_height), color = 'lightgray')
        # plt.text(Zend_buff  , -edge_height,   'end', color = 'black')
        plt.plot((arrival_time1,arrival_time1), (-edge_height, edge_height), color = 'black')
        plt.text(arrival_time1, 0.8*edge_height, phase1, color = 'black')
        if phase2 != 'no':
            plt.plot((arrival_time2, arrival_time2), (-edge_height, edge_height), color = 'gray')
            plt.text(arrival_time2, edge_height, phase2, color = 'black')
        if phase3 != 'no':
            plt.text(arrival_time3, edge_height, phase3, color = 'black')
            plt.plot((arrival_time3, arrival_time3), (-edge_height, edge_height), color = 'gray')
        if phase4 != 'no':
            plt.text(arrival_time4, edge_height, phase4, color = 'black')
            plt.plot((arrival_time4, arrival_time4), (-edge_height, edge_height), color = 'gray')
        if phasePKP_single or phasePKP_double:
            plt.text(arrival_timePKP1, 0.9*edge_height, 'PKP1', color = 'black')
            plt.plot((arrival_timePKP1, arrival_timePKP1), (-edge_height, edge_height), color = 'gray')
        if phasePKP_double:
            plt.text(arrival_timePKP2, 0.9*edge_height, 'PKP2', color = 'black')
            plt.plot((arrival_timePKP2, arrival_timePKP2), (-edge_height, edge_height), color = 'gray')

        plt.xlabel('Time (s)')
        plt.ylabel('Amplitude')

        plt.title(f'{pair_name}  {arrayname}  {phase1} {date_label1} {date_label2}  stack at predicted slowness \n in events {eq_num1} and {eq_num2}  Dist {ref_dist:.2f} Rslow  {Rslow_pred:.4f}  Tslow {Tslow_pred:.4f}')
        os.chdir('/Users/vidale/Documents/Research/IC/Plots_hold')
        plt.savefig(save_name + '_' + 'pred_wig' + str(fig_index) + '.png')

#  Save processed files
#    fname = 'HD' + date_label + '_slice.mseed'
#    stack.write(fname,format = 'MSEED')
    print(f'1 {phase1} has slowness {pred_Eslo:.3f} {pred_Nslo:.3f} pred arrival time {arrival_time:.3f} at dist {ref_dist:.3f}')
    if phase2 != 'no':
        print(f'2 {phase2} has slowness {pred_Eslo2:.3f} {pred_Nslo2:.3f} pred arrival time {arrival_time2:.3f}')
    if phase3 != 'no':
        print(f'3 {phase3} has slowness {pred_Eslo3:.3f} {pred_Nslo3:.3f} pred arrival time {arrival_time3:.3f}')
    if phase4 != 'no':
        print(f'4 {phase4} has slowness {pred_Eslo4:.3f} {pred_Nslo4:.3f} pred arrival time {arrival_time4:.3f}')
    if phasePKP_single or phasePKP_double:
        print(f'PKP1 has slowness {pred_EsloPKP1:.3f} {pred_NsloPKP1:.3f} pred arrival time {arrival_timePKP1:.3f}')
    if phasePKP_double:
        print(f'PKP2 has slowness {pred_EsloPKP2:.3f} {pred_NsloPKP2:.3f} pred arrival time {arrival_timePKP2:.3f}')

    elapsed_time_wc = time.time() - start_time_wc
    print(f'Pro7 took {elapsed_time_wc:.1f} seconds')
    os.system('say "seven done"')
