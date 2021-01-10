#!/usr/bin/env python
# looks at envelope stack differential parameters from pro6
# Reads in tdiff, ave_amp, cc computed from a pair of events
# window by signal quality
# eleven_slice == T produces sections at 4 radial slownesses (0.0, 0.05, 0.01, 0.015)
# plus plots sections at 7 transv slownesses (-0.015 to 0.015)
# eleven_slice == F plots one radial and one transverse plot, plus snaps
# John Vidale 2/2019, overhauled 1/2021

def pro7_pair_scan_junk(eq_file1, eq_file2, slow_delta = 0.0005,
              slowR_lo = -0.1, slowR_hi = 0.1, slowT_lo = -0.1, slowT_hi = 0.1,
              start_buff = 50, end_buff = 50,fig_index = 401, skip_T = 0, skip_R = 0,
              ZslowR_lo = -0.1, ZslowR_hi = 0.1, ZslowT_lo = -0.1, ZslowT_hi = 0.1,
              Zstart_buff = 50, Zend_buff = 50, zoom = 0, tdiff_clip = 1,
              ref_phase = 'blank', cc_thres = 0.8, min_amp = 0.2, eleven_slice = True,
              R_slow_plot = 0.06, T_slow_plot = 0.0, snaptime = 8, snaps = 10, skip_snaps = 0):

    from obspy import read
    import numpy as np
    import os
    import matplotlib.pyplot as plt
    import time
    from obspy import UTCDateTime
    from obspy import Stream
    from termcolor import colored

    print(colored('Running pro7b_plot_stack', 'cyan'))
    start_time_wc = time.time()

    if zoom == 1:
        if Zstart_buff  < start_buff:
            print(f'Zstart_buff of {Zstart_buff:.1f} cannot be < start_buff of {start_buff:.1f}')
            Zstart_buff = start_buff
        if Zend_buff    > end_buff:
            print(f'Zend_buff of {Zend_buff:.1f} cannot be > end_buff of {end_buff:.1f}')
            Zend_buff   = end_buff

    #%% Input parameters and computed files
    folder_name = '/Users/vidale/Documents/Research/IC/'
    file1 = open(folder_name + 'EvLocs/' + eq_file1, 'r')
    file2 = open(folder_name + 'EvLocs/' + eq_file2, 'r')
    lines1=file1.readlines()
    lines2=file2.readlines()
    split_line1 = lines1[0].split()
    split_line2 = lines2[0].split()
    t1   = UTCDateTime(split_line1[1])
    # t2 = UTCDateTime(split_line2[1])
    date_label1  = split_line1[1][0:10]
    date_label2  = split_line2[1][0:10]

    # date_label = '2018-04-02' # dates in filename
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

#%% -- compute timing time series
    ttt = (np.arange(len(tdiff[0].data)) * tdiff[0].stats.delta + start_buff) # in units of seconds
    #%% -- take 7 R slices at fixed T
    lowest_Tslow = 1000000
    for slow_i in range(slowT_n):
        if abs(stack_Tslows[slow_i] - 0.005) < lowest_Tslow:
            lowest_TindexM05 = slow_i
            lowest_Tslow = abs(stack_Tslows[slow_i] - 0.005)

    print(f'{lowest_TindexM05:4d} is T slow nearest -0.005, difference is {lowest_Tslow:.3f}')

    # Select only stacks with that slowness for Transverse plot
    centralR05_st = Stream()
    for slowR_i in range(slowR_n):
        centralR05_st += tdiff[slowR_i*slowT_n + lowest_TindexM05]


    #%% -- take 4 T slices at fixed R
    # Slice near radial slowness of 0.005
    lowest_Rslow = 1000000
    for slow_i in range(slowR_n):
        if abs(stack_Rslows[slow_i] - 0.005) < lowest_Rslow:
            lowest_Rindex05 = slow_i
            lowest_Rslow = abs(stack_Rslows[slow_i] - 0.005)

    print(f'{lowest_Rindex05:4d} is R slow nearest 0.005, difference is {lowest_Rslow:.3f}')

    # Select only stacks with that slowness for Radial plot
    centralT05_st = Stream()
    for slowT_i in range(slowT_n):
        centralT05_st += tdiff[lowest_Rindex05*slowT_n + slowT_i]

#%% -- Plot R-time stacks
    stack_array = np.zeros((slowR_n,stack_nt))

    fig_index += 1
    stack_array = np.zeros((slowR_n,stack_nt))

    for it in range(stack_nt):  # check points one at a time
        for slowR_i in range(slowR_n):  # for this station, loop over slownesses
            num_val = centralR05_st[slowR_i].data[it]
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
    plt.title(ref_phase + ' Time lag at 0.005 s/km transverse slowness, ' + fname1[48:58] + '   ' + fname1[59:69] + ' min amp ' + str(min_amp) + ' cc_thres ' + str(cc_thres))
    plt.show()

#%%  -- Plot T-time stacks
    stack_array = np.zeros((slowT_n,stack_nt))

    for it in range(stack_nt):  # check points one at a time
        for slowT_i in range(slowT_n):  # for this station, loop over slownesses
            num_val = centralT05_st[slowT_i].data[it]
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
    plt.title(ref_phase + ' Time lag at 0.005 s/km radial slowness, ' + fname1[48:58] + '   ' + fname1[59:69] + ' min amp ' + str(min_amp) + ' cc_thres ' + str(cc_thres))
    plt.show()