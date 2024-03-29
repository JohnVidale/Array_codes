#!/usr/bin/env python
# 1D Slant stack for a single event
# Input is set of selected traces "*sel.mseed"
# traces have already been aligned and corrected for near-vertical statics
#   to have specified phase start at the earthquake origin time
# plots either traces or envelopes of traces
# saves 1D stack "_1Dstack.mseed"
# John Vidale 2/2019

def pro5stack1d(eq_num, plot_scale_fac = 0.05, slowR_lo = -0.1, slowR_hi = 0.1, stack_option = 1,
            slow_delta = 0.0005, start_buff = -50, end_buff = 50,
            ref_lat = 36.3, ref_lon = 138.5, ref_loc = 0, envelope = 1, plot_dyn_range = 1000,
            log_plot = 1, norm = 1, global_norm_plot = 1, color_plot = 1, fig_index = 401, ARRAY = 0):

#%% Import functions
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
    from matplotlib.colors import LogNorm
    model = TauPyModel(model='iasp91')
    from scipy.signal import hilbert
    import math
    import time
    from termcolor import colored

    env_stack = 0  # flag to stack envelopes instead of oscillating seismograms

#    import sys # don't show any warnings
#    import warnings

    print(colored('Running pro5a_stack', 'cyan'))

#%% Get saved event info, also used to name files
    start_time_wc = time.time()

    fname = '/Users/vidale/Documents/Research/IC/EvLocs/event' + str(eq_num) + '.txt'
    file = open(fname, 'r')
    lines=file.readlines()

    split_line = lines[0].split()
#            ids.append(split_line[0])  ignore label for now
    t           = UTCDateTime(split_line[1])
    date_label  = split_line[1][0:10]
    ev_lat      = float(      split_line[2])
    ev_lon      = float(      split_line[3])
    ev_depth    = float(      split_line[4])

    #if not sys.warnoptions:
    #    warnings.simplefilter("ignore")

#%% Get station location file
    if ARRAY == 0: # Hinet set and center
        sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/sta_hinet.txt'
        if ref_loc == 0:
            ref_lat = 36.3
            ref_lon = 138.5
    elif ARRAY == 1: # LASA set and center
        sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/sta_LASA.txt'
        if ref_loc == 0:
            ref_lat = 46.69
            ref_lon = -106.22
    elif ARRAY == 2: # China set and center
        sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/sta_ch.txt'
        if ref_loc == 0:
            ref_lat = 38      # °N
            ref_lon = 104.5   # °E
    else:         # NORSAR set and center
        sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/sta_NORSAR.txt'
        if ref_loc == 0:
            ref_lat = 61
            ref_lon = 11
    with open(sta_file, 'r') as file:
        lines = file.readlines()
    print('    ' + str(len(lines)) + ' stations of metadata read from ' + sta_file)
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

#%% Name file, read data
    # date_label = '2018-04-02' # date for filename
    fname = 'HD' + date_label + 'sel.mseed'
    goto = '/Users/vidale/Documents/Research/IC/Pro_Files'
    os.chdir(goto)

    # fname = '/Users/vidale/Documents/PyCode/Pro_Files/HD' + date_label + 'sel.mseed'

    st = Stream()
    print('        reading ' + fname)
    print('        Stack option is ' + str(stack_option))
    st = read(fname)
    print('    ' + str(len(st)) + ' traces read in')
    nt = len(st[0].data)
    dt = st[0].stats.delta
    print(f'        First trace has {nt} time pts, time sampling of {dt:.2f} and thus duration of {(nt-1)*dt:.0f} and max amp of {max(abs(st[0].data)):.1f}')
    print(f'st[0].stats.starttime-t {(st[0].stats.starttime-t):.2f} start_buff {start_buff:.2f}')

#%% Build Stack arrays
    stack = Stream()
    tr = Trace()
    tr.stats.delta = dt
    tr.stats.network = 'stack'
    tr.stats.channel = 'BHZ'
    slow_n = int(1 + (slowR_hi - slowR_lo)/slow_delta)  # number of slownesses
    stack_nt = int(1 + ((end_buff - start_buff)/dt))  # number of time points
    # In English, stack_slows = range(slow_n) * slow_delta - slowR_lo
    a1 = range(slow_n)
    stack_slows = [(x * slow_delta + slowR_lo) for x in a1]
    print('        ' + str(slow_n) + ' slownesses.')
    tr.stats.starttime = t + start_buff
    # print(f'tr.stats.starttime-t {(tr.stats.starttime-t):.2f} start_buff {start_buff:.2f}')
    tr.data = np.zeros(stack_nt)
    done = 0
    for stack_one in stack_slows:
        tr1 = tr.copy()
        tr1.stats.station = str(int(done))
        stack.extend([tr1])
        done += 1
    #    stack.append([tr])
    #    stack += tr

    #  Only need to compute ref location to event distance once
    ref_distance = gps2dist_azimuth(ev_lat,ev_lon,ref_lat,ref_lon)

#%% Select traces by distance, window and adjust start time to align picked times
    done = 0
    if env_stack == 1: #convert oscillating seismograms to envelopes
        for tr in st:
            tr.data = np.abs(hilbert(tr.data))

    for tr in st: # traces one by one
        if tr.stats.station in st_names:  # find station in station list
            ii = st_names.index(tr.stats.station)
            if norm == 1:
                tr.normalize()
            stalat = float(st_lats[ii])
            stalon = float(st_lons[ii]) # look up lat & lon again to find distance
            distance = gps2dist_azimuth(stalat,stalon,ev_lat,ev_lon) # Get traveltimes again, hard to store
            tr.stats.distance=distance[0] # distance in m
            del_dist = (ref_distance[0] - distance[0])/(1000) # in km
            rel_start_buff = tr.stats.starttime - (t + start_buff)
            print(f'{tr.stats.station} del_dist {del_dist:.2f} ref_dist {ref_distance[0]/1000.:.2f} distance {distance[0]/1000.:.2f} rel_start_buff {rel_start_buff:.2f} tr.stats.starttime-t {(tr.stats.starttime-t):.2f} start_buff {start_buff:.2f}')

            for slow_i in range(slow_n):  # for this station, loop over slownesses
                time_lag = -del_dist * stack_slows[slow_i]  # time shift due to slowness, flipped to match 2D
                time_correction = (rel_start_buff + time_lag)/dt
                # print(f'{slow_i} time_lag {time_lag:.1f} time correction {time_correction:.1f}')

                if stack_option == 0:
                    for it in range(stack_nt):  # check points one at a time
                        it_in = int(it + time_correction)
                        if it_in >= 0 and it_in < nt - 1: # does data lie within seismogram?
                            stack[slow_i].data[it] += tr[it_in]

                if stack_option == 1:
                    arr = tr.data
                    nshift = int(time_correction)
                    if time_correction < 0:
                        nshift = nshift-1
                    if nshift <= 0:
                        nbeg1 = -nshift
                        nend1 = stack_nt
                        nbeg2 = 0
                        nend2 = stack_nt + nshift;
                    elif nshift > 0:
                        nbeg1 = 0
                        nend1 = stack_nt - nshift
                        nbeg2 = nshift
                        nend2 = stack_nt
                    if nend1 >= 0 and nbeg1 <= stack_nt:
                        stack[slow_i].data[nbeg1:nend1] += arr[nbeg2:nend2]

            done += 1
            if done % 50 == 0:
                print('        Done stacking ' + str(done) + ' out of ' + str(len(st)) + ' stations.')
        else:
            print(tr.stats.station + ' not found in station list')

#%% Plot traces
    global_max = 0
    for slow_i in range(slow_n): # find global max, and if requested, take envelope
        if len(stack[slow_i].data) == 0:
                print('%d data has zero length ' % (slow_i))
        if envelope == 1 or color_plot == 1:
            stack[slow_i].data = np.abs(hilbert(stack[slow_i].data))
        local_max = max(abs(stack[slow_i].data))
        if local_max > global_max:
            global_max = local_max
    if global_max <= 0:
        print('        global_max ' + str(global_max) + ' slow_n ' + str(slow_n))

    # create time axis (x-axis), use of slow_i here is arbitrary, oops
    ttt = (np.arange(len(stack[slow_i].data)) * stack[slow_i].stats.delta +
         (stack[slow_i].stats.starttime - t)) # in units of seconds

    # Plotting
    if color_plot == 1: # 2D color plot
        stack_array = np.zeros((slow_n,stack_nt))

    #    stack_array = np.random.rand(int(slow_n),int(stack_nt))  # test with random numbers
        min_allowed = global_max/plot_dyn_range
        if log_plot == 1:
            for it in range(stack_nt):  # check points one at a time
                for slow_i in range(slow_n):  # for this station, loop over slownesses
                    num_val = stack[slow_i].data[it]
                    if num_val < min_allowed:
                        num_val = min_allowed
                    stack_array[slow_i, it] = math.log10(num_val) - math.log10(min_allowed)
        else:
            for it in range(stack_nt):  # check points one at a time
                for slow_i in range(slow_n):  # for this station, loop over slownesses
                    stack_array[slow_i, it] = stack[slow_i].data[it]/global_max
        y, x = np.mgrid[slice(stack_slows[0], stack_slows[-1] + slow_delta, slow_delta),
                     slice(ttt[0], ttt[-1] + dt, dt)]  # make underlying x-y grid for plot
    #    y, x = np.mgrid[ stack_slows , time ]  # make underlying x-y grid for plot
        plt.close(fig_index)

        fig, ax = plt.subplots(1, figsize=(9,9))
        fig.subplots_adjust(bottom=0.3)
        c = ax.pcolormesh(x, y, stack_array, cmap=plt.cm.gist_rainbow_r)
        # c = ax.pcolormesh(x, y, stack_array, cmap=plt.cm.gist_yarg)
        # c = ax.pcolormesh(x, y, stack_array, cmap=plt.cm.binary)
        ax.axis([x.min(), x.max(), y.min(), y.max()])
        if log_plot == 1:
            fig.colorbar(c, ax=ax, label='log amplitude')
        else:
            fig.colorbar(c, ax=ax, label='linear amplitude')
        plt.figure(fig_index,figsize=(6,8))
        plt.close(fig_index)
    else: # line plot
        for slow_i in range(slow_n):
            dist_offset = stack_slows[slow_i] # in units of slowness
            if global_norm_plot != 1:
                plt.plot(ttt, stack[slow_i].data*plot_scale_fac / (stack[slow_i].data.max()
            - stack[slow_i].data.min()) + dist_offset, color = 'black')
            else:
                plt.plot(ttt, stack[slow_i].data*plot_scale_fac / (global_max
            - stack[slow_i].data.min()) + dist_offset, color = 'black')
        plt.ylim(slowR_lo,slowR_hi)
        plt.xlim(start_buff,end_buff)
    plt.xlabel('Time (s)')
    plt.ylabel('Slowness (s/km)')
    plt.title('1Dstack   ' + str(eq_num) + '  ' + date_label)
    # os.chdir('/Users/vidale/Documents/PyCode/Plots')
    # plt.savefig(date_label + '_' + str(start_buff) + '_' + str(end_buff) + '_1D.png')
    plt.show()

#%% Save processed files
    print('        Stack has ' + str(len(stack)) + ' slownesses')
#
#    if ARRAY == 0:
#        goto = '/Users/vidale/Documents/PyCode/Hinet'
#    if ARRAY == 1:
#        goto = '/Users/vidale/Documents/PyCode/LASA/Pro_Files'
#    os.chdir(goto)
#    fname = 'HD' + date_label + '_1dstack.mseed'
#    stack.write(fname,format = 'MSEED')

    elapsed_time_wc = time.time() - start_time_wc
    print(f'    This job took   {elapsed_time_wc:.1f}   seconds')
    os.system('say "Done"')