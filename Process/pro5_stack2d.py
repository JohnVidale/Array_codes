#!/usr/bin/env python
# 2D Slant stack for a single event
# Input is set of selected traces "*sel.mseed"
# traces have already been aligned and corrected for near-vertical statics
#   to have specified phase start at the earthquake origin time
# doesn't plot
# saves 2D stack "_2Dstack.mseed" and envelope of 2D stack "_2Dstack_env.mseed"
# John Vidale 2/2019

def pro5stack2d(eq_num, slow_delta = 0.0005, slowR_lo = -0.1, slowR_hi = 0.1, slowT_lo = -0.1, slowT_hi = 0.1,
              start_buff = -50, end_buff = 50, norm = True, ARRAY = 0, NS = False, decimate_fac = 0,
              ref_loc = False, ref_lat = 36.3, ref_lon = 138.5, stack_option = 1, min_dist = 0, max_dist = 180):

    from obspy import UTCDateTime
    from obspy import Stream, Trace
    from obspy import read
    from obspy.geodetics import gps2dist_azimuth
    import numpy as np
    import os
    from scipy.signal import hilbert
    import math
    import time

    import sys # don't show any warnings
    import warnings
    from termcolor import colored
    print(colored('Running pro5b_stack2d', 'cyan'))

    env_stack = 0  # flag to stack envelopes instead of oscillating seismograms
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
#    ev_depth    = float(      split_line[4])

    if not sys.warnoptions:
        warnings.simplefilter("ignore")

#%% Get location file
    if ARRAY == 0: # Hinet set
        sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/sta_hinet.txt'
    elif ARRAY == 1:         # LASA set
        sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/sta_LASA.txt'
    elif ARRAY == 2: # China set and center
        sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/sta_ch.txt'
    elif ARRAY == 3: #         NORSAR set
        sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/sta_NORSAR.txt'
    elif ARRAY == 4: #         Warramunga set
        sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/sta_AU_WR.txt'
    elif ARRAY == 5: #         Yellowknife set
        sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/sta_CN_YK.txt'

#%% Set array reference location if not input
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
        elif ARRAY == 4:
            ref_lat = -19.90  # °N Warramunga
            ref_lon = 134.42  # °E
        elif ARRAY == 5:
            ref_lat =  62.49  # °N Yellowknife
            ref_lon = -114.6  # °E
    with open(sta_file, 'r') as file:
        lines = file.readlines()
    print(str(len(lines)) + ' stations read from ' + sta_file)
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

#%% Input parameters
    # date_label = '2018-04-02' # date for filename
    fname = 'HD' + date_label + 'sel.mseed'
    goto = '/Users/vidale/Documents/Research/IC/Pro_Files'
    os.chdir(goto)

    st = Stream()
    st = read(fname)
    print('Read in: ' + str(len(st)) + ' traces')
    nt = len(st[0].data)
    dt = st[0].stats.delta
    print('First trace has : ' + str(nt) + ' time pts, time sampling of '
          + str(dt) + ' and thus duration of ' + str((nt-1)*dt))

    #%% Make grid of slownesses
    slowR_n = int(round(1 + (slowR_hi - slowR_lo)/slow_delta))  # number of slownesses
    slowT_n = int(round(1 + (slowT_hi - slowT_lo)/slow_delta))  # number of slownesses
    stack_nt = int(round(1 + ((end_buff - start_buff)/dt)))  # number of time points

    # In English, stack_slows = range(slow_n) * slow_delta - slow_lo
    a1R = range(slowR_n)
    a1T = range(slowT_n)
    stack_Rslows = [(x * slow_delta + slowR_lo) for x in a1R]
    stack_Tslows = [(x * slow_delta + slowT_lo) for x in a1T]

    # testing slownesses in indexing
    print(str(slowR_n) + ' radial slownesses, ' + str(slowT_n) + ' trans slownesses, ')
    print('Radial     slownesses 0' + ' ' + str(stack_Rslows[0]) + '   ' 'end' + ' ' + str(stack_Rslows[-1]))
    print('Transverse slownesses 1' + ' ' + str(stack_Tslows[0]) + '   ' 'end' + ' ' + str(stack_Tslows[-1]))


#%% Build empty Stack array
    stack = Stream()
    tr = Trace()
    tr.stats.delta = dt
    tr.stats.starttime = t + start_buff
    tr.stats.npts = stack_nt
    tr.stats.network = 'stack'
    tr.stats.channel = 'BHZ'
    tr.data = np.zeros(stack_nt)
    done = 0
    for stackR_one in stack_Rslows:
        for stackT_one in stack_Tslows:
            tr1 = tr.copy()
            tr1.stats.station = str(int(round(done)))
            stack.extend([tr1])
            done += 1

    #  Only need to compute ref location to event distance once
    ref_dist_az = gps2dist_azimuth(ev_lat,ev_lon,ref_lat,ref_lon)
    ref_back_az = ref_dist_az[2]

#%% select by window, norm, and adjust start time to align picked times
    if env_stack == 1:
        for tr in st: #  #convert oscillating seismograms to envelopes
            tr.data = np.abs(hilbert(tr.data))

    done = 0  # kludge to get centroid of array
    ave_lat = 0
    ave_lon = 0
    cnt = 0
    bad_trace = False
    for tr in st: #  #convert oscillating seismograms to envelopes
        cnt = cnt + 1
        ave_lat = ave_lat + float(st_lats[ii])
        ave_lon = ave_lon + float(st_lons[ii])
    ave_lat = ave_lat / cnt
    ave_lon = ave_lon / cnt
    print(f'Average latitude  {ave_lat:.3f}   Average longitude  {ave_lon:.3f}')

    for tr in st: # traces one by one, find lat-lon by searching entire inventory.  Inefficient but cheap
        if tr.stats.station in st_names:  # find station in station list
            ii = st_names.index(tr.stats.station)
            if norm == True:
                tr.normalize() # trace divided abs(max of trace)
            stalat = float(st_lats[ii])
            stalon = float(st_lons[ii]) # use lat & lon to find distance and back-az
            # rel_dist_az = gps2dist_azimuth(stalat,stalon,ref_lat,ref_lon)
            rel_dist_az = gps2dist_azimuth(stalat,stalon,ave_lat,ave_lon)
            rel_dist    = rel_dist_az[0]/1000  # km
            rel_back_az = rel_dist_az[1]       # radians

            if NS == False:
                del_distR = rel_dist * math.cos((rel_back_az - ref_back_az)* math.pi/180)
                del_distT = rel_dist * math.sin((rel_back_az - ref_back_az)* math.pi/180)
            # North and east
            else:
                del_distR = rel_dist * math.cos( rel_back_az * math.pi/180)
                del_distT = rel_dist * math.sin( rel_back_az * math.pi/180)
            for slowR_i in range(slowR_n):  # for this station, loop over radial slownesses
                for slowT_i in range(slowT_n):  # loop over transverse slownesses
                    time_lag  = del_distR * stack_Rslows[slowR_i]  # time shift due to radial slowness
                    time_lag += del_distT * stack_Tslows[slowT_i]  # time shift due to transverse slowness
                    time_correction = ((t-tr.stats.starttime) + (time_lag + start_buff))/dt
                    indx = int(round(slowR_i*slowT_n + slowT_i))
                    # could do a little better by sampling finer, 20 sps?, before applying statics in pro3

                    if stack_option == 0:  # my old inefficient method
                        for it in range(stack_nt):  # check points one at a time
                            it_in = int(round(it + time_correction))
                            if it_in >= 0 and it_in < nt - 1: # does data lie within seismogram?
                                stack[indx].data[it] += tr[it_in]

                    if stack_option == 1:  #  Wei's much faster method
                        arr = tr.data
                        nshift = round(time_correction)
                        if time_correction < 0:
                            nshift = nshift - 1
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
                        if nend1 - nbeg1 != nend2 - nbeg2:
                            print('nend1 - nbeg1 != nend2 - nbeg2:  nbeg1 ' + str(nbeg1) + ', nend1 '+ str(nend1) + ', nbeg2 ' + str(nbeg2) + ', nend2 ' + str(nend2))
                        # print('str(len(arr)) < nend1 or len(arr) < nend2:  len(arr) ' + str(len(arr)) + ', nend1 ' + str(nend1) + ', nend2 ' + str(nend2))
                        if len(arr) < nend1 or len(arr) < nend2:
                            # print('str(len(arr)) < nend1 or len(arr) < nend2:  len(arr) ' + str(len(arr)) + ', nend1 ' + str(nend1) + ', nend2 ' + str(nend2))
                            # print('Sorry, fast code cannot handle running into end of trace')
                            # print('A trace from ' + tr.stats.station + ' is rejected')
                            bad_trace = True
                            continue
                            # print(tr.stats.station + ' being rejected, 1st break')
                            # sys.exit(-1)
                        if nend1 >= 0 and nbeg1 <= stack_nt:
                            stack[indx].data[nbeg1 : nend1] += arr[nbeg2 : nend2]
            done += 1
            if done == 1:
                print('Done stacking ' + str(done) + ' out of ' + str(len(st)) + ' stations.')
                elapsed_time_wc = time.time() - start_time_wc
                print(f'So far it has taken   {elapsed_time_wc:.1f}   seconds')
            elif done == 10:
                print('Done stacking ' + str(done) + ' out of ' + str(len(st)) + ' stations.')
                elapsed_time_wc = time.time() - start_time_wc
                print(f'So far it has taken   {elapsed_time_wc:.1f}   seconds')
            elif done == 100:
                print('Done stacking ' + str(done) + ' out of ' + str(len(st)) + ' stations.')
                elapsed_time_wc = time.time() - start_time_wc
                print(f'So far it has taken   {elapsed_time_wc:.1f}   seconds')
            elif done%200 == 0:
                print('Done stacking ' + str(done) + ' out of ' + str(len(st)) + ' stations.')
                elapsed_time_wc = time.time() - start_time_wc
                print(f'So far it has taken   {elapsed_time_wc:.1f}   seconds')

        else:
            print(tr.stats.station + ' not found in station list')
    if bad_trace == True:
        print(colored('There was a trace not long enough for stacking, check start and stop times', 'red'))

#%% take envelope, decimate envelope
    stack_raw = stack.copy()
    for slowR_i in range(slowR_n):  # loop over radial slownesses
        for slowT_i in range(slowT_n):  # loop over transverse slownesses
            indx = slowR_i*slowT_n + slowT_i
            stack[indx].data = np.abs(hilbert(stack[indx].data))
            if decimate_fac != 0:
                stack[indx].decimate(decimate_fac, no_filter=True)

#%%  Save processed files
    fname = 'HD' + date_label + '_2dstack_env.mseed'
    stack.write(fname,format = 'MSEED')

    fname = 'HD' + date_label + '_2dstack.mseed'
    stack_raw.write(fname,format = 'MSEED')

    elapsed_time_wc = time.time() - start_time_wc
    print(f'This job took   {elapsed_time_wc:.1f}   seconds')
    os.system('say "Five"')
