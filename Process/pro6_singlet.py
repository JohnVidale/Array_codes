#!/usr/bin/env python
# Read in 2D stacks for two events
# Compute and save tdiff, cc, ave_amp
# John Vidale 1/2021

def pro6_singlet(eq_num, slow_delta = 0.0005,
              slowR_lo = -0.1, slowR_hi = 0.1, slowT_lo = -0.1, slowT_hi = 0.1,
              start_buff = 1040, end_buff = 1180, cc_delta = -1):

    # import obspy
    # from obspy import UTCDateTime
    from obspy import read
    from obspy import Stream
    from scipy.signal import hilbert
    import numpy as np
    from termcolor import colored
    import os
    import time
    import sys
    # from obspy import Trace

#%% Get info
    # get locations
    print(colored('Running pro6_singlet', 'cyan'))

    start_time_wc = time.time()

        #%% Input parameters and computed files
    # folder_name = '/Users/vidale/Documents/Research/IC/'
    # file = open(folder_name + 'EvLocs/' + eq_file, 'r')
    fname = '/Users/vidale/Documents/Research/IC/EvLocs/event' + str(eq_num) + '.txt'
    file = open(fname, 'r')
    lines=file.readlines()
    split_line = lines[0].split()
    # t1          = UTCDateTime(split_line1[1])
    date_label  = split_line[1][0:10]
    # date_label = '2018-04-02' # dates in filename

    #%% -- read files
    # Get saved event info, also used to name files
    # date_label = '2018-04-02' # date for filename
    goto = '/Users/vidale/Documents/Research/IC/Pro_files'
    os.chdir(goto)
    fname = 'HD' + date_label + '_2dstack.mseed'
    st = Stream()
    st = read(fname)

    amp_ave       = st.copy()  # make array for average amplitude

    # decimate several arrays to sampling of time shift
    # if cc_delta == -1:
    #     print('cc_delta not set, needed to determine decimation')
    #     sys.exit(-1)
    # dec_fac = int(round(cc_delta/st[0].stats.delta))
    # if (dec_fac - cc_delta/st[0].stats.delta)/cc_delta/st[0].stats.delta > 0.00001:
    #     print('dec_fac must be an integer, pick more suitable parameters')
    #     sys.quit()
    # print(f'decimation factor {dec_fac:.3f} original sampling  {st[0].stats.delta:.3f} correlation sampling  {cc_delta:.3f}.')
    #%% Decimate some files from input dt to output dt

    print('Beam read: event: ' + str(len(st)))
    nt = len(st[0].data)
    dt = st[0].stats.delta
    print('Event: 1st trace has ' + str(nt) + ' time pts, time sampling of '
          + str(dt) + ' thus duration ' + str((nt-1)*dt))

#%% -- Count slownesses
    # count rows and columns
    slowR_n = int(round(1 + (slowR_hi - slowR_lo)/slow_delta))  # number of slownesses
    slowT_n = int(round(1 + (slowT_hi - slowT_lo)/slow_delta))  # number of slownesses

    #  Loop over slowness
    total_slows = slowR_n * slowT_n
    print(f'{slowR_n} radial slownesses, low {slowR_lo:.4f} high {slowR_hi:.4f}')
    print(f'{slowT_n} transv slownesses, low {slowT_lo:.4f} high {slowT_hi:.4f}  total slows: {total_slows}')

#%% Find envelope and tshift
    for slow_i in range(total_slows):

        if slow_i % 1000 == 0:
            print('Envelope: ' + str(slow_i) + ' finished slownesses out of ' + str(total_slows))

        if len(st[slow_i].data) == 0: # test for zero-length traces, indexing errors
            print('Slowness ' + str(slow_i) + ' trace has zero length, problem!')

        seismogram = hilbert(st[slow_i].data)  # make analytic seismograms
        env = np.abs(seismogram) # amplitude and amp ratio
        amp_ave[slow_i].data    = env

    #%% Decimate amp files
    # if dec_fac > 1:
    #     amp_ave.decimate(dec_fac, no_filter=True)

#%%  Save processed file
    goto = '/Users/vidale/Documents/Research/IC/Pro_Files'
    os.chdir(goto)

    fname = 'HD' + date_label + '_amp_ave.mseed'
    amp_ave.write(    fname,format = 'MSEED')

    print(f'len(amp_ave)  {len(amp_ave)}  len(amp_ave[0].data)  {len(amp_ave[0].data)} dt {amp_ave[0].stats.delta}')

    elapsed_time_wc = time.time() - start_time_wc
    print(f'This job took   {elapsed_time_wc:.1f}   seconds')
    os.system('say "Six"')