#!/usr/bin/env python
# Read in 2D stacks for one event
# Compute ave_amp
# Write out tdiff, ave_amp results
# John Vidale 3/2019

def pro6_cc_pair(eq_num1, eq_num2, repeater = '0', slow_delta = 0.0005, Spyder = True,
              slowR_lo = -0.1, slowR_hi = 0.1, slowT_lo = -0.1, slowT_hi = 0.1,
              start_buff = 1040, end_buff = 1180,
              cc_twin = 2, cc_len = 0.5, cc_interp1d = 5, cc_delta = 0.1, cc_thres = 0.8):

#%% Import functions
    # import obspy
    # from obspy import UTCDateTime
    from obspy import read
    from obspy import Stream
    # from obspy import Trace
    from scipy.signal import hilbert
    import numpy as np
    import os
    import time
    import sys
    import statistics
    import math
    import pandas as pd
    from termcolor import colored
    file_directory = '/Users/vidale/Documents/GitHub/Array_codes/Process'
    os.chdir(file_directory)
    from pro_proceed_john import cc_measure_tshift

#%% Get info
    # get locations
    print(colored('Running pro6_pair_cc', 'cyan'))

    start_time_wc = time.time()

    def search_df(df, column, value, partial_match=True):
        df = df.astype({column:'string'})
        if partial_match:
            return df.loc[df[column].str.contains(value, na=False)]
        else:
            return df.loc[df[column] == value]

    # look up pair of earthquakes and time shifts in pairs
    df = pd.read_excel('/Users/vidale/Documents/GitHub/Array_codes/Files/ICevents_full.xlsx', sheet_name='pairs')
    lines0       = search_df(df,'label'      ,repeater,partial_match=True)
    eq_num1      = lines0.index1.iloc[0]
    eq_num2      = lines0.index2.iloc[0]

    # read origin times for that pair in events
    df = pd.read_excel('/Users/vidale/Documents/GitHub/Array_codes/Files/ICevents_full.xlsx', sheet_name='events')
    lines1 = search_df(df,'INDEX',str(eq_num1),partial_match=True)
    lines2 = search_df(df,'INDEX',str(eq_num2),partial_match=True)

    time1 = lines1.TIME.iloc[0]
    time2 = lines2.TIME.iloc[0]

    #  new lines to match more specific naming
    date_label1  = time1[0:10]
    date_label2  = time2[0:10]
    # date_label = '2018-04-02' # dates in filename

    #%% -- read files
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

    test1 = statistics.stdev(st1[0].data)  # check for NANs, quit if bad
    test2 = statistics.stdev(st2[0].data)
    if math.isnan(test1) or math.isnan(test2):
        print('Stdev first row in st1 and st2 are ' + str(test1) + ' and ' + str(test2))
        if math.isnan(test1):
            print(colored('st1 has NANs', 'yellow'))
        if math.isnan(test2):
            print(colored('st2 has NANs', 'yellow'))
        sys.exit(-1)

    tshift        = st1.copy()  # make array for time shift
    tshift_cc     = st1.copy()  # make array for new time shift
    cc            = st1.copy()  # make array for cc coefficient
    amp_ave       = st1.copy()  # make array for average amplitude

    # decimate several arrays to sampling of time shift
    dec_fac = int(round(cc_delta/st1[0].stats.delta))
    if (dec_fac - cc_delta/st1[0].stats.delta)/cc_delta/st1[0].stats.delta > 0.00001:
        print('dec_fac must be an integer, pick more suitable parameters')
        sys.quit()
    print(f'decimation factor {dec_fac:.3f} original sampling  {st1[0].stats.delta:.3f} correlation sampling  {cc_delta:.3f}.')
    #%% Decimate some files from input dt to output dt
    print(f'0 Before len(tshift)  {len(tshift)}  len(tshift[0].data)  {len(tshift[0].data)} dt {tshift[0].stats.delta}')
    if dec_fac > 1:
        tshift.decimate(   dec_fac, no_filter=True)
        tshift_cc.decimate(dec_fac, no_filter=True)
        cc.decimate(       dec_fac, no_filter=True)
    print(f'0 After len(tshift)  {len(tshift)}  len(tshift[0].data)  {len(tshift[0].data)} dt {tshift[0].stats.delta}')

    print('Beams read: event1: ' + str(len(st1)) + ' event2: ' + str(len(st2)))
    nt1 = len(st1[0].data)
    nt2 = len(st2[0].data)
    dt1 = st1[0].stats.delta
    dt2 = st2[0].stats.delta
    print('Event1: 1st trace has ' + str(nt1) + ' time pts, time sampling of '
          + str(dt1) + ' thus duration ' + str((nt1-1)*dt1))
    print('Event2: 1st trace has ' + str(nt2) + ' time pts, time sampling of '
          + str(dt2) + ' thus duration ' + str((nt2-1)*dt2))
    if nt1 != nt2 or dt1 != dt2:
        print('Trouble, nt or dt not does not match')
        sys.exit(-1)

#%% -- Count slownesses
    # count rows and columns
    slowR_n = int(round(1 + (slowR_hi - slowR_lo)/slow_delta))  # number of slownesses
    slowT_n = int(round(1 + (slowT_hi - slowT_lo)/slow_delta))  # number of slownesses

    #  Loop over slowness
    total_slows = slowR_n * slowT_n
    print(f'{slowR_n} radial slownesses, low {slowR_lo:.4f} high {slowR_hi:.4f}')
    print(f'{slowT_n} transv slownesses, low {slowT_lo:.4f} high {slowT_hi:.4f}  total slows: {total_slows}')

    ttt             = start_buff + (np.arange( st1[0].stats.npts) * st1[0].stats.delta) # time array for inputs
    cc_ttt_full     = start_buff + (np.arange((st1[0].stats.npts-1)/dec_fac + 1) * st1[0].stats.delta * dec_fac) # time array for outputs
    print(f'input npts {len(ttt)} end points {ttt[0]:.3f} {ttt[-1]:.3f} ')
    print(f'output npts {len(cc_ttt_full)}  end points {cc_ttt_full[0]:.3f} {cc_ttt_full[-1]:.3f} ')

#%% Find envelope and tshift
    for slow_i in range(total_slows):

        if slow_i % 100 == 0:
            print('Measuring time shifts, finished ' + str(slow_i) + ' of ' + str(total_slows) + ' slownesses')

        if len(st1[slow_i].data) == 0: # test for zero-length traces, indexing errors
            print('Slowness ' + str(slow_i) + ' trace has zero length, problem!')

        seismogram1 = hilbert(st1[slow_i].data)  # make analytic seismograms
        seismogram2 = hilbert(st2[slow_i].data)

        env1 = np.abs(seismogram1) # amplitude and amp ratio
        env2 = np.abs(seismogram2)
        amp_ave[slow_i].data = 0.5 * (env1 + env2)

        tr1 = st1[slow_i]
        tr2 = st2[slow_i]
        cc_ttt, cc_coef, tshift_new = cc_measure_tshift( tr1=tr1, tr2=tr2, tarr_beg = start_buff,
                          cc_twin=cc_twin, cc_len=cc_len, cc_delta=cc_delta, cc_interp1d=cc_interp1d)

        # restore arrays to full length by adding zeroes to ends
        zero_fill1 = int(round((len(tshift_cc[0].data) - len(tshift_new)) / 2))  # zeroes or NAN to be added to each end
        misfit = len(cc[0].data) - (2*zero_fill1 + len(cc_coef))
        if misfit == 0:  # if total zeroes to fill is odd, different fill is required front vs back end
            zero_fill2 = zero_fill1
        elif misfit == 1:
            zero_fill2 = zero_fill1 + 1
        elif misfit == -1:
            zero_fill2 = zero_fill1 - 1
        else:
            print(f'zero_fill1: {zero_fill1}, len(tshift_cc[0].data {len(tshift_cc[0].data)}, len(tshift_new): {len(tshift_new)}')
            print(f'len(cc[0].data): {len(cc[0].data)}, len(cc_coef) {len(cc_coef)}')
            print(f'misfit == {misfit}, this is not going to work, fix to make cc_coef + zero_fill1 + zero_fill2 = data length')
            sys.exit()
        cc_coef_full = np.concatenate([np.zeros(zero_fill1),    cc_coef, np.zeros(zero_fill2)])
        tshift_full  = np.concatenate([np.zeros(zero_fill1), tshift_new, np.zeros(zero_fill2)])

        tshift[slow_i].data =  tshift_full  # tshift just gets  NaNs for plotting, full is file that is saved
        cc[    slow_i].data = cc_coef_full

    #%% Decimate amp files
    if dec_fac > 1:
        amp_ave.decimate(dec_fac, no_filter=True)

#%%  Save processed file
    goto = '/Users/vidale/Documents/Research/IC/Pro_Files'
    os.chdir(goto)

    fname = 'HD' + date_label1 + '_' + date_label2 + '_tshift.mseed'
    tshift.write(fname,format = 'MSEED')

    fname = 'HD' + date_label1 + '_' + date_label2 + '_amp_ave.mseed'
    amp_ave.write(    fname,format = 'MSEED')

    fname = 'HD' + date_label1 + '_' + date_label2 + '_cc.mseed'
    cc.write(         fname,format = 'MSEED')

    print(f'len(tshift)  {len(tshift)}  len(tshift[0].data)  {len(tshift[0].data)} dt {tshift[0].stats.delta}')
    print(f'len(amp_ave)  {len(amp_ave)}  len(amp_ave[0].data)  {len(amp_ave[0].data)} dt {amp_ave[0].stats.delta}')
    print(f'len(cc)  {len(cc)}  cc[0].data)  {len(cc[0].data)} dt {cc[0].stats.delta}')

    elapsed_time_wc = time.time() - start_time_wc
    print(f'This job took   {elapsed_time_wc:.1f}   seconds')
    os.system('say "six done"')
