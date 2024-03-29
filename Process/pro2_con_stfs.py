#!/usr/bin/env python
# this program convolves a time function with an mseed file
# John Vidale 6/2019

def pro2_convstf(eq_num, conv_file):

    from obspy import UTCDateTime
    from obspy import Stream, Trace
    from obspy import read
#    from obspy.signal import correlate_template
    import os
    import time
    import numpy as np

    import sys
    import warnings

    # if not sys.warnoptions: # don't show any warnings
    #     warnings.simplefilter('ignore')

    print('Running pro2_con_stfs')
    start_time_wc = time.time()

    #%%
    taper_frac = .05      #Fraction of window tapered on both ends
#    conv_file = 'HD1971-11-06_stf'

    #%% input event data with 1-line file of format
    #  event 2016-05-28T09:47:00.000 -56.241 -26.935 78
    # folder_name = '/Users/vidale/Documents/Research/IC/'
    # file = open(folder_name + 'EvLocs/' + eq_file, 'r')
    fname = '/Users/vidale/Documents/Research/IC/EvLocs/event' + str(eq_num) + '.txt'
    file = open(fname, 'r')

    lines=file.readlines()
    split_line = lines[0].split()
#            ids.append(split_line[0])  ignore label for now
    t           = UTCDateTime(split_line[1])
    date_label  = split_line[1][0:10]
    print('date_label ' + date_label + ' time ' + str(t))

    #%% Load waveforms and convolution trace
    st        = Stream()
    con_trace = Stream()
    st_out    = Stream()
    tr = Trace()

    fname_sel     = '/Users/vidale/Documents/Research/IC/Pro_Files/HD' + date_label + 'sel.mseed'

    st        = read(fname_sel)
    fname     = conv_file
    con_trace = read(fname)
    con_trace.taper(0.5)  # added June 10, 2019 to shorten stf


    nt = len(st[0].data)
    dt = st[0].stats.delta
    print('Read in:\n' + str(len(st)) + ' traces' + ' from file ' + fname +
       ', \n' + str(nt) + ' time pts, time sampling of '
          + str(dt) + ' and thus duration of ' + str((nt-1)*dt))

        #%%  detrend, taper
    st.detrend(type='simple')
    st.taper(taper_frac)

#    st3 = Stream()
#    st3 = correlate_template(st, con_trace, normalize='none')

    done = 0

    for tr in st: # traces one by one, find lat-lon by searching entire inventory.  Inefficient but cheap
#        print('con_trace data has length ' + str(len(con_trace[0].data)))
#        print('Tr data has length ' + str(len(tr.data)) + ' con_trace data has length ' + str(len(con_trace[0].data)))
        tr.data = np.convolve(tr.data, con_trace[0].data)
#        print('Now, Tr data has length ' + str(len(tr.data)))
        tr.stats.starttime = tr.stats.starttime - 9 # shift timing to reflect convolution delay
        st_out += tr
        done += 1
        if done%50 == 0:
            print('Done stacking ' + str(done) + ' out of ' + str(len(st)) + ' stations.')

    nt = len(st_out[0].data)
    dt = st_out[0].stats.delta
    print('After detrend and taper:\n' + str(len(st_out)) + ' traces written to file ' + fname_sel +
       ', ' + str(nt) + ' time pts, time sampling of '
          + str(dt) + ' and thus duration of ' + str((nt-1)*dt))

    #  Save processed files
    st_out.write(fname_sel,format = 'MSEED')

    elapsed_time_wc = time.time() - start_time_wc
    print('This job took ' + str(elapsed_time_wc) + ' seconds')

    os.system('say "Done"')