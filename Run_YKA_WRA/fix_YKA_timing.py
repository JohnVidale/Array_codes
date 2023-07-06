#!/usr/bin/env python
# adjusts timing for 0.125s jumps and transition statics
# John Vidale 11/2022

def fix_YKA_timing(eq_num = 737):

    import os
    pro_directory = '/Users/vidale/Documents/GitHub/Array_codes/Run_YKA_WRA'
    os.chdir(pro_directory)

    #%% Import
    #%% -- Functions
    from obspy import UTCDateTime
    from obspy import Stream
    from obspy import read
    import os
    import time
    from termcolor import colored

    # eq_num = 737

    # print(colored('Running fix_YKA_timing', 'cyan'))
    start_time_wc = time.time()

    #%% -- Event info
    #  input event data with 1-line file of format
    #  event 2016-05-28T09:47:00.000 -56.241 -26.935 78
    fname = '/Users/vidale/Documents/Research/IC/EvLocs/event' + str(eq_num) + '.txt'
    # print('Opening event file ' + fname)
    file = open(fname, 'r')
    lines=file.readlines()
    split_line = lines[0].split()
    #            ids.append(split_line[0])  ignore label for now
    t           = UTCDateTime(split_line[1])
    date_label  = split_line[1][0:10]
    year_label  = split_line[1][0:4]
    year_short_label  = split_line[1][2:4]
    month_label   = split_line[1][5:7]
    day_label     = split_line[1][8:10]
    hour_label    = split_line[1][11:13]
    minute_label  = split_line[1][14:16]

    #%% -- Station locations and statics
    sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/sta_CN_YK.txt'
    with open(sta_file, 'r') as file:
        lines = file.readlines()
    # print('    ' + str(len(lines)) + ' stations read from ' + sta_file)
    # Load station coords into arrays
    station_index = range(len(lines))
    st_names = []
    st_lats  = []
    st_lons  = []
    for ii in station_index:
        line = lines[ii]
        split_line = line.split()
        # print('input station ' + split_line[0])
        st_names.append(split_line[0])
        st_lats.append( split_line[1])
        st_lons.append( split_line[2])

    # print('Event: date_label ' + date_label)

    #%% Load waveforms
    st = Stream()
        # fname     = '/Users/vidale/Documents/GitHub/LASA_data/HD' + date_label + '.mseed'
    mseed_name = year_label + month_label + day_label + '_' + hour_label + minute_label
    fname     = '/Users/vidale/Documents/Research/IC/Mseed/YKA/before_shift/' + mseed_name + '.mseed'

    # print('Opening data file: ' + fname)
    st=read(fname)
    # print('data len ' + str(len(st[0].data)) + ' # traces ' + str(len(st)) + ' delta ' + str(st[0].stats.delta))
    #%% Select and process traces
    stgood = Stream()

    tra_sta_found = 0
    nodata = 0
    for tr in st: # traces one by one, find lat-lon
        if tr.stats.station in st_names:  # find station in station list
            tra_sta_found += 1

    #%% -- Apply static
            st_name = tr.stats.station
            # correction for shift to new instrumentation 2013, no idea what it's intermittent for 6 years
            shifter = 0

            if eq_num <= 746 or (eq_num >= 800 and eq_num <= 827):

                shifter = 0.1

            # correction for missed data packets, 2014 to early 2020, picked by hand
            if eq_num == 747:
                if st_name == 'YKAR4' or st_name == 'YKAB6' or st_name == 'YKAR1':
                    shifter = shifter -0.125
            if eq_num == 748:
                if st_name == 'YKAR1':
                    shifter = shifter +0.125
                if st_name == 'YKAR2' or st_name == 'YKAR3' or st_name == 'YKAR4' or st_name == 'YKAB7':
                    shifter = shifter -0.125
            if eq_num == 749:
                if st_name == 'YKAR2' or st_name == 'YKAR3' or st_name == 'YKAR4' or st_name == 'YKAR5' or st_name == 'YKAB7':
                    shifter = shifter -0.125
                if st_name == 'YKAB6':
                    shifter = shifter -0.250
            if eq_num == 750:
                if st_name == 'YKAR3' or st_name == 'YKAR4' or st_name == 'YKAR5' or st_name == 'YKAB7':
                    shifter = shifter -0.125
                if st_name == 'YKAB6':
                    shifter = shifter -0.250
            if eq_num == 751:
                if st_name == 'YKAR4' or st_name == 'YKAR5' or st_name == 'YKAR3' or st_name == 'YKAB7':
                    shifter = shifter -0.125
                if st_name == 'YKAB6':
                    shifter = shifter -0.250
                if st_name == 'YKAR8':
                    shifter = shifter + 0.125
            if eq_num == 753:
                if st_name == 'YKAB2' or st_name == 'YKAR5' or st_name == 'YKAR6' or st_name == 'YKAB9':
                    shifter = shifter -0.125
            if eq_num == 754:
                if st_name == 'YKAB2' or st_name == 'YKAB3' or st_name == 'YKAB7' or st_name == 'YKAB9' or st_name == 'YKAR5':
                    shifter = shifter -0.125
            if eq_num == 755:
                if st_name == 'YKAB2' or st_name == 'YKAR3' or st_name == 'YKAB7' or st_name == 'YKAB9' or st_name == 'YKAR5':
                    shifter = shifter -0.125
            if eq_num == 756:
                if st_name == 'YKAB2' or st_name == 'YKAB7' or st_name == 'YKAR3' or st_name == 'YKAR5':
                    shifter = shifter -0.125
            if eq_num == 757:
                if st_name == 'YKAB1' or st_name == 'YKAB2' or st_name == 'YKAB9' or st_name == 'YKAR3' or st_name == 'YKAR5' or st_name == 'YKAR6':
                    shifter = shifter -0.125
                if st_name == 'YKAB7':
                    shifter = shifter -0.250
            if eq_num == 758:
                if st_name == 'YKAB1' or st_name == 'YKAB2' or st_name == 'YKAB6' or st_name == 'YKAB9' or st_name == 'YKAR6':
                    shifter = shifter -0.125
                if st_name == 'YKAR3' or st_name == 'YKAB7':
                    shifter = shifter -0.250
            if eq_num == 832:
                if st_name == 'YKAB7' or st_name == 'YKAB8' or  st_name == 'YKAR2' or st_name == 'YKAR4':
                    shifter = shifter -0.125
                if st_name == 'YKAB1':
                    shifter = shifter -0.250
                if st_name == 'YKAR1':
                    shifter = shifter + 0.125
            if eq_num == 833:
                if st_name == 'YKAR3' or st_name == 'YKAR4' or st_name == 'YKAR5' or st_name == 'YKAB7':
                    shifter = shifter -0.125
                if st_name == 'YKAR1' or st_name == 'YKAB6':
                    shifter = shifter -0.250
                if st_name == 'YKAR8':
                    shifter = shifter + 0.125
            if eq_num == 836:
                if st_name == 'YKAB1':
                    shifter = shifter -0.125
                if st_name == 'YKAB2' or st_name == 'YKAR5' or st_name == 'YKAR6':
                    shifter = shifter -0.250
            if eq_num == 837:
                if st_name == 'YKAB2' or st_name == 'YKAB9' or st_name == 'YKAR5' or st_name == 'YKAR6' or st_name == 'YKAR9':
                    shifter = shifter -0.125
            if eq_num == 838:
                if st_name == 'YKAB2' or st_name == 'YKAB3' or st_name == 'YKAB7' or st_name == 'YKAB9' or st_name == 'YKAR5':
                    shifter = shifter -0.125
            if eq_num == 839:
                if st_name == 'YKAB2' or st_name == 'YKAR3' or st_name == 'YKAB7' or st_name == 'YKAB9' or st_name == 'YKAR5':
                    shifter = shifter -0.125
            if eq_num == 840:
                if st_name == 'YKAB9' or st_name == 'YKAR1' or st_name == 'YKAR3' :
                    shifter = shifter -0.125
            if eq_num == 842:
                if st_name == 'YKAB2' or st_name == 'YKAB4' or st_name == 'YKAB6' or st_name == 'YKAB7' or st_name == 'YKAB9' or st_name == 'YKAR5' or st_name == 'YKAR6':
                    shifter = shifter -0.125
                if st_name == 'YKAR3':
                    shifter = shifter -0.250
            if eq_num == 843:
                if st_name == 'YKAB2' or st_name == 'YKAB4' or st_name == 'YKAB9':
                    shifter = shifter -0.125
                if st_name == 'YKAR6' or st_name == 'YKAR5' or st_name == 'YKAB7' or st_name == 'YKAR3':
                    shifter = shifter -0.250
            if eq_num == 844:
                if st_name == 'YKAB2' or st_name == 'YKAB3' or st_name == 'YKAB9' or st_name == 'YKAR5' or st_name == 'YKAR8':
                    shifter = shifter -0.125
                if st_name == 'YKAB7' or st_name == 'YKAR3':
                    shifter = shifter -0.250
            if eq_num == 845:
                if st_name == 'YKAB1' or st_name == 'YKAB2' or st_name == 'YKAB6' or st_name == 'YKAR5' or st_name == 'YKAR6':
                    shifter = shifter -0.125
                if st_name == 'YKAR8':
                    shifter = shifter -0.250
                if st_name == 'YKAR4':
                    shifter = shifter + 0.125 # oddity
                if st_name == 'YKAB4' or st_name == 'YKAB7':
                    shifter = shifter -0.375
            tr.stats.starttime += shifter
    #         print('static correction is ' + str(shifter) + ' for station ' + tr.stats.station)
    # #%% -- Reject if not in time window
            if len(tr.data) > 0:
                stgood += tr
            else:
                nodata += 1
    # #%% -- Reject if not in station (static) list
    #     else:
    #         print(tr.stats.station + ' not found in station list with statics')
    # print('    ' + str(tra_sta_found) + '  traces')
    # print('        ' + str(nodata) + '  traces with no data,  ')
    # print('    ' + str(len(stgood)) + '  traces after check for data, tt calc')

    #%%  Save processed seismograms
    fname3 = '/Users/vidale/Documents/Research/IC/Mseed/YKA/' + mseed_name + '.mseed'

    stgood.write(fname3,format = 'MSEED')

    elapsed_time_wc = time.time() - start_time_wc
    # print(f'    This job took   {elapsed_time_wc:.1f}   seconds')
    os.system('say "three"')
