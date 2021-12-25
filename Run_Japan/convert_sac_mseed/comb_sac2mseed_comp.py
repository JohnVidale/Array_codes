#!/usr/bin/env python
# This programs reads all SAC files in a directory matching the search string,
#   then writes out a single mseed file.
# John Vidale 3/2021

def comb_sac2mseed_comp(d_in, d_out, ev_num, labl):
    import os
    from obspy import read, Stream
    import glob

    os.chdir(d_in)
    search_field = '*' + labl + '.SAC'   # for syn files
    st = Stream()

    comb_sac2mseed_comp
    file_list = glob.glob(search_field)
#    print('Number of files ' + str(len(file_list)) + '; first -  ' + file_list[0])
    for sgrams in file_list:
        new_one = read(sgrams)
        st += new_one
        # print('Channel is ' + str(new_one[0].stats.network))

    cnt = len(st)
    print(labl + ' traces output ' + str(cnt))
    if cnt > 0:
        os.chdir(d_out)
        st.write(ev_num + labl + '.mseed', format='MSEED')
