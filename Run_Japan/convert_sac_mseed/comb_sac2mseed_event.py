#!/usr/bin/env python
# this program runs python code to combine sac files into an mseed file for each component
# John Vidale 9/2020

def comb_sac2mseed_event(d_in, d_out, ev_num):
    import os
    from comb_sac2mseed_comp import comb_sac2mseed_comp

    # dir_in  = '/Users/vidale/Documents/Research/IC/Japan/Tanaka_events/' + d_in
    # dir_out = '/Users/vidale/Documents/Research/IC/Japan/Tanaka_events/'

    print(ev_num + ' in comb_sac2mseed_comp')
    comb_sac2mseed_comp(d_in = d_in, d_out = d_out, ev_num = ev_num, labl = 'U')
    comb_sac2mseed_comp(d_in = d_in, d_out = d_out, ev_num = ev_num, labl = 'N')
    comb_sac2mseed_comp(d_in = d_in, d_out = d_out, ev_num = ev_num, labl = 'E')