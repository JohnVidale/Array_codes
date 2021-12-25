#!/usr/bin/env python
# this program runs python code to combine sac files into a single mseed file
# John Vidale 9/2020

import os
from comb_sac2mseed_event import comb_sac2mseed_event

d_out = '/Users/vidale/Documents/Research/IC/Japan/Tanaka_events/'

ev_num = '20040725_'
print(ev_num + ' in comb_run')
d_in  = '/Users/vidale/Documents/Research/IC/Japan/Tanaka_events/200407252335'
comb_sac2mseed_event(d_in = d_in, d_out = d_out, ev_num = ev_num)

ev_num = '20060128_'
print(ev_num + ' in comb_run')
d_in  = '/Users/vidale/Documents/Research/IC/Japan/Tanaka_events/200601280158'
comb_sac2mseed_event(d_in = d_in, d_out = d_out, ev_num = ev_num)

ev_num = '20070928_'
print(ev_num + ' in comb_run')
d_in  = '/Users/vidale/Documents/Research/IC/Japan/Tanaka_events/200709282238'
comb_sac2mseed_event(d_in = d_in, d_out = d_out, ev_num = ev_num)

ev_num = '20080323_'
print(ev_num + ' in comb_run')
d_in  = '/Users/vidale/Documents/Research/IC/Japan/Tanaka_events/200803230624'
comb_sac2mseed_event(d_in = d_in, d_out = d_out, ev_num = ev_num)

ev_num = '20080705_'
print(ev_num + ' in comb_run')
d_in  = '/Users/vidale/Documents/Research/IC/Japan/Tanaka_events/200807051112'
comb_sac2mseed_event(d_in = d_in, d_out = d_out, ev_num = ev_num)

ev_num = '20081124_'
print(ev_num + ' in comb_run')
d_in  = '/Users/vidale/Documents/Research/IC/Japan/Tanaka_events/200811241802'
comb_sac2mseed_event(d_in = d_in, d_out = d_out, ev_num = ev_num)

ev_num = '20090828_'
print(ev_num + ' in comb_run')
d_in  = '/Users/vidale/Documents/Research/IC/Japan/Tanaka_events/200908281051'
comb_sac2mseed_event(d_in = d_in, d_out = d_out, ev_num = ev_num)

ev_num = '20090930_'
print(ev_num + ' in comb_run')
d_in  = '/Users/vidale/Documents/Research/IC/Japan/Tanaka_events/200909301916'
comb_sac2mseed_event(d_in = d_in, d_out = d_out, ev_num = ev_num)

ev_num = '20100724_'
print(ev_num + ' in comb_run')
d_in  = '/Users/vidale/Documents/Research/IC/Japan/Tanaka_events/201007240815'
comb_sac2mseed_event(d_in = d_in, d_out = d_out, ev_num = ev_num)