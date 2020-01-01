#!/usr/bin/env python
# input is set of hinet traces
# this program tapers, filters, selects range and SNR
# plots against traveltime curves, either raw or reduced against traveltimes
# This programs deals with a single event.
# John Vidale 2/2019

import os
os.chdir('/Users/vidale/Documents/GitHub/Array_codes')
from pro3b_sort_plot_singlet import pro3singlet
from pro4_get_shifts import pro4statics
os.chdir('/Users/vidale/Documents/PyCode/LASA')

os.environ['PATH'] += os.pathsep + '/usr/local/bin'
ARRAY      = 1
start_buff = 10
end_buff   = 20
eq_file = 'event10.txt'
ref_phase = 'PKIKP'
freq_min = 0.4
freq_max = 1.
min_dist = 163.9
max_dist = 165.5
simple_taper = 1

#%% --LASA Cull seismic section event 1
pro3singlet(ARRAY = ARRAY, stat_corr = 1, rel_time = 1, eq_file = eq_file, simple_taper = simple_taper,
			start_buff = start_buff, end_buff = end_buff,
			plot_scale_fac = 0.1, skip_SNR = 1,
			dphase = ref_phase, dphase2 = 'PKiKP', dphase3 = 'PKIKP', dphase4 = 'PKiKP',
			freq_min = freq_min, freq_max = freq_max,
			min_dist = min_dist, max_dist = max_dist, ref_loc = 0, fig_index = 101)

#%% -- LASA make a set of statics
os.chdir('/Users/vidale/Documents/PyCode/LASA/EvLocs/')
pro4statics(eq_file = eq_file, ref_trace = 'A041z',
			dphase = ref_phase, dphase2 = 'PKiKP', dphase3 = 'PKIKP', dphase4 = 'PKiKP',
			start_corr_win = 0, end_corr_win = 23, plot_scale_fac = 0.05,
			qual_threshold = 0, corr_threshold = 0, max_time_shift = 2, ARRAY = 1,
			min_dist = min_dist, max_dist = max_dist)

#%% -- HiNet version
#os.chdir('/Users/vidale/Documents/PyCode/Hinet/Biggies/SSI_2016-05-28')
#pro4statics(eq_file = 'event.txt', stat_corr = 1, ref_trace = 'N.SZW',
#			dphase = 'PKIKP', dphase2 = 'PKiKP', dphase3 = 'PKIKP', dphase4 = 'PKiKP',
#			start_corr_win = -1, end_corr_win = 3, plot_scale_fac = 0.05,
#			qual_threshold = 0, corr_threshold = 0, max_time_shift = 2)

#%% --LASA PKIKP relative to origin time rather than predicted PKIKP time
#start_buff = -1200
#end_buff   = 1230
#pro3singlet(ARRAY = ARRAY, stat_corr = 1, rel_time = 0, eq_file = eq_file, simple_taper = simple_taper,
#			start_buff = start_buff, end_buff = end_buff,
#			plot_scale_fac = 0.1, skip_SNR = 1,
#			dphase = ref_phase, dphase2 = 'PKiKP', dphase3 = 'PKIKP', dphase4 = 'PKiKP',
#			freq_min = freq_min, freq_max = freq_max,
#			min_dist = min_dist, max_dist = max_dist, ref_loc = 0, fig_index = 101)
