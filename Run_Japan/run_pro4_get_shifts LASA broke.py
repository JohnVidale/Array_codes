#!/usr/bin/env python
# input is set of hinet traces
# this program tapers, filters, selects range and SNR
# plots against traveltime curves, either raw or reduced against traveltimes
# This programs deals with a single event.
# John Vidale 2/2019
# re-visited to align repeating events on Hi-Net 4/2022

import os
import time
import matplotlib.pyplot as plt

#%% Import functions
from pro3_sort_plot_singlet import pro3singlet
from pro4_get_shifts import pro4statics

start_time_wc = time.time()

# os.chdir('/Users/vidale/Documents/PyCode/LASA')

# os.environ['PATH'] += os.pathsep + '/usr/local/bin'

ARRAY      = 0
start_buff = 10
end_buff   = 20
eq_file = 'event161.txt'
eq_num = 161
# ref_phase = 'PKIKP'
freq_min = 0.4
freq_max = 1.
min_dist = 163.9
max_dist = 165.5
simple_taper = 1
stat_corr = True

#%% -- HiNet version
# pro3singlet(ARRAY = ARRAY, stat_corr = stat_corr, rel_time = 1, eq_file = eq_file, simple_taper = simple_taper,
# 			start_buff = start_buff, end_buff = end_buff,
# 			plot_scale_fac = 0.1, skip_SNR = 1,
# 			dphase = ref_phase, dphase2 = 'PKiKP', dphase3 = 'PKIKP', dphase4 = 'PKiKP',
# 			freq_min = freq_min, freq_max = freq_max,
# 			min_dist = min_dist, max_dist = max_dist, ref_loc = 0, fig_index = 101)

pro4statics(eq_num = eq_num, stat_corr = stat_corr, ref_trace = 'N.SZW',
			dphase = 'PKiKP', dphase2 = 'PKP', dphase3 = 'PKIKP', dphase4 = 'iPKiKP',
			start_corr_win = -1, end_corr_win = 3, plot_scale_fac = 0.05,
			qual_threshold = 0, corr_threshold = 0, max_time_shift = 2)

#%% --LASA Cull seismic section event 1
# pro3singlet(ARRAY = ARRAY, stat_corr = stat_corr, rel_time = 1, eq_file = eq_file, simple_taper = simple_taper,
# 			start_buff = start_buff, end_buff = end_buff,
# 			plot_scale_fac = 0.1, skip_SNR = 1,
# 			dphase = ref_phase, dphase2 = 'PKiKP', dphase3 = 'PKIKP', dphase4 = 'PKiKP',
# 			freq_min = freq_min, freq_max = freq_max,
# 			min_dist = min_dist, max_dist = max_dist, ref_loc = 0, fig_index = 101)

#%% -- LASA make a set of statics
# os.chdir('/Users/vidale/Documents/PyCode/LASA/EvLocs/')
# pro4statics(eq_file = eq_file, ref_trace = 'A041z',
# 			dphase = ref_phase, dphase2 = 'PKiKP', dphase3 = 'PKIKP', dphase4 = 'PKiKP',
# 			start_corr_win = 0, end_corr_win = 23, plot_scale_fac = 0.05,
# 			qual_threshold = 0, corr_threshold = 0, max_time_shift = 2, ARRAY = 1,
# 			min_dist = min_dist, max_dist = max_dist)

#%% --LASA PKIKP relative to origin time rather than predicted PKIKP time
#start_buff = -1200
#end_buff   = 1230
#pro3singlet(ARRAY = ARRAY, stat_corr = stat_corr, rel_time = 0, eq_file = eq_file, simple_taper = simple_taper,
#			start_buff = start_buff, end_buff = end_buff,
#			plot_scale_fac = 0.1, skip_SNR = 1,
#			dphase = ref_phase, dphase2 = 'PKiKP', dphase3 = 'PKIKP', dphase4 = 'PKiKP',
#			freq_min = freq_min, freq_max = freq_max,
#			min_dist = min_dist, max_dist = max_dist, ref_loc = 0, fig_index = 101)
elapsed_time_wc = time.time() - start_time_wc
print(f'    Get shift code took   {elapsed_time_wc:.1f}   seconds')
os.system('say "done"')
