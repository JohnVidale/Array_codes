#!/usr/bin/env python
# input is set of hinet traces
# this program tapers, filters, selects range and SNR
# plots against traveltime curves, either raw or reduced against traveltimes
# This programs deals with a single event.
# John Vidale 2/2019

import os

os.environ['PATH'] += os.pathsep + '/usr/local/bin'
os.chdir('/Users/vidale/Documents/GitHub/Array_codes')

#%% Import functions
from pro3a_sort_plot_pair    import pro3pair
from pro3b_sort_plot_singlet import pro3singlet
from pro5a_stack             import pro5stack
from pro5b_stack2d           import pro5stack2d
from pro6_plot_stacked_seis  import pro6stacked_seis
from pro7a_plot_envstack     import pro7plotstack
from pro7b_plot_stack        import pro7plotstack2
from pro7c_plot_stack        import pro7plotstack3
from pro7b_dec               import pro7dec
os.chdir('/Users/vidale/Documents/PyCode/LASA')

#%% Common parameters
ARRAY      = 1
eq_file1   = 'event1.txt'
eq_file2   = 'event2.txt'

# PcP
#start_buff = 640
#end_buff   = 680

# PKiKP
#start_buff = 1030
#end_buff   = 1160

# PKKP
#start_buff = 1875
#end_buff   = 1890

# P'P'
#start_buff = 2370
#end_buff   = 2395

# Full array
min_dist = 58.5
max_dist = 60.5

# Just inner rings
#min_dist = 59.2
#max_dist = 59.75

# HF
freq_min = 1
freq_max = 3

# LF
#freq_min = 0.5
#freq_max = 2

slowR_lo   =  0.02
slowR_hi   =  0.12
slowT_lo   = -0.05
slowT_hi   =  0.05
slow_delta =  0.0025
R_slow_plot = 0.04
T_slow_plot = 0.00

max_rat = 10
min_amp = 0.2
tdiff_clip = 0.2

slowR_lo_1D = -0.04
slowR_hi_1D = 0.1
slow_delta_1D = 0.001

decimate_fac   =  5
simple_taper   =  1
skip_SNR       =  1
snaptime = 660
snaps = 1
freq_corr = 1.2
ref_phase = 'PKiKP'
stat_corr = 1

#%% Comparing events
#%% --Cull seismic section for common stations
pro3pair(ARRAY = ARRAY, eq_file1 = eq_file1, eq_file2 = eq_file2, simple_taper = simple_taper, skip_SNR = skip_SNR,
			rel_time = 0, start_buff = start_buff, end_buff = end_buff,
			freq_min = freq_min, freq_max = freq_max,
			plot_scale_fac = 0.025, stat_corr = stat_corr,
			dphase = ref_phase, dphase2 = 'PKKP', dphase3 = 'PP', dphase4 = 'S',
			min_dist = min_dist, max_dist = max_dist, ref_loc = 0)

##%%  --2D stacks
pro5stack2d(eq_file = eq_file1, plot_scale_fac = 0.2,
			slowR_lo = slowR_lo, slowR_hi = slowR_hi, slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
			start_buff = start_buff, end_buff = end_buff,
			norm = 1, global_norm_plot = 1,
			ARRAY = ARRAY, decimate_fac = decimate_fac, NS = 0)

pro5stack2d(eq_file = eq_file2, plot_scale_fac = 0.2,
			slowR_lo = slowR_lo, slowR_hi = slowR_hi, slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
			start_buff = start_buff, end_buff = end_buff,
			norm = 1, global_norm_plot = 1,
			ARRAY = ARRAY, decimate_fac = decimate_fac, NS = 0)

##%% --Compare pair of 2D stack results
pro6stacked_seis(eq_file1 = eq_file1, eq_file2 = eq_file2, plot_scale_fac = 0.002,
			slowR_lo = slowR_lo, slowR_hi = slowR_hi, slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
			start_buff = start_buff, end_buff = end_buff, freq_corr = freq_corr, ref_phase = ref_phase,
			fig_index = 301, plot_dyn_range = 100, ARRAY = ARRAY,
			turn_off_black = 0, R_slow_plot = R_slow_plot, T_slow_plot = T_slow_plot,
			max_rat = max_rat, min_amp = min_amp, tdiff_clip = tdiff_clip)

##%% --decimate stacking files to shorten processing time
#pro7dec(eq_file1 = eq_file1, eq_file2 = eq_file2, decimate_fac = decimate_fac, ARRAY = ARRAY)

#%% --slowness slices of time shift
#pro7plotstack2(eq_file1 = eq_file1, eq_file2 = eq_file2, plot_scale_fac = 0.05,
#			slowR_lo = slowR_lo, slowR_hi = slowR_hi, slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
#			zoom = 0, ZslowR_lo = -0.03, ZslowR_hi = 0.03, ZslowT_lo = -0.03, ZslowT_hi = 0.03, Zstart_buff = 5, Zend_buff = 15,
#			start_buff = start_buff, end_buff = end_buff, skip_T = 0, skip_R = 0, skip_snaps = 1, tdiff_clip = tdiff_clip,
#			fig_index = 301, plot_dyn_range = 100, snaptime = snaptime, snaps=snaps, decimate_fac = 1, in_dec = 1,
#			ref_phase = ref_phase, ARRAY = ARRAY,
#			max_rat = max_rat, min_amp = min_amp)

#%% --single slowness slices of time shift
#pro7plotstack3(eq_file1 = eq_file1, eq_file2 = eq_file2, plot_scale_fac = 0.02,
#			slowR_lo = slowR_lo, slowR_hi = slowR_hi, slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
#			start_buff = start_buff, end_buff = end_buff, skip_T = 0, skip_R = 0, skip_snaps = 1, tdiff_clip = tdiff_clip,
#			fig_index = 301, plot_dyn_range = 100, snaptime = snaptime, snaps=snaps, decimate_fac = 1, in_dec = 1,
#			ref_phase = ref_phase, ARRAY = ARRAY, R_slow_plot = R_slow_plot, T_slow_plot = T_slow_plot,
#			max_rat = max_rat, min_amp = min_amp)

#%% Individual events
#%% --Cull seismic section event 1
#pro3singlet(ARRAY = ARRAY, stat_corr = stat_corr, eq_file = eq_file1, simple_taper = simple_taper,
#			rel_time = 0, start_buff = start_buff, end_buff = end_buff,
#			plot_scale_fac = 0.1, skip_SNR = 1,
#			dphase = ref_phase, dphase2 = 'SKKP', dphase3 = 'PKPPcP', dphase4 = 'pPKIKKIKP',
#			freq_min = freq_min, freq_max = freq_max,
#			min_dist = min_dist, max_dist = max_dist, ref_loc = 0, fig_index = 101)

#%% --Cull seismic section event 2
#pro3singlet(ARRAY = ARRAY, stat_corr = stat_corr, eq_file = eq_file2, simple_taper = simple_taper,
#			rel_time = 0, start_buff = start_buff, end_buff = end_buff,
#			plot_scale_fac = 0.1, skip_SNR = 1,
#			dphase = ref_phase, dphase2 = 'PcP', dphase3 = 'PP', dphase4 = 'P',
#			freq_min = freq_min, freq_max = freq_max,
#			min_dist = min_dist, max_dist = max_dist, ref_loc = 0, fig_index = 102)

#%% --1D stack
#pro5stack(ARRAY = ARRAY, eq_file = eq_file1, plot_scale_fac = 0.05,
#			slowR_lo = slowR_lo_1D, slowR_hi = slowR_hi_1D, slow_delta = slow_delta_1D,
#			start_buff = start_buff, end_buff = end_buff,
#			log_plot = 0, envelope = 1, plot_dyn_range = 50,
#			norm = 1, global_norm_plot = 1, color_plot = 1, fig_index = 301)

#pro5stack(ARRAY = ARRAY, eq_file = eq_file2, plot_scale_fac = 0.05,
#			slowR_lo = slowR_lo_1D, slowR_hi = slowR_hi_1D, slow_delta = slow_delta_1D,
#			start_buff = start_buff, end_buff = end_buff,
#			log_plot = 0, envelope = 1, plot_dyn_range = 50,
#			norm = 1, global_norm_plot = 1, color_plot = 1, fig_index = 302)

#%%  --2D stacks
#pro5stack2d(eq_file = eq_file1, plot_scale_fac = 0.05,
#			slowR_lo = slowR_lo, slowR_hi = slowR_hi, slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
#			start_buff = start_buff, end_buff = end_buff,
#			norm = 1, global_norm_plot = 1,
#			ARRAY = ARRAY, decimate_fac = decimate_fac, NS = 0)

#pro5stack2d(eq_file = eq_file2, plot_scale_fac = 0.05,
#			slowR_lo = slowR_lo, slowR_hi = slowR_hi, slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
#			start_buff = start_buff, end_buff = end_buff,
#			norm = 1, global_norm_plot = 1,
#			ARRAY = ARRAY, decimate_fac = decimate_fac, NS = 0)

#%% --2D envelop stack results for individual events
#pro7plotstack(eq_file = eq_file1, plot_scale_fac = 0.05,
#			slowR_lo = slowR_lo, slowR_hi = slowR_hi, slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
#			start_buff = start_buff, end_buff = end_buff, skip_T = 0, skip_R = 0,
#			zoom = 0, ZslowR_lo = -0.03, ZslowR_hi = 0.03, ZslowT_lo = -0.03, ZslowT_hi = 0.03, Zstart_buff = 0, Zend_buff = 200,
#			fig_index = 401, plot_dyn_range = 50, snaptime = snaptime, snaps=0, ARRAY = ARRAY)
#
#pro7plotstack(eq_file = eq_file2, plot_scale_fac = 0.05,
#			slowR_lo = slowR_lo, slowR_hi = slowR_hi, slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
#			start_buff = start_buff, end_buff = end_buff, skip_T = 0, skip_R = 0,
#			zoom = 0, ZslowR_lo = -0.03, ZslowR_hi = 0.03, ZslowT_lo = -0.03, ZslowT_hi = 0.03, Zstart_buff = 0, Zend_buff = 200,
#			fig_index = 402, plot_dyn_range = 50, snaptime = snaptime, snaps=0, ARRAY = ARRAY)
