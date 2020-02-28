#!/usr/bin/env python3
# this program runs a suite of programs on Hi-net data
# John Vidale 8/2019

import os
os.environ['PATH'] += os.pathsep + '/usr/local/bin'
os.chdir('/Users/vidale/Documents/GitHub/Array_codes')

from pro1_get                import pro1get
from pro2_dec                import pro2_decimate
from pro3a_sort_plot_pair    import pro3pair
from pro3b_sort_plot_singlet import pro3singlet
from pro5a_stack             import pro5stack
from pro5b_stack2d           import pro5stack2d
from pro6_plot_stacked_seis  import pro6stacked_seis
from pro7a_plot_envstack     import pro7plotstack
from pro7b_plot_stack        import pro7plotstack2
from pro7b_dec               import pro7dec

#%% Common parameters
ev_directory = '/Users/vidale/Documents/PyCode/Hinet/Reps/Aleutians/Pair1'
os.chdir(ev_directory)

ARRAY      = 0
eq_file   = 'event.txt'
eq_file1   = 'event1.txt'
eq_file2   = 'event2.txt'

start_buff = -10
end_buff   = 30
freq_min = 0.2
freq_max = 1
corr_threshold = 0
qual_threshold = 1.5

#min_dist = 26
#max_dist = 29

min_dist = 25
max_dist = 47

slowR_lo   = -0.04
slowR_hi   =  0.08
slowT_lo   = -0.04
slowT_hi   =  0.04
slow_delta =  0.005

slowR_lo_1D = -0.04
slowR_hi_1D = 0.1
slow_delta_1D = 0.001

decimate_fac   =  5
simple_taper   =  1
skip_SNR       =  1
snaptime = 30
freq_corr = 1.2
ref_phase = 'P'

#%% Singlets

# get data from Japan
#pro1get(eq_file)

# decimate, in 100 sps, out 20 sps
#pro2_decimate(eq_file, decimate_fac = decimate_fac)

#pro3singlet(eq_file = eq_file, start_buff = start_buff,
#			end_buff = end_buff, plot_scale_fac = 0.1, qual_threshold = qual_threshold,
#			dphase = ref_phase, dphase2 = 'pP', dphase3 = 'PcP', dphase4 = 'sP',
#			corr_threshold = corr_threshold, freq_min = freq_min, freq_max = freq_max,
#			min_dist = min_dist, max_dist = max_dist)

#pro5stack(eq_file = eq_file, plot_scale_fac = 0.05, slowR_lo = slowR_lo_1D, slowR_hi = slowR_hi_1D,
#			  slow_delta = slow_delta_1D, start_buff = start_buff, end_buff = end_buff, envelope = 1,
#			  norm = 1, global_norm_plot = 1, color_plot = 1, fig_index = 301)

#%% Pairs

#pro1get(eq_file1)
#pro1get(eq_file2)

#pro2_decimate(eq_file1, decimate_fac = decimate_fac)
#pro2_decimate(eq_file2, decimate_fac = decimate_fac)

#pro3singlet(eq_file = eq_file1, start_buff = start_buff,
#			end_buff = end_buff, plot_scale_fac = 0.1, qual_threshold = qual_threshold,
#			dphase = ref_phase, dphase2 = 'pP', dphase3 = 'PcP', dphase4 = 'sP',
#			corr_threshold = corr_threshold, freq_min = freq_min, freq_max = freq_max,
#			min_dist = min_dist, max_dist = max_dist)
#pro3singlet(eq_file = eq_file2, start_buff = start_buff,
#			end_buff = end_buff, plot_scale_fac = 0.1, qual_threshold = qual_threshold,
#			dphase = ref_phase, dphase2 = 'pP', dphase3 = 'PcP', dphase4 = 'sP',
#			corr_threshold = corr_threshold, freq_min = freq_min, freq_max = freq_max,
#			min_dist = min_dist, max_dist = max_dist)

pro3pair(eq_file1 = eq_file2, eq_file2 = eq_file1, start_buff = start_buff,
			end_buff = end_buff, plot_scale_fac = 0.2, qual_threshold = qual_threshold,
			dphase = ref_phase, dphase2 = 'pP', dphase3 = 'PcP', dphase4 = 'sP',
			corr_threshold = corr_threshold, freq_min = freq_min, freq_max = freq_max,
			min_dist = min_dist, max_dist = max_dist)

pro5stack(eq_file = eq_file1, plot_scale_fac = 0.05, slowR_lo = slowR_lo_1D, slowR_hi = slowR_hi_1D,
			  slow_delta = slow_delta_1D, start_buff = start_buff, end_buff = end_buff, envelope = 1,
			  norm = 1, global_norm_plot = 1, color_plot = 1, fig_index = 301)
pro5stack(eq_file = eq_file2, plot_scale_fac = 0.05, slowR_lo = slowR_lo_1D, slowR_hi = slowR_hi_1D,
			  slow_delta = slow_delta_1D, start_buff = start_buff, end_buff = end_buff, envelope = 1,
			  norm = 1, global_norm_plot = 1, color_plot = 1, fig_index = 301)

#%%  --2D stacks
#pro5stack2d(eq_file = eq_file1, plot_scale_fac = 0.05,
#			slowR_lo = slowR_lo, slowR_hi = slowR_hi, slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
#			start_buff = start_buff, end_buff = end_buff,
#			norm = 1, global_norm_plot = 1,
#			ARRAY = ARRAY, decimate_fac = decimate_fac, NS = 0)
#
#pro5stack2d(eq_file = eq_file2, plot_scale_fac = 0.05,
#			slowR_lo = slowR_lo, slowR_hi = slowR_hi, slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
#			start_buff = start_buff, end_buff = end_buff,
#			norm = 1, global_norm_plot = 1,
#			ARRAY = ARRAY, decimate_fac = decimate_fac, NS = 0)

#%% --Compare pair of 2D stack results
#pro6stacked_seis(eq_file1 = eq_file1, eq_file2 = eq_file2, plot_scale_fac = 0.003,
#			slowR_lo = slowR_lo, slowR_hi = slowR_hi, slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
#			start_buff = start_buff, end_buff = end_buff, freq_corr = freq_corr, ref_phase = ref_phase,
#			fig_index = 301, plot_dyn_range = 100, ARRAY = ARRAY)

#%% just copied from LASA workflow

#%% --decimate stacking files to shorten processing time
#pro7dec(eq_file1 = eq_file1, eq_file2 = eq_file2, decimate_fac = decimate_fac, ARRAY = ARRAY)

#%% --slowness slices of time shift
#pro7plotstack2(eq_file1 = eq_file1, eq_file2 = eq_file2, plot_scale_fac = 0.05,
#			slowR_lo = slowR_lo, slowR_hi = slowR_hi, slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
#			zoom = 0, ZslowR_lo = -0.03, ZslowR_hi = 0.03, ZslowT_lo = -0.03, ZslowT_hi = 0.03, Zstart_buff = 5, Zend_buff = 15,
#			start_buff = start_buff, end_buff = end_buff, skip_T = 0, skip_R = 1, skip_snaps = 0, tdiff_clip = 0.2,
#			fig_index = 301, plot_dyn_range = 100, snaptime = snaptime, snaps=1, decimate_fac = 1, in_dec = 1,
#			ref_phase = ref_phase, ARRAY = ARRAY)

#%% Individual events
#%% --Cull seismic section event 1
#pro3singlet(ARRAY = ARRAY, stat_corr = 1, eq_file = eq_file1, simple_taper = simple_taper,
#			start_buff = start_buff, end_buff = end_buff,
#			plot_scale_fac = 0.1, skip_SNR = 1,
#			dphase = 'PKiKP', dphase2 = 'PcP', dphase3 = 'PP', dphase4 = 'P',
#			freq_min = freq_min, freq_max = freq_max,
#			min_dist = 58.5, max_dist = 60.5, ref_loc = 0)

#%% --Cull seismic section event 2
#pro3singlet(ARRAY = ARRAY, stat_corr = 1, eq_file = eq_file2, simple_taper = simple_taper,
#			start_buff = start_buff, end_buff = end_buff,
#			plot_scale_fac = 0.1, skip_SNR = 1,
#			dphase = 'PKiKP', dphase2 = 'PcP', dphase3 = 'PP', dphase4 = 'P',
#			freq_min = freq_min, freq_max = freq_max,
#			min_dist = 58.5, max_dist = 60.5, ref_loc = 0)

#%% --1D stack
#pro5stack(ARRAY = ARRAY, eq_file = eq_file1, plot_scale_fac = 0.05,
#			slowR_lo = slowR_lo_1D, slowR_hi = slowR_hi_1D, slow_delta = slow_delta_1D,
#			start_buff = start_buff, end_buff = end_buff,
#			log_plot = 0, envelope = 1, plot_dyn_range = 50,
#			norm = 1, global_norm_plot = 1, color_plot = 1, fig_index = 301)
#
#pro5stack(ARRAY = ARRAY, eq_file = eq_file2, plot_scale_fac = 0.05,
#			slowR_lo = slowR_lo_1D, slowR_hi = slowR_hi_1D, slow_delta = slow_delta_1D,
#			start_buff = start_buff, end_buff = end_buff,
#			log_plot = 0, envelope = 1, plot_dyn_range = 50,
#			norm = 1, global_norm_plot = 1, color_plot = 1, fig_index = 302)

#%% --2D envelop stack results for individual events
#pro7plotstack(eq_file = eq_file1, plot_scale_fac = 0.05,
#			slowR_lo = slowR_lo, slowR_hi = slowR_hi, slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
#			start_buff = start_buff, end_buff = end_buff, skip_T = 1, skip_R = 1,
#			zoom = 0, ZslowR_lo = -0.03, ZslowR_hi = 0.03, ZslowT_lo = -0.03, ZslowT_hi = 0.03, Zstart_buff = 0, Zend_buff = 200,
#			fig_index = 401, plot_dyn_range = 50, snaptime = snaptime, snaps=1, ARRAY = ARRAY)

#pro7plotstack(eq_file = eq_file2, plot_scale_fac = 0.05,
#			slowR_lo = slowR_lo, slowR_hi = slowR_hi, slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
#			start_buff = start_buff, end_buff = end_buff, skip_T = 1, skip_R = 0,
#			zoom = 0, ZslowR_lo = -0.03, ZslowR_hi = 0.03, ZslowT_lo = -0.03, ZslowT_hi = 0.03, Zstart_buff = 0, Zend_buff = 200,
#			fig_index = 402, plot_dyn_range = 50, snaptime = snaptime, snaps=1, ARRAY = ARRAY)
