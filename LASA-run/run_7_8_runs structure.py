#!/usr/bin/env python
# For the pair of explosions under Amchitka Island.
# John Vidale 2/2019

import os

os.environ['PATH'] += os.pathsep + '/usr/local/bin'
os.chdir('/Users/vidale/Documents/GitHub/Array_codes')

#%% Import functions
from pro2_con_stfs           import pro2_convstf
#from pro2_plot_stfs          import pro2_test
from pro2_plot_conv          import pro2_test
from pro3a_sort_plot_pair    import pro3pair
from pro3b_sort_plot_singlet import pro3singlet
from pro5a_stack             import pro5stack
from pro5b_stack2d           import pro5stack2d
from pro6_plot_stacked_seis  import pro6stacked_seis
from pro7a_plot_envstack     import pro7plotstack
from pro7b_plot_stack        import pro7plotstack2
from pro7b_dec               import pro7dec
os.chdir('/Users/vidale/Documents/PyCode/LASA')

#%% Common parameters
ARRAY      = 1
eq_file1   = 'event7.txt'
eq_file2   = 'event8.txt'

# PKiKP
start_buff = -1000
end_buff   = 1100 # middle rather than end of ICS
#end_buff   = 1200

conv_file1 = 'Pro_Files/HD1971-11-06_stf.mseed'
conv_file2 = 'Pro_Files/HD1969-10-02_stf.mseed'

freq_min = 0.5
freq_max = 2

# full array
min_dist = 46.2
max_dist = 48.2

# only central array
#min_dist = 46.85
#max_dist = 47.6

slowR_lo   = -0.08
slowR_hi   =  0.08
slowT_lo   = -0.08
slowT_hi   =  0.08
slow_delta =  0.005

Zstart_buff =  0
Zend_buff   = 20
ZslowR_lo = -0.03
ZslowR_hi =  0.03
ZslowT_lo = -0.03
ZslowT_hi =  0.03

min_rat = 0.6
max_rat = 1.8
min_amp = 0.2

slowR_lo_1D = -0.04
slowR_hi_1D = 0.1
slow_delta_1D = 0.001

decimate_fac   =  5
simple_taper   =  1
skip_SNR       =  1
skip_snaps = 0
snaptime = 28
snaps = 8
freq_corr = 1.2
tdiff_clip = 0.3
ref_phase = 'PKiKP'
dphase = 'PKiKP' # phase for start_buff and end_buff for initial trace selection

#%% Comparing events
#%% --Cull seismic section for common stations
pro3pair(ARRAY = ARRAY, eq_file1 = eq_file1, eq_file2 = eq_file2, simple_taper = simple_taper, skip_SNR = skip_SNR,
			rel_time = 0, start_buff = start_buff, end_buff = end_buff,
			freq_min = freq_min, freq_max = freq_max,
			plot_scale_fac = 0.1, stat_corr = 1,
			dphase = dphase, dphase2 = 'PKKP', dphase3 = 'PP', dphase4 = 'S',
			min_dist = min_dist, max_dist = max_dist, ref_loc = 0)

#%%  --Cross_convolve time functions
#pro2_convstf(eq_file = eq_file1, conv_file = conv_file1)
#pro2_convstf(eq_file = eq_file2, conv_file = conv_file2)
#pro2_test(eq_file1 = eq_file1, conv_file1 = conv_file1, eq_file2 = eq_file2, conv_file2 = conv_file2)

#%%  --2D stacks
pro5stack2d(eq_file = eq_file1, plot_scale_fac = 0.05,
			slowR_lo = slowR_lo, slowR_hi = slowR_hi, slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
			start_buff = start_buff, end_buff = end_buff,
			norm = 1, global_norm_plot = 1,
			ARRAY = ARRAY, decimate_fac = decimate_fac, NS = 0)

pro5stack2d(eq_file = eq_file2, plot_scale_fac = 0.05,
			slowR_lo = slowR_lo, slowR_hi = slowR_hi, slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
			start_buff = start_buff, end_buff = end_buff,
			norm = 1, global_norm_plot = 1,
			ARRAY = ARRAY, decimate_fac = decimate_fac, NS = 0)

#%% --Compare pair of 2D stack results
pro6stacked_seis(eq_file1 = eq_file1, eq_file2 = eq_file2, plot_scale_fac = 0.003,
			slowR_lo = slowR_lo, slowR_hi = slowR_hi, slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
			start_buff = start_buff, end_buff = end_buff, freq_corr = freq_corr, ref_phase = ref_phase,
			fig_index = 301, plot_dyn_range = 100, ARRAY = ARRAY, min_rat = min_rat, max_rat = max_rat, min_amp = min_amp)

#%% --decimate stacking files to shorten processing time
#pro7dec(eq_file1 = eq_file1, eq_file2 = eq_file2, decimate_fac = decimate_fac, ARRAY = ARRAY)

#%% --slowness slices of time shift
#pro7plotstack2(eq_file1 = eq_file1, eq_file2 = eq_file2, plot_scale_fac = 0.05,
#			slowR_lo = slowR_lo, slowR_hi = slowR_hi, slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
#			zoom = 0, ZslowR_lo = ZslowR_lo, ZslowR_hi = ZslowR_hi, ZslowT_lo = ZslowT_lo, ZslowT_hi = ZslowT_hi, Zstart_buff = Zstart_buff, Zend_buff = Zend_buff,
#			start_buff = start_buff, end_buff = end_buff, skip_T = 1, skip_R = 0, tdiff_clip = tdiff_clip,
#			fig_index = 301, plot_dyn_range = 100, skip_snaps = skip_snaps, snaptime = snaptime, snaps=snaps,
#			decimate_fac = 1, in_dec = 1, ref_phase = ref_phase, ARRAY = ARRAY, min_rat = min_rat, max_rat = max_rat, min_amp = min_amp)

#%% Individual events
#%% --Cull seismic section event 1
#pro3singlet(ARRAY = ARRAY, stat_corr = 1, eq_file = eq_file1, simple_taper = simple_taper,
#			start_buff = start_buff, end_buff = end_buff,
#			plot_scale_fac = 0.1, skip_SNR = 1,
#			dphase = 'PKiKP', dphase2 = 'PcP', dphase3 = 'PP', dphase4 = 'P',
#			freq_min = freq_min, freq_max = freq_max,
#			min_dist = 46.2, max_dist = 48.2, ref_loc = 0)

#%% --Cull seismic section event 2
#pro3singlet(ARRAY = ARRAY, stat_corr = 1, eq_file = eq_file2,
#			start_buff = start_buff, end_buff = end_buff,
#			plot_scale_fac = 0.1,
#			dphase = 'PKiKP', dphase2 = 'PcP', dphase3 = 'PP', dphase4 = 'P',
#			freq_min = freq_min, freq_max = freq_max,
#			min_dist = 46.2, max_dist = 48.2, ref_loc = 0)

#%% --1D stack
#pro5stack(ARRAY = ARRAY, eq_file = eq_file1, plot_scale_fac = 0.05,
#			slowR_lo = slowR_lo_1D, slowR_hi = slowR_hi_1D, slow_delta = slow_delta_1D,
#			start_buff = start_buff, end_buff = end_buff,
#			envelope = 1, plot_dyn_range = 50,
#			norm = 1, global_norm_plot = 1, color_plot = 1, fig_index = 301)

#pro5stack(ARRAY = ARRAY, eq_file = eq_file2, plot_scale_fac = 0.05,
#			slowR_lo = slowR_lo_1D, slowR_hi = slowR_hi_1D, slow_delta = slow_delta_1D,
#			start_buff = start_buff, end_buff = end_buff,
#			envelope = 1, plot_dyn_range = 50,
#			norm = 1, global_norm_plot = 1, color_plot = 1, fig_index = 302)

#%% --2D envelop stack results for individual events
#pro7plotstack(eq_file = eq_file1, plot_scale_fac = 0.05,
#			slowR_lo = slowR_lo, slowR_hi = slowR_hi, slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
#			start_buff = start_buff, end_buff = end_buff, skip_T = 1, skip_R = 1,
#			zoom = 0, ZslowR_lo = -0.03, ZslowR_hi = 0.03, ZslowT_lo = -0.03, ZslowT_hi = 0.03, Zstart_buff = 0, Zend_buff = 200,
#			fig_index = 301, plot_dyn_range = 50, snaptime = snaptime, snaps=1, ARRAY = ARRAY)
#
#pro7plotstack(eq_file = eq_file2, plot_scale_fac = 0.05,
#			slowR_lo = slowR_lo, slowR_hi = slowR_hi, slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
#			start_buff = start_buff, end_buff = end_buff,
#			fig_index = 401, plot_dyn_range = 100, snaptime = snaptime, snaps=1, ARRAY = ARRAY)