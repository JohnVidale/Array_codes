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
from pro3b_sort_plot_singlet import pro3singlet
from pro5a_stack             import pro5stack
from pro5b_stack2d           import pro5stack2d
from pro6_plot_stacked_seis  import pro6stacked_seis
from pro7a_plot_envstack     import pro7plotstack
os.chdir('/Users/vidale/Documents/PyCode/LASA')

#%% Common parameters
ARRAY      = 1
eq_file   = 'event10.txt'

start_buff = -1200
end_buff   = 1220

# fine
#freq_min = 1
#freq_max = 3

# low-freq
freq_min = 0.5
freq_max = 2

# all stations, with outer rings
min_dist = 163.9
max_dist = 166.0

# omit outer rings
#min_dist = 164.4
#max_dist = 165.0

slowR_lo   = -0.04
slowR_hi   =  0.04
slowT_lo   = -0.04
slowT_hi   =  0.04
slow_delta =  0.0025

slowR_lo_1D = -0.1
slowR_hi_1D = 0.1
slow_delta_1D = 0.0025

decimate_fac   =  5
simple_taper   =  1
skip_SNR       =  1
snaptime = 10 - start_buff
freq_corr = 1.2
ref_phase = 'PKIKP'

#%% Individual event
#%% --Cull seismic section event 1
pro3singlet(ARRAY = ARRAY, stat_corr = 1, eq_file = eq_file, simple_taper = simple_taper,
			rel_time = 0, start_buff = start_buff, end_buff = end_buff,
			plot_scale_fac = 0.1, skip_SNR = 1,
			dphase = ref_phase, dphase2 = 'PP', dphase3 = 'PKP', dphase4 = 'PPP',
			freq_min = freq_min, freq_max = freq_max,
			min_dist = min_dist, max_dist = max_dist, ref_loc = 0, fig_index = 101)

#%% --1D stack
#pro5stack(ARRAY = ARRAY, eq_file = eq_file, plot_scale_fac = 0.05,
#			slowR_lo = slowR_lo_1D, slowR_hi = slowR_hi_1D, slow_delta = slow_delta_1D,
#			start_buff = start_buff, end_buff = end_buff,
#			log_plot = 1, envelope = 1, plot_dyn_range = 200,
#			norm = 1, global_norm_plot = 1, color_plot = 1, fig_index = 301)

#%%  --2D stacks
pro5stack2d(eq_file = eq_file, plot_scale_fac = 0.05,
			slowR_lo = slowR_lo, slowR_hi = slowR_hi, slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
			start_buff = start_buff, end_buff = end_buff,
			norm = 1, global_norm_plot = 1,
			ARRAY = ARRAY, decimate_fac = decimate_fac, NS = 0)

##%% --Compare pair of 2D stack results
pro6stacked_seis(eq_file1 = eq_file, eq_file2 = eq_file, plot_scale_fac = 0.003,
			slowR_lo = slowR_lo, slowR_hi = slowR_hi, slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
			start_buff = start_buff, end_buff = end_buff, freq_corr = freq_corr, ref_phase = ref_phase,
			fig_index = 301, plot_dyn_range = 100, ARRAY = ARRAY)

#%% --2D envelop stack results for individual events
#pro7plotstack(eq_file = eq_file, plot_scale_fac = 0.05,
#			slowR_lo = slowR_lo, slowR_hi = slowR_hi, slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
#			start_buff = start_buff, end_buff = end_buff, skip_T = 1, skip_R = 1,
#			zoom = 0, ZslowR_lo = -0.03, ZslowR_hi = 0.03, ZslowT_lo = -0.03, ZslowT_hi = 0.03, Zstart_buff = 0 - start_buff, Zend_buff = 200 - start_buff,
#			fig_index = 401, plot_dyn_range = 50, snaptime = snaptime, snaps=1, ARRAY = ARRAY)