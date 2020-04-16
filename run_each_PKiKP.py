#!/usr/bin/env python3
# John Vidale 4/2020

def run_each(start_buff = 980, end_buff = 1180, event_no = 35, min_dist = 0,
			  max_dist = 180, freq_min = 1, freq_max = 3, slow_delta = 0.0025):

	import os
	os.environ['PATH'] += os.pathsep + '/usr/local/bin'
	os.chdir('/Users/vidale/Documents/GitHub/Array_codes')

	#%% Import functions
	from pro3b_sort_plot_singlet import pro3singlet
	from pro5a_stack             import pro5stack
	from pro5b_stack2d           import pro5stack2d
	from pro6_plot_stacked_seis  import pro6stacked_seis
	#from junk     import pro7plotstack
	from pro7a_plot_envstack     import pro7plotstack
	from pro7b_plot_stack        import pro7plotstack2
	from pro7b_dec               import pro7dec
	import matplotlib.pyplot as plt
	os.chdir('/Users/vidale/Documents/PyCode/LASA')

	#%% Common parameters
	ARRAY      = 1
	eq_file   = 'event' + str(event_no) +'.txt'
	#eq_file   = 'event35.txt'

	# window
	#start_buff = 980
	#end_buff   = 1180
	NS = 1
	slowR_lo   = -0.04
	slowR_hi   =  0.04
	slowT_lo   = -0.04
	slowT_hi   =  0.04
#	slow_delta =  0.0025

	dphase = 'PKiKP'
	dphase2 = 'PP'
	dphase3 = 'PKIKKIKP'
	dphase4 = 'PKIKPPKIKP'

	# Full array
	min_dist = min_dist
	max_dist = max_dist
	auto_dist = 1

	ref_loc = 0 # all stations
	#ref_loc = 1 # only rings A-D

	freq_min = freq_min
	freq_max = freq_max

	slowR_lo_1D = -0.04
	slowR_hi_1D = 0.1
	slow_delta_1D = 0.005

	decimate_fac   = 10
	simple_taper   =  1
	skip_SNR       =  1
	snaptime = 995.5
	snaps = 0
	freq_corr = 1.2
	stat_corr = 1

	#%% --Cull seismic section
	# plot lines are blue, orange, yellow, purple for phases 1 through 4
	pro3singlet(ARRAY = ARRAY, stat_corr = stat_corr, eq_file = eq_file, simple_taper = simple_taper,
				rel_time = 1, rel_slow = 1, start_buff = start_buff, end_buff = end_buff,
				plot_scale_fac = 0.1, skip_SNR = 1,
				dphase = dphase, dphase2 = dphase2, dphase3 = dphase3, dphase4 = dphase4,
				freq_min = freq_min, freq_max = freq_max,
				min_dist = min_dist, max_dist = max_dist, auto_dist = auto_dist,
				ref_loc = ref_loc, fig_index = 102)

	#%% --1D stack
	pro5stack(ARRAY = ARRAY, eq_file = eq_file, plot_scale_fac = 0.05,
				slowR_lo = slowR_lo_1D, slowR_hi = slowR_hi_1D, slow_delta = slow_delta_1D,
				start_buff = start_buff, end_buff = end_buff,
				log_plot = 0, envelope = 1, plot_dyn_range = 50,
				norm = 1, global_norm_plot = 1, color_plot = 1, fig_index = 302)

	##%%  --2D stack
	pro5stack2d(eq_file = eq_file, plot_scale_fac = 0.05,
				slowR_lo = slowR_lo, slowR_hi = slowR_hi, slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
				start_buff = start_buff, end_buff = end_buff,
				norm = 1, global_norm_plot = 1,
				ARRAY = ARRAY, decimate_fac = decimate_fac, NS = NS)

	#%% --Compare 2D stack results with themselves
	pro6stacked_seis(eq_file1 = eq_file, eq_file2 = eq_file, plot_scale_fac = 0.003,
				slowR_lo = slowR_lo, slowR_hi = slowR_hi, slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
				start_buff = start_buff, end_buff = end_buff, freq_corr = freq_corr, ref_phase = dphase,
				fig_index = 301, plot_dyn_range = 100, ARRAY = ARRAY, event_no = event_no)

	#%% --2D envelope stack
	pro7plotstack(eq_file = eq_file, plot_scale_fac = 0.05,
				slowR_lo = slowR_lo, slowR_hi = slowR_hi, slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
				start_buff = start_buff, end_buff = end_buff, skip_T = 0, skip_R = 0,
				zoom = 0, ZslowR_lo = -0.03, ZslowR_hi = 0.03, ZslowT_lo = -0.03, ZslowT_hi = 0.03, Zstart_buff = 0, Zend_buff = 200,
				fig_index = 402, plot_dyn_range = 50, snaptime = snaptime, snaps=snaps, ARRAY = ARRAY)

#	plt.close('all')

	code_directory = '/Users/vidale/Documents/GitHub/Array_codes'
	os.chdir(code_directory)