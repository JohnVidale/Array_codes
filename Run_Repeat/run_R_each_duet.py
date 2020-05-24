#!/usr/bin/env python3
# John Vidale 4/2020

def runR_each_duet(start_buff = 980, end_buff = 1180, min_dist = 0,
			  max_dist = 180, slow_delta = 0.0025,
			  start_beam = 0, end_beam = 0, dphase  = 'PKiKP',
			  event1_no = 0, event2_no = 0):

	import os
	os.environ['PATH'] += os.pathsep + '/usr/local/bin'
	os.chdir('/Users/vidale/Documents/GitHub/Array_codes')

	#%% Import functions
	from pro3b_sort_plot_singlet import pro3singlet
	from pro3a_sort_plot_pair    import pro3pair
	from pro5a_stack             import pro5stack
	from pro5b_stack2d           import pro5stack2d
	from pro6_plot_singlet       import pro6stacked_singlet
	from pro6_plot_pair          import pro6stacked_pair
	from pro7a_plot_envstack     import pro7plotstack
	from pro7b_plot_stack        import pro7plotstack2
	from pro7b_dec               import pro7dec
	import matplotlib.pyplot as plt

	#%% Common parameters
	ev_directory = '/Users/vidale/Documents/PyCode/EvLocs'
	os.chdir(ev_directory)

	eq1_file   = 'event' + str(event1_no) + '.txt'
	eq2_file   = 'event' + str(event2_no) + '.txt'

	ARRAY     = 1  # 0 Japan, 1 LASA, 2 China

	freq_min = 1
	freq_max = 3

	rel_time = 1   # phase alignment details
#	rel_time == 0  window in absolute time after origin time
#	rel_time == 1  each window has a shift proportional to (dist - ref_dist) at phase slowness at ref_dist
#	rel_time == 2  each window has a distinct phase-chose shift, but time offset is common to all stations
#	rel_time == 3  each station has an individual, chosen-phase shift, phase arrival set to common time
#	rel_time == 4  use same window around chosen phase for all stations, using ref distance

	ref_loc  = 0   # 0 select stations by distance from epicenter, 1 select stations by distance from ref location
	ref_rad  = 0.4 # radius of stations around ref_loc chosen
#	ref_lat  =  36 # Japan
#	ref_lon  = 138
	ref_lat = 46.7      # LASA °N keep only inner rings A-D if radius is 0.4°
	ref_lon = -106.22   # °E
	NS       = 1   # 0 plot slowness R-T, 1 plot slowness N-S

	auto_dist = 1  #  automatically plot only real distance range
	min_dist = 0
	max_dist = 180

	slowR_lo   = -0.03
	slowR_hi   =  0.03
	slowT_lo   = -0.03
	slowT_hi   =  0.03
	slow_delta =  0.001

	slowR_lo_1D = -0.05
	slowR_hi_1D =  0.1
	slow_delta_1D = 0.001

	decimate_fac   =  10
	simple_taper   =  1
	snaptime  = -0.5
	freq_corr = 1.2

	stat_corr      = 1
	corr_threshold = 0.6
	skip_SNR       = 1
	qual_threshold = 1.5

#	dphase  = 'PcP'
	dphase2 = 'PP'
	dphase3 = 'PcP'
	dphase4 = 'sP'

	#%% Pairs
	# get data from Japan
	#os.chdir('/Users/vidale/Documents/PyCode/Hinet/Reps/Aleutians/Pair1')
	#pro1get(eq_file1)
	#pro1get(eq_file2)

	#pro2_decimate(eq_file1, decimate_fac = decimate_fac)
	#pro2_decimate(eq_file2, decimate_fac = decimate_fac)

	pro3pair(ARRAY = ARRAY, eq_file1 = eq1_file, eq_file2 = eq2_file, start_buff = start_buff, skip_SNR = skip_SNR,
				end_buff = end_buff, plot_scale_fac = 0.01, qual_threshold = qual_threshold, simple_taper = 1,
				dphase = dphase, dphase2 = dphase2, dphase3 = dphase3, dphase4 = dphase4,
				corr_threshold = corr_threshold, freq_min = freq_min, freq_max = freq_max,
				min_dist = min_dist, max_dist = max_dist, ref_loc = ref_loc, ref_rad = ref_rad, ref_lat = ref_lat, ref_lon = ref_lon)

	pro5stack(ARRAY = ARRAY, eq_file = eq1_file, plot_scale_fac = 0.05, event_no = event1_no,
				slowR_lo = slowR_lo_1D, slowR_hi = slowR_hi_1D, slow_delta = slow_delta_1D,
				start_buff = start_buff, end_buff = end_buff,
				ref_loc = ref_loc, ref_lat = ref_lat, ref_lon = ref_lon,
				log_plot = 0, envelope = 1, plot_dyn_range = 50,
				norm = 1, global_norm_plot = 1, color_plot = 1, fig_index = 301)

	pro5stack(ARRAY = ARRAY, eq_file = eq2_file, plot_scale_fac = 0.05, event_no = event2_no,
				slowR_lo = slowR_lo_1D, slowR_hi = slowR_hi_1D, slow_delta = slow_delta_1D,
				start_buff = start_buff, end_buff = end_buff,
				ref_loc = ref_loc, ref_lat = ref_lat, ref_lon = ref_lon,
				log_plot = 0, envelope = 1, plot_dyn_range = 50,
				norm = 1, global_norm_plot = 1, color_plot = 1, fig_index = 302)

	#%%  --2D stacks
#	pro5stack2d(eq_file = eq1_file, plot_scale_fac = 0.05,
#				slowR_lo = slowR_lo, slowR_hi = slowR_hi, slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
#				start_buff = start_buff, end_buff = end_buff,
#				norm = 1, ref_loc = ref_loc, ref_lat = ref_lat, ref_lon = ref_lon,
#				ARRAY = ARRAY, decimate_fac = decimate_fac, NS = NS)
#	pro5stack2d(eq_file = eq2_file, plot_scale_fac = 0.05,
#				slowR_lo = slowR_lo, slowR_hi = slowR_hi, slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
#				start_buff = start_buff, end_buff = end_buff,
#				norm = 1, ref_loc = ref_loc, ref_lat = ref_lat, ref_lon = ref_lon,
#				ARRAY = ARRAY, decimate_fac = decimate_fac, NS = NS)

	#%% --Compare pair of 2D stack results
#	pro6stacked_pair(eq_file1 = eq1_file, eq_file2 = eq2_file, plot_scale_fac = 0.003, tdiff_clip = 0.1,
#				slowR_lo = slowR_lo, slowR_hi = slowR_hi, slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
#				start_buff = start_buff, end_buff = end_buff, R_slow_plot = 0, T_slow_plot = 0, dphase = dphase, freq_corr = freq_corr,
#				fig_index = 301, plot_dyn_range = 100, ARRAY = ARRAY, event_no = event1_no, start_beam = start_beam, end_beam = end_beam)

	#%% --Examine each one of a pair of 2D stack results
#	pro6stacked_singlet(eq_file = eq1_file, plot_scale_fac = 0.003,
#				slowR_lo = slowR_lo, slowR_hi = slowR_hi, slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
#				start_buff = start_buff, end_buff = end_buff, R_slow_plot = 0, T_slow_plot = 0, dphase = dphase,
#				fig_index = 301, plot_dyn_range = 100, ARRAY = ARRAY, event_no = event1_no, start_beam = start_beam, end_beam = end_beam)
#
#	pro6stacked_singlet(eq_file = eq2_file, plot_scale_fac = 0.003,
#				slowR_lo = slowR_lo, slowR_hi = slowR_hi, slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
#				start_buff = start_buff, end_buff = end_buff, R_slow_plot = -0.015, T_slow_plot = 0, dphase = dphase,
#				fig_index = 301, plot_dyn_range = 100, ARRAY = ARRAY, event_no = event2_no, start_beam = start_beam, end_beam = end_beam)

	#%% --decimate stacking files to shorten processing time
	#pro7dec(eq_file1 = eq_file1, eq_file2 = eq_file2, decimate_fac = decimate_fac, ARRAY = ARRAY)

	#%% --slowness slices of time shift
	#pro7plotstack2(eq_file1 = eq_file1, eq_file2 = eq_file2, plot_scale_fac = 0.05,
	#			slowR_lo = slowR_lo, slowR_hi = slowR_hi, slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
	#			zoom = 0, ZslowR_lo = -0.03, ZslowR_hi = 0.03, ZslowT_lo = -0.03, ZslowT_hi = 0.03, Zstart_buff = 5, Zend_buff = 15,
	#			start_buff = start_buff, end_buff = end_buff, skip_T = 0, skip_R = 1, skip_snaps = 0, tdiff_clip = 0.2,
	#			fig_index = 301, plot_dyn_range = 100, snaptime = snaptime, snaps=1, decimate_fac = 1, in_dec = 1,
	#			ref_phase = ref_phase, ARRAY = ARRAY)

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

	code_directory = '/Users/vidale/Documents/GitHub/Array_codes'
	os.chdir(code_directory)
