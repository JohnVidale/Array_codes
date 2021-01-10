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
from pro3a_sort_plot_pair      import pro3pair
from pro3b_sort_plot_singlet   import pro3singlet
from pro5a_stack               import pro5stack
from pro5b_stack2d             import pro5stack2d
from pro6_pair_cc              import pro6_cc_pair      # new experimental code
from pro7a_plot_envstack       import pro7plotstack
from pro7_pair_scan            import pro7_pair_scan

#%% Common parameters
ARRAY      = 1
eq_file1   = 'event1.txt'
eq_file2   = 'event2.txt'

# PcP
#start_buff = 640
#end_buff   = 680

# PKiKP best
# start_buff = 1070
# end_buff   = 1170
start_buff = 1000
end_buff   = 1300

# PKiKP whole
#start_buff = 1000
#end_buff   = 1240

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

# Pro5 stacking
stat_corr = 1
decimate_fac   =    5 # set for pro5stack2d for single event envelopes, set to 0 for other codes
simple_taper   =    1
skip_SNR       =    1
ref_phase      = 'PKiKP'
slowR_lo       = -0.08
slowR_hi       =  0.08
slowT_lo       = -0.08
slowT_hi       =  0.08
slow_delta     =  0.0025
NS = False  # 1 for N-S co=ords, 0 for R-T

# Pro5 1D plot options
slowR_lo_1D   = -0.04
slowR_hi_1D   =  0.1
slow_delta_1D =  0.001

# Pro6 options: to restrict time range
start_beam = 0
end_beam =   0
tdiff_clip   =  0.15
cc_twin      =  4      # time window for cross-correlation (s)
cc_len       =  0.0625 # max time window shift to compute CC (fraction of whole time window)
cc_delta     =  0.4    # time interval for cc (s)
cc_interp1d  =  5      # interpolation factor
cc_thres     =  0.7    # threshold beam correlation to use in stack
min_amp      =  0.0    # threshold amp to use in stack
freq_corr    =  1.2    # only used in old time shift estimation


# Pro 7 options
skip_T = 0       # do not do T plots
skip_R = 1       # do not do R plots
zoom =  True     # to restrict time range in pro7_pair_scan
ZslowR_lo = -0.03
ZslowR_hi =  0.03
ZslowT_lo = -0.03
ZslowT_hi =  0.03
Zstart_buff = 1020
Zend_buff =   1150

# Pro 7 auto_slice == True options
auto_slice      = False  # slices span wide range of R and T slownesses
two_slice_plots = False  # makes R-T pair and snap through time span
beam_sums       = False  # sums tdiff and amp over time
wiggly_plots    = True  # shows wiggly plots

# Pro7 auto-plot options
nR_plots  = 0     # number of plots along the radial axis, makes (2 x nR_plots - 1) total
nT_plots  = 0     # number of plots along the trans  axis
slow_incr = 0.01  # increment at which amp and tdiff are plotted

# Pro7 rest_of_plot options
R_slow_plot    =  0.010
T_slow_plot    =  0.000
skip_snaps     =    1
snaptime       = 1100
snaps          =    2

#%% Comparing events
#%% -- Cull seismic section for common stations
# pro3pair(ARRAY = ARRAY, eq_file1 = eq_file1, eq_file2 = eq_file2, simple_taper = simple_taper, skip_SNR = skip_SNR,
#             rel_time = 0, start_buff = start_buff, end_buff = end_buff,
#             freq_min = freq_min, freq_max = freq_max,
#             plot_scale_fac = 0.025, stat_corr = stat_corr,
#             dphase = ref_phase, dphase2 = 'PKKP', dphase3 = 'PP', dphase4 = 'S',
#             min_dist = min_dist, max_dist = max_dist, ref_loc = 0)

#%%  -- 2D stacks
# pro5stack2d(eq_file = eq_file1,
#             slowR_lo = slowR_lo, slowR_hi = slowR_hi, slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
#             start_buff = start_buff, end_buff = end_buff, norm = 1, ARRAY = ARRAY, decimate_fac = decimate_fac, NS = NS)

# pro5stack2d(eq_file = eq_file2,
#             slowR_lo = slowR_lo, slowR_hi = slowR_hi, slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
#             start_buff = start_buff, end_buff = end_buff, norm = 1, ARRAY = ARRAY, decimate_fac = decimate_fac, NS = NS)

#%% -- Compare pair of 2D stack results to find shift, amp, amp ratio, uses cc rather than instant phase
# pro6_cc_pair(eq_file1 = eq_file1, eq_file2 = eq_file2, plot_scale_fac = 0.002,
#             slowR_lo = slowR_lo, slowR_hi = slowR_hi, slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
#             start_buff = start_buff, end_buff = end_buff, freq_corr = freq_corr, dphase = ref_phase,
#             ARRAY = ARRAY, turn_off_black = 0, R_slow_plot = R_slow_plot, T_slow_plot = T_slow_plot,
#             min_amp = min_amp, tdiff_clip = tdiff_clip, cc_twin = cc_twin, cc_len = cc_len,
#             cc_interp1d = cc_interp1d, cc_delta = cc_delta, cc_thres = cc_thres, NS = NS)

#%% -- Scan of 11 slices of time shift and amplitude
pro7_pair_scan(eq_file1 = eq_file1, eq_file2 = eq_file2, plot_scale_fac = 0.002,
            slowR_lo = slowR_lo, slowR_hi = slowR_hi, slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
            zoom = zoom, ZslowR_lo = ZslowR_lo, ZslowR_hi = ZslowR_hi, ZslowT_lo = ZslowT_lo,
            ZslowT_hi = ZslowT_hi, Zstart_buff = Zstart_buff, Zend_buff = Zend_buff,
            start_buff = start_buff, end_buff = end_buff, skip_T = skip_T, skip_R = skip_R, tdiff_clip = tdiff_clip,
            min_amp = min_amp, ref_phase = ref_phase, cc_thres = cc_thres,
            R_slow_plot = R_slow_plot, T_slow_plot = T_slow_plot, snaptime = snaptime, snaps = snaps, skip_snaps = skip_snaps,
            nR_plots  = nR_plots, nT_plots = nT_plots, slow_incr = slow_incr, NS = NS, dphase = ref_phase,
            ARRAY = ARRAY, auto_slice = auto_slice, two_slice_plots = two_slice_plots, beam_sums = beam_sums,
            wiggly_plots = wiggly_plots)

#%% Individual events
#%% -- Cull seismic section event 1
#pro3singlet(ARRAY = ARRAY, stat_corr = stat_corr, eq_file = eq_file1, simple_taper = simple_taper,
#            rel_time = 0, start_buff = start_buff, end_buff = end_buff,
#            plot_scale_fac = 0.1, skip_SNR = 1,
#            dphase = ref_phase, dphase2 = 'SKKP', dphase3 = 'PKPPcP', dphase4 = 'pPKIKKIKP',
#            freq_min = freq_min, freq_max = freq_max,
#            min_dist = min_dist, max_dist = max_dist, ref_loc = 0, fig_index = 101)

#%% -- Cull seismic section event 2
#pro3singlet(ARRAY = ARRAY, stat_corr = stat_corr, eq_file = eq_file2, simple_taper = simple_taper,
#            rel_time = 0, start_buff = start_buff, end_buff = end_buff,
#            plot_scale_fac = 0.1, skip_SNR = 1,
#            dphase = ref_phase, dphase2 = 'PcP', dphase3 = 'PP', dphase4 = 'P',
#            freq_min = freq_min, freq_max = freq_max,
#            min_dist = min_dist, max_dist = max_dist, ref_loc = 0, fig_index = 102)

#%% -- 1D stack
#pro5stack(ARRAY = ARRAY, eq_file = eq_file1, plot_scale_fac = 0.05,
#            slowR_lo = slowR_lo_1D, slowR_hi = slowR_hi_1D, slow_delta = slow_delta_1D,
#            start_buff = start_buff, end_buff = end_buff,
#            log_plot = 0, envelope = 1, plot_dyn_range = 50,
#            norm = 1, global_norm_plot = 1, color_plot = 1, fig_index = 301)

#pro5stack(ARRAY = ARRAY, eq_file = eq_file2, plot_scale_fac = 0.05,
#            slowR_lo = slowR_lo_1D, slowR_hi = slowR_hi_1D, slow_delta = slow_delta_1D,
#            start_buff = start_buff, end_buff = end_buff,
#            log_plot = 0, envelope = 1, plot_dyn_range = 50,
#            norm = 1, global_norm_plot = 1, color_plot = 1, fig_index = 302)

#%%  -- 2D stacks
#pro5stack2d(eq_file = eq_file1, plot_scale_fac = 0.05,
#            slowR_lo = slowR_lo, slowR_hi = slowR_hi, slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
#            start_buff = start_buff, end_buff = end_buff,
#            norm = 1, global_norm_plot = 1,
#            ARRAY = ARRAY, decimate_fac = decimate_fac, NS = False)

#pro5stack2d(eq_file = eq_file2, plot_scale_fac = 0.05,
#            slowR_lo = slowR_lo, slowR_hi = slowR_hi, slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
#            start_buff = start_buff, end_buff = end_buff,
#            norm = 1, global_norm_plot = 1,
#            ARRAY = ARRAY, decimate_fac = decimate_fac, NS = False)

#%% -- 2D envelop stack results for individual events
#pro7plotstack(eq_file = eq_file1, plot_scale_fac = 0.05,
#            slowR_lo = slowR_lo, slowR_hi = slowR_hi, slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
#            start_buff = start_buff, end_buff = end_buff, skip_T = 0, skip_R = 0,
#            zoom = 0, ZslowR_lo = -0.03, ZslowR_hi = 0.03, ZslowT_lo = -0.03, ZslowT_hi = 0.03, Zstart_buff = 0, Zend_buff = 200,
#            fig_index = 401, plot_dyn_range = 50, snaptime = snaptime, snaps=0, ARRAY = ARRAY)
#
#pro7plotstack(eq_file = eq_file2, plot_scale_fac = 0.05,
#            slowR_lo = slowR_lo, slowR_hi = slowR_hi, slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
#            start_buff = start_buff, end_buff = end_buff, skip_T = 0, skip_R = 0,
#            zoom = 0, ZslowR_lo = -0.03, ZslowR_hi = 0.03, ZslowT_lo = -0.03, ZslowT_hi = 0.03, Zstart_buff = 0, Zend_buff = 200,
#            fig_index = 402, plot_dyn_range = 50, snaptime = snaptime, snaps=0, ARRAY = ARRAY)
