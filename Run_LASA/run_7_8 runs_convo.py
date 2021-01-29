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
from pro2_con_stfs             import pro2_convstf
from pro2_plot_conv            import pro2_test
from pro3a_sort_plot_pair      import pro3pair
from pro3b_sort_plot_singlet   import pro3singlet
from pro5a_stack               import pro5stack
from pro5b_stack2d             import pro5stack2d
from pro6_pair_cc              import pro6_cc_pair
from pro6_singlet              import pro6_singlet
from pro7a_plot_envstack       import pro7plotstack
from pro7_pair_scan            import pro7_pair_scan
from pro7_singlet              import pro7_singlet

#%% Common parameters

do_3  = True
do_2a = True
do_2b = True
do_5  = True

do_6 = False # treats pairs of events with time diffs
do_7 = False

do_6a = True # treats single events, no time shifts calculated or plotted
do_7a = True

ARRAY      = 1
eq_file1 = 'event7.txt'  # pair
eq_file2 = 'event8.txt'
eq_file  = 'event7.txt'  # singlet

# P
# start_buff =  600
# end_buff   = 2500
start_buff = 950
end_buff   = 1300

# Full array
min_dist = 46.2
max_dist = 48.2
conv_file1 = '/Users/vidale/Documents/GitHub/Array_codes/Files/HD1971-11-06_stf.mseed'
conv_file2 = '/Users/vidale/Documents/GitHub/Array_codes/Files/HD1969-10-02_stf.mseed'

# HF
freq_min = 1
freq_max = 3

# Pro5 stacking
stat_corr = 1
decimate_fac   =    5 # set for pro5stack2d for single event envelopes, set to 0 for other codes
simple_taper   =    1
max_taper_length = 5. # taper is minimum of taper_frac (0.05) and this number of seconds
skip_SNR       =    1
ref_phase      = 'PKiKP'
slowR_lo       = -0.03
slowR_hi       =  0.03
slowT_lo       = -0.03
slowT_hi       =  0.03
slow_delta     =  0.0025
NS = False  # 1 for N-S co=ords, 0 for R-T

# Pro5 1D plot options
slowR_lo_1D   = -0.04
slowR_hi_1D   =  0.1
slow_delta_1D =  0.001

# Pro6 options: to restrict time range
tdiff_clip   =  0.15
no_tdiff_plot = False  # also to speed plots of only amplitude
log_plot      = True
plot_scale_fac = 1
wig_scale_fac = 0.5
tdiff_scale_fac = 1
log_plot_range = 1.5

cc_twin      =  10     # time window for cross-correlation (s)
cc_len       =  0.05 # max time window shift to compute CC (fraction of whole time window)
cc_delta     =  0.4    # time interval for cc (s)
cc_interp1d  =  5      # interpolation factor
cc_thres     =  0.25    # threshold beam correlation to use in stack
min_amp      =  0.0    # threshold amp to use in stack

# Pro 7 options
do_T = True       # present T plots
do_R = True       # present R plots
zoom = False      # to restrict time range and slowness range in pro7_pair_scan
ZslowR_lo = -0.06
ZslowR_hi =  0.06
ZslowT_lo = -0.06
ZslowT_hi =  0.06
Zstart_buff = 1000
Zend_buff =   1250
start_beam = 0  # Limit time window for summary slowness beam in beam sums
end_beam   = 0  # better be within Zstart and Zend, if zoom is set

# Pro 7 auto_slice == True options
auto_slice      = True  # slices span wide range of R and T slownesses
two_slice_plots = True  # makes R-T pair and snap through time span
beam_sums       = True  # sums tdiff and amp over time
wiggly_plots    = True  # shows wiggly plots

# Pro7 auto-plot options
nR_plots  = 2     # number of plots along the radial axis, makes (2 x nR_plots - 1) total
nT_plots  = 2     # number of plots along the transv axis
slow_incr = 0.01  # increment at which amp and tdiff are plotted

# Pro7 two_slice and snap options
R_slow_plot    =    0.010
T_slow_plot    =    0.000
snaptime       = 1100
snaps          =    0

#%% Comparing events
#%% -- Cull seismic section for common stations
if do_3 == True:
    pro3pair(ARRAY = ARRAY, eq_file1 = eq_file1, eq_file2 = eq_file2, skip_SNR = skip_SNR,
                rel_time = 0, start_buff = start_buff, end_buff = end_buff,
                freq_min = freq_min, freq_max = freq_max,
                max_taper_length = max_taper_length, simple_taper = simple_taper,
                plot_scale_fac = 0.025, stat_corr = stat_corr,
                dphase = ref_phase, dphase2 = 'PKKP', dphase3 = 'PP', dphase4 = 'S',
                min_dist = min_dist, max_dist = max_dist, ref_loc = 0)

#%%  --Cross_convolve time functions
if do_2a == True:
    pro2_convstf(eq_file = eq_file1, conv_file = conv_file1)
    pro2_convstf(eq_file = eq_file2, conv_file = conv_file2)
if do_2b == True:
    print('made it to 3')
    pro2_test(eq_file1 = eq_file1, conv_file1 = conv_file1, eq_file2 = eq_file2, conv_file2 = conv_file2)
    print('made it to 4')

#%%  -- 2D stacks
if do_5 == True:
    pro5stack2d(eq_file = eq_file1,
                slowR_lo = slowR_lo, slowR_hi = slowR_hi, slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
                start_buff = start_buff, end_buff = end_buff, norm = 1, ARRAY = ARRAY, decimate_fac = decimate_fac, NS = NS)

    pro5stack2d(eq_file = eq_file2,
                slowR_lo = slowR_lo, slowR_hi = slowR_hi, slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
                start_buff = start_buff, end_buff = end_buff, norm = 1, ARRAY = ARRAY, decimate_fac = decimate_fac, NS = NS)

# %% -- Compare pair of 2D stack results to find shift, amp, amp ratio, uses cc rather than instant phase
if do_6 == True:
    pro6_cc_pair(eq_file1 = eq_file1, eq_file2 = eq_file2,
                slowR_lo = slowR_lo, slowR_hi = slowR_hi, slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
                start_buff = start_buff, end_buff = end_buff,
                cc_twin = cc_twin, cc_len = cc_len, cc_interp1d = cc_interp1d, cc_delta = cc_delta, cc_thres = cc_thres)

#%% -- Make a variety of plots
if do_7 == True:
    pro7_pair_scan(eq_file1 = eq_file1, eq_file2 = eq_file2, wig_scale_fac = wig_scale_fac, tdiff_scale_fac = tdiff_scale_fac,
                slowR_lo = slowR_lo, slowR_hi = slowR_hi, slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
                zoom = zoom, ZslowR_lo = ZslowR_lo, ZslowR_hi = ZslowR_hi, ZslowT_lo = ZslowT_lo,
                ZslowT_hi = ZslowT_hi, Zstart_buff = Zstart_buff, Zend_buff = Zend_buff,
                start_buff = start_buff, end_buff = end_buff, do_T = do_T, do_R = do_R, tdiff_clip = tdiff_clip,
                min_amp = min_amp, ref_phase = ref_phase, cc_thres = cc_thres,
                R_slow_plot = R_slow_plot, T_slow_plot = T_slow_plot, snaptime = snaptime, snaps = snaps,
                nR_plots  = nR_plots, nT_plots = nT_plots, slow_incr = slow_incr, NS = NS,
                ARRAY = ARRAY, auto_slice = auto_slice, two_slice_plots = two_slice_plots, beam_sums = beam_sums,
                wiggly_plots = wiggly_plots, log_plot = log_plot, log_plot_range = log_plot_range,
                no_tdiff_plot = no_tdiff_plot, start_beam = start_beam, end_beam = end_beam)

#%% just amp, no time shifts estimates
if do_6a == True:
    pro6_singlet(eq_file = eq_file,
                slowR_lo = slowR_lo, slowR_hi = slowR_hi, slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
                start_buff = start_buff, end_buff = end_buff, cc_delta = cc_delta)

#%% -- Make a variety of plots
if do_7a == True:
    pro7_singlet(eq_file = eq_file, wig_scale_fac = wig_scale_fac,
                slowR_lo = slowR_lo, slowR_hi = slowR_hi, slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
                zoom = zoom, ZslowR_lo = ZslowR_lo, ZslowR_hi = ZslowR_hi, ZslowT_lo = ZslowT_lo,
                ZslowT_hi = ZslowT_hi, Zstart_buff = Zstart_buff, Zend_buff = Zend_buff,
                start_buff = start_buff, end_buff = end_buff, do_T = do_T, do_R = do_R,
                min_amp = min_amp, ref_phase = ref_phase,
                R_slow_plot = R_slow_plot, T_slow_plot = T_slow_plot, snaptime = snaptime, snaps = snaps,
                nR_plots  = nR_plots, nT_plots = nT_plots, slow_incr = slow_incr, NS = NS,
                ARRAY = ARRAY, auto_slice = auto_slice, two_slice_plots = two_slice_plots, beam_sums = beam_sums,
                wiggly_plots = wiggly_plots, log_plot = log_plot, log_plot_range = log_plot_range,
                start_beam = start_beam, end_beam = end_beam)

#%% Individual events
#%% -- Cull seismic section event 1
# if do_3 == True:
    # pro3singlet(ARRAY = ARRAY, stat_corr = stat_corr, eq_file = eq_file1,
    #             max_taper_length = max_taper_length, simple_taper = simple_taper,
    #             rel_time = 0, start_buff = start_buff, end_buff = end_buff,
    #             plot_scale_fac = 0.1, skip_SNR = 1,
    #             dphase = ref_phase, dphase2 = 'SKKP', dphase3 = 'PKPPcP', dphase4 = 'pPKIKKIKP',
    #             freq_min = freq_min, freq_max = freq_max,
    #             min_dist = min_dist, max_dist = max_dist, ref_loc = 0, fig_index = 101)

#%% -- Cull seismic section event 2
# if do_3 == True:
    # pro3singlet(ARRAY = ARRAY, stat_corr = stat_corr, eq_file = eq_file2,
    #             max_taper_length = max_taper_length, simple_taper = simple_taper,
    #             rel_time = 0, start_buff = start_buff, end_buff = end_buff,
    #             plot_scale_fac = 0.1, skip_SNR = 1,
    #             dphase = ref_phase, dphase2 = 'PcP', dphase3 = 'PP', dphase4 = 'P',
    #             freq_min = freq_min, freq_max = freq_max,
    #             min_dist = min_dist, max_dist = max_dist, ref_loc = 0, fig_index = 102)

#%% -- 1D stack
# if do_5 == True:
# pro5stack(ARRAY = ARRAY, eq_file = eq_file1, plot_scale_fac = 0.05,
#             slowR_lo = slowR_lo_1D, slowR_hi = slowR_hi_1D, slow_delta = slow_delta_1D,
#             start_buff = start_buff, end_buff = end_buff,
#             log_plot = 0, envelope = 1, plot_dyn_range = 50,
#             norm = 1, global_norm_plot = 1, color_plot = 1, fig_index = 301)

#pro5stack(ARRAY = ARRAY, eq_file = eq_file2, plot_scale_fac = 0.05,
#            slowR_lo = slowR_lo_1D, slowR_hi = slowR_hi_1D, slow_delta = slow_delta_1D,
#            start_buff = start_buff, end_buff = end_buff,
#            log_plot = 0, envelope = 1, plot_dyn_range = 50,
#            norm = 1, global_norm_plot = 1, color_plot = 1, fig_index = 302)

#%%  -- 2D stacks
# if do_5 == True:
#     pro5stack2d(eq_file = eq_file1, slowR_lo = slowR_lo, slowR_hi = slowR_hi,
#                 slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
#                 start_buff = start_buff, end_buff = end_buff,
#                 norm = 1, ARRAY = ARRAY, NS = False)

# pro5stack2d(eq_file = eq_file2, slowR_lo = slowR_lo, slowR_hi = slowR_hi,
#             slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
#             start_buff = start_buff, end_buff = end_buff,
#             norm = 1, ARRAY = ARRAY, NS = False)

#%% -- 2D envelop stack results for individual events
# if do_7 == True:
#     pro7plotstack(eq_file = eq_file1, plot_scale_fac = 0.05,
#                 slowR_lo = slowR_lo, slowR_hi = slowR_hi, slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
#                 start_buff = start_buff, end_buff = end_buff, do_T = True, do_R = True,
#                 zoom = 0, ZslowR_lo = -0.03, ZslowR_hi = 0.03, ZslowT_lo = -0.03, ZslowT_hi = 0.03, Zstart_buff = 0, Zend_buff = 200,
#                 fig_index = 401, plot_dyn_range = 50, snaptime = snaptime, snaps=0, ARRAY = ARRAY)
#
#pro7plotstack(eq_file = eq_file2, plot_scale_fac = 0.05,
#            slowR_lo = slowR_lo, slowR_hi = slowR_hi, slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
#            start_buff = start_buff, end_buff = end_bufff, do_T = True, do_R = Tru,
#            zoom = 0, ZslowR_lo = -0.03, ZslowR_hi = 0.03, ZslowT_lo = -0.03, ZslowT_hi = 0.03, Zstart_buff = 0, Zend_buff = 200,
#            fig_index = 402, plot_dyn_range = 50, snaptime = snaptime, snaps=0, ARRAY = ARRAY)
