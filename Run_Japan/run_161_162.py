#!/usr/bin/env python
# input is set of hinet traces
# this program tapers, filters, selects range and SNR
# plots against traveltime curves, either raw or reduced against traveltimes
# This programs deals with a single event.
# John Vidale 2/2019

import os
import matplotlib.pyplot as plt
#%% close plots
plt.close('all')

os.environ['PATH'] += os.pathsep + '/usr/local/bin'
os.chdir('/Users/vidale/Documents/GitHub/Array_codes')

#%% Import functions
ev_directory = '/Users/vidale/Documents/GitHub/Array_codes/Process'
os.chdir(ev_directory)
from pro3_sort_plot_pair       import pro3pair
from pro3_sort_plot_singlet    import pro3singlet
from pro5_stack1d              import pro5stack1d
from pro5_stack2d              import pro5stack2d
from pro6_pair_cc              import pro6_cc_pair
from pro6_singlet              import pro6_singlet
from pro7_pair_scan            import pro7_pair_scan
from pro7_singlet              import pro7_singlet

#%% Workflow selection
do_3  = False  # pair of events
do_5  = False
do_6  = False
do_7  = True
# eq_file1 = 'event1.txt'  # pair
# eq_file2 = 'event2.txt'
eq_num1 = '161'  # pair
eq_num2 = '162'

do_3a = False # single event
do_5a = False # stack
do_6a = False # treats single events, no time shifts calculated or plotted
do_7a = False
eq_num  = '161'  # singlet

#%% Common parameters
ARRAY      = 0
auto_dist = True
min_dist = 0
max_dist = 180

# HF
freq_min = 2
freq_max = 4

# Pro5 stacking
stat_corr = 1         # 0 no corr, 1 775 by Wei, 2 219 from SSI, 3 300 from Kawakatsu
rel_time  = 3         # time alignment in shift
decimate_fac     =  5 # set for pro5stack2d for single event envelopes, set to 0 for other codes
simple_taper     =  1
max_taper_length =  5 # taper is minimum of taper_frac (0.05) and this number of seconds
apply_SNR        =  False
SNR_thres        =  1.2
corr_threshold   =  0.5
precursor_shift  = -2
signal_dur       =  4

ref_phase      = 'PKIKP'

# Window
zoom = True                                     # to restrict time range and slowness range in pro7_pair_scan

# time
wind_buff   = 30                                # buffer before time window of analysis
Zstart_buff =  0                                # analysis window start relative to phase arrival
wind_len    = 20                                # analysis window length
Zend_buff   = Zstart_buff + wind_len            # analysis window end
start_buff  = Zstart_buff - wind_buff            # trace start relative to analysis window (if Zoom is True)
end_buff    = Zstart_buff + wind_len + wind_buff # trace end

# slowness
beam_offset = 0.02
beam_width  = 0.02
slow_delta  = 0.0004
slowR_lo    = -beam_width + beam_offset
slowR_hi    =  beam_width + beam_offset
slowT_lo    = -beam_width
slowT_hi    =  beam_width
ZslowR_lo   = -beam_width + beam_offset
ZslowR_hi   =  beam_width + beam_offset
ZslowT_lo   = -beam_width
ZslowT_hi   =  beam_width

NS = False  # True for N-S co=ords, False for R-T

# Pro5 1D plot options
slowR_lo_1D   = -0.04
slowR_hi_1D   =  0.1
slow_delta_1D =  0.001

# Pro6 options: mostly time shift measurement
cc_twin      =  3      # time window for cross-correlation (s)
cc_len       =  0.033  # max time window shift to compute CC (fraction of whole time window)
cc_delta     =  0.2    # temporal frequency of cc (s)
cc_interp1d  =  5      # interpolation factor
cc_thres     =  0.75   # threshold beam correlation to use in stack

# Pro 7 range selection options
start_beam = 0  # Limit time window for summary slowness beam in beam sums
end_beam   = 0  # better be within Zstart and Zend, if zoom is set
min_amp    = 0.0    # threshold amp to use in stack

# Pro 7 auto_slice == True options
auto_slice      = False  # slices span wide range of R and T slownesses
two_slice_plots = True  # makes R-T pair and snap through time span
beam_sums       = True  # sums tdiff and amp over time
wiggly_plots    = True  # shows wiggly plots

# Pro7 auto-plot options
nR_plots  = 2     # number of plots along the radial axis, makes (2 x nR_plots - 1) total
nT_plots  = 2     # number of plots along the transv axis
slow_incr = 0.01  # increment at which amp and tdiff are plotted

# Pro7 two_slice and snap options
R_slow_plot    =    0.019
T_slow_plot    =    0.000
snaptime       =    0  # relative to start_buff
snaps          =    0
snap_depth     =   20  # time window over which snap is integrated (s)

# Pro 7 more plotting options
do_T = False       # present T plots
do_R = True       # present R plots
tdiff_plots_too = False  # only applies to auto_slice - to speed plots of only amplitude
log_plot      = True
# tdiff_clip   =  0.15
tdiff_clip   =  cc_twin * cc_len * 0.99
wig_scale_fac = 0.5
tdiff_scale_fac = 1
log_plot_range = 1.5
plot_scale_fac = 0.1

#%% Comparing events
#%% -- Cull seismic section for common stations
if do_3 == True:
    pro3pair(ARRAY = ARRAY, eq_num1 = eq_num1, eq_num2 = eq_num2, apply_SNR = apply_SNR,
                rel_time = rel_time, start_buff = start_buff, end_buff = end_buff,
                freq_min = freq_min, freq_max = freq_max,
                SNR_thres = SNR_thres, corr_threshold = corr_threshold,
                precursor_shift = precursor_shift, signal_dur = signal_dur,
                max_taper_length = max_taper_length, simple_taper = simple_taper,
                plot_scale_fac = plot_scale_fac, stat_corr = stat_corr,
                dphase = ref_phase, dphase2 = 'PKiKP', dphase3 = 'PKP', dphase4 = 'pPKIKP',
                min_dist = min_dist, max_dist = max_dist, auto_dist = auto_dist, ref_loc = 0)

#%%  -- 2D stacks
if do_5 == True:
    pro5stack2d(eq_num = eq_num1,
                slowR_lo = slowR_lo, slowR_hi = slowR_hi, slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
                start_buff = start_buff, end_buff = end_buff, norm = 1, ARRAY = ARRAY, decimate_fac = decimate_fac, NS = NS)

    pro5stack2d(eq_num = eq_num2,
                slowR_lo = slowR_lo, slowR_hi = slowR_hi, slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
                start_buff = start_buff, end_buff = end_buff, norm = 1, ARRAY = ARRAY, decimate_fac = decimate_fac, NS = NS)

# %% -- Compare pair of 2D stack results to find shift, amp, amp ratio, uses cc rather than instant phase
if do_6 == True:
    pro6_cc_pair(eq_num1 = eq_num1, eq_num2 = eq_num2,
                slowR_lo = slowR_lo, slowR_hi = slowR_hi, slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
                start_buff = start_buff, end_buff = end_buff,
                cc_twin = cc_twin, cc_len = cc_len, cc_interp1d = cc_interp1d, cc_delta = cc_delta, cc_thres = cc_thres)

#%% -- Make a variety of plots
if do_7 == True:
    pro7_pair_scan(eq_num1 = eq_num1, eq_num2 = eq_num2, wig_scale_fac = wig_scale_fac, tdiff_scale_fac = tdiff_scale_fac,
                slowR_lo = slowR_lo, slowR_hi = slowR_hi, slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
                zoom = zoom, ZslowR_lo = ZslowR_lo, ZslowR_hi = ZslowR_hi, ZslowT_lo = ZslowT_lo,
                ZslowT_hi = ZslowT_hi, Zstart_buff = Zstart_buff, Zend_buff = Zend_buff,
                start_buff = start_buff, end_buff = end_buff, do_T = do_T, do_R = do_R, tdiff_clip = tdiff_clip,
                min_amp = min_amp, ref_phase = ref_phase, cc_thres = cc_thres,
                R_slow_plot = R_slow_plot, T_slow_plot = T_slow_plot,
                snaptime = snaptime, snaps = snaps, snap_depth = snap_depth,
                nR_plots  = nR_plots, nT_plots = nT_plots, slow_incr = slow_incr, NS = NS,
                ARRAY = ARRAY, auto_slice = auto_slice, two_slice_plots = two_slice_plots, beam_sums = beam_sums,
                wiggly_plots = wiggly_plots, log_plot = log_plot, log_plot_range = log_plot_range,
                tdiff_plots_too = tdiff_plots_too, start_beam = start_beam, end_beam = end_beam)

#%% Individual event
#%% -- Cull seismic section event
# if do_3a == True:
#     pro3singlet(ARRAY = ARRAY, stat_corr = stat_corr, eq_num = eq_num,
#                 max_taper_length = max_taper_length, simple_taper = simple_taper,
#                 rel_time = rel_time, start_buff = start_buff, end_buff = end_buff,
#                 plot_scale_fac = 0.1, apply_SNR = False,
#                 dphase = ref_phase, dphase2 = 'SKKP', dphase3 = 'PKPPcP', dphase4 = 'pPKIKKIKP',
#                 freq_min = freq_min, freq_max = freq_max,
#                 min_dist = min_dist, max_dist = max_dist, ref_loc = 0, fig_index = 101)
#%% -- 1D stack
# if do_5 == True:
# pro5stack(ARRAY = ARRAY, eq_num = eq_num, plot_scale_fac = 0.05,
#             slowR_lo = slowR_lo_1D, slowR_hi = slowR_hi_1D, slow_delta = slow_delta_1D,
#             start_buff = start_buff, end_buff = end_buff,
#             log_plot = 0, envelope = 1, plot_dyn_range = 50,
#             norm = 1, global_norm_plot = 1, color_plot = 1, fig_index = 301)

#%%  -- 2D stack
# if do_5a == True:
#     pro5stack2d(eq_num = eq_num, slowR_lo = slowR_lo, slowR_hi = slowR_hi, slowT_lo = slowT_lo,
#                 slowT_hi = slowT_hi, slow_delta = slow_delta,
#                 start_buff = start_buff, end_buff = end_buff, norm = 1,
#                 ARRAY = ARRAY, decimate_fac = decimate_fac, NS = NS)

#%% just amp, no time shifts estimates
# if do_6a == True:
#     pro6_singlet(eq_num = eq_num,
#                 slowR_lo = slowR_lo, slowR_hi = slowR_hi, slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
#                 start_buff = start_buff, end_buff = end_buff, cc_delta = cc_delta)

#%% -- Make a variety of plots
# if do_7a == True:
#     pro7_singlet(eq_num = eq_num, wig_scale_fac = wig_scale_fac,
#                 slowR_lo = slowR_lo, slowR_hi = slowR_hi, slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
#                 zoom = zoom, ZslowR_lo = ZslowR_lo, ZslowR_hi = ZslowR_hi, ZslowT_lo = ZslowT_lo,
#                 ZslowT_hi = ZslowT_hi, Zstart_buff = Zstart_buff, Zend_buff = Zend_buff,
#                 start_buff = start_buff, end_buff = end_buff, do_T = do_T, do_R = do_R,
#                 min_amp = min_amp, ref_phase = ref_phase,
#                 R_slow_plot = R_slow_plot, T_slow_plot = T_slow_plot,
#                 snaptime = snaptime, snaps = snaps, snap_depth = snap_depth,
#                 nR_plots  = nR_plots, nT_plots = nT_plots, slow_incr = slow_incr, NS = NS,
#                 ARRAY = ARRAY, auto_slice = auto_slice, two_slice_plots = two_slice_plots, beam_sums = beam_sums,
#                 wiggly_plots = wiggly_plots, log_plot = log_plot, log_plot_range = log_plot_range,
#                 start_beam = start_beam, end_beam = end_beam)
