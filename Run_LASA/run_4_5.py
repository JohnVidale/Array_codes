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
from pro7_pair_scan            import pro7_pair_scan
from pro7_singlet              import pro7_singlet

#%% Workflow selection
do_3  = True  # pair of events
do_5  = True
do_6  = True
do_7  = True
eq_num1 = '4'  # pair
eq_num2 = '5'

do_3a = False # single event
do_5a = False # stack
do_6a = False # treats single events, no time shifts calculated or plotted
do_7a = False
eq_num  = '4'  # singlet

#%% Common parameters
ARRAY      = 1
ref_loc = False   # if true,  use ref_rad + distance to filter station distance
                  # if false, use earthquake distance to filter station distance
ref_rad = 0.4 # radius of stations around the array center included
freq_min = 1
freq_max = 3

# Time and slowness window for stacking
# start_buff =  990
# end_buff   = 1320
# Zstart_buff = 1025
# Zend_buff =   1250
start_buff =  980
end_buff   = 1250
# Pro 7 range selection options, necessary in thin window stacks
zoom = True      # to restrict time range and slowness range in pro7_pair_scan
Zstart_buff = 1030
Zend_buff =   1200
slowR_lo       = -0.04
slowR_hi       =  0.04
slowT_lo       = -0.04
slowT_hi       =  0.04
slow_delta     =  0.002
ZslowR_lo = -0.04
ZslowR_hi =  0.04
ZslowT_lo = -0.04
ZslowT_hi =  0.04

# Station distance range
# min_dist = 58.5 # events 1&2
# max_dist = 60.5
min_dist = 60 # events 4&5
max_dist = 64
# min_dist = 46.2 # events 7&8
# max_dist = 48.2

# Pro5 stacking
stat_corr        = 1 # applied statics from reference event
decimate_fac     = 5 # set for pro5stack2d for single event envelopes, set to 0 for other codes
simple_taper     = 1
max_taper_length = 5. # taper is minimum of taper_frac (0.05) and this number of seconds
skip_SNR         = 1  # only useful when there is a sharp phase onset
rel_time         = 0  # sets relation of window to reference phase and its slowness
ref_phase        = 'PKiKP'  # only used in some rel_time options
NS               = False  # True for N-S co=ords, False for R-T

# Pro5 1D plot options
slowR_lo_1D   = -0.04
slowR_hi_1D   =  0.1
slow_delta_1D =  0.001

# Pro6 options: mostly time shift measurement
cc_twin      =  10   # time window for cross-correlation (s)
cc_len       =  0.02 # max time window shift to compute CC (fraction of whole time window)
cc_delta     =  2    # temporal frequency of cc (s)
cc_interp1d  =  5    # interpolation factor
cc_thres     =  0.2  # threshold beam correlation to use in stack

# more Pro 7 range selection options
start_beam = 0  # Limit time window for summary slowness beam in beam sums
end_beam   = 0  # better be within Zstart and Zend, if zoom is set
min_amp    = 0.1    # threshold ratio of ave amp to peak amp to use in stack

# Pro 7 auto_slice == True options
auto_slice      = True  # slices span wide range of R and T slownesses
two_slice_plots = True  # makes R-T pair and snap through time span
beam_sums       = True  # sums tdiff and amp over time
wiggly_plots    = False  # shows wiggly plots

# Pro7 auto-plot options
nR_plots  = 2     # number of plots along the radial axis, makes (2 x nR_plots - 1) total
nT_plots  = 2     # number of plots along the transv axis
slow_incr = 0.01  # increment at which amp and tdiff are plotted

# Pro7 two_slice and snap options
R_slow_plot    =    0.010
T_slow_plot    =    0.000
snaptime       = 1100
snaps          =    0

# Pro 7 more plotting options
do_T = False      # present T plots
do_R = True       # present R plots
no_tdiff_plot = False  # also to speed plots of only amplitude
log_plot      = False
tdiff_clip   =  0.2
wig_scale_fac = 0.5
tdiff_scale_fac = 1
log_plot_range = 1
plot_scale_fac = 1

#%% Comparing events
#%% -- Cull seismic section for common stations
if do_3 == True:
    pro3pair(ARRAY = ARRAY, eq_num1 = eq_num1, eq_num2 = eq_num2, skip_SNR = skip_SNR,
                rel_time = rel_time, start_buff = start_buff, end_buff = end_buff,
                freq_min = freq_min, freq_max = freq_max,
                max_taper_length = max_taper_length, simple_taper = simple_taper,
                plot_scale_fac = 0.025, stat_corr = stat_corr,
                dphase = ref_phase, dphase2 = 'PKKP', dphase3 = 'PP', dphase4 = 'S',
                min_dist = min_dist, max_dist = max_dist, ref_loc = ref_loc, ref_rad = ref_rad)

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
                R_slow_plot = R_slow_plot, T_slow_plot = T_slow_plot, snaptime = snaptime, snaps = snaps,
                nR_plots  = nR_plots, nT_plots = nT_plots, slow_incr = slow_incr, NS = NS,
                ARRAY = ARRAY, auto_slice = auto_slice, two_slice_plots = two_slice_plots, beam_sums = beam_sums,
                wiggly_plots = wiggly_plots, log_plot = log_plot, log_plot_range = log_plot_range,
                no_tdiff_plot = no_tdiff_plot, start_beam = start_beam, end_beam = end_beam)

#%% Individual event
#%% -- Cull seismic section event
if do_3a == True:
    pro3singlet(ARRAY = ARRAY, stat_corr = stat_corr, eq_num = eq_num,
                max_taper_length = max_taper_length, simple_taper = simple_taper,
                rel_time = rel_time, start_buff = start_buff, end_buff = end_buff,
                plot_scale_fac = 0.1, skip_SNR = 1,
                dphase = ref_phase, dphase2 = 'SKKP', dphase3 = 'PKPPcP', dphase4 = 'pPKIKKIKP',
                freq_min = freq_min, freq_max = freq_max,
                min_dist = min_dist, max_dist = max_dist, ref_loc = 0, fig_index = 101)
#%% -- 1D stack
# if do_5 == True:
# pro5stack(ARRAY = ARRAY, eq_num = eq_num, plot_scale_fac = 0.05,
#             slowR_lo = slowR_lo_1D, slowR_hi = slowR_hi_1D, slow_delta = slow_delta_1D,
#             start_buff = start_buff, end_buff = end_buff,
#             log_plot = 0, envelope = 1, plot_dyn_range = 50,
#             norm = 1, global_norm_plot = 1, color_plot = 1, fig_index = 301)

#%%  -- 2D stack
if do_5a == True:
    pro5stack2d(eq_num = eq_num, slowR_lo = slowR_lo, slowR_hi = slowR_hi, slowT_lo = slowT_lo,
                slowT_hi = slowT_hi, slow_delta = slow_delta,
                start_buff = start_buff, end_buff = end_buff, norm = 1,
                ARRAY = ARRAY, decimate_fac = decimate_fac, NS = NS)

#%% just amp, no time shifts estimates
if do_6a == True:
    pro6_singlet(eq_num = eq_num,
                slowR_lo = slowR_lo, slowR_hi = slowR_hi, slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
                start_buff = start_buff, end_buff = end_buff, cc_delta = cc_delta)

#%% -- Make a variety of plots
if do_7a == True:
    pro7_singlet(eq_num = eq_num, wig_scale_fac = wig_scale_fac,
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