#!/usr/bin/env python
# This programs deals with a single event.
# John Vidale 2/2019, still modifying 2/2021

import os

os.environ['PATH'] += os.pathsep + '/usr/local/bin'
os.chdir('/Users/vidale/Documents/GitHub/Array_codes')

#%% Import functions
from pro3b_sort_plot_singlet   import pro3singlet
from pro5a_stack               import pro5stack
from pro5b_stack2d             import pro5stack2d
from pro6_singlet              import pro6_singlet
from pro7_singlet              import pro7_singlet

#%% Workflow selection

do_3a = False # single event
do_5a = False # stack
do_6a = False # treats single events, no time shifts calculated or plotted
do_7a = True
eq_num = 91  # index on location file

#%% Common parameters
ARRAY      = 1
auto_dist = True
min_dist = 0
max_dist = 180

# P
# start_buff =  600
# end_buff   = 2500
start_buff = 1350
end_buff   = 1450
# start_buff = 1550
# end_buff   = 1750

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
slowR_lo       = -0.1
slowR_hi       =  0.1
slowT_lo       = -0.1
slowT_hi       =  0.1
slow_delta     =  0.0025
NS = True  # 1 for N-S co=ords, 0 for R-T

# Pro5 1D plot options
slowR_lo_1D   = -0.04
slowR_hi_1D   =  0.1
slow_delta_1D =  0.001

# Pro6 options: mostly time shift measurement
cc_twin      =  3     # time window for cross-correlation (s)
cc_len       =  0.2 # max time window shift to compute CC (fraction of whole time window)
cc_delta     =  0.4    # time interval for cc (s)
cc_interp1d  =  5      # interpolation factor
cc_thres     =  0.5    # threshold beam correlation to use in stack

# Pro 7 range selection options
zoom = False      # to restrict time range and slowness range in pro7_pair_scan
ZslowR_lo = -0.06
ZslowR_hi =  0.06
ZslowT_lo = -0.06
ZslowT_hi =  0.06
Zstart_buff = 1000
Zend_buff =   1250
start_beam = 0  # Limit time window for summary slowness beam in beam sums
end_beam   = 0  # better be within Zstart and Zend, if zoom is set
min_amp      =  0.0    # threshold amp to use in stack

# Pro 7 auto_slice == True options
auto_slice      = True  # slices span wide range of R and T slownesses
two_slice_plots = True  # makes R-T pair and snap through time span
beam_sums       = True  # sums tdiff and amp over time
wiggly_plots    = False  # shows wiggly plots

# Pro7 auto-plot options
nR_plots  = 1     # number of plots along the radial axis, makes (2 x nR_plots - 1) total
nT_plots  = 1     # number of plots along the transv axis
slow_incr = 0.06  # increment at which amp and tdiff are plotted

# Pro7 two_slice and snap options
R_slow_plot    =    0.010
T_slow_plot    =    0.000
snaptime       = 1350
snaps          =   10
snap_depth     =   10  # time window over which snap is integrated, 0 is one time point

# Pro 7 more plotting options
do_T = False       # present T plots
do_R = True       # present R plots
no_tdiff_plot = True  # also to speed plots of only amplitude
log_plot      = False
tdiff_clip   =  0.15
wig_scale_fac = 0.5
tdiff_scale_fac = 1
log_plot_range = 1.5
plot_scale_fac = 1

#%% Individual event
#%% -- Cull seismic section event
if do_3a == True:
    pro3singlet(ARRAY = ARRAY, stat_corr = stat_corr, eq_num = eq_num,
                max_taper_length = max_taper_length, simple_taper = simple_taper,
                rel_time = 0, start_buff = start_buff, end_buff = end_buff,
                plot_scale_fac = 0.1, skip_SNR = 1,
                dphase = ref_phase, dphase2 = 'SKKP', dphase3 = 'PKPPcP', dphase4 = 'pPKIKKIKP',
                freq_min = freq_min, freq_max = freq_max,
                ref_loc = 0, fig_index = 101, min_dist = min_dist, max_dist = max_dist, auto_dist = auto_dist)
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
    pro6_singlet(eq_num = eq_num, slowR_lo = slowR_lo, slowR_hi = slowR_hi,
                slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
                start_buff = start_buff, end_buff = end_buff, cc_delta = cc_delta)

#%% -- Make a variety of plots
if do_7a == True:
    pro7_singlet(eq_num = eq_num, wig_scale_fac = wig_scale_fac,
                slowR_lo = slowR_lo, slowR_hi = slowR_hi, slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
                zoom = zoom, ZslowR_lo = ZslowR_lo, ZslowR_hi = ZslowR_hi, ZslowT_lo = ZslowT_lo,
                ZslowT_hi = ZslowT_hi, Zstart_buff = Zstart_buff, Zend_buff = Zend_buff,
                start_buff = start_buff, end_buff = end_buff, do_T = do_T, do_R = do_R,
                min_amp = min_amp, ref_phase = ref_phase,
                R_slow_plot = R_slow_plot, T_slow_plot = T_slow_plot,
                snaptime = snaptime, snaps = snaps, snap_depth =snap_depth,
                nR_plots  = nR_plots, nT_plots = nT_plots, slow_incr = slow_incr, NS = NS,
                ARRAY = ARRAY, auto_slice = auto_slice, two_slice_plots = two_slice_plots, beam_sums = beam_sums,
                wiggly_plots = wiggly_plots, log_plot = log_plot, log_plot_range = log_plot_range,
                start_beam = start_beam, end_beam = end_beam)

os.system('say "All Done"')