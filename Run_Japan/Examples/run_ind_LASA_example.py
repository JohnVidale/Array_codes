#!/usr/bin/env python
# purpose
# John Vidale 6/2021

def run_ind_LASA_example(eq_num = '1', start_time = 0, wind_len = 100, ref_phase = 'PKiKP',
                   slow_limit = 0.04, slow_delta = 0.0025, wind_buff = 50):

    import os
    import matplotlib.pyplot as plt
    #%% close plots
    plt.close('all')

    #%% Import functions
    from pro3b_sort_plot_singlet   import pro3singlet
    from pro5a_stack               import pro5stack
    from pro5b_stack2d             import pro5stack2d
    from pro6_singlet              import pro6_singlet
    from pro7_singlet              import pro7_singlet

    do_3a = True # single event
    do_5a = True # stack
    do_6a = True # treats single events, no time shifts calculated or plotted
    do_7a = True
    # eq_num  = '12'  # singlet

    #%% Common parameters
    ref_loc = True   # if true,  use ref_rad + distance to filter station distance
                      # if false, use earthquake distance to filter station distance
    ref_rad = 0.4 # radius of stations around the array center included
    ARRAY      = 1
    auto_dist = True
    min_dist = 0
    max_dist = 180

    # Window
    zoom = True      # to restrict time range and slowness range in pro7_pair_scan
    Zstart_buff = start_time
    wind_buff = wind_buff
    wind_len = wind_len
    Zend_buff =   Zstart_buff + wind_len

    start_buff = Zstart_buff - wind_buff
    end_buff   = Zstart_buff + wind_len + wind_buff

    # HF
    freq_min = 0.5
    freq_max = 1.5

    # Pro5 stacking
    stat_corr = 1
    decimate_fac   =    5 # set for pro5stack2d for single event envelopes, set to 0 for other codes
    simple_taper   =    1
    max_taper_length =  5 # taper is minimum of taper_frac (0.05) and this number of seconds
    skip_SNR       =    1
    # ref_phase      = 'PKiKP'
    slowR_lo       = -slow_limit
    slowR_hi       =  slow_limit
    slowT_lo       = -slow_limit
    slowT_hi       =  slow_limit
    # slow_delta     =  0.005
    NS = False  # True for N-S co=ords, False for R-T

    # Pro5 1D plot options
    slowR_lo_1D   = -0.04
    slowR_hi_1D   =  0.1
    slow_delta_1D =  0.001

    # Pro6 decimation
    cc_delta     =  0.1    # temporal frequency of output (s)

    # Pro 7 range selection options
    ZslowR_lo       = -slow_limit
    ZslowR_hi       =  slow_limit
    ZslowT_lo       = -slow_limit
    ZslowT_hi       =  slow_limit
    start_beam = 0  # Limit time window for summary slowness beam in beam sums
    end_beam   = 0  # better be within Zstart and Zend, if zoom is set
    min_amp      =  0.0    # threshold amp to use in stack

    # Pro 7 auto_slice == True options
    auto_slice      = False  # slices span wide range of R and T slownesses
    two_slice_plots = True  # makes R-T pair and snap through time span
    beam_sums       = True  # sum amp over time
    wiggly_plots    = False  # shows wiggly plots

    # Pro7 auto-plot options
    nR_plots  = 2     # number of plots along the radial axis, makes (2 x nR_plots - 1) total
    nT_plots  = 2     # number of plots along the transv axis
    slow_incr = 0.01  # increment at which amp and tdiff are plotted

    # Pro7 two_slice and snap options
    R_slow_plot    =    0.012
    T_slow_plot    =    0.00
    snaptime       =    0  # relative to start_buff
    snaps          =    0
    snap_depth     =    5  # time window over which snap is integrated (s)

    # Pro 7 more plotting options
    do_T = True       # present T plots
    do_R = True       # present R plots
    log_plot      = True
    wig_scale_fac = 0.5
    log_plot_range = 2.0
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
                    min_dist = min_dist, max_dist = max_dist, ref_loc = ref_loc, ref_rad = ref_rad, fig_index = 101)
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
                    R_slow_plot = R_slow_plot, T_slow_plot = T_slow_plot,
                    snaptime = snaptime, snaps = snaps, snap_depth = snap_depth,
                    nR_plots  = nR_plots, nT_plots = nT_plots, slow_incr = slow_incr, NS = NS,
                    ARRAY = ARRAY, auto_slice = auto_slice, two_slice_plots = two_slice_plots, beam_sums = beam_sums,
                    wiggly_plots = wiggly_plots, log_plot = log_plot, log_plot_range = log_plot_range,
                    start_beam = start_beam, end_beam = end_beam)