#!/usr/bin/env python
# input is set of hinet traces
# this program tapers, filters, selects range and SNR
# plots against traveltime curves, either raw or reduced against traveltimes
# This programs deals with a single event.
# John Vidale 2/2019

def runpair(Tstart, Tend, eq_num1, eq_num2, decon78 = False):

    import os
    import matplotlib.pyplot as plt
    #%% close plots
    plt.close('all')

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
    eq_num1 = str(eq_num1)  # pair
    eq_num2 = str(eq_num2)

    do_3a = False # single event
    do_5a = False # stack
    do_6a = False # treats single events, no time shifts calculated or plotted
    do_7a = False
    eq_num  = '2'  # singlet

    #%% Common parameters
    ARRAY      = 1
    no_plots = True
    auto_dist = True
    min_dist = 0
    max_dist = 180

    # Window
    start_buff = Tstart - 20
    end_buff   = Tend + 20
    zoom = True      # to restrict time range and slowness range in pro7_pair_scan
    Zstart_buff = Tstart
    Zend_buff =   Tend

    # HF
    freq_min = 1
    freq_max = 3

    # Pro5 stacking
    stat_corr      =    1
    decimate_fac   =    5 # set for pro5stack2d for single event envelopes, set to 0 for other codes
    simple_taper   =    1
    max_taper_length =  5 # taper is minimum of taper_frac (0.05) and this number of seconds
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

    # Pro6 options: mostly time shift measurement
    cc_twin      =  7     # time window for cross-correlation (s)
    cc_len       =  0.05 # max time window shift to compute CC (fraction of whole time window)
    cc_delta     =  0.4    # temporal frequency of cc (s)
    cc_interp1d  =  5      # interpolation factor
    cc_thres     =  0.7    # threshold beam correlation to use in stack
    min_amp      =  0.2    # threshold amp to use in stack

    # Pro 7 range selection options
    ZslowR_lo       = -0.03
    ZslowR_hi       =  0.03
    ZslowT_lo       = -0.03
    ZslowT_hi       =  0.03
    start_beam = 0  # Limit time window for summary slowness beam in beam sums
    end_beam   = 0  # better be within Zstart and Zend, if zoom is set

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
    R_slow_plot    =    0.010
    T_slow_plot    =    0.000
    snaptime       =    0  # relative to start_buff
    snaps          =   10
    snap_depth     =   30  # time window over which snap is integrated (s)

    # Pro 7 more plotting options
    do_T = False           # present T plots
    do_R = True            # present R plots
    no_tdiff_plot = False  # also to speed plots of only amplitude, only applies to auto_slice
    turn_off_black = False # controls whether wiggle plot also has time shift plotted
    log_plot      = True
    tdiff_clip   =  0.2
    wig_scale_fac = 0.5
    tdiff_scale_fac = 3
    log_plot_range = 1.5
    plot_scale_fac = 1

    #%% Comparing events
    #%% -- Cull seismic section for common stations
    if do_3 == True:

        pro3pair(ARRAY = ARRAY, eq_num1 = eq_num1, eq_num2 = eq_num2, skip_SNR = skip_SNR,
                    rel_time = 0, start_buff = start_buff, end_buff = end_buff,
                    freq_min = freq_min, freq_max = freq_max, no_plots = no_plots,
                    max_taper_length = max_taper_length, simple_taper = simple_taper,
                    plot_scale_fac = 0.025, stat_corr = stat_corr,
                    dphase = ref_phase, dphase2 = 'PKKP', dphase3 = 'PP', dphase4 = 'S',
                    min_dist = min_dist, max_dist = max_dist, auto_dist = auto_dist, ref_loc = 0)

        if decon78 == True:
            # Cross_convolve time functions only set up for Amchitka events 7 & 8
            conv_file1 = '/Users/vidale/Documents/GitHub/Array_codes/Files/HD1971-11-06_stf.mseed'
            conv_file2 = '/Users/vidale/Documents/GitHub/Array_codes/Files/HD1969-10-02_stf.mseed'
            pro2_convstf(eq_num = eq_num1, conv_file = conv_file1)
            pro2_convstf(eq_num = eq_num2, conv_file = conv_file2)
            # pro2_test(eq_num1 = eq_num1, conv_file1 = conv_file1, eq_num2 = eq_num2, conv_file2 = conv_file2)

    #%%  -- 2D stacks
    if do_5 == True:
        pro5stack2d(eq_num = eq_num1,
                    slowR_lo = slowR_lo, slowR_hi = slowR_hi, slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
                    start_buff = start_buff, end_buff = end_buff, norm = 1, ARRAY = ARRAY, decimate_fac = decimate_fac, NS = NS)

        pro5stack2d(eq_num = eq_num2,
                    slowR_lo = slowR_lo, slowR_hi = slowR_hi, slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
                    start_buff = start_buff, end_buff = end_buff, norm = 1, ARRAY = ARRAY, decimate_fac = decimate_fac, NS = NS)

    #%% -- Compare pair of 2D stack results to find shift, amp, amp ratio, uses cc rather than instant phase
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
                    min_amp = min_amp, ref_phase = ref_phase, cc_thres = cc_thres, turn_off_black = turn_off_black,
                    R_slow_plot = R_slow_plot, T_slow_plot = T_slow_plot, no_plots = no_plots,
                    snaptime = snaptime, snaps = snaps, snap_depth = snap_depth,
                    nR_plots  = nR_plots, nT_plots = nT_plots, slow_incr = slow_incr, NS = NS,
                    ARRAY = ARRAY, auto_slice = auto_slice, two_slice_plots = two_slice_plots, beam_sums = beam_sums,
                    wiggly_plots = wiggly_plots, log_plot = log_plot, log_plot_range = log_plot_range,
                    no_tdiff_plot = no_tdiff_plot, start_beam = start_beam, end_beam = end_beam)

    #%% Individual event
    #%% -- Cull seismic section event
    if do_3a == True:
        pro3singlet(ARRAY = ARRAY, stat_corr = stat_corr, eq_num = eq_num,
                    max_taper_length = max_taper_length, simple_taper = simple_taper,
                    rel_time = 0, start_buff = start_buff, end_buff = end_buff,
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
                    min_amp = min_amp, ref_phase = ref_phase, turn_off_black = turn_off_black,
                    R_slow_plot = R_slow_plot, T_slow_plot = T_slow_plot,
                    snaptime = snaptime, snaps = snaps, snap_depth = snap_depth,
                    nR_plots  = nR_plots, nT_plots = nT_plots, slow_incr = slow_incr, NS = NS,
                    ARRAY = ARRAY, auto_slice = auto_slice, two_slice_plots = two_slice_plots, beam_sums = beam_sums,
                    wiggly_plots = wiggly_plots, log_plot = log_plot, log_plot_range = log_plot_range,
                    start_beam = start_beam, end_beam = end_beam)