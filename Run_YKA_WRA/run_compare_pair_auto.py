#!/usr/bin/env python
# John Vidale 7/2022
# input is pair of traces
# this program tapers, filters, selects range and SNR
# plots against traveltime curves, either raw or reduced against traveltimes

def run_compare_pair(eq_num1 = 401, eq_num2 = 402, Zstart_buff =  0, wind_len    = 20, wind_buff = 30,
                     precursor_shift  = -1000, signal_dur     = -1000, ARRAY = 5, pair_name = '',
                     beam_width = 0.04, beam_offset = 0.00, slow_delta = 0.0025, flip = False, plot_peak = 1,
                     zerophase = True, freq_min = 1, freq_max = 3, temp_shift_both = 0,
                     min_dist = 0, max_dist = 180, dphase  = 'PKiKP', fig_index = 100,
                     JST = False, R_slow_plot = 0.017, T_slow_plot = 0, temp_shift2 = 0,
                     stat_corr = 1, apply_SNR = False, beam_stack_rad = 0.01, do_interpolate = False):

    import os
    #%% load functions

    os.environ['PATH'] += os.pathsep + '/usr/local/bin'
    ev_directory = '/Users/vidale/Documents/GitHub/Array_codes/Process'
    os.chdir(ev_directory)
    #%% Import functions
    # from pro2_dec                import pro2decimate
    from pro3_sort_plot_pair     import pro3pair
    from pro3_sort_plot_singlet  import pro3singlet
    from pro4_get_shifts         import pro4_get_shifts
    from pro5_stack1d            import pro5stack1d
    from pro5_stack2d            import pro5stack2d
    from pro6_singlet            import pro6_singlet
    from pro6_pair_cc            import pro6_cc_pair
    from pro7_pair_scan          import pro7_pair_scan
    from pro7_singlet            import pro7_singlet
    import matplotlib.pyplot as plt
    from termcolor import colored

    # plt.close('all')
    print(colored('Running events ' + str(eq_num1) + ' and ' + str(eq_num2) + ' phase ' +  dphase, 'green'))
    #%% Workflow selection
    do_3  = True  # pair of events
    do_5  = True
    do_6  = True
    do_7  = True

    # eq_file1 = 'event1.txt'  # pair
    # eq_file2 = 'event2.txt'
    # eq_num1 = '633'  # pair
    # eq_num2 = '649'

    do_3a = False # single event
    do_5a = False # stack
    do_6a = False # treats single events, no time shifts calculated or plotted
    do_7a = False
    eq_num  = '161'  # singlet

    #%% Common parameters
    # ARRAY      = 6
    auto_dist = True
    min_dist = 0
    max_dist = 180

    # HF
    # freq_min = 1
    # freq_max = 2

    # Pro5 stacking
    # stat_corr = 0         # 0 no corr, 1 775 by Wei, 2 219 from SSI, 3 300 from Kawakatsu
    rel_time  = 1         # time alignment in shift
    decimate_fac     = 40 # set for pro5stack2d for single event envelopes, set to 0 for other codes
    simple_taper     =  1
    max_taper_length =  5 # taper is minimum of taper_frac (0.05) and this number of seconds
    apply_SNR        =  False
    SNR_thres        =  1.2
    corr_threshold   =  0.5
    precursor_shift  = -2
    signal_dur       =  4

    phase1 = 'PP'
    phase2 = 'PP'
    phase3 = 'PP'
    phase4 = 'PP'
    phase1 = 'PKiKP'
    phase2 = 'PKIKP'
    phase3 = 'PKP'
    phase4 = 'pPKIKP'

    # Window
    zoom = True                                     # to restrict time range and slowness range in pro7_pair_scan

    # time
    # wind_buff   = 30                                # buffer before time window of analysis
    # Zstart_buff =  0                                # analysis window start relative to phase arrival
    # wind_len    = 20                                # analysis window length
    Zend_buff   = Zstart_buff + wind_len            # analysis window end
    start_buff  = Zstart_buff - wind_buff            # trace start relative to analysis window (if Zoom is True)
    end_buff    = Zstart_buff + wind_len + wind_buff # trace end

    # slowness
    # beam_offset = 0.01
    # beam_width  = 0.02
    # slow_delta  = 0.002
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
    cc_twin      =  2      # time window for cross-correlation (s)
    cc_len       =  0.05   # max time window shift to compute CC (fraction of cc_twin time window)
    cc_delta     =  0.25   # temporal frequency (time spacing) of amp and shift estimates (s)
    cc_interp1d  =  5      # interpolation factor
    cc_thres     =  0.7    # threshold beam correlation to use in stack

    # Pro 7 range selection options
    start_beam = 0  # Limit time window for summary slowness beam in beam sums
    end_beam   = 0  # better be within Zstart and Zend, if zoom is set
    min_amp    = 0.0    # threshold amp to use in stack

    # Pro 7 auto_slice == True options
    auto_slice      = False  # slices span wide range of R and T slownesses
    two_slice_plots = True  # makes R-T pair and snap through time span
    beam_sums       = True  # sums tdiff and amp over time
    wiggly_plots    = False  # shows wiggly plots

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
    do_R = False       # present R plots
    tdiff_plots_too = False  # only applies to auto_slice - to speed plots of only amplitude
    log_plot      = False
    # tdiff_clip   =  0.15
    tdiff_clip   =  cc_twin * cc_len * 0.99
    wig_scale_fac = 0.5
    tdiff_scale_fac = 1
    log_plot_range = 2
    plot_scale_fac = 0.015

    #%% Comparing events
    #%% -- Cull seismic section for common stations
    if do_3 == True:
        pro3pair(ARRAY = ARRAY, eq_num1 = eq_num1, eq_num2 = eq_num2, apply_SNR = apply_SNR,
                    rel_time = rel_time, start_buff = start_buff, end_buff = end_buff,
                    freq_min = freq_min, freq_max = freq_max, do_interpolate = do_interpolate,
                    zoom = zoom, Zstart_buff = Zstart_buff, Zend_buff = Zend_buff, flip = flip,
                    SNR_thres = SNR_thres, corr_threshold = corr_threshold, temp_shift_both = temp_shift_both,
                    precursor_shift = precursor_shift, signal_dur = signal_dur, fig_index = fig_index,
                    max_taper_length = max_taper_length, simple_taper = simple_taper,
                    plot_scale_fac = plot_scale_fac, stat_corr = stat_corr, temp_shift2 = temp_shift2,
                    phase1 = phase1, phase2 = phase2, phase3 = phase3, phase4 = phase4,
                    min_dist = min_dist, max_dist = max_dist, auto_dist = auto_dist, ref_loc = 0)

    #%%  -- 2D stacks
    if do_5 == True:
        pro5stack2d(eq_num = eq_num1, slowR_lo = slowR_lo, slowR_hi = slowR_hi, slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta, start_buff = start_buff, end_buff = end_buff, norm = 1, ARRAY = ARRAY, decimate_fac = decimate_fac, NS = NS)

        pro5stack2d(eq_num = eq_num2, slowR_lo = slowR_lo, slowR_hi = slowR_hi, slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,start_buff = start_buff, end_buff = end_buff, norm = 1, ARRAY = ARRAY, decimate_fac = decimate_fac, NS = NS)

    # %% -- Compare pair of 2D stack results to find shift, amp, amp ratio, uses cc rather than instant phase
    if do_6 == True:
        pro6_cc_pair(eq_num1 = eq_num1, eq_num2 = eq_num2, slowR_lo = slowR_lo, slowR_hi = slowR_hi, slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,start_buff = start_buff, end_buff = end_buff, cc_twin = cc_twin, cc_len = cc_len, cc_interp1d = cc_interp1d, cc_delta = cc_delta, cc_thres = cc_thres)

    #%% -- Make a variety of plots
    if do_7 == True:
        pro7_pair_scan(eq_num1 = eq_num1, eq_num2 = eq_num2, wig_scale_fac = wig_scale_fac, pair_name = pair_name,
                    tdiff_scale_fac = tdiff_scale_fac, slowR_lo = slowR_lo, slowR_hi = slowR_hi,
                    slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
                    zoom = zoom, ZslowR_lo = ZslowR_lo, ZslowR_hi = ZslowR_hi, ZslowT_lo = ZslowT_lo,
                    ZslowT_hi = ZslowT_hi, Zstart_buff = Zstart_buff, Zend_buff = Zend_buff,
                    start_buff = start_buff, end_buff = end_buff, do_T = do_T, do_R = do_R, tdiff_clip = tdiff_clip,
                    phase1 = phase1, phase2 = phase2, phase3 = phase3, phase4 = phase4,
                    min_amp = min_amp, cc_thres = cc_thres, cc_twin = cc_twin, cc_len = cc_len,
                    R_slow_plot = R_slow_plot, T_slow_plot = T_slow_plot, freq_min = freq_min, freq_max = freq_max,
                    snaptime = snaptime, snaps = snaps, snap_depth = snap_depth,
                    nR_plots  = nR_plots, nT_plots = nT_plots, slow_incr = slow_incr, NS = NS,
                    ARRAY = ARRAY, auto_slice = auto_slice, two_slice_plots = two_slice_plots, beam_sums = beam_sums,
                    wiggly_plots = wiggly_plots, log_plot = log_plot, log_plot_range = log_plot_range, plot_peak = plot_peak,
                    tdiff_plots_too = tdiff_plots_too, start_beam = start_beam, end_beam = end_beam)
