#!/usr/bin/env python
# John Vidale 7/2020

def run_get_shift_PKiKP(start_buff = 980, end_buff = 1180, event_no = 1, min_dist = 0,
              max_dist = 180, freq_min = 1, freq_max = 3, slow_delta = 0.0025,
              start_SNRbeam = 0, end_SNRbeam = 0, start_beam = 0, end_beam = 0, dphase  = 'PKiKP'):

    import os
    os.environ['PATH'] += os.pathsep + '/usr/local/bin'
    os.chdir('/Users/vidale/Documents/GitHub/Array_codes')

    #%% Import functions
    from pro2_dec                import pro2decimate
    from pro3b_sort_plot_singlet import pro3singlet
    from pro5a_stack             import pro5stack
    from pro5b_stack2d           import pro5stack2d
    from pro6_plot_singlet       import pro6stacked_singlet
    from pro7a_plot_envstack     import pro7plotstack
    from pro7b_plot_stack        import pro7plotstack2
    from pro7b_dec               import pro7dec
    from pro4_get_shifts_Hinet   import pro4statics
    import matplotlib.pyplot as plt

    #%% Common parameters
    ev_directory = '/Users/vidale/Documents/PyCode/Hinet/Tian_events'
    os.chdir(ev_directory)

    eq_file   = 'event' + str(event_no) + '.txt'
    ARRAY     = 0

    ref_loc = 1   # 0 select stations by distance from epicenter, 1 select stations by distance from ref location
    ref_rad = 0.25 # radius of stations around ref_loc chosen
    ref_lat = 36
    ref_lon = 138

    auto_dist = 1  #  automatically plot only real distance range
    min_dist = 18.6
    max_dist = 19.2

    # dphase  = 'PKiKP'
    dphase2 = 'PcP'
    dphase3 = 'PKiKP'
    dphase4 = 'pPcP'

    # decimate_fac = 5
    # decimate, in 100 sps, out 20 sps
    # pro2decimate(eq_file, decimate_fac = decimate_fac)

#%% Get static corrections
    freq_min     = 0.7
    freq_max     = 5.0
    start_corr_win = 4
    end_corr_win = 2  # before pick is POSITIVE
    corr_threshold = 0.6
    skip_SNR       = 1
    qual_threshold = 3
    simple_taper   =  1

    pro3singlet(ARRAY = ARRAY, stat_corr = 1, eq_file = eq_file, simple_taper = 1,
        rel_time = 3, start_buff = start_buff, end_buff = end_buff,
        start_beam = start_SNRbeam, end_beam = end_SNRbeam,
        plot_scale_fac = 0.03, skip_SNR = skip_SNR, event_no = event_no,
        dphase = dphase, dphase2 = dphase2, dphase3 = dphase3, dphase4 = dphase4,
        freq_min = freq_min, freq_max = freq_max,
        min_dist = min_dist, max_dist = max_dist, auto_dist = auto_dist,
        qual_threshold = qual_threshold, corr_threshold = corr_threshold,
        ref_loc = ref_loc, ref_lat = ref_lat, ref_lon = ref_lon, ref_rad = ref_rad,
        fig_index = 102, JST = 1)

    pro4statics(eq_file, use_ref_trace = 0, ref_trace = 'nothing',
        dphase = dphase, dphase2 = dphase2, dphase3 = dphase3, dphase4 = dphase4,
        start_corr_win = start_corr_win, end_corr_win = end_corr_win, plot_scale_fac = 0.05,
        qual_threshold = 0, corr_threshold = corr_threshold,
        max_time_shift = 0.25, min_dist = min_dist, max_dist = max_dist, ARRAY = 0)

    #%% Filtering, sorting, and beaming
    freq_min     = 0.7
    freq_max     = 5.0

    slowR_lo   = -0.1
    slowR_hi   =  0.1
    slowT_lo   = -0.1
    slowT_hi   =  0.1
    slow_delta =  0.005
    NS         =  1   # 0 plot slowness R-T, 1 plot slowness N-S
    plot_scale_fac = 0.03

    snaptime  = 10

    stat_corr      = 1

    slowR_lo_1D = -0.1
    slowR_hi_1D =  0.1
    slow_delta_1D = 0.0005
    dec_fac = 10

    rel_time = 1   # phase alignment details
    # rel_time == 0  window in absolute time after origin time
    # rel_time == 1  each window has a shift proportional to (dist - ref_dist) at phase slowness at ref_dist
    # rel_time == 2  each window has a distinct phase-chose shift, but time offset is common to all stations
    # rel_time == 3  each station has an individual, chosen-phase shift, phase arrival set to common time
    # rel_time == 4  use same window around chosen phase for all stations, using ref distance

    pro3singlet(ARRAY = ARRAY, stat_corr = stat_corr, fine_stats = 1, eq_file = eq_file, simple_taper = simple_taper,
        rel_time = rel_time, start_buff = start_buff, end_buff = end_buff,
        start_beam = start_SNRbeam, end_beam = end_SNRbeam,
        plot_scale_fac = plot_scale_fac, skip_SNR = skip_SNR, event_no = event_no,
        dphase = dphase, dphase2 = dphase2, dphase3 = dphase3, dphase4 = dphase4,
        freq_min = freq_min, freq_max = freq_max,
        min_dist = min_dist, max_dist = max_dist, auto_dist = auto_dist,
        qual_threshold = qual_threshold, corr_threshold = corr_threshold,
        ref_loc = ref_loc, ref_lat = ref_lat, ref_lon = ref_lon, ref_rad = ref_rad,
        fig_index = 102, JST = 1)

    # pro5stack(ARRAY = ARRAY, eq_file = eq_file, plot_scale_fac = 0.05,
    #     slowR_lo = slowR_lo_1D, slowR_hi = slowR_hi_1D, slow_delta = slow_delta_1D,
    #     start_buff = start_buff, end_buff = end_buff,
    #     ref_loc = ref_loc, ref_lat = ref_lat, ref_lon = ref_lon,
    #     log_plot = 0, envelope = 1, plot_dyn_range = 50,
    #     norm = 1, global_norm_plot = 1, color_plot = 1, fig_index = 302)

    #%%  --2D stacks
    # pro5stack2d(eq_file = eq_file, plot_scale_fac = 0.05,
    #     slowR_lo = slowR_lo, slowR_hi = slowR_hi, slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
    #     	start_buff = start_buff, end_buff = end_buff, norm = 1,
    #     ref_loc = ref_loc, ref_lat = ref_lat, ref_lon = ref_lon, ARRAY = ARRAY, decimate_fac = dec_fac, NS = NS)

    #%% --Sum beams over time windows
    # pro6stacked_singlet(eq_file = eq_file, plot_scale_fac = 0.003,
    #     slowR_lo = slowR_lo, slowR_hi = slowR_hi, slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
    #     start_buff = start_buff, end_buff = end_buff, R_slow_plot = 0, T_slow_plot = 0, dphase = dphase,
    #     fig_index = 301, plot_dyn_range = 100, ARRAY = ARRAY, event_no = event_no, start_beam = start_beam, end_beam = end_beam)

    #%% --2D envelop stack results for individual events
    # pro7plotstack(eq_file = eq_file, plot_scale_fac = 0.05,
    #     slowR_lo = slowR_lo, slowR_hi = slowR_hi, slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
    #     start_buff = start_buff, end_buff = end_buff, skip_T = 1, skip_R = 0,
    #     zoom = 0, ZslowR_lo = -0.03, ZslowR_hi = 0.03, ZslowT_lo = -0.03, ZslowT_hi = 0.03, Zstart_buff = 0, Zend_buff = 200,
    #     fig_index = 401, plot_dyn_range = 50, snaptime = snaptime, snaps=1, ARRAY = ARRAY)
