#!/usr/bin/env python
# John Vidale 7/2020

def run_individual_little(eq_num = 401, start_buff = 980, end_buff = 1180,
    precursor_shift  = -1000, signal_dur     = -1000,
    start_buff_stack =   -30, end_buff_stack =    40,
    start_beam_stack =     0, end_beam_stack =     0,
    beam_width = 0.04, beam_step = 0.005, beam_offset = 0.00, freq_min = 1, freq_max = 3, slow_delta = 0.0025,
    min_dist = 0, max_dist = 180, dphase  = 'PKiKP', JST = False, R_slow_plot = 0, T_slow_plot = 0,
    fig_index = 401, stat_corr = 1, apply_SNR = False):

    import os
    #%% close plots

    #%% Import functions
    # from pro2_dec                import pro2decimate
    from pro3_sort_plot_singlet  import pro3singlet
    from pro4_get_shifts         import pro4_get_shifts
    from pro5_stack1d            import pro5stack1d
    from pro5_stack2d            import pro5stack2d
    from pro6_singlet            import pro6_singlet
    from pro7_singlet            import pro7_singlet
    import matplotlib.pyplot as plt
    from termcolor import colored

    # plt.close('all')
    print(colored('Running event ' + str(eq_num) + ' phase ' +  dphase, 'green'))
    #%% Workflow selection
    do_3  = True  # pair of events
    do_5  = True
    do_6  = True
    do_7  = True

#%%  Parameters - basic
    ev_directory = '/Users/vidale/Documents/Research/IC/EvLocs'
    os.chdir(ev_directory)

    if eq_num < 100 or (eq_num < 400 and eq_num >= 300):
        ARRAY = 1 # LASA
    elif (eq_num >= 100 and eq_num < 200):
        ARRAY = 0 # HiNet
    elif (eq_num >= 200 and eq_num < 300):
        ARRAY = 2 # China
    elif (eq_num >= 400 and eq_num < 500):
        ARRAY = 4 # WRA
    elif (eq_num >= 500 and eq_num < 600):
        ARRAY = 5 # YKA

    ref_loc = False   # Override default array center with selection
    ref_rad = 180     # selected radius of stations around ref_loc
    ref_lat = 0       # ref_loc choice if ref_loc == True
    ref_lon = 0
    # ref_lat = 38  # middle 1
    # ref_lon = 140
    # ref_lat = 34  # middle 2
    # ref_lon = 133
    # ref_rad = 1 # selected radius of stations around ref_loc
    # ref_lat = 34  # south
    # ref_lon = 133
    # ref_lat = 42.5  # north
    # ref_lon = 143.5

    # min_dist = 12
    # max_dist = 25

    # dphase  = 'PKiKP'
    dphase2 = 'pPKiKP'
    dphase3 = 'sPKiKP'
    dphase4 = 'S'

    # decimate_fac = 5
    # decimate, in 100 sps, out 20 sps
    # eq_file   = 'event' + str(event_no) + '.txt'
    # pro2decimate(eq_file, decimate_fac = decimate_fac)

    # start_beam_align = 4  # times before pick are now normal, negative
    # end_beam_align   = 2
    # precursor_shift = 0
    # signal_dur      = 5
    corr_threshold  = 0.0
    # apply_SNR       = True
    SNR_thres       = 1.3

    simple_taper    = True
    snaptime        = 0
    snaps           = 0

#%%  Parameters for filtering, sorting, and beaming
    # freq_min     = 1
    # freq_max     = 3

    # beam_offset = 0.0
    # beam_width = 0.02
    zoom = 1
    slowR_lo   = -beam_width + beam_offset
    slowR_hi   =  beam_width + beam_offset
    slowT_lo   = -beam_width
    slowT_hi   =  beam_width
    ZslowR_lo  = -beam_width + beam_offset
    ZslowR_hi  =  beam_width + beam_offset
    ZslowT_lo  = -beam_width
    ZslowT_hi  =  beam_width
    slow_delta =  beam_step
    NS         =  True   # 0 plot slowness R-T, 1 plot slowness N-S
    plot_scale_fac = 0.02
    log_plot = 0

    # stat_corr = 3

    slowR_lo_1D = -0.04
    slowR_hi_1D =  0.04
    slow_delta_1D = 0.0005
    dec_fac = 1 # decimation factor

    rel_time = 1   # phase alignment details
    # rel_time == 0  no shift - window in absolute time after origin time
    # rel_time == 1  chosen arrival is flattened, then local slowness of chosen arrival is added back to moveout
    # rel_time == 2  each window has a distinct phase-chosen shift, but time offset is common to all stations
    # rel_time == 3  each station has an individual, chosen-phase shift, phase arrival set to common time
    # rel_time == 4  use same window around chosen phase for all stations, baseline set at ref distance

#%%  pro3singlet -- selects data
    if do_3 == True:
        pro3singlet(ARRAY = ARRAY, stat_corr = stat_corr,
            eq_num = eq_num, simple_taper = simple_taper, rel_time = rel_time,
            start_buff = start_buff_stack, end_buff = end_buff_stack,
            zoom = zoom, Zstart_buff = start_beam_stack, Zend_buff = end_beam_stack,
            start_beam = start_beam_stack, end_beam = end_beam_stack,
            precursor_shift = precursor_shift, signal_dur = signal_dur,
            plot_scale_fac = plot_scale_fac, fig_index = fig_index, apply_SNR = apply_SNR,
            dphase = dphase, dphase2 = dphase2, dphase3 = dphase3, dphase4 = dphase4,
            freq_min = freq_min, freq_max = freq_max,
            min_dist = min_dist, max_dist = max_dist,
            SNR_thres = SNR_thres, corr_threshold = corr_threshold,
            ref_loc = ref_loc, ref_rad = ref_rad, JST = JST, decimate_fac = dec_fac, )

#%%  pro5stack -- 1D stack mostly good for long time window surveys
    # pro5stack1d(ARRAY = ARRAY, eq_num = event_no, plot_scale_fac = 0.05,
    #     slowR_lo = slowR_lo_1D, slowR_hi = slowR_hi_1D, slow_delta = slow_delta_1D,
    #     start_buff = start_buff_stack, end_buff = end_buff_stack,
    #     ref_loc = ref_loc, ref_lat = ref_lat, ref_lon = ref_lon,
    #     log_plot = log_plot, envelope = 1, plot_dyn_range = 50,
    #     norm = 1, global_norm_plot = 1, color_plot = 1, fig_index = 302)

#%%  pro5stack2d -- 2D stacks, generates no plots
    if do_5 == True:
        pro5stack2d(eq_num = eq_num, slowR_lo = slowR_lo, slowR_hi = slowR_hi,
            slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
            start_buff = start_buff_stack, end_buff = end_buff_stack,
            norm = True, min_dist = min_dist, max_dist = max_dist,
            ref_loc = ref_loc, ref_lat = ref_lat, ref_lon = ref_lon, ARRAY = ARRAY,
            decimate_fac = dec_fac, NS = NS)

#%%  pro6stacked_singlet -- summed beams over 2D stacks
    if do_6 == True:
        pro6_singlet(eq_num = eq_num,
            slowR_lo = slowR_lo, slowR_hi = slowR_hi,
            slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
            start_buff = start_buff_stack, end_buff = end_buff_stack)

#%%  pro7plotstack --  snapshots of 2D stack results for individual events
    if do_7 == True:
        pro7_singlet(eq_num = eq_num,
            slowR_lo = slowR_lo, slowR_hi = slowR_hi, slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
            ref_phase = dphase, start_buff = start_buff_stack, end_buff = end_buff_stack,
            zoom = False, ZslowR_lo = ZslowR_lo, ZslowR_hi = ZslowR_hi, ZslowT_lo = ZslowT_lo, ZslowT_hi = ZslowT_hi,
            Zstart_buff = start_beam_stack, Zend_buff = end_beam_stack,
            two_slice_plots = True, wiggly_plots = False, freq_min = freq_min, freq_max = freq_max,
            fig_index = fig_index + 10, snaptime = snaptime, snaps=snaps, ARRAY = ARRAY, NS = NS,
            ref_loc = ref_loc, ref_lat = ref_lat, ref_lon = ref_lon, R_slow_plot = R_slow_plot, T_slow_plot = T_slow_plot)

#%%  pro7plotstack --  snapshots of 2D stack results for individual events
    if do_7 == True:
        pro7_singlet(eq_num = eq_num,
            slowR_lo = slowR_lo, slowR_hi = slowR_hi, slowT_lo = slowT_lo, slowT_hi = slowT_hi, slow_delta = slow_delta,
            ref_phase = dphase, start_buff = start_buff_stack, end_buff = end_buff_stack,
            zoom = True, ZslowR_lo = ZslowR_lo, ZslowR_hi = ZslowR_hi, ZslowT_lo = ZslowT_lo, ZslowT_hi = ZslowT_hi,
            Zstart_buff = start_beam_stack, Zend_buff = end_beam_stack,
            two_slice_plots = True, wiggly_plots = False,
            fig_index = fig_index + 20, snaptime = snaptime, snaps=snaps, ARRAY = ARRAY, NS = NS,
            ref_loc = ref_loc, ref_lat = ref_lat, ref_lon = ref_lon, R_slow_plot = R_slow_plot, T_slow_plot = T_slow_plot)
