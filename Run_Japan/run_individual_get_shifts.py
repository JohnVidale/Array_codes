#!/usr/bin/env python
# John Vidale 7/2020

def run_individual_get_shifts(eq_num = 1, min_dist = 0,
                    max_dist = 180, freq_min = 1, freq_max = 2.5,
                    start_beam = 0, end_beam = 0, start_buff = 0, end_buff = 0,
                    precursor_shift = 0, signal_dur = 0, dphase  = 'PKiKP'):

    import os

    #%% Import functions
    from pro3_sort_plot_singlet import pro3singlet
    from pro4_get_shifts        import pro4_get_shifts

#%%  Parameters - basic
    ev_directory = '/Users/vidale/Documents/Research/IC/EvLocs'
    os.chdir(ev_directory)

    ARRAY = 0 # 0 is Hi-Net
    JST   = False # UTC or JST

    ref_loc = False   # 0 select stations by distance from epicenter, 1 select stations by distance from ref location
    ref_rad = 180 # radius of stations around ref_loc chosen
    ref_lat = 38  # middle
    ref_lon = 140
    # ref_lat = 34  # south
    # ref_lon = 133
    # ref_lat = 42.5  # north
    # ref_lon = 143.5

    plot_auto_dist = True #  automatically plot only real distance range
    # min_dist = 12
    # max_dist = 25

    dphase  = 'PKiKP'
    dphase2 = 'PKiKP'
    dphase3 = 'PKIKP'
    dphase4 = 'PKP'

#%% Parameters for static calculation
    freq_min     = 2
    freq_max     = 8
    # start_beam_align = 4  # times before pick are now normal, negative
    # end_beam_align   = 2
    corr_threshold = 0.5
    max_time_shift = 1.2
    apply_SNR      = False
    SNR_thres      = 1.5
    simple_taper   = True
    stat_corr = 0

#%%  pro3singlet -- selects data for statics
    pro3singlet(ARRAY = ARRAY, stat_corr = stat_corr, simple_taper = simple_taper, rel_time = 3,
        start_buff = start_buff, end_buff = end_buff,
        precursor_shift = precursor_shift, signal_dur = signal_dur,
        plot_scale_fac = 0.05, apply_SNR = apply_SNR, eq_num = eq_num,
        dphase = dphase, dphase2 = dphase2, dphase3 = dphase3, dphase4 = dphase4,
        freq_min = freq_min, freq_max = freq_max,
        min_dist = min_dist, max_dist = max_dist, plot_auto_dist = plot_auto_dist,
        SNR_thres = SNR_thres, corr_threshold = corr_threshold,
        ref_loc = ref_loc, ref_rad = ref_rad,
        fig_index = 102, JST = JST)

#%%  pro4_get_shifts -- measures statics
    # pro4_get_shifts(eq_num = eq_num, use_ref_trace = True, ref_trace = 'N.TOW',
    #     dphase = dphase, dphase2 = dphase2, dphase3 = dphase3, dphase4 = dphase4,
    #     start_beam = precursor_shift, end_beam = signal_dur + precursor_shift, start_buff = start_buff, end_buff = end_buff,
    #     plot_scale_fac = 0.05, corr_threshold = corr_threshold,
    #     max_time_shift = max_time_shift, min_dist = min_dist, max_dist = max_dist, ARRAY = ARRAY)
