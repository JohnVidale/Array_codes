#!/usr/bin/env python
# John Vidale 7/2020

def run_individual_align(start_buff = 980, end_buff = 1180, event_no = 1, min_dist = 0,
                    max_dist = 180, freq_min = 1, freq_max = 3, slow_delta = 0.0025,
                    start_beam_align = 0, end_beam_align = 0, start_buff_align = 0, end_buff_align = 0,
                    start_beam_stack = 0, end_beam_stack = 0, start_buff_stack = 0, end_buff_stack = 0,
                    dphase  = 'PKiKP'):

    import os

    #%% Import functions
    from pro3_sort_plot_singlet import pro3singlet
    from pro4_get_shifts        import pro4statics

#%%  Parameters - basic
    ev_directory = '/Users/vidale/Documents/Research/IC/EvLocs'
    os.chdir(ev_directory)

    eq_file   = 'event' + str(event_no) + '.txt'
    ARRAY     = 0

    ref_loc = 1   # 0 select stations by distance from epicenter, 1 select stations by distance from ref location
    ref_rad = 3 # radius of stations around ref_loc chosen
    ref_lat = 38  # middle
    ref_lon = 140
    # ref_lat = 34  # south
    # ref_lon = 133
    # ref_lat = 42.5  # north
    # ref_lon = 143.5

    auto_dist = 1  #  automatically plot only real distance range
    # min_dist = 12
    # max_dist = 25

    # dphase  = 'PKiKP'
    dphase2 = 'PcP'
    dphase3 = 'PKiKP'
    dphase4 = 'pPcP'

#%% Parameters for static calculation
    freq_min     = 0.5
    freq_max     = 2
    # start_beam_align = 4  # times before pick are now normal, negative
    # end_beam_align   = 2
    corr_threshold = 0.4
    max_time_shift = 1.5
    skip_SNR       = 1
    qual_threshold = 1.3
    simple_taper   = 1

#%%  pro3singlet -- selects data for statics
    pro3singlet(ARRAY = ARRAY, stat_corr = 1, eq_file = eq_file, simple_taper = 1, rel_time = 3,
        start_buff = start_buff_align, end_buff = end_buff_align,
        start_beam = start_beam_align, end_beam = end_beam_align,
        plot_scale_fac = 0.03, skip_SNR = skip_SNR, event_no = event_no,
        dphase = dphase, dphase2 = dphase2, dphase3 = dphase3, dphase4 = dphase4,
        freq_min = freq_min, freq_max = freq_max,
        min_dist = min_dist, max_dist = max_dist, auto_dist = auto_dist,
        qual_threshold = qual_threshold, corr_threshold = corr_threshold,
        ref_loc = ref_loc, ref_lat = ref_lat, ref_lon = ref_lon, ref_rad = ref_rad,
        fig_index = 102, JST = 1)

#%%  pro4statics -- measures statics
    pro4statics(eq_file, use_ref_trace = 0, ref_trace = 'N.SUY',
        dphase = dphase, dphase2 = dphase2, dphase3 = dphase3, dphase4 = dphase4,
        start_beam = start_beam_align, end_beam = end_beam_align,
        start_buff = start_buff_align, end_buff = end_buff_align,
        plot_scale_fac = 0.05, event_no = event_no,
        qual_threshold = 0, corr_threshold = corr_threshold,
        max_time_shift = max_time_shift, min_dist = min_dist, max_dist = max_dist, ARRAY = ARRAY)
