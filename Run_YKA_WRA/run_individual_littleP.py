#!/usr/bin/env python
# John Vidale 7/2020

def run_individual_littleP(eq_num = 401, start_buff = 980, end_buff = 1180,
    precursor_shift  = -1000, signal_dur     = -1000,
    start_buff_stack =   -30, end_buff_stack =    40,
    start_beam_stack =     0, end_beam_stack =     0,
    freq_min = 1, freq_max = 3, slow_delta = 0.0025,
    min_dist = 0, max_dist = 180, dphase  = 'PKiKP', JST = False, R_slow_plot = 0, T_slow_plot = 0,
    fig_index = 401, stat_corr = 1, apply_SNR = False, shift_tt = 0, zerophase = True):

    import os

    #%% Import functions
    # from pro2_dec                import pro2decimate
    from pro3_sort_plot_singlet  import pro3singlet
    from termcolor import colored

    # plt.close('all')
    print(colored('Running event ' + str(eq_num) + ' phase ' +  dphase, 'green'))

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

    # dphase  = 'PKiKP'
    dphase2 = 'pP'
    dphase3 = 'sP'
    dphase4 = 'PP'

    # decimate_fac = 5
    # decimate, in 100 sps, out 20 sps

    # precursor_shift = 0
    # signal_dur      = 5
    corr_threshold  = 0.0
    # apply_SNR       = True
    SNR_thres       = 1.3

    simple_taper    = True

#%%  Parameters for filtering, sorting, and beaming
    # freq_min     = 1
    # freq_max     = 3

    zoom = 0
    plot_scale_fac = 0.3

    # stat_corr = 3

    dec_fac = 1 # decimation factor

    rel_time = 3   # phase alignment details

#%%  pro3singlet -- selects data
    pro3singlet(ARRAY = ARRAY, stat_corr = stat_corr, shift_tt = shift_tt,
        eq_num = eq_num, simple_taper = simple_taper, rel_time = rel_time,
        start_buff = start_buff_stack, end_buff = end_buff_stack,
        zoom = zoom, Zstart_buff = start_beam_stack, Zend_buff = end_beam_stack,
        start_beam = start_beam_stack, end_beam = end_beam_stack,
        precursor_shift = precursor_shift, signal_dur = signal_dur,
        plot_scale_fac = plot_scale_fac, fig_index = fig_index, apply_SNR = apply_SNR,
        dphase = dphase, dphase2 = dphase2, dphase3 = dphase3, dphase4 = dphase4,
        freq_min = freq_min, freq_max = freq_max, zerophase = zerophase,
        min_dist = min_dist, max_dist = max_dist,
        SNR_thres = SNR_thres, corr_threshold = corr_threshold,
        ref_loc = ref_loc, ref_rad = ref_rad, JST = JST, decimate_fac = dec_fac, )
