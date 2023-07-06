#!/usr/bin/env python
# John Vidale 7/2022
# input is pair of traces
# this program tapers, filters, selects range and SNR
# plots against traveltime curves, either raw or reduced against traveltimes

def run_compare_global(repeater = 'NoName', eq_num1 = 401, eq_num2 = 402, start_buff =  0, wind_len    = 20,
                     precursor_shift  = -1000, signal_dur = -1000, ARRAY = 5,
                     freq_min = 1, freq_max = 2, shift_both = 0, tshift = 0,
                     min_dist = 0, max_dist = 180,
                     apply_SNR = False):

    import os
    #%% load functions

    os.environ['PATH'] += os.pathsep + '/usr/local/bin'
    #%% Import functions
    pro_directory = '/Users/vidale/Documents/GitHub/Array_codes/Run_YKA_WRA'
    os.chdir(pro_directory)
    from get_global_stas         import get_global_stas
    from get_test_stas         import get_test_stas

    pro_directory = '/Users/vidale/Documents/GitHub/Array_codes/Process'
    os.chdir(pro_directory)
    from pro3_sort_plot_pair import pro3pair

    import matplotlib.pyplot as plt
    from termcolor import colored

    # plt.close('all')
    print(colored('Running global events ' + str(eq_num1) + ' and ' + str(eq_num2), 'green'))
    #%% Workflow selection
    get_small_1st  = False # get waveform data from tight arrays for 1st event
    get_small_2nd  = False # get waveform data from tight arrays for 2nd event
    get_big_set_1st  = False # get waveform data from subset of best stations for 1st event
    get_big_set_2nd  = False  # get waveform data from subset of best stations for 2nd event
    make_plots  = True  # plot pair of events

    #%% Common parameters
    apply_SNR        =  False
    SNR_thres        =  1.2
    precursor_shift  = -2
    signal_dur       =  4

    phase1 = 'P'
    phase2 = 'pP'
    phase3 = 'sP'
    phase4 = 'PP'
    auto_dist = True
    min_dist = 0
    max_dist = 180

    # time
    end_buff = start_buff + wind_len # trace end

    plot_scale_fac = 3

    #%% Extract seismograms from subset of best stations
    if get_small_1st == True:
        get_test_stas(eq_num = eq_num1, fig_index = 1)
    if get_small_2nd == True:
        get_test_stas(eq_num = eq_num2, fig_index = 2)

    #%% Extract seismograms from all stations
    if get_big_set_1st == True:
        get_global_stas(eq_num = eq_num1, fig_index = 1)
    if get_big_set_2nd == True:
        get_global_stas(eq_num = eq_num2, fig_index = 2)

    #%% Cull seismic section for common stations
    if make_plots == True:
        pro3pair(ARRAY = ARRAY, repeater = repeater, eq_num1 = eq_num1, eq_num2 = eq_num2, stat_corr = 0,
                    rel_time = 3, start_buff = start_buff, end_buff = end_buff,
                    freq_min = freq_min, freq_max = freq_max, do_interpolate = False,
                    apply_SNR = apply_SNR, SNR_thres = SNR_thres,
                    precursor_shift = precursor_shift, signal_dur = signal_dur, fig_index = 3,
                    max_taper_length = 5, simple_taper = True,
                    plot_scale_fac = plot_scale_fac, tshift = tshift, shift_both = shift_both,
                    phase1 = phase1, phase2 = phase2, phase3 = phase3, phase4 = phase4,
                    min_dist = min_dist, max_dist = max_dist, auto_dist = auto_dist, ref_loc = False)
