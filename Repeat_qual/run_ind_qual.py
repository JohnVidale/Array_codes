#!/usr/bin/env python
# John Vidale 7/2022
# input is pair of traces
# this program tapers, filters, selects range and SNR
# plots against traveltime curves, either raw or reduced against traveltimes

def run_ind_qual(repeater = 'NoName', eq_num1 = 401, eq_num2 = 402, start_buff =  0, wind_len    = 20,
                     precursor_shift  = -1000, signal_dur = -1000, ARRAY = 5,
                     freq_min = 1, freq_max = 2, temp_shift_both = 0, temp_shift2 = 0,
                     min_dist = 0, max_dist = 180, dphase  = 'PKiKP',
                     apply_SNR = False):

    import os
    #%% load functions

    os.environ['PATH'] += os.pathsep + '/usr/local/bin'
    #%% Import functions
    pro_directory = '/Users/vidale/Documents/GitHub/Array_codes/Repeat_qual'
    os.chdir(pro_directory)
    from get_13_stas         import get_13_stas

    pro_directory = '/Users/vidale/Documents/GitHub/Array_codes/Process'
    os.chdir(pro_directory)
    from pro3_sort_plot_pair import pro3pair

    import matplotlib.pyplot as plt
    from termcolor import colored

    plt.close('all')
    print(colored('Running events ' + str(eq_num1) + ' and ' + str(eq_num2) + ' phase ' +  dphase, 'green'))
    #%% Workflow selection
    do_2  = False  # get waveform data
    do_3  = True  # plot pair of events

    #%% Common parameters
    apply_SNR        =  False
    SNR_thres        =  1.2
    precursor_shift  = -2
    signal_dur       =  4

    phase1 = 'P'
    phase2 = 'pP'
    phase3 = 'sP'
    phase4 = 'PP'
    # phase1 = 'PKIKP'
    # phase2 = 'PKiKP'
    # phase3 = 'PKP'
    # phase4 = 'pPKIKP'
    auto_dist = True
    min_dist = 0
    max_dist = 180

    # time
    end_buff = start_buff + wind_len # trace end

    plot_scale_fac = 3

    #%% Extract seismograms
    if do_2 == True:
        get_13_stas(eq_num = eq_num1, fig_index = 1)
        get_13_stas(eq_num = eq_num2, fig_index = 2)

    #%% Cull seismic section for common stations
    if do_3 == True:
        pro3pair(ARRAY = ARRAY, repeater = repeater, eq_num1 = eq_num1, eq_num2 = eq_num2, stat_corr = 0,
                    rel_time = 3, start_buff = start_buff, end_buff = end_buff,
                    freq_min = freq_min, freq_max = freq_max, do_interpolate = False,
                    apply_SNR = apply_SNR, SNR_thres = SNR_thres,
                    precursor_shift = precursor_shift, signal_dur = signal_dur, fig_index = 3,
                    max_taper_length = 5, simple_taper = True,
                    plot_scale_fac = plot_scale_fac, temp_shift2 = temp_shift2, temp_shift_both = temp_shift_both,
                    phase1 = phase1, phase2 = phase2, phase3 = phase3, phase4 = phase4,
                    min_dist = min_dist, max_dist = max_dist, auto_dist = auto_dist, ref_loc = False)
