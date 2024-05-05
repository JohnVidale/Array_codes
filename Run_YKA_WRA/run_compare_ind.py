#!/usr/bin/env python3
# John Vidale 4/2020
# last modified 8/2022

def run_compare_ind(repeater = 'NoName',do_global = False, do_YKA = False, do_ILAR = False):

    import os
    import sys
    import time
    import matplotlib.pyplot as plt
    from termcolor import colored
    import pandas as pd

    #%% Import functions
    # pro_directory = '/Users/vidale/Documents/GitHub/Array_codes/Run_YKA_WRA'
    # os.chdir(pro_directory)
    # # Print the current working directory (CWD)
    # cwd = os.getcwd()
    # print("Run_compare current working directory: ", cwd)

    from run_compare_global     import run_compare_global
    from run_individual_df      import run_individual_df
    from run_compare_pair       import run_compare_pair

    start_time_wc = time.time()

    plt.close('all')

#
#%% Parameters

    def search_df(df, column, value, partial_match=True):
        df = df.astype({column:'string'})
        if partial_match:
            return df.loc[df[column].str.contains(value, na=False)]
        else:
            return df.loc[df[column] == value]

    # look up pair of earthquakes and time shifts in pairs
    df = pd.read_excel('/Users/vidale/Documents/GitHub/Array_codes/Files/ICevents_full.xlsx', sheet_name='pairs')
    lines0       = search_df(df,'label'      ,repeater,partial_match=True)
    eq_num1      = lines0.index1.iloc[0]
    eq_num2      = lines0.index2.iloc[0]
    tshift       = lines0.tshift.iloc[0]
    shift_both   = lines0.shift_both.iloc[0]
    Y_shift      = lines0.Y_shift.iloc[0]
    shift_bothY  = lines0.shift_bothY.iloc[0]
    do_YKA       = True
    do_ILAR      = lines0.ILAR.iloc[0]
    do_global    = lines0.global_sta.iloc[0]

    # wretched Excel quirk
    if do_ILAR   == 'TRUE':  do_ILAR = 'True'
    if do_ILAR   == 'FALSE': do_ILAR = 'False'
    if do_global == 'TRUE':  do_global = 'True'
    if do_global == 'FALSE': do_global = 'False'

    # freq_min = 0.6; freq_max = 1.5
    freq_min = 1; freq_max = 2
    # freq_min = 2; freq_max = 4
    do_ILAR_pre = do_ILAR

    # Skip some arrays?
    # do_YKA      = False
    # do_ILAR     = False
    # do_ILAR_pre = False
    # do_global   = False

    print(colored('do_global ' + str(do_global) + ' and do_ILAR ' + str(do_ILAR), 'green'))

    if do_global:
        start_buff = -10 # analysis window start relative to phase arrival
        wind_len    = 30 # analysis window length
        run_compare_global(repeater = repeater, freq_min = freq_min, freq_max = freq_max, ARRAY = 7,
                start_buff = start_buff, wind_len = wind_len)

    wind_buff = 30 # buffer before and after time window of analysis
    plot_peak = 1

    beam_offset = 0.02
    beam_width  = 0.010
    slow_delta  = 0.002

#%% YKA PKIKP
    if do_YKA:
        Zstart_buff = -20 # analysis window start relative to phase arrival
        wind_len    =  40 # analysis window length
        plot_peak = 1.0
        run_compare_pair(repeater = repeater, dphase = 'PKIKP',
                beam_width = beam_width, slow_delta = slow_delta, beam_offset = beam_offset,
                freq_min = freq_min, freq_max = freq_max, stat_corr = 0, ARRAY = 5,
                Zstart_buff = Zstart_buff, wind_len = wind_len, wind_buff = wind_buff,
                fig_index = 200, do_interpolate =  True, pair_name = repeater, plot_peak = plot_peak)

#%% YKA PKIKP - long window
    # if do_YKA:
    #     Zstart_buff = -20 # analysis window start relative to phase arrival
    #     wind_len    =  40 # analysis window length
    #     plot_peak = 1
    #     run_compare_pair(repeater = repeater, dphase = 'PKIKP',
    #             beam_width = beam_width, slow_delta = slow_delta, beam_offset = beam_offset,
    #             freq_min = freq_min, freq_max = freq_max, stat_corr = 0, ARRAY = 5,
    #             Zstart_buff = Zstart_buff, wind_len = wind_len, wind_buff = wind_buff,
    #             tshift = tshift + Y_shift, fig_index = 500, do_interpolate =  True, pair_name = repeater, plot_peak = plot_peak)

#%% ILAR PKP
    if do_ILAR:
        Zstart_buff = -10 # analysis window start relative to phase arrival
        wind_len    =  30 # analysis window length
        plot_peak = 1.0
        run_compare_pair(repeater = repeater, dphase = 'PKIKP',
                beam_width = beam_width, slow_delta = slow_delta, beam_offset = beam_offset,
                freq_min = freq_min, freq_max = freq_max, stat_corr = 0, ARRAY = 6,
                Zstart_buff = Zstart_buff, wind_len = wind_len, wind_buff = wind_buff,
                fig_index = 300, do_interpolate =  True, pair_name = repeater, plot_peak = plot_peak)

#%% ILAR PKIKP
    win_norm =  True
    trace_norm = False
    if do_ILAR_pre:
        Zstart_buff = -15 # analysis window start relative to phase arrival
        wind_len    =  16 # analysis window length
        plot_peak = 0.05
        run_compare_pair(repeater = repeater, dphase = 'PKIKP',
                beam_width = beam_width, slow_delta = slow_delta, beam_offset = beam_offset, win_norm = win_norm,trace_norm = trace_norm,
                freq_min = freq_min, freq_max = freq_max, stat_corr = 0, ARRAY = 6, trace_amp = 0.5,
                Zstart_buff = Zstart_buff, wind_len = wind_len, wind_buff = wind_buff,
                fig_index = 400, do_interpolate =  True, pair_name = repeater, plot_peak = plot_peak)

    elapsed_time_wc = time.time() - start_time_wc
    print(f'This job took {elapsed_time_wc:.1f} seconds')
    os.system('say "All done"')
