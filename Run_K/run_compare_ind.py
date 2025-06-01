#!/usr/bin/env python3
# John Vidale 4/2020
# last modified 8/2022, modified for KK 4/2025

def run_compare_ind(repeater = 'NoName',do_global = False, do_KK = False, do_KUR = False):

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
    do_global    = lines0.global_sta.iloc[0] # NOT WORKING
    print(f'eq_num1 {eq_num1} eq_num2 {eq_num2} tshift {tshift} do_global {do_global}')

    # wretched Excel quirk
    if   do_global == 'TRUE':  do_global = 'True'
    elif do_global == 'FALSE': do_global = 'False'

    freq_min = 2; freq_max = 5

    print(colored('do_global ' + str(do_global) + ' and do_KK ' + str(do_KK) + ' and do_KUR ' + str(do_KUR), 'green'))

    if do_global:
        start_buff = -10 # analysis window start relative to phase arrival
        wind_len    = 30 # analysis window length
        run_compare_global(repeater = repeater, freq_min = freq_min, freq_max = freq_max, ARRAY = 7,
                start_buff = start_buff, wind_len = wind_len)

    wind_buff = 30 # buffer before and after time window of analysis
    plot_peak = 1

    beam_offset = 0.02 # tuned for speed and resolution
    beam_width  = 0.010
    slow_delta  = 0.002

    # beam_offsetW = 0.02 # tuned for wide angle view
    # beam_widthW  = 0.050
    # slow_deltaW  = 0.003

#%% Array KUR
    if do_KUR:
        Zstart_buff = -80 # analysis window start relative to phase arrival
        wind_len    = 200 # analysis window length
        plot_peak = 1.0
        run_compare_pair(repeater = repeater, dphase = 'PcP',
                beam_width = beam_width, slow_delta = slow_delta, beam_offset = beam_offset,
                freq_min = freq_min, freq_max = freq_max, stat_corr = 0, ARRAY = 11, win_norm = True,
                Zstart_buff = Zstart_buff, wind_len = wind_len, wind_buff = wind_buff,
                fig_index = 200, do_interpolate =  True, pair_name = repeater, plot_peak = plot_peak)

#%% Array KK
    if do_KK:
        Zstart_buff =  -80 # analysis window start relative to phase arrival
        wind_len    =  200 # analysis window length
        run_compare_pair(ARRAY=12, fig_index=500, repeater=repeater, dphase='PcP',
                beam_width=beam_width, slow_delta=slow_delta, beam_offset=beam_offset,
                freq_min=freq_min, freq_max=freq_max, stat_corr=0, win_norm = True, 
                Zstart_buff=Zstart_buff, wind_len=wind_len, wind_buff=wind_buff,
                do_interpolate=True, pair_name=repeater, plot_peak=plot_peak)

    elapsed_time_wc = time.time() - start_time_wc
    print(f'This job took {elapsed_time_wc:.1f} seconds')
    os.system('say "All done"')
