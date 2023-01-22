#!/usr/bin/env python3
# John Vidale 4/2020
# last modified 8/2022

def run_compare_ind(repeater = 'NoName',do_global = False, do_YKA = False, do_ILAR = False):

    import os
    import sys
    import time
    import matplotlib.pyplot as plt
    from termcolor import colored

    #%% Import functions
    pro_directory = '/Users/vidale/Documents/GitHub/Array_codes/Run_YKA_WRA'
    os.chdir(pro_directory)
    from run_compare_global     import run_compare_global
    from run_individual_df      import run_individual_df
    from run_compare_pair       import run_compare_pair

    start_time_wc = time.time()

    plt.close('all')

    if   repeater == 'P01':
        eq_num1 = 701; eq_num2 = 726; do_global = False; do_ILAR = False; do_YKA =  True; tshift_g = -5.75; tshift_I = -5.75; tshift_Y = -5.75; shift_both =  3.40; flip =  True
    elif repeater == 'P02':
        eq_num1 = 701; eq_num2 = 753; do_global = False; do_ILAR = False; do_YKA =  True; tshift_g = -4.85; tshift_I = -4.85; tshift_Y = -4.85; shift_both =  3.40; flip =  True
    elif repeater == 'P03':
        eq_num1 = 702; eq_num2 = 708; do_global = False; do_ILAR = False; do_YKA =  True; tshift_g = -8.08; tshift_I = -8.08; tshift_Y = -8.08; shift_both =  2.00; flip =  True
    elif repeater == 'P04':
        eq_num1 = 702; eq_num2 = 723; do_global = False; do_ILAR = False; do_YKA =  True; tshift_g = -6.80; tshift_I = -6.80; tshift_Y = -6.80; shift_both =  2.50; flip =  True
    elif repeater == 'P05':
        eq_num1 = 702; eq_num2 = 739; do_global = False; do_ILAR = False; do_YKA =  True; tshift_g = -4.72; tshift_I = -4.72; tshift_Y = -4.72; shift_both =  3.00; flip =  True
    elif repeater == 'P06':
        eq_num1 = 702; eq_num2 = 757; do_global = False; do_ILAR = False; do_YKA =  True; tshift_g = -4.81; tshift_I = -4.81; tshift_Y = -4.81; shift_both =  3.00; flip =  True
    elif repeater == 'P07':
        eq_num1 = 703; eq_num2 = 721; do_global =  True; do_ILAR = False; do_YKA =  True; tshift_g =  0.00; tshift_I =  0.00; tshift_Y =  0.00; shift_both = -1.30; flip =  True
    elif repeater == 'P08':
        eq_num1 = 704; eq_num2 = 706; do_global =  True; do_ILAR = False; do_YKA =  True; tshift_g = -5.35; tshift_I = -5.35; tshift_Y = -5.35; shift_both =  0.00; flip = False
    elif repeater == 'P09':
        eq_num1 = 704; eq_num2 = 713; do_global =  True; do_ILAR = False; do_YKA =  True; tshift_g = -3.40; tshift_I = -3.40; tshift_Y = -3.40; shift_both =  5.20; flip = False
    elif repeater == 'P10':
        eq_num1 = 705; eq_num2 = 711; do_global =  True; do_ILAR = False; do_YKA = False; tshift_g =  0.05; tshift_I =  0.05; tshift_Y =  0.05; shift_both =  0.00; flip = False
    elif repeater == 'P11':
        eq_num1 = 705; eq_num2 = 729; do_global =  True; do_ILAR = False; do_YKA =  True; tshift_g = -3.90; tshift_I = -3.90; tshift_Y = -3.90; shift_both = -0.50; flip = False
    elif repeater == 'P12':
        eq_num1 = 706; eq_num2 = 713; do_global =  True; do_ILAR = False; do_YKA =  True; tshift_g =  1.93; tshift_I =  1.93; tshift_Y =  1.93; shift_both =  0.00; flip = False
    elif repeater == 'P13':
        eq_num1 = 707; eq_num2 = 756; do_global =  True; do_ILAR = False; do_YKA =  True; tshift_g =  1.70; tshift_I =  1.70; tshift_Y =  1.70; shift_both =  0.00; flip = False
    elif repeater == 'P14':
        eq_num1 = 709; eq_num2 = 730; do_global =  True; do_ILAR = False; do_YKA =  True; tshift_g =  2.25; tshift_I =  2.25; tshift_Y =  2.25; shift_both =  0.00; flip = False
    elif repeater == 'P15':
        eq_num1 = 709; eq_num2 = 744; do_global =  True; do_ILAR = False; do_YKA =  True; tshift_g = -2.60; tshift_I = -2.60; tshift_Y = -2.60; shift_both =  0.00; flip = False
    elif repeater == 'P16':
        eq_num1 = 708; eq_num2 = 723; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift_g =  1.30; tshift_I =  1.30; tshift_Y =  1.30; shift_both = -6.00; flip =  True
    elif repeater == 'P17':
        eq_num1 = 708; eq_num2 = 739; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift_g =  3.22; tshift_I =  3.22; tshift_Y =  3.22; shift_both = -6.00; flip =  True
    elif repeater == 'P18':
        eq_num1 = 708; eq_num2 = 757; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift_g =  3.07; tshift_I =  3.07; tshift_Y =  3.07; shift_both = -6.00; flip =  True
    elif repeater == 'P19':
        eq_num1 = 710; eq_num2 = 728; do_global =  True; do_ILAR = False; do_YKA =  True; tshift_g = -0.60; tshift_I = -0.60; tshift_Y = -0.60; shift_both =  0.00; flip = False
    elif repeater == 'P20':
        eq_num1 = 711; eq_num2 = 729; do_global =  True; do_ILAR = False; do_YKA =  True; tshift_g = -3.87; tshift_I = -3.87; tshift_Y = -3.87; shift_both =  0.00; flip = False
    elif repeater == 'P21':
        eq_num1 = 712; eq_num2 = 718; do_global =  True; do_ILAR = False; do_YKA =  True; tshift_g =  0.35; tshift_I =  0.35; tshift_Y =  0.35; shift_both =  0.00; flip = False
    elif repeater == 'P22':
        eq_num1 = 714; eq_num2 = 741; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift_g = -2.90; tshift_I = -2.90; tshift_Y = -2.90; shift_both =  1.80; flip = False
    elif repeater == 'P23':
        eq_num1 = 714; eq_num2 = 748; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift_g =  1.95; tshift_I =  2.00; tshift_Y =  2.00; shift_both =  1.80; flip = False
    elif repeater == 'P24':
        eq_num1 = 715; eq_num2 = 758; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift_g = -3.60; tshift_I = -3.60; tshift_Y = -3.60; shift_both = -1.18; flip =  True
    elif repeater == 'P25':
        eq_num1 = 716; eq_num2 = 737; do_global =  True; do_ILAR =  True; do_YKA = False; tshift_g = -2.25; tshift_I = -2.25; tshift_Y = -2.25; shift_both =  0.00; flip =  True
    elif repeater == 'P26':
        eq_num1 = 717; eq_num2 = 743; do_global =  True; do_ILAR =  True; do_YKA = False; tshift_g = -4.98; tshift_I = -4.98; tshift_Y = -4.98; shift_both =  0.28; flip = False
    elif repeater == 'P27':
        eq_num1 = 719; eq_num2 = 754; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift_g =  0.02; tshift_I =  0.02; tshift_Y =  0.02; shift_both = -0.10; flip =  True
    elif repeater == 'P28':
        eq_num1 = 720; eq_num2 = 736; do_global =  True; do_ILAR =  True; do_YKA = False; tshift_g = -3.95; tshift_I = -3.95; tshift_Y = -3.95; shift_both =  0.37; flip = False
    elif repeater == 'P29':
        eq_num1 = 722; eq_num2 = 740; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift_g = -2.30; tshift_I = -2.30; tshift_Y = -2.30; shift_both = -0.12; flip =  True
    elif repeater == 'P30':
        eq_num1 = 723; eq_num2 = 739; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift_g =  1.95; tshift_I =  1.95; tshift_Y =  1.95; shift_both = -4.72; flip =  True
    elif repeater == 'P31':
        eq_num1 = 723; eq_num2 = 757; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift_g =  1.78; tshift_I =  1.78; tshift_Y =  1.78; shift_both = -4.72; flip = False
    elif repeater == 'P32':
        eq_num1 = 724; eq_num2 = 746; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift_g =  0.80; tshift_I =  0.80; tshift_Y =  0.80; shift_both = -1.14; flip = False
    elif repeater == 'P33':   # weird that global disagrees, but ILAR is very similar
        eq_num1 = 725; eq_num2 = 749; do_global =  True; do_ILAR =  True; do_YKA = False; tshift_g =  0.56; tshift_I =  0.56; tshift_Y =  0.56; shift_both = -2.15; flip = False
    elif repeater == 'P34':
        eq_num1 = 726; eq_num2 = 753; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift_g =  0.74; tshift_I =  0.74; tshift_Y =  0.74; shift_both = -2.09; flip =  True
    elif repeater == 'P35':
        eq_num1 = 727; eq_num2 = 752; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift_g =  3.00; tshift_I =  3.00; tshift_Y =  3.00; shift_both = -2.36; flip =  True
    elif repeater == 'P36':
        eq_num1 = 730; eq_num2 = 744; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift_g = -4.82; tshift_I = -4.82; tshift_Y = -4.82; shift_both = -1.81; flip =  True
    elif repeater == 'P37':
        eq_num1 = 731; eq_num2 = 747; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift_g = -1.45; tshift_I = -1.45; tshift_Y = -1.45; shift_both = -2.41; flip =  True
    elif repeater == 'P38':  # big changes only in YKA
        eq_num1 = 732; eq_num2 = 750; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift_g =  0.60; tshift_I =  0.60; tshift_Y =  0.60; shift_both = -2.27; flip =  True
    elif repeater == 'P39':  # some surprising global shifts
        eq_num1 = 733; eq_num2 = 745; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift_g = -2.12; tshift_I = -2.12; tshift_Y = -2.12; shift_both = -1.62; flip =  True
    elif repeater == 'P40':
        eq_num1 = 734; eq_num2 = 742; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift_g =  7.70; tshift_I =  7.70; tshift_Y =  7.70; shift_both = -1.18; flip = False
    elif repeater == 'P41':
        eq_num1 = 735; eq_num2 = 751; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift_g =  2.48; tshift_I =  2.58; tshift_Y =  2.68; shift_both = -1.85; flip = False
    elif repeater == 'P42':
        eq_num1 = 738; eq_num2 = 755; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift_g = -0.15; tshift_I = -0.15; tshift_Y = -0.15; shift_both = -2.56; flip =  True
    elif repeater == 'P43':  # big changes only in YKA
        eq_num1 = 739; eq_num2 = 757; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift_g = -0.15; tshift_I = -0.15; tshift_Y = -0.15; shift_both = -2.77; flip =  True
    elif repeater == 'P44':  # big changes only in YKA
        eq_num1 = 741; eq_num2 = 748; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift_g =  4.86; tshift_I =  4.86; tshift_Y =  4.86; shift_both = -6.61; flip = False
    elif repeater == 'P45':  # Sumatra, IV2 is too close
        eq_num1 = 759; eq_num2 = 760; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift_g = -0.67; tshift_I = -0.67; tshift_Y = -0.67; shift_both = -0.45; flip =  True
    else:
        print(colored('No such event pair ' + repeater, 'yellow'))
        sys.exit(-1)

    freq_min = 1; freq_max = 2

    if do_global:
        start_buff = -20 # analysis window start relative to phase arrival
        wind_len    = 100 # analysis window length

        run_compare_global(repeater = repeater, eq_num1 = eq_num1, eq_num2 = eq_num2,
                freq_min = freq_min, freq_max = freq_max, ARRAY = 7, tshift = tshift_g,
                start_buff = start_buff, wind_len = wind_len)

    wind_buff   =  10 # buffer before and after time window of analysis

    Zstart_buff = -5 # analysis window start relative to phase arrival
    wind_len    = 15 # analysis window length
    plot_peak = 1

    beam_offset = 0.01
    beam_width  = 0.03
    slow_delta  = 0.004

    if do_YKA:
        run_compare_pair(eq_num1 = eq_num1, eq_num2 = eq_num2, dphase = 'PKIKP', shift_both = shift_both,
                beam_width = beam_width, slow_delta = slow_delta, beam_offset = beam_offset,
                freq_min = freq_min, freq_max = freq_max, stat_corr = 0, ARRAY = 5,
                Zstart_buff = Zstart_buff, wind_len = wind_len, wind_buff = wind_buff, flip = flip,
                tshift = tshift_Y, fig_index = 200, do_interpolate =  True, pair_name = repeater, plot_peak = plot_peak)

    if do_ILAR:
        run_compare_pair(eq_num1 = eq_num1, eq_num2 = eq_num2, dphase = 'PKIKP', shift_both = shift_both,
                beam_width = beam_width, slow_delta = slow_delta, beam_offset = beam_offset,
                freq_min = freq_min, freq_max = freq_max, stat_corr = 0, ARRAY = 6,
                Zstart_buff = Zstart_buff, wind_len = wind_len, wind_buff = wind_buff, flip = flip,
                tshift = tshift_I, fig_index = 300, do_interpolate =  True, pair_name = repeater, plot_peak = plot_peak)

    elapsed_time_wc = time.time() - start_time_wc
    print(f'This job took {elapsed_time_wc:.1f} seconds')
    os.system('say "All done"')
