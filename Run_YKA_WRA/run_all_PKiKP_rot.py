#!/usr/bin/env python3
# John Vidale 4/2020
# last modified 8/2022

import os
import time
import matplotlib.pyplot as plt

#%% Import functions
ev_directory = '/Users/vidale/Documents/GitHub/Array_codes/Run_YKA_WRA'
os.chdir(ev_directory)
from run_individual_df      import run_individual_df
from run_compare_pair_auto       import run_compare_pair

start_time_wc = time.time()

plt.close('all')

# eq_num1 = 635; eq_num2 = 642  # compare both YKA and ILAR for both repeating events in individual plots

# run_individual_df(eq_num = eq_num1, dphase = 'PKiKP', start_buff_stack = -50, end_buff_stack = 75, start_beam_stack = -15, end_beam_stack = 40, beam_width = 0.03, beam_step = 0.0025, freq_min = 0.5, freq_max = 1, fig_index = 100, stat_corr = 0, ARRAY = 5)

# run_individual_df(eq_num = eq_num2, dphase = 'PKiKP', start_buff_stack = -50, end_buff_stack = 75, start_beam_stack = -15, end_beam_stack = 40, beam_width = 0.03, beam_step = 0.0025, freq_min = 0.5, freq_max = 1, fig_index = 200, stat_corr = 0, ARRAY = 5)

# run_individual_df(eq_num = eq_num1, dphase = 'PKiKP', start_buff_stack = -50, end_buff_stack = 75, start_beam_stack = -15, end_beam_stack = 40, beam_width = 0.03, beam_step = 0.0025, freq_min = 0.5, freq_max = 1, fig_index = 300, stat_corr = 0, ARRAY = 6)

# run_individual_df(eq_num = eq_num2, dphase = 'PKiKP', start_buff_stack = -50, end_buff_stack = 75, start_beam_stack = -15, end_beam_stack = 40, beam_width = 0.03, beam_step = 0.0025, freq_min = 0.5, freq_max = 1, fig_index = 400, stat_corr = 0, ARRAY = 6)

# look at individual event
eq_num = 650
ARRAY = 6

# inspect 1 event, 1 frequency
# run_individual_df(eq_num = eq_num, dphase = 'PKIKP', start_buff_stack = -30, end_buff_stack = 5, start_beam_stack = -12, end_beam_stack = -1.5, beam_width = 0.05, beam_step = 0.0025, freq_min = 1, freq_max = 4, fig_index = 100, stat_corr = 1, ARRAY = ARRAY)

# do 4 frequencies
# run_individual_df(eq_num = eq_num, dphase = 'PKIKP', start_buff_stack = -50, end_buff_stack = 75, start_beam_stack = -15, end_beam_stack = 40, beam_width = 0.03, beam_step = 0.0025, freq_min = 0.5, freq_max = 1, fig_index = 100, stat_corr = 0, ARRAY = ARRAY)
# run_individual_df(eq_num = eq_num, dphase = 'PKIKP', start_buff_stack = -50, end_buff_stack = 75, start_beam_stack = -15, end_beam_stack = 40, beam_width = 0.03, beam_step = 0.0025, freq_min = 1, freq_max = 3, fig_index = 100, stat_corr = 0, ARRAY = ARRAY)
# run_individual_df(eq_num = eq_num, dphase = 'PKiKP', start_buff_stack = -50, end_buff_stack = 75, start_beam_stack = -15, end_beam_stack = 40, beam_width = 0.03, beam_step = 0.0025, freq_min = 2, freq_max = 4, fig_index = 200, stat_corr = 0, ARRAY = ARRAY)
# run_individual_df(eq_num = eq_num, dphase = 'PKiKP', start_buff_stack = -50, end_buff_stack = 75, start_beam_stack = -15, end_beam_stack = 40, beam_width = 0.03, beam_step = 0.0025, freq_min = 3, freq_max = 5, fig_index = 300, stat_corr = 0, ARRAY = ARRAY)

# Compare a pair of repeating events
repeater = 'YV25'

# Time shifts to align
#           YKA events 136°-139°, plus one at 115°
#           ILAR events 150°-153°

# ILAR pair is in Pang and Koper
if repeater == 'ILAR01' or repeater == 'YV21':
    eq_num1 = 606; eq_num2 = 634; temp_shift2_ILAR =  3.22; temp_shift2_YKA =  2.62; temp_shift_both = -6.00; flip = True
if repeater == 'ILAR02' or repeater == 'YKA08':
    eq_num1 = 606; eq_num2 = 650; temp_shift2_ILAR =  3.07; temp_shift2_YKA =  3.27; temp_shift_both = -6.00; flip = True
if repeater == 'ILAR03' or repeater == 'YKA07':
    eq_num1 = 606; eq_num2 = 618; temp_shift2_ILAR =  1.30; temp_shift2_YKA =  1.30; temp_shift_both = -6.00; flip = True
if repeater == 'ILAR04' or repeater == 'YKA13':
    eq_num1 = 611; eq_num2 = 651; temp_shift2_ILAR = -3.60; temp_shift2_YKA = -3.40; temp_shift_both = -1.18; flip = True
if repeater == 'ILAR05' or repeater == 'NoName':
    eq_num1 = 612; eq_num2 = 632; temp_shift2_ILAR = -2.25; temp_shift2_YKA = -2.25; temp_shift_both =  0.02; flip = True
if repeater == 'ILAR06' or repeater == 'NoName':
    eq_num1 = 613; eq_num2 = 637; temp_shift2_ILAR = -4.98; temp_shift2_YKA = -5.00; temp_shift_both =  0.28; flip = False
if repeater == 'ILAR07' or repeater == 'YKA14':
    eq_num1 = 615; eq_num2 = 648; temp_shift2_ILAR =  0.02; temp_shift2_YKA =  0.22; temp_shift_both = -0.10; flip = True
if repeater == 'ILAR08' or repeater == 'NoName':
    eq_num1 = 616; eq_num2 = 631; temp_shift2_ILAR = -3.98; temp_shift2_YKA = -3.95; temp_shift_both =  0.37; flip = False
if repeater == 'ILAR09' or repeater == 'YV22':
    eq_num1 = 618; eq_num2 = 634; temp_shift2_ILAR =  1.95; temp_shift2_YKA =  1.85; temp_shift_both = -4.72; flip = True
if repeater == 'ILAR10' or repeater == 'YV23':
    eq_num1 = 618; eq_num2 = 650; temp_shift2_ILAR =  1.78; temp_shift2_YKA =  1.78; temp_shift_both = -4.72; flip = False
if repeater == 'ILAR11' or repeater == 'NoName':
    eq_num1 = 621; eq_num2 = 643; temp_shift2_ILAR =  0.56; temp_shift2_YKA =  0.61; temp_shift_both = -2.15; flip = False
if repeater == 'ILAR12' or repeater == 'YV24':
    eq_num1 = 619; eq_num2 = 647; temp_shift2_ILAR =  0.74; temp_shift2_YKA =  0.74; temp_shift_both = -2.09; flip = True
if repeater == 'ILAR13' or repeater == 'YKA12':
    eq_num1 = 624; eq_num2 = 638; temp_shift2_ILAR = -4.82; temp_shift2_YKA = -4.72; temp_shift_both = -1.81; flip = True
if repeater == 'ILAR14' or repeater == 'YV1':
    eq_num1 = 625; eq_num2 = 641; temp_shift2_ILAR = -1.45; temp_shift2_YKA = -1.30; temp_shift_both = -2.41; flip = True
if repeater == 'ILAR15' or repeater == 'YV2':
    eq_num1 = 626; eq_num2 = 644; temp_shift2_ILAR =  0.62; temp_shift2_YKA =  0.70; temp_shift_both = -2.27; flip = True
if repeater == 'ILAR16' or repeater == 'YV7':
    eq_num1 = 634; eq_num2 = 650; temp_shift2_ILAR = -0.15; temp_shift2_YKA = -0.05; temp_shift_both = -2.77; flip = True
if repeater == 'ILAR17' or repeater == 'YV3':
    eq_num1 = 627; eq_num2 = 640; temp_shift2_ILAR = -2.12; temp_shift2_YKA = -2.10; temp_shift_both = -1.62; flip = True
if repeater == 'ILAR18' or repeater == 'YV4':
    eq_num1 = 628; eq_num2 = 636; temp_shift2_ILAR =  7.70; temp_shift2_YKA =  7.70; temp_shift_both = -1.18; flip = False
if repeater == 'ILAR19' or repeater == 'YV5':
    eq_num1 = 630; eq_num2 = 645; temp_shift2_ILAR =  2.58; temp_shift2_YKA =  2.68; temp_shift_both = -1.85; flip = False

# ILAR pair not in Pang and Koper
if repeater == 'IV1' or repeater == 'YKA11':  # YKA11 not retrieved yet
    eq_num1 = 620; eq_num2 = 646; temp_shift2_ILAR = 3.00; temp_shift2_YKA = 3.20; temp_shift_both = -2.36; flip = True
if repeater == 'IV2' or repeater == 'YKA17':  # IV2 is too close
    eq_num1 = 629; eq_num2 = 639; temp_shift2_ILAR =     0; temp_shift2_YKA = -0.67; temp_shift_both = -0.45; flip = True
if repeater == 'IV4' or repeater == 'YKA15':
    eq_num1 = 633; eq_num2 = 649; temp_shift2_ILAR = -0.15; temp_shift2_YKA = -0.11; temp_shift_both = -2.56; flip = True
if repeater == 'IV5' or repeater == 'YKA16':
    eq_num1 = 635; eq_num2 = 642; temp_shift2_ILAR =  4.86; temp_shift2_YKA =  4.93; temp_shift_both = -6.61; flip = False
if repeater == 'IV6'    or repeater == 'YV8':
    eq_num1 = 652; eq_num2 = 606; temp_shift2_ILAR =  0.00; temp_shift2_YKA = -8.08; temp_shift_both =  2.00; flip = True
if repeater == 'IV7'    or repeater == 'YV9':
    eq_num1 = 652; eq_num2 = 618; temp_shift2_ILAR =  0.00; temp_shift2_YKA = -6.80; temp_shift_both =  2.50; flip = True
if repeater == 'IV8'    or repeater == 'YV10':
    eq_num1 = 652; eq_num2 = 634; temp_shift2_ILAR =  0.00; temp_shift2_YKA = -4.72; temp_shift_both =  3.00; flip = True
if repeater == 'IV9'    or repeater == 'YV11':
    eq_num1 = 652; eq_num2 = 650; temp_shift2_ILAR =  0.00; temp_shift2_YKA = -4.81; temp_shift_both =  3.00; flip = True
if repeater == 'IV10'   or repeater == 'YV12':
    eq_num1 = 653; eq_num2 = 660; temp_shift2_ILAR =  0.00; temp_shift2_YKA =  1.85; temp_shift_both =  0.00; flip = False
if repeater == 'IV11'   or repeater == 'YV13':
    eq_num1 = 654; eq_num2 = 624; temp_shift2_ILAR =  0.00; temp_shift2_YKA =  2.30; temp_shift_both =  0.00; flip = False
if repeater == 'IV12'   or repeater == 'YV14':
    eq_num1 = 654; eq_num2 = 638; temp_shift2_ILAR =  0.00; temp_shift2_YKA = -2.50; temp_shift_both =  0.00; flip = False
if repeater == 'IV13'   or repeater == 'YV15':
    eq_num1 = 655; eq_num2 = 635; temp_shift2_ILAR = -2.90; temp_shift2_YKA = -2.65; temp_shift_both =  1.80; flip = False
if repeater == 'IV14'   or repeater == 'YV16':
    eq_num1 = 655; eq_num2 = 642; temp_shift2_ILAR =  2.00; temp_shift2_YKA =  2.37; temp_shift_both =  1.80; flip = False
if repeater == 'IV15'   or repeater == 'YV17':
    eq_num1 = 656; eq_num2 = 658; temp_shift2_ILAR = -2.30; temp_shift2_YKA = -2.16; temp_shift_both = -0.12; flip = True
if repeater == 'IV16'   or repeater == 'YV18':  # YV18 not retrieved yet
    eq_num1 = 657; eq_num2 = 659; temp_shift2_ILAR =  0.80; temp_shift2_YKA =  0.90; temp_shift_both = -1.14; flip = False

# Too old for our ILAR dataset
if repeater == 'NoName'   or repeater == 'YKA01':
    eq_num1 = 601; eq_num2 = 619; temp_shift2_ILAR =  0.00; temp_shift2_YKA = -5.75; temp_shift_both =  3.40; flip = True
if repeater == 'NoName'   or repeater == 'YKA02':
    eq_num1 = 602; eq_num2 = 617; temp_shift2_ILAR =  0.00; temp_shift2_YKA =  0.00; temp_shift_both = -1.30; flip = True
if repeater == 'NoName'   or repeater == 'YKA03':
    eq_num1 = 603; eq_num2 = 605; temp_shift2_ILAR =  0.00; temp_shift2_YKA = -5.40; temp_shift_both =  5.20; flip = False
if repeater == 'NoName'   or repeater == 'YKA04':
    eq_num1 = 604; eq_num2 = 608; temp_shift2_ILAR =  0.00; temp_shift2_YKA =  1.95; temp_shift_both =  0.00; flip = False
if repeater == 'NoName'   or repeater == 'YKA05':
    eq_num1 = 604; eq_num2 = 623; temp_shift2_ILAR =  0.00; temp_shift2_YKA = -3.90; temp_shift_both = -0.50; flip = False
if repeater == 'NoName'   or repeater == 'YKA06':
    eq_num1 = 603; eq_num2 = 610; temp_shift2_ILAR =  0.00; temp_shift2_YKA = -3.42; temp_shift_both =  4.50; flip = False
if repeater == 'NoName'   or repeater == 'YV19':
    eq_num1 = 605; eq_num2 = 610; temp_shift2_ILAR =  0.00; temp_shift2_YKA =  1.95; temp_shift_both =  0.00; flip = False
if repeater == 'NoName'   or repeater == 'YV20':
    eq_num1 = 608; eq_num2 = 623; temp_shift2_ILAR =  0.00; temp_shift2_YKA =  1.95; temp_shift_both =  0.00; flip = False
if repeater == 'NoName'   or repeater == 'YKA09':
    eq_num1 = 607; eq_num2 = 622; temp_shift2_ILAR =  0.00; temp_shift2_YKA = -0.60; temp_shift_both =  0.00; flip = False
if repeater == 'NoName'   or repeater == 'YKA10':
    eq_num1 = 609; eq_num2 = 614; temp_shift2_ILAR =  0.00; temp_shift2_YKA =  0.35; temp_shift_both =  0.00; flip = False
if repeater == 'NoName'   or repeater == 'YV25':
    eq_num1 = 601; eq_num2 = 647; temp_shift2_ILAR =  0.00; temp_shift2_YKA = -5.00; temp_shift_both =  4.00; flip = False

if repeater[0] == 'Y':
    ARRAY = 5; temp_shift2 = temp_shift2_YKA
elif repeater[0] == 'I':
    ARRAY = 6; temp_shift2 = temp_shift2_ILAR

freq_min = 1; freq_max = 4

wind_buff   =  10 # buffer before and after time window of analysis
#  ILAR
# Zstart_buff = -10 # analysis window start relative to phase arrival
# wind_len    =  12 # analysis window length
# plot_peak = 0.1
# Zstart_buff = -20 # analysis window start relative to phase arrival
# wind_len    =  30 # analysis window length
# plot_peak = 1
# Zstart_buff = -4 # analysis window start relative to phase arrival
# wind_len    =  14 # analysis window length
# plot_peak = 1

# YKA
# Zstart_buff = -10 # analysis window start relative to phase arrival
# wind_len    =  40 # analysis window length
# plot_peak = 1
Zstart_buff = -5 # analysis window start relative to phase arrival
wind_len    = 15 # analysis window length
plot_peak = 1

beam_offset = 0.01
beam_width  = 0.03
slow_delta  = 0.004

run_compare_pair(eq_num1 = eq_num1, eq_num2 = eq_num2, dphase = 'PKIKP', temp_shift_both = temp_shift_both,
                  beam_width = beam_width, slow_delta = slow_delta, beam_offset = beam_offset,
                  freq_min = freq_min, freq_max = freq_max, stat_corr = 0, ARRAY = ARRAY,
                  Zstart_buff = Zstart_buff, wind_len = wind_len, wind_buff = wind_buff, flip = flip,
                  temp_shift2 = temp_shift2, fig_index = 200, do_interpolate = True, pair_name = repeater, plot_peak = plot_peak)

elapsed_time_wc = time.time() - start_time_wc
print(f'This job took {elapsed_time_wc:.1f} seconds')
os.system('say "All done"')
