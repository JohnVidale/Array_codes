#!/usr/bin/env python3
# John Vidale 4/2020
# last modified 8/2022

import os
import time
import matplotlib.pyplot as plt
from termcolor import colored

#%% Import functions
ev_directory = '/Users/vidale/Documents/GitHub/Array_codes/Repeat_qual'
os.chdir(ev_directory)
from run_ind_qual      import run_ind_qual

start_time_wc = time.time()

plt.close('all')


# Compare a pair of repeating events
repeater = 'YKA07'

if repeater ==    'YV21':
    eq_num1 = 606; eq_num2 = 634; temp_shift2 =  3.22; temp_shift_both = -2.30;
elif repeater == 'YKA08':
    eq_num1 = 606; eq_num2 = 650; temp_shift2 =  3.07; temp_shift_both = -3.00;
elif repeater == 'YKA07':
    eq_num1 = 606; eq_num2 = 618; temp_shift2 =  1.30; temp_shift_both = -2.50;
elif repeater == 'YKA13':
    eq_num1 = 611; eq_num2 = 651; temp_shift2 = -3.60; temp_shift_both =  0.30;
elif repeater == 'YV25':
    eq_num1 = 612; eq_num2 = 632; temp_shift2 = -2.25; temp_shift_both =  2.50;
elif repeater == 'YV26':
    eq_num1 = 613; eq_num2 = 637; temp_shift2 = -5.00; temp_shift_both =  0.28;
elif repeater == 'YKA14':
    eq_num1 = 615; eq_num2 = 648; temp_shift2 =  0.02; temp_shift_both =  1.80;
elif repeater == 'YV27':
    eq_num1 = 616; eq_num2 = 631; temp_shift2 = -0.35; temp_shift_both =  2.00;
elif repeater == 'YV22':
    eq_num1 = 618; eq_num2 = 634; temp_shift2 =  1.95; temp_shift_both = -2.62;
elif repeater == 'YV23':
    eq_num1 = 618; eq_num2 = 650; temp_shift2 =  1.78; temp_shift_both = -4.72;
elif repeater == 'YV28':
    eq_num1 = 621; eq_num2 = 643; temp_shift2 =  0.61; temp_shift_both = -2.15;
elif repeater == 'YV24':
    eq_num1 = 619; eq_num2 = 647; temp_shift2 =  0.74; temp_shift_both = -2.09;
elif repeater == 'YKA12':
    eq_num1 = 624; eq_num2 = 638; temp_shift2 = -4.82; temp_shift_both = -1.81;
elif repeater == 'YV1':
    eq_num1 = 625; eq_num2 = 641; temp_shift2 = -1.43; temp_shift_both =  1.59;
elif repeater == 'YV2':
    eq_num1 = 626; eq_num2 = 644; temp_shift2 =  0.60; temp_shift_both = -2.27;
elif repeater == 'YV7':
    eq_num1 = 634; eq_num2 = 650; temp_shift2 = -0.18; temp_shift_both = -2.77;
elif repeater == 'YV3':
    eq_num1 = 627; eq_num2 = 640; temp_shift2 = -2.10; temp_shift_both = -1.62;
elif repeater == 'YV4':
    eq_num1 = 628; eq_num2 = 636; temp_shift2 =  7.70; temp_shift_both = -1.18;
elif repeater == 'YV5':
    eq_num1 = 630; eq_num2 = 645; temp_shift2 =  2.48; temp_shift_both = -1.85;

# ILAR pair not in Pang and Koper
elif repeater == 'YKA11':  # YKA11 not retrieved yet
    eq_num1 = 620; eq_num2 = 646; temp_shift2 = 3.00; temp_shift_both = -2.36;
elif repeater == 'YKA17':  # IV2 is too close
    eq_num1 = 629; eq_num2 = 639; temp_shift2 = -0.67; temp_shift_both = -0.45;
elif repeater == 'YKA15':
    eq_num1 = 633; eq_num2 = 649; temp_shift2 = -0.18; temp_shift_both = -2.56;
elif repeater == 'YKA16':
    eq_num1 = 635; eq_num2 = 642; temp_shift2 =  4.85; temp_shift_both = -6.61;
elif repeater == 'YV8':
    eq_num1 = 652; eq_num2 = 606; temp_shift2 = -8.08; temp_shift_both =  2.00;
elif repeater == 'YV9':
    eq_num1 = 652; eq_num2 = 618; temp_shift2 = -6.80; temp_shift_both =  2.50;
elif repeater == 'YV10':
    eq_num1 = 652; eq_num2 = 634; temp_shift2 = -4.72; temp_shift_both =  3.00;
elif repeater == 'YV11':
    eq_num1 = 652; eq_num2 = 650; temp_shift2 = -4.81; temp_shift_both =  3.00;
elif repeater == 'YV12':
    eq_num1 = 653; eq_num2 = 660; temp_shift2 =  1.85; temp_shift_both =  0.00;
elif repeater == 'YV13':
    eq_num1 = 654; eq_num2 = 624; temp_shift2 =  2.30; temp_shift_both =  0.00;
elif repeater == 'YV14':
    eq_num1 = 654; eq_num2 = 638; temp_shift2 = -2.50; temp_shift_both =  0.00;
elif repeater == 'YV15':
    eq_num1 = 655; eq_num2 = 635; temp_shift2 = -2.65; temp_shift_both =  1.80;
elif repeater == 'YV16':
    eq_num1 = 655; eq_num2 = 642; temp_shift2 =  2.37; temp_shift_both =  1.80;
elif repeater == 'YV17':
    eq_num1 = 656; eq_num2 = 658; temp_shift2 = -2.16; temp_shift_both = -0.12;
elif repeater == 'YV18':  # YV18 not retrieved yet
    eq_num1 = 657; eq_num2 = 659; temp_shift2 =  0.90; temp_shift_both = -1.14;

# Too old for our ILAR dataset
elif repeater == 'YKA01':
    eq_num1 = 601; eq_num2 = 619; temp_shift2 = -5.75; temp_shift_both =  3.40;
elif repeater == 'YKA02':
    eq_num1 = 602; eq_num2 = 617; temp_shift2 =  0.00; temp_shift_both = -1.30;
elif repeater == 'YKA03':
    eq_num1 = 603; eq_num2 = 605; temp_shift2 = -5.40; temp_shift_both =  5.20;
elif repeater == 'YKA04':
    eq_num1 = 604; eq_num2 = 608; temp_shift2 =  1.95; temp_shift_both =  0.00;
elif repeater == 'YKA05':
    eq_num1 = 604; eq_num2 = 623; temp_shift2 = -3.90; temp_shift_both = -0.50;
elif repeater == 'YKA06':
    eq_num1 = 603; eq_num2 = 610; temp_shift2 = -3.42; temp_shift_both =  4.50;
elif repeater == 'YV19':
    eq_num1 = 605; eq_num2 = 610; temp_shift2 =  1.95; temp_shift_both =  0.00;
elif repeater == 'YV20':
    eq_num1 = 608; eq_num2 = 623; temp_shift2 =  1.95; temp_shift_both =  0.00;
elif repeater == 'YKA09':
    eq_num1 = 607; eq_num2 = 622; temp_shift2 = -0.60; temp_shift_both =  0.00;
elif repeater == 'YKA10':
    eq_num1 = 609; eq_num2 = 614; temp_shift2 =  0.35; temp_shift_both =  0.00;
elif repeater == 'YV25':
    eq_num1 = 601; eq_num2 = 647; temp_shift2 = -5.00; temp_shift_both =  4.00;
else:
    print(colored('No such event pair ' + repeater, 'yellow'))


freq_min = 1; freq_max = 2

start_buff = -10 # analysis window start relative to phase arrival
wind_len    = 40 # analysis window length

run_ind_qual(repeater = repeater, eq_num1 = eq_num1, eq_num2 = eq_num2, dphase = 'PKIKP', temp_shift_both = temp_shift_both,
                  freq_min = freq_min, freq_max = freq_max, ARRAY = 7, temp_shift2 = temp_shift2,
                  start_buff = start_buff, wind_len = wind_len)

elapsed_time_wc = time.time() - start_time_wc
print(f'This job took {elapsed_time_wc:.1f} seconds')
os.system('say "All done"')
