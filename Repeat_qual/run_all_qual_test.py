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
start_buff = -20 # analysis window start relative to phase arrival
wind_len    = 1800 # analysis window length

run_ind_qual(eq_num1 = 606, eq_num2 = 650, dphase = 'PKIKP', temp_shift_both = 0,
                  freq_min = 0.1, freq_max = 10, ARRAY = 7,
                  start_buff = start_buff, wind_len = wind_len)

elapsed_time_wc = time.time() - start_time_wc
print(f'This job took {elapsed_time_wc:.1f} seconds')
os.system('say "All done"')
