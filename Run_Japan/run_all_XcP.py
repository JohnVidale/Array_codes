#!/usr/bin/env python3
# John Vidale 4/2020

import os
import time
import matplotlib.pyplot as plt

#%% Import functions
from run_individual      import run_individual
from run_dual            import run_dual
from run_bin_one      import run_bin_one
from run_individual_get_shifts      import run_individual_get_shifts

start_time_wc = time.time()

plt.close('all')

# 101
# run_individual(eq_num = 101, dphase = 'PcP'  , start_beam_stack =    2, end_beam_stack =    7, R_slow_plot = -0.0098, T_slow_plot = 0.0072, fig_index = 100, stat_corr = 1, apply_SNR = True)
# run_individual(eq_num = 101, dphase = 'ScP'  , start_beam_stack =    3, end_beam_stack =    8, R_slow_plot = 0.013, T_slow_plot = 0.000, fig_index = 200, stat_corr = 1, apply_SNR = True)
# run_individual(eq_num = 101, dphase = 'PKiKP', start_beam_stack =    1, end_beam_stack =    6, R_slow_plot = 0.013, T_slow_plot = 0.000, fig_index = 300, stat_corr = 1, apply_SNR = True)

# 102
# run_individual(eq_num = 102, dphase = 'PcP'  , start_beam_stack =    0, end_beam_stack =    3, R_slow_plot = 0.013, T_slow_plot = 0.000, fig_index = 100, stat_corr = 1, apply_SNR = True)
# run_individual(eq_num = 102, dphase = 'ScP'  , start_beam_stack = -3.5, end_beam_stack =  0.5, R_slow_plot = 0.014, T_slow_plot = 0.000, fig_index = 200, stat_corr = 1, apply_SNR = True)
# run_individual(eq_num = 102, dphase = 'PKiKP', start_beam_stack =    0, end_beam_stack =    3, R_slow_plot = 0.003, T_slow_plot = 0.000, fig_index = 300, stat_corr = 1, apply_SNR = True)

# 103
# run_individual(eq_num = 103, dphase = 'PcP'  , start_beam_stack =   -3, end_beam_stack =    0, R_slow_plot = 0.013, T_slow_plot = 0.000, fig_index = 100, stat_corr = 1, apply_SNR = True)
# run_individual(eq_num = 103, dphase = 'ScP'  , start_beam_stack =   -2, end_beam_stack =    1, R_slow_plot = 0.013, T_slow_plot = 0.000, fig_index = 200, stat_corr = 1, apply_SNR = True)
# run_individual(eq_num = 103, dphase = 'PKiKP', start_beam_stack =   -3, end_beam_stack =    0, R_slow_plot = 0.013, T_slow_plot = 0.000, fig_index = 300, stat_corr = 1, apply_SNR = True)

# 104
# run_individual(eq_num = 104, dphase = 'PcP'  , start_beam_stack =   -1, end_beam_stack =    2, R_slow_plot = 0.013, T_slow_plot = 0.000, fig_index = 100, stat_corr = 1, apply_SNR = True)
# run_individual(eq_num = 104, dphase = 'ScP'  , start_beam_stack = -0.5, end_beam_stack =    1, R_slow_plot = 0.014, T_slow_plot = 0.000, fig_index = 200, stat_corr = 1, apply_SNR = True)
# run_individual(eq_num = 104, dphase = 'PKiKP', start_beam_stack =   -2, end_beam_stack =    1, R_slow_plot = 0.003, T_slow_plot = 0.000, fig_index = 300, stat_corr = 1, apply_SNR = True)

# 105
# run_individual(eq_num = 105, dphase = 'PcP'  , start_beam_stack =   -2, end_beam_stack =  0.5, R_slow_plot = 0.013, T_slow_plot = 0.000, fig_index = 100, stat_corr = 1, apply_SNR = True)
# run_individual(eq_num = 105, dphase = 'ScP'  , start_beam_stack =  0.5, end_beam_stack =    2, R_slow_plot = 0.014, T_slow_plot = 0.000, fig_index = 200, stat_corr = 1, apply_SNR = True)
# run_individual(eq_num = 105, dphase = 'PKiKP', start_beam_stack =   -1, end_beam_stack =  1.5, R_slow_plot = 0.003, T_slow_plot = 0.000, fig_index = 300, stat_corr = 1, apply_SNR = True)

# 106
# run_individual(eq_num = 106, dphase = 'PcP'  , start_beam_stack =    5, end_beam_stack =    7, R_slow_plot = 0.013, T_slow_plot = 0.000, fig_index = 100, stat_corr = 1, apply_SNR = True)
# run_individual(eq_num = 106, dphase = 'ScP'  , start_beam_stack =  6.5, end_beam_stack =    9, R_slow_plot = 0.014, T_slow_plot = 0.000, fig_index = 200, stat_corr = 1, apply_SNR = True)
# run_individual(eq_num = 106, dphase = 'PKiKP', start_beam_stack =    2, end_beam_stack =    4, R_slow_plot = 0.003, T_slow_plot = 0.000, fig_index = 300, stat_corr = 1, apply_SNR = True)

# investigate multipathing
# run_individual(eq_num = 106, dphase = 'ScP'  , start_beam_stack = 6, end_beam_stack = 10, precursor_shift = 0, signal_dur = 10,R_slow_plot = 0.0075, T_slow_plot = 0.000,
#                fig_index = 200, stat_corr = 1, apply_SNR = True)
# run_individual(eq_num = 106, dphase = 'PKiKP'  , start_beam_stack = 2, end_beam_stack = 4, precursor_shift = 2, signal_dur = 4,R_slow_plot = 0.0075, T_slow_plot = 0.000,
#                fig_index = 200, stat_corr = 1, apply_SNR = True)

# 107
# run_individual(eq_num = 107, dphase = 'PcP'  , start_beam_stack = -0.5, end_beam_stack =    2, R_slow_plot = -0.0116, T_slow_plot = 0.0084, fig_index = 100, stat_corr = 1, apply_SNR = True)
# run_individual(eq_num = 107, dphase = 'ScP'  , start_beam_stack = -3.5, end_beam_stack =   -1, R_slow_plot = 0.014, T_slow_plot = 0.000, fig_index = 200, stat_corr = 1, apply_SNR = True)
# run_individual(eq_num = 107, dphase = 'PKiKP', start_beam_stack =   -2, end_beam_stack =    1, R_slow_plot = 0.003, T_slow_plot = 0.000, fig_index = 300, stat_corr = 1, apply_SNR = True)

# 108
# run_individual(eq_num = 108, dphase = 'PcP'  , start_beam_stack = -2, end_beam_stack =    1, R_slow_plot = 0.013, T_slow_plot = 0.000, fig_index = 100, stat_corr = 3, apply_SNR = True)
# run_individual(eq_num = 108, dphase = 'ScP'  , start_beam_stack =   -3, end_beam_stack =  0.5, R_slow_plot = 0.014, T_slow_plot = 0.000, fig_index = 200, stat_corr = 1, apply_SNR = True)
# run_individual(eq_num = 108, dphase = 'PKiKP', start_beam_stack =   -1, end_beam_stack =  1.5, R_slow_plot = 0.003, T_slow_plot = 0.000, fig_index = 300, stat_corr = 3, apply_SNR = True)

# 109
# run_individual(eq_num = 109, dphase = 'PcP'  , start_beam_stack =   -1, end_beam_stack =    2, R_slow_plot = -0.0089, T_slow_plot = 0.0021, fig_index = 100, stat_corr = 1, apply_SNR = True)
# run_individual(eq_num = 109, dphase = 'ScP'  , start_beam_stack = -1.5, end_beam_stack =  1.5, R_slow_plot = -0.009, T_slow_plot = 0.0045, fig_index = 200, stat_corr = 1, apply_SNR = True)
# run_individual(eq_num = 109, dphase = 'PKiKP', start_beam_stack = -1.5, end_beam_stack =    1, R_slow_plot = 0.003, T_slow_plot = 0.000, fig_index = 300, stat_corr = 1, apply_SNR = True)

# 110
# run_individual(eq_num = 110, dphase = 'PcP'  , start_beam_stack =  0.5, end_beam_stack =  2.5, R_slow_plot = -0.0055, T_slow_plot = 0.0011, fig_index = 100, stat_corr = 1, apply_SNR = True)
# run_individual(eq_num = 110, dphase = 'ScP'  , start_beam_stack =    0, end_beam_stack =    2, R_slow_plot = 0.014, T_slow_plot = 0.000, fig_index = 200, stat_corr = 1, apply_SNR = True)
# run_individual(eq_num = 110, dphase = 'PKiKP', start_beam_stack =   -1, end_beam_stack =  1.5, R_slow_plot = 0.003, T_slow_plot = 0.000, fig_index = 300, stat_corr = 1, apply_SNR = True)

# 111
# run_individual(eq_num = 111, dphase = 'PcP'  , start_beam_stack =  1.5, end_beam_stack =  3.5, R_slow_plot = 0.013, T_slow_plot = 0.000, fig_index = 100, stat_corr = 1, apply_SNR = True)
# run_individual(eq_num = 111, dphase = 'ScP'  , start_beam_stack =    1, end_beam_stack =    5, R_slow_plot = 0.014, T_slow_plot = 0.000, fig_index = 200, stat_corr = 1, apply_SNR = True)
# run_individual(eq_num = 111, dphase = 'PKiKP', start_beam_stack = -2.5, end_beam_stack =  0.5, R_slow_plot = 0.003, T_slow_plot = 0.000, fig_index = 300, stat_corr = 1, apply_SNR = True)

# 112
# run_individual(eq_num = 112, dphase = 'PcP'  , start_beam_stack =  1.5, end_beam_stack =  3.5, R_slow_plot = 0.013, T_slow_plot = 0.000, fig_index = 100, stat_corr = 1, apply_SNR = True)
# run_individual(eq_num = 112, dphase = 'ScP'  , start_beam_stack =   -1, end_beam_stack =    2, R_slow_plot = 0.014, T_slow_plot = 0.000, fig_index = 200, stat_corr = 1, apply_SNR = True)
# run_individual(eq_num = 112, dphase = 'PKiKP', start_beam_stack =  2.5, end_beam_stack =  5.5, R_slow_plot = 0.003, T_slow_plot = 0.000, fig_index = 300, stat_corr = 1, apply_SNR = True)

# 113 - Long window
# run_individual(eq_num = 113, dphase = 'PcP'  , start_beam_stack =   -5, end_beam_stack =   12, R_slow_plot =  0.03, T_slow_plot = 0.000, fig_index = 100, stat_corr = 1, apply_SNR = True)
# run_individual(eq_num = 113, dphase = 'ScP'  , start_beam_stack =   -3, end_beam_stack =   15, R_slow_plot = 0.035, T_slow_plot = 0.000, fig_index = 200, stat_corr = 1, apply_SNR = True)
# run_individual(eq_num = 113, dphase = 'PKiKP', start_beam_stack =   -3, end_beam_stack =   15, R_slow_plot = 0.003, T_slow_plot = 0.000, fig_index = 300, stat_corr = 1, apply_SNR = True)

# 113
# run_individual(eq_num = 113, dphase = 'PcP'  , start_beam_stack =   -2, end_beam_stack =    2, R_slow_plot = 0.013, T_slow_plot = 0.000, fig_index = 100, stat_corr = 1, apply_SNR = True)
# run_individual(eq_num = 113, dphase = 'ScP'  , start_beam_stack =   -1, end_beam_stack =    1, R_slow_plot = 0.014, T_slow_plot = 0.000, fig_index = 200, stat_corr = 1, apply_SNR = True)
# run_individual(eq_num = 113, dphase = 'PKiKP', start_beam_stack =    0, end_beam_stack =    2, R_slow_plot = 0.003, T_slow_plot = 0.000, fig_index = 300, stat_corr = 1, apply_SNR = True)

# 114
run_individual(eq_num = 114, dphase = 'PcP'  , start_beam_stack =    0, end_beam_stack =    4, R_slow_plot = -0.0105, T_slow_plot = 0.0070, fig_index = 100, stat_corr = 1, apply_SNR = True)
# run_individual(eq_num = 114, dphase = 'ScP'  , start_beam_stack =   -2, end_beam_stack =    5, R_slow_plot = 0.014, T_slow_plot = 0.000, fig_index = 200, stat_corr = 1, apply_SNR = True)
# run_individual(eq_num = 114, dphase = 'PKiKP', start_beam_stack =  0.5, end_beam_stack =    9, R_slow_plot = 0.003, T_slow_plot = 0.000, fig_index = 300, stat_corr = 1, apply_SNR = True)

# 115
# run_individual(eq_num = 115, dphase = 'PcP'  , start_beam_stack =    2, end_beam_stack =    8, R_slow_plot = 0.013, T_slow_plot = 0.000, fig_index = 100, stat_corr = 1, apply_SNR = True)
# run_individual(eq_num = 115, dphase = 'ScP'  , start_beam_stack =    2, end_beam_stack =    7, R_slow_plot = 0.014, T_slow_plot = 0.000, fig_index = 200, stat_corr = 1, apply_SNR = True)
# run_individual(eq_num = 115, dphase = 'PKiKP', start_beam_stack =   -2, end_beam_stack =    3, R_slow_plot = 0.003, T_slow_plot = 0.000, fig_index = 300, stat_corr = 1, apply_SNR = True)

# 116
# run_individual(eq_num = 116, dphase = 'PcP'  , start_beam_stack =    0, end_beam_stack =   10, R_slow_plot = 0.013, T_slow_plot = 0.000, fig_index = 100, stat_corr = 1, apply_SNR = True)
# run_individual(eq_num = 116, dphase = 'ScP'  , start_beam_stack =    0, end_beam_stack =   10, R_slow_plot = 0.014, T_slow_plot = 0.000, fig_index = 200, stat_corr = 1, apply_SNR = True)
# run_individual(eq_num = 116, dphase = 'PKiKP', start_beam_stack =   -2, end_beam_stack =   10, R_slow_plot = 0.003, T_slow_plot = 0.000, fig_index = 300, stat_corr = 1, apply_SNR = True)

# 117
# run_individual(eq_num = 117, dphase = 'PcP'  , start_beam_stack =    8, end_beam_stack =   15, R_slow_plot = 0.013, T_slow_plot = 0.000, fig_index = 100, stat_corr = 1, apply_SNR = True)
# run_individual(eq_num = 117, dphase = 'ScP'  , start_beam_stack =    0, end_beam_stack =   10, R_slow_plot = 0.014, T_slow_plot = 0.000, fig_index = 200, stat_corr = 1, apply_SNR = True)
# run_individual(eq_num = 117, dphase = 'PKiKP', start_beam_stack =    5, end_beam_stack =   10, R_slow_plot = 0.003, T_slow_plot = 0.000, fig_index = 300, stat_corr = 1, apply_SNR = True)

# 118
# run_individual(eq_num = 118, dphase = 'PcP'  , start_beam_stack =   -4, end_beam_stack =    3, R_slow_plot = 0.013, T_slow_plot = 0.000, fig_index = 100, stat_corr = 1, apply_SNR = True)
# run_individual(eq_num = 118, dphase = 'ScP'  , start_beam_stack =   -4, end_beam_stack =    2, R_slow_plot = 0.014, T_slow_plot = 0.000, fig_index = 200, stat_corr = 1, apply_SNR = True)
# run_individual(eq_num = 118, dphase = 'PKiKP', start_beam_stack =   -2, end_beam_stack =    5, R_slow_plot = 0.003, T_slow_plot = 0.000, fig_index = 300, stat_corr = 1, apply_SNR = True)

# 119
# run_individual(eq_num = 119, dphase = 'PcP'  , start_beam_stack =    0, end_beam_stack =    5, R_slow_plot = 0.013, T_slow_plot = 0.000, fig_index = 100, stat_corr = 1, apply_SNR = True)
# run_individual(eq_num = 119, dphase = 'ScP'  , start_beam_stack =    2, end_beam_stack =    8, R_slow_plot = 0.014, T_slow_plot = 0.000, fig_index = 200, stat_corr = 1, apply_SNR = True)
# run_individual(eq_num = 119, dphase = 'PKiKP', start_beam_stack =    2, end_beam_stack =    8, R_slow_plot = 0.003, T_slow_plot = 0.000, fig_index = 300, stat_corr = 1, apply_SNR = True)

# 120
# run_individual(eq_num = 120, dphase = 'PcP'  , start_beam_stack =    2, end_beam_stack =   10, R_slow_plot = 0.013, T_slow_plot = 0.000, fig_index = 100, stat_corr = 1, apply_SNR = True)
# run_individual(eq_num = 120, dphase = 'ScP'  , start_beam_stack =   -2, end_beam_stack =   10, R_slow_plot = 0.014, T_slow_plot = 0.000, fig_index = 200, stat_corr = 1, apply_SNR = True)
# run_individual(eq_num = 120, dphase = 'PKiKP', start_beam_stack =    2, end_beam_stack =    6, R_slow_plot = 0.003, T_slow_plot = 0.000, fig_index = 300, stat_corr = 1, apply_SNR = True)


elapsed_time_wc = time.time() - start_time_wc
print(f'This job took {elapsed_time_wc:.1f} seconds')
os.system('say "All done"')