#!/usr/bin/env python3
# John Vidale 4/2020

import os
import time
import matplotlib.pyplot as plt

#%% Import functions
from run_individual_little      import run_individual_little
from run_individual_littleP      import run_individual_littleP

start_time_wc = time.time()

plt.close('all')

# eq_num = 412
eq_num = 746
shift_tt = 0

run_individual_little(eq_num = eq_num, dphase = 'PKiKP', start_buff_stack = -6, end_buff_stack = 2, start_beam_stack = -6, end_beam_stack = 2, beam_width = 0.03, beam_step = 0.0025, freq_min = 1, freq_max = 3, fig_index = 100, stat_corr = 0)
# run_individual_little(eq_num = eq_num, dphase = 'PKiKP', start_buff_stack = -50, end_buff_stack = 200, start_beam_stack = 10, end_beam_stack = 60, beam_width = 0.03, beam_step = 0.0025, freq_min = 2, freq_max = 4, fig_index = 200, stat_corr = 0)
# run_individual_little(eq_num = eq_num, dphase = 'PKiKP', start_buff_stack = -50, end_buff_stack = 200, start_beam_stack = 10, end_beam_stack = 60, beam_width = 0.03, beam_step = 0.0025, freq_min = 3, freq_max = 5, fig_index = 300, stat_corr = 0)
# plt.close('all')

# P run
# run_individual_littleP(eq_num = eq_num, dphase = 'P', start_buff_stack = -10, end_buff_stack = 10, start_beam_stack = -0.7, end_beam_stack = -0.4, freq_min = 3, freq_max = 5, stat_corr = 0, shift_tt = shift_tt, zerophase = False)

# for ii in range(546,560):
#     eq_num = ii

    # run_individual_little(eq_num = eq_num, dphase = 'PKiKP', start_buff_stack = -50, end_buff_stack = 200, start_beam_stack = 10, end_beam_stack = 60, beam_width = 0.03, beam_step = 0.0025, freq_min = 1, freq_max = 3, fig_index = 100, stat_corr = 0)
    # run_individual_little(eq_num = eq_num, dphase = 'PKiKP', start_buff_stack = -50, end_buff_stack = 200, start_beam_stack = 10, end_beam_stack = 60, beam_width = 0.03, beam_step = 0.0025, freq_min = 2, freq_max = 4, fig_index = 200, stat_corr = 0)
    # run_individual_little(eq_num = eq_num, dphase = 'PKiKP', start_buff_stack = -50, end_buff_stack = 200, start_beam_stack = 10, end_beam_stack = 60, beam_width = 0.03, beam_step = 0.0025, freq_min = 3, freq_max = 5, fig_index = 300, stat_corr = 0)
    # plt.close('all')

# fine tuning PKiKP run
# run_individual_little(eq_num = eq_num, dphase = 'PKiKP', start_buff_stack = -500, end_buff_stack = 0, start_beam_stack = -60, end_beam_stack = -40, beam_width = 0.1, beam_step = 0.05, beam_offset = 0.0, freq_min = 2, freq_max = 4, stat_corr = 0)
# run_individual_little(eq_num = eq_num, dphase = 'PKiKP', start_buff_stack = -20, end_buff_stack = 50, start_beam_stack = -1.2, end_beam_stack = -0.2, beam_width = 0.04, beam_step = 0.0025, beam_offset = 0.0, freq_min = 2, freq_max = 4, stat_corr = 0)
# run_individual_little(eq_num = eq_num, dphase = 'PKiKP', start_buff_stack = -150, end_buff_stack = 0, start_beam_stack = -42, end_beam_stack = -38, beam_width = 0.1, beam_step = 0.005, freq_min = 2, freq_max = 4, fig_index = 100, stat_corr = 0)

# funny western pattern
# run_individual_little(eq_num = 533, dphase = 'PKiKP', start_buff_stack = -50, end_buff_stack = 150, start_beam_stack = 0, end_beam_stack = 100, beam_width = 0.1, beam_step = 0.0025, freq_min = 1, freq_max = 5, fig_index = 100, stat_corr = 0, R_slow_plot = 0.03, T_slow_plot = 0.03)
# run_individual_little(eq_num = 411, dphase = 'PKiKP', start_buff_stack = -150, end_buff_stack = 50, start_beam_stack = -100, end_beam_stack = 0, beam_width = 0.6, beam_step = 0.01, freq_min = 2, freq_max = 3, fig_index = 100, stat_corr = 0, R_slow_plot = 0.03, T_slow_plot = 0.03)

# PKiKP runs WRA
# run_individual_little(eq_num = 401, dphase = 'PKiKP', start_buff_stack = -10, end_buff_stack = 40, start_beam_stack = -5, end_beam_stack = -2, beam_width = 0.03, beam_step = 0.0025, beam_offset = 0.0, freq_min = 1, freq_max = 3, stat_corr = 0)
# run_individual_little(eq_num = 402, dphase = 'PKiKP', start_buff_stack = -10, end_buff_stack = 40, start_beam_stack = -2, end_beam_stack = 1, beam_width = 0.03, beam_step = 0.0025, beam_offset = 0.0, freq_min = 2, freq_max = 4, stat_corr = 0)
# run_individual_little(eq_num = 403, dphase = 'PKiKP', start_buff_stack = -10, end_buff_stack = 40, start_beam_stack = 0, end_beam_stack = 3, beam_width = 0.03, beam_step = 0.0025, beam_offset = 0.0, freq_min = 1, freq_max = 3, stat_corr = 0)
# run_individual_little(eq_num = 404, dphase = 'PKiKP', start_buff_stack = -10, end_buff_stack = 40, start_beam_stack = 0, end_beam_stack = 3, beam_width = 0.03, beam_step = 0.0025, beam_offset = 0.0, freq_min = 1, freq_max = 3, stat_corr = 0)
# run_individual_little(eq_num = 405, dphase = 'PKiKP', start_buff_stack = -10, end_buff_stack = 40, start_beam_stack = 20, end_beam_stack = 30, beam_width = 0.03, beam_step = 0.0025, beam_offset = 0.0, freq_min = 3, freq_max = 5, stat_corr = 0)
# run_individual_little(eq_num = 406, dphase = 'PKiKP', start_buff_stack = -10, end_buff_stack = 40, start_beam_stack = 1, end_beam_stack = 2, beam_width = 0.02, beam_step = 0.001, beam_offset = 0.0, freq_min = 3, freq_max = 5, stat_corr = 0)
# run_individual_little(eq_num = 407, dphase = 'PKiKP', start_buff_stack = -10, end_buff_stack = 40, start_beam_stack = 0, end_beam_stack = 4, beam_width = 0.03, beam_step = 0.0025, beam_offset = 0.0, freq_min = 3, freq_max = 5, stat_corr = 0)
# run_individual_little(eq_num = 408, dphase = 'PKiKP', start_buff_stack = -10, end_buff_stack = 40, start_beam_stack = 0.5, end_beam_stack = 2.5, beam_width = 0.04, beam_step = 0.0025, beam_offset = 0.0, freq_min = 3, freq_max = 5, stat_corr = 0)
# run_individual_little(eq_num = 409, dphase = 'PKiKP', start_buff_stack = -10, end_buff_stack = 40, start_beam_stack = 4.5, end_beam_stack = 5.5, beam_width = 0.04, beam_step = 0.0025, beam_offset = 0.0, freq_min = 2, freq_max = 4, stat_corr = 0)
# run_individual_little(eq_num = 410, dphase = 'PKiKP', start_buff_stack = -10, end_buff_stack = 40, start_beam_stack = 1.2, end_beam_stack = 4.2, beam_width = 0.04, beam_step = 0.0025, beam_offset = 0.0, freq_min = 2, freq_max = 4, stat_corr = 0)
# run_individual_little(eq_num = 411, dphase = 'PKiKP', start_buff_stack = -20, end_buff_stack = 50, start_beam_stack = -1.2, end_beam_stack = -0.2, beam_width = 0.04, beam_step = 0.0025, beam_offset = 0.0, freq_min = 2, freq_max = 4, stat_corr = 0)
# run_individual_little(eq_num = 412, dphase = 'PKiKP', start_buff_stack = -10, end_buff_stack = 40, start_beam_stack = -0.5, end_beam_stack = 1, beam_width = 0.02, beam_step = 0.001, beam_offset = 0.0, freq_min = 1, freq_max = 3, stat_corr = 0)

# PKiKP runs YKA
# run_individual_little(eq_num = 528, dphase = 'PKiKP', start_buff_stack = -20, end_buff_stack = 30, start_beam_stack = 10, end_beam_stack = 11, beam_width = 0.04, beam_step = 0.0025, freq_min = 1, freq_max = 2, fig_index = 100, stat_corr = 0)
# run_individual_little(eq_num = 529, dphase = 'PKiKP', start_buff_stack = -10, end_buff_stack = 40, start_beam_stack = -1, end_beam_stack = 1, beam_width = 0.04, beam_step = 0.0025, beam_offset = 0.0, freq_min = 2, freq_max = 4, stat_corr = 0)
# run_individual_little(eq_num = 530, dphase = 'PKiKP', start_buff_stack = -10, end_buff_stack = 40, start_beam_stack = 0.5, end_beam_stack = 4, beam_width = 0.035, beam_step = 0.0025, beam_offset = 0.0, freq_min = 2, freq_max = 4, stat_corr = 0)
# run_individual_little(eq_num = 531, dphase = 'PKiKP', start_buff_stack = -10, end_buff_stack = 40, start_beam_stack = 13, end_beam_stack = 14, beam_width = 0.04, beam_step = 0.0025, beam_offset = 0.0, freq_min = 1, freq_max = 3, stat_corr = 0)
# run_individual_little(eq_num = 532, dphase = 'PKiKP', start_buff_stack = -10, end_buff_stack = 40, start_beam_stack = -3, end_beam_stack = 2, beam_width = 0.04, beam_step = 0.0025, beam_offset = 0.0, freq_min = 1, freq_max = 3, stat_corr = 0)
# run_individual_little(eq_num = 533, dphase = 'PKiKP', start_buff_stack = -10, end_buff_stack = 40, start_beam_stack = 0, end_beam_stack = 3, beam_width = 0.04, beam_step = 0.0025, beam_offset = 0.0, freq_min = 2, freq_max = 4, stat_corr = 0)
# run_individual_little(eq_num = 534, dphase = 'PKiKP', start_buff_stack = -10, end_buff_stack = 40, start_beam_stack = 7, end_beam_stack = 8, beam_width = 0.04, beam_step = 0.0025, beam_offset = 0.0, freq_min = 2, freq_max = 4, stat_corr = 0)
# run_individual_little(eq_num = 535, dphase = 'PKiKP', start_buff_stack = -10, end_buff_stack = 40, start_beam_stack = 12.5, end_beam_stack = 16, beam_width = 0.04, beam_step = 0.0025, beam_offset = 0.0, freq_min = 2, freq_max = 4, stat_corr = 0)
# BAD run_individual_little(eq_num = 536, dphase = 'PKiKP', start_buff_stack = -10, end_buff_stack = 40, start_beam_stack = 0, end_beam_stack = 3, beam_width = 0.04, beam_step = 0.0025, beam_offset = 0.0, freq_min = 1, freq_max = 3, stat_corr = 0)
# run_individual_little(eq_num = 537, dphase = 'PKiKP', start_buff_stack = -10, end_buff_stack = 40, start_beam_stack = -2, end_beam_stack = 3.5, beam_width = 0.04, beam_step = 0.0025, beam_offset = 0.0, freq_min = 3, freq_max = 5, stat_corr = 0)
# BAD run_individual_little(eq_num = 538, dphase = 'PKiKP', start_buff_stack = -10, end_buff_stack = 40, start_beam_stack = -2, end_beam_stack = 3.5, beam_width = 0.04, beam_step = 0.0025, beam_offset = 0.0, freq_min = 1, freq_max = 3, stat_corr = 0)
# run_individual_little(eq_num = 539, dphase = 'PKiKP', start_buff_stack = -10, end_buff_stack = 40, start_beam_stack = 4.5, end_beam_stack = 5.5, beam_width = 0.04, beam_step = 0.0025, beam_offset = 0.0, freq_min = 2, freq_max = 4, stat_corr = 0)
# run_individual_little(eq_num = 540, dphase = 'PKiKP', start_buff_stack = -10, end_buff_stack = 40, start_beam_stack = 1.5, end_beam_stack = 4.5, beam_width = 0.04, beam_step = 0.0025, beam_offset = 0.0, freq_min = 3, freq_max = 5, stat_corr = 0)
# run_individual_little(eq_num = 541, dphase = 'PKiKP', start_buff_stack = -10, end_buff_stack = 40, start_beam_stack = 13, end_beam_stack = 15, beam_width = 0.04, beam_step = 0.0025, beam_offset = 0.0, freq_min = 3, freq_max = 5, stat_corr = 0)
# run_individual_little(eq_num = 542, dphase = 'PKiKP', start_buff_stack = -10, end_buff_stack = 40, start_beam_stack = 1, end_beam_stack = 3, beam_width = 0.04, beam_step = 0.0025, beam_offset = 0.0, freq_min = 3, freq_max = 5, stat_corr = 0)
# run_individual_little(eq_num = 543, dphase = 'PKiKP', start_buff_stack = -10, end_buff_stack = 40, start_beam_stack = 4, end_beam_stack = 7, beam_width = 0.04, beam_step = 0.0025, beam_offset = 0.0, freq_min = 2, freq_max = 4, stat_corr = 0)
# run_individual_little(eq_num = 544, dphase = 'PKiKP', start_buff_stack = -10, end_buff_stack = 40, start_beam_stack = 3.5, end_beam_stack = 6.5, beam_width = 0.04, beam_step = 0.0025, beam_offset = 0.0, freq_min = 3, freq_max = 5, stat_corr = 0)
# run_individual_little(eq_num = 545, dphase = 'PKiKP', start_buff_stack = -10, end_buff_stack = 40, start_beam_stack = 8, end_beam_stack = 11, beam_width = 0.04, beam_step = 0.0025, beam_offset = 0.0, freq_min = 1, freq_max = 3, stat_corr = 0)
# run_individual_little(eq_num = 546, dphase = 'PKiKP', start_buff_stack = -10, end_buff_stack = 40, start_beam_stack = 1, end_beam_stack = 4, beam_width = 0.04, beam_step = 0.0025, beam_offset = 0.0, freq_min = 2, freq_max = 4, stat_corr = 0)
# run_individual_little(eq_num = 547, dphase = 'PKiKP', start_buff_stack = -10, end_buff_stack = 40, start_beam_stack = 9.5, end_beam_stack = 11.5, beam_width = 0.04, beam_step = 0.0025, beam_offset = 0.0, freq_min = 2, freq_max = 4, stat_corr = 0)
# run_individual_little(eq_num = 548, dphase = 'PKiKP', start_buff_stack = -10, end_buff_stack = 40, start_beam_stack = 12, end_beam_stack = 13, beam_width = 0.04, beam_step = 0.0025, beam_offset = 0.0, freq_min = 2, freq_max = 4, stat_corr = 0)
# run_individual_little(eq_num = 549, dphase = 'PKiKP', start_buff_stack = -10, end_buff_stack = 40, start_beam_stack = -1.5, end_beam_stack = 0.5, beam_width = 0.04, beam_step = 0.0025, beam_offset = 0.0, freq_min = 3, freq_max = 5, stat_corr = 0)
# run_individual_little(eq_num = 550, dphase = 'PKiKP', start_buff_stack = -10, end_buff_stack = 40, start_beam_stack = 3.5, end_beam_stack = 6.5, beam_width = 0.04, beam_step = 0.0025, beam_offset = 0.0, freq_min = 2, freq_max = 4, stat_corr = 0)
# run_individual_little(eq_num = 551, dphase = 'PKiKP', start_buff_stack = -10, end_buff_stack = 40, start_beam_stack = 9.5, end_beam_stack = 11, beam_width = 0.04, beam_step = 0.0025, beam_offset = 0.0, freq_min = 3, freq_max = 5, stat_corr = 0)
# run_individual_little(eq_num = 552, dphase = 'PKiKP', start_buff_stack = -10, end_buff_stack = 40, start_beam_stack = 5.5, end_beam_stack = 6.5, beam_width = 0.04, beam_step = 0.0025, beam_offset = 0.0, freq_min = 3, freq_max = 5, stat_corr = 0)
# run_individual_little(eq_num = 553, dphase = 'PKiKP', start_buff_stack = -10, end_buff_stack = 40, start_beam_stack = 1, end_beam_stack = 3, beam_width = 0.04, beam_step = 0.0025, beam_offset = 0.0, freq_min = 1, freq_max = 3, stat_corr = 0)
# run_individual_little(eq_num = 554, dphase = 'PKiKP', start_buff_stack = -10, end_buff_stack = 40, start_beam_stack = -2, end_beam_stack = 0, beam_width = 0.04, beam_step = 0.0025, beam_offset = 0.0, freq_min = 3, freq_max = 5, stat_corr = 0)
# run_individual_little(eq_num = 555, dphase = 'PKiKP', start_buff_stack = -10, end_buff_stack = 40, start_beam_stack = -1, end_beam_stack = 0, beam_width = 0.04, beam_step = 0.0025, beam_offset = 0.0, freq_min = 3, freq_max = 5, stat_corr = 0)
# run_individual_little(eq_num = 556, dphase = 'PKiKP', start_buff_stack = -10, end_buff_stack = 40, start_beam_stack = -1, end_beam_stack = 3, beam_width = 0.04, beam_step = 0.0025, beam_offset = 0.0, freq_min = 3, freq_max = 5, stat_corr = 0)
# run_individual_little(eq_num = 557, dphase = 'PKiKP', start_buff_stack = -10, end_buff_stack = 40, start_beam_stack = -2, end_beam_stack = 1.5, beam_width = 0.04, beam_step = 0.0025, beam_offset = 0.0, freq_min = 1, freq_max = 3, stat_corr = 0)
# run_individual_little(eq_num = 558, dphase = 'PKiKP', start_buff_stack = -10, end_buff_stack = 40, start_beam_stack = 4, end_beam_stack = 6.5, beam_width = 0.03, beam_step = 0.0025, beam_offset = 0.0, freq_min = 2, freq_max = 4, stat_corr = 0)
# run_individual_little(eq_num = 559, dphase = 'PKiKP', start_buff_stack = -10, end_buff_stack = 40, start_beam_stack = 3.5, end_beam_stack = 6.5, beam_width = 0.03, beam_step = 0.0025, beam_offset = 0.0, freq_min = 1, freq_max = 3, stat_corr = 0)

elapsed_time_wc = time.time() - start_time_wc
print(f'This job took {elapsed_time_wc:.1f} seconds')
os.system('say "All done"')
