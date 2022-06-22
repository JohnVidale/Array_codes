#!/usr/bin/env python3
# John Vidale 4/2020

import os
import time
import matplotlib.pyplot as plt

#%% Import functions
from run_individual_little      import run_individual_little

start_time_wc = time.time()

plt.close('all')

eq_num = 529

# run_individual_little(eq_num = eq_num, dphase = 'PKiKP', start_buff_stack = -100, end_buff_stack = 200, start_beam_stack = 0, end_beam_stack = 100, beam_width = 0.04, beam_step = 0.0025, freq_min = 1, freq_max = 3, fig_index = 100, stat_corr = 0)
# run_individual_little(eq_num = eq_num, dphase = 'PKiKP', start_buff_stack = -100, end_buff_stack = 200, start_beam_stack = 0, end_beam_stack = 100, beam_width = 0.04, beam_step = 0.0025, freq_min = 2, freq_max = 4, fig_index = 200, stat_corr = 0)
# run_individual_little(eq_num = eq_num, dphase = 'PKiKP', start_buff_stack = -100, end_buff_stack = 200, start_beam_stack = 0, end_beam_stack = 100, beam_width = 0.04, beam_step = 0.0025, freq_min = 3, freq_max = 5, fig_index = 300, stat_corr = 0)

# PKiKP run
# run_individual_little(eq_num = eq_num, dphase = 'PKiKP', start_buff_stack = -20, end_buff_stack = 30, start_beam_stack = 10, end_beam_stack = 11, beam_width = 0.04, beam_step = 0.0025, freq_min = 1, freq_max = 3, fig_index = 100, stat_corr = 0)
# run_individual_little(eq_num = 555, dphase = 'PKiKP', start_buff_stack = -10, end_buff_stack = 40, start_beam_stack = -2, end_beam_stack = 3, beam_width = 0.04, beam_step = 0.0025, beam_offset = 0.0, freq_min = 2, freq_max = 4, stat_corr = 0)
run_individual_little(eq_num = 529, dphase = 'PKiKP', start_buff_stack = -10, end_buff_stack = 40, start_beam_stack = -1, end_beam_stack = 1, beam_width = 0.04, beam_step = 0.0025, beam_offset = 0.0, freq_min = 2, freq_max = 4, stat_corr = 0)

# PKiKP runs WRA
# run_individual_little(eq_num = 401, dphase = 'PKiKP', start_buff_stack = -50, end_buff_stack = 50, start_beam_stack = 1, end_beam_stack = 2, beam_width = 0.04, beam_step = 0.0025, freq_min = 2, freq_max = 4)
# run_individual_little(eq_num = 401, dphase = 'PKiKP', start_buff_stack = -50, end_buff_stack = 50, start_beam_stack = 4, end_beam_stack = 5.5, beam_width = 0.04, beam_step = 0.0025, freq_min = 2, freq_max = 4)
# run_individual_little(eq_num = 402, dphase = 'PKiKP', start_buff_stack = -50, end_buff_stack = 50, start_beam_stack = -2, end_beam_stack = 1, beam_width = 0.04, beam_step = 0.0025, freq_min = 2, freq_max = 4)
# run_individual_little(eq_num = 403, dphase = 'PKiKP', start_buff_stack = -10, end_buff_stack = 40, start_beam_stack = 3, end_beam_stack = 5, beam_width = 0.04, beam_step = 0.0025, beam_offset = 0.0, freq_min = 2, freq_max = 4, stat_corr = 0)
# 404 - no
# 405 - no
# run_individual_little(eq_num = 406, dphase = 'PKiKP', start_buff_stack = -50, end_buff_stack = 50, start_beam_stack = 0, end_beam_stack = 3, beam_width = 0.04, beam_step = 0.0025, freq_min = 4, freq_max = 6)
# 407 - no
# run_individual_little(eq_num = 408, dphase = 'PKiKP', start_buff_stack = -10, end_buff_stack = 40, start_beam_stack = 0, end_beam_stack = 2.5, beam_width = 0.04, beam_step = 0.0025, beam_offset = 0.0, freq_min = 3, freq_max = 5, stat_corr = 0)
# run_individual_little(eq_num = 409, dphase = 'PKiKP', start_buff_stack = -10, end_buff_stack = 40, start_beam_stack = 4.5, end_beam_stack = 5.5, beam_width = 0.04, beam_step = 0.0025, beam_offset = 0.0, freq_min = 2, freq_max = 4, stat_corr = 0)
# run_individual_little(eq_num = 410, dphase = 'PKiKP', start_buff_stack = -10, end_buff_stack = 40, start_beam_stack = 5, end_beam_stack = 10, beam_width = 0.04, beam_step = 0.0025, beam_offset = 0.0, freq_min = 1, freq_max = 3, stat_corr = 0)
# 411 - no
# run_individual_little(eq_num = 412, dphase = 'PKiKP', start_buff_stack = -10, end_buff_stack = 40, start_beam_stack = 0, end_beam_stack = 1, beam_width = 0.04, beam_step = 0.0025, beam_offset = 0.0, freq_min = 3, freq_max = 5, stat_corr = 0)

# PKiKP runs YKA
# BAD run_individual_little(eq_num = eq_num, dphase = 'PKiKP', start_buff_stack = -20, end_buff_stack = 30, start_beam_stack = 10, end_beam_stack = 11, beam_width = 0.04, beam_step = 0.0025, freq_min = 1, freq_max = 3, fig_index = 100, stat_corr = 0)
# run_individual_little(eq_num = 529, dphase = 'PKiKP', start_buff_stack = -10, end_buff_stack = 40, start_beam_stack = 0, end_beam_stack = 5, beam_width = 0.04, beam_step = 0.0025, beam_offset = 0.0, freq_min = 2, freq_max = 4, stat_corr = 0)
# 530 - no
# 531 - no
# run_individual_little(eq_num = 532, dphase = 'PKiKP', start_buff_stack = -10, end_buff_stack = 40, start_beam_stack = -2, end_beam_stack = 2, beam_width = 0.04, beam_step = 0.0025, beam_offset = 0.0, freq_min = 1, freq_max = 3, stat_corr = 0)
# run_individual_little(eq_num = 533, dphase = 'PKiKP', start_buff_stack = -10, end_buff_stack = 40, start_beam_stack = 0, end_beam_stack = 3, beam_width = 0.04, beam_step = 0.0025, beam_offset = 0.0, freq_min = 2, freq_max = 4, stat_corr = 0)
# 534 - no
# 535 - no
# 536 - no
# run_individual_little(eq_num = 537, dphase = 'PKiKP', start_buff_stack = -10, end_buff_stack = 40, start_beam_stack = -2, end_beam_stack = 3.5, beam_width = 0.04, beam_step = 0.0025, beam_offset = 0.0, freq_min = 3, freq_max = 5, stat_corr = 0)
# 538 - no
# 539 - no
# run_individual_little(eq_num = 540, dphase = 'PKiKP', start_buff_stack = -10, end_buff_stack = 40, start_beam_stack = 1, end_beam_stack = 3, beam_width = 0.04, beam_step = 0.0025, beam_offset = 0.0, freq_min = 3, freq_max = 5, stat_corr = 0)
# run_individual_little(eq_num = 541, dphase = 'PKiKP', start_buff_stack = -10, end_buff_stack = 40, start_beam_stack = 12, end_beam_stack = 16.5, beam_width = 0.04, beam_step = 0.0025, beam_offset = 0.0, freq_min = 3, freq_max = 5, stat_corr = 0)
# run_individual_little(eq_num = 542, dphase = 'PKiKP', start_buff_stack = -10, end_buff_stack = 40, start_beam_stack = 0, end_beam_stack = 15, beam_width = 0.04, beam_step = 0.0025, beam_offset = 0.0, freq_min = 3, freq_max = 5, stat_corr = 0)
# run_individual_little(eq_num = 543, dphase = 'PKiKP', start_buff_stack = -10, end_buff_stack = 40, start_beam_stack = 4, end_beam_stack = 6, beam_width = 0.04, beam_step = 0.0025, beam_offset = 0.0, freq_min = 2, freq_max = 4, stat_corr = 0)
# run_individual_little(eq_num = 544, dphase = 'PKiKP', start_buff_stack = -10, end_buff_stack = 40, start_beam_stack = 4, end_beam_stack = 6, beam_width = 0.04, beam_step = 0.0025, beam_offset = 0.0, freq_min = 3, freq_max = 5, stat_corr = 0)
# 545 - no
# run_individual_little(eq_num = 546, dphase = 'PKiKP', start_buff_stack = -10, end_buff_stack = 40, start_beam_stack = 1, end_beam_stack = 4, beam_width = 0.04, beam_step = 0.0025, beam_offset = 0.0, freq_min = 2, freq_max = 4, stat_corr = 0)
# run_individual_little(eq_num = 547, dphase = 'PKiKP', start_buff_stack = -10, end_buff_stack = 40, start_beam_stack = 4, end_beam_stack = 5, beam_width = 0.04, beam_step = 0.0025, beam_offset = 0.0, freq_min = 1, freq_max = 3, stat_corr = 0)
# 548 - no
# run_individual_little(eq_num = 549, dphase = 'PKiKP', start_buff_stack = -10, end_buff_stack = 40, start_beam_stack = -1, end_beam_stack = 1, beam_width = 0.04, beam_step = 0.0025, beam_offset = 0.0, freq_min = 2, freq_max = 4, stat_corr = 0)
# run_individual_little(eq_num = 551, dphase = 'PKiKP', start_buff_stack = -10, end_buff_stack = 40, start_beam_stack = 12, end_beam_stack = 15, beam_width = 0.04, beam_step = 0.0025, beam_offset = 0.0, freq_min = 1, freq_max = 3, stat_corr = 0)
# run_individual_little(eq_num = 552, dphase = 'PKiKP', start_buff_stack = -10, end_buff_stack = 40, start_beam_stack = 1, end_beam_stack = 3, beam_width = 0.04, beam_step = 0.0025, beam_offset = 0.0, freq_min = 1, freq_max = 3, stat_corr = 0)
# 553 - no sign of PKiKP
# run_individual_little(eq_num = 554, dphase = 'PKiKP', start_buff_stack = -10, end_buff_stack = 40, start_beam_stack = -2, end_beam_stack = 1, beam_width = 0.04, beam_step = 0.0025, beam_offset = 0.0, freq_min = 3, freq_max = 5, stat_corr = 0)
# run_individual_little(eq_num = 555, dphase = 'PKiKP', start_buff_stack = -10, end_buff_stack = 40, start_beam_stack = -2, end_beam_stack = 3, beam_width = 0.04, beam_step = 0.0025, beam_offset = 0.0, freq_min = 2, freq_max = 4, stat_corr = 0)
# run_individual_little(eq_num = 556, dphase = 'PKiKP', start_buff_stack = -100, end_buff_stack = 200, start_beam_stack = 0, end_beam_stack = 3, beam_width = 0.04, beam_step = 0.0025, freq_min = 3, freq_max = 5)
# run_individual_little(eq_num = 557, dphase = 'PKiKP', start_buff_stack = -100, end_buff_stack = 200, start_beam_stack = -1.5, end_beam_stack = 1, beam_width = 0.04, beam_step = 0.0025, freq_min = 1, freq_max = 3)
# 558 - no sign of PKiKP
# run_individual_little(eq_num = 559, dphase = 'PKiKP', start_buff_stack = -100, end_buff_stack = 200, start_beam_stack = 0, end_beam_stack = 100, beam_width = 0.04, beam_step = 0.0025, freq_min = 1, freq_max = 3)
# run_individual_little(eq_num = 559, dphase = 'PKiKP', start_buff_stack = -50, end_buff_stack = 50, start_beam_stack = 1, end_beam_stack = 3, beam_width = 0.04, beam_step = 0.0025, freq_min = 2, freq_max = 4)

elapsed_time_wc = time.time() - start_time_wc
print(f'This job took {elapsed_time_wc:.1f} seconds')
os.system('say "All done"')
