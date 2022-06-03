# !/usr/bin/env python3
# John Vidale 4/2020
# revisited   9/2021

import os
# os.environ['PATH'] += os.pathsep + '/usr/local/bin'
# os.chdir('/Users/vidale/Documents/GitHub/Array_codes/Run_ICS')

#%% Import functions
from run_individual_get_shifts      import run_individual_get_shifts

#%% PKiKP alignment
# run_individual_get_shifts(eq_num = 161, start_buff = -10, end_buff = 10, precursor_shift = -3, signal_dur = 5, min_dist = 148, max_dist = 155.43)
run_individual_get_shifts(eq_num = 101, start_buff = -10, end_buff = 10, precursor_shift = 1, signal_dur = 7, min_dist = 13.5, max_dist = 24.5)

# run_individual_get_shifts(
#     start_buff = -10, end_buff = 10, start_beam_align = -4, end_beam_align = 1, min_dist = 148, max_dist = 155,
#     eq_num = 161, dphase = 'PKiKP')

#%% PcP
# run_individual_get_shifts(
#     start_buff_align = -10, end_buff_align = 10, start_beam_align = -9, end_beam_align = 9,
#     start_buff_stack = -20, end_buff_stack = 20, start_beam_stack = -10, end_beam_stack = 10, eq_num = 111, dphase = 'PcP')