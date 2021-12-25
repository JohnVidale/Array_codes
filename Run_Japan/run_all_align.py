# !/usr/bin/env python3
# John Vidale 4/2020
# revisited   9/2021

import os
# os.environ['PATH'] += os.pathsep + '/usr/local/bin'
# os.chdir('/Users/vidale/Documents/GitHub/Array_codes/Run_ICS')

#%% Import functions
from run_individual_align      import run_individual_align

#%% PKiKP alignment (more refined than in the PKiKP/ISC powerpoint)
# run_individual_align(
#     start_buff_align = -10, end_buff_align = 10, start_beam_align = -4, end_beam_align = 1,
#     start_buff_stack =   0, end_buff_stack = 80, start_beam_stack = 20, end_beam_stack = 60, event_no = 111, dphase = 'PKiKP')

#%% PcP
run_individual_align(
    start_buff_align = -10, end_buff_align = 10, start_beam_align = -9, end_beam_align = 9,
    start_buff_stack = -20, end_buff_stack = 20, start_beam_stack = -10, end_beam_stack = 10, event_no = 111, dphase = 'PcP')