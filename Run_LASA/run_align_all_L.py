#!/usr/bin/env python
# John Vidale 4/2020

import os
os.environ['PATH'] += os.pathsep + '/usr/local/bin'
os.chdir('/Users/vidale/Documents/GitHub/Array_codes/Run_ICS')

#%% Import functions
from run_align_L      import run_get_shift_L

#%% explosion PcP precursors
# run_get_shift_L(
#     start_buff_align = -40, end_buff_align = 40, start_beam_align =  -4, end_beam_align = 3,
#     start_buff_stack = -30, end_buff_stack =  0, start_beam_stack = -30, end_beam_stack = -4, event_no = 1, dphase = 'PcP')
# run_get_shift_L(
#     start_buff_align = -40, end_buff_align = 40, start_beam_align =   0, end_beam_align = 5,
#     start_buff_stack = -30, end_buff_stack =  0, start_beam_stack = -30, end_beam_stack = -2, event_no = 4, dphase = 'PcP')
# run_get_shift_L(
#     start_buff_align = -40, end_buff_align = 40, start_beam_align =  -4, end_beam_align = 3,
#     start_buff_stack = -15, end_buff_stack =  0, start_beam_stack = -7, end_beam_stack = -4, event_no = 7, dphase = 'PcP')
# run_get_shift_L(
#     start_buff_align = -40, end_buff_align = 40, start_beam_align =  -4, end_beam_align = 15,
#     start_buff_stack = -15, end_buff_stack =  0, start_beam_stack = -7, end_beam_stack = -4, event_no = 8, dphase = 'PcP')

#%% earthquake PcP precursors
run_get_shift_L(
    start_buff_align = -40, end_buff_align = 40, start_beam_align =  -4, end_beam_align = 3,
    start_buff_stack = -30, end_buff_stack =  0, start_beam_stack = -30, end_beam_stack = -4, event_no = 16, dphase = 'PcP')
