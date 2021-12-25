#!/usr/bin/env python3
# John Vidale 4/2020
# revisited 9/2021

import os
# os.environ['PATH'] += os.pathsep + '/usr/local/bin'

#%% Import functions
from run_align_J      import run_get_shift_J

run_get_shift_J(
    start_buff_align = -10, end_buff_align = 10, start_beam_align = -4, end_beam_align = 1,
    start_buff_stack =   0, end_buff_stack = 80, start_beam_stack = 20, end_beam_stack = 60, event_no = 111, dphase = 'PKiKP')
