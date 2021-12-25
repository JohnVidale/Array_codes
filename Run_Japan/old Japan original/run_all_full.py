#!/usr/bin/env python3
# John Vidale 4/2020

import os
# os.environ['PATH'] += os.pathsep + '/usr/local/bin'
# os.chdir('/Users/vidale/Documents/GitHub/Array_codes/Run_ICS')

#%% Import functions
from run_individual      import run_individual

#%% PKiKP alignment (more refined than in the PKiKP/ISC powerpoint)
# run_individual(event_no =102, dphase = 'PKiKP',
#     start_buff_stack = -10.0, end_buff_stack =  10.0,
#     start_beam_stack =  -0.5, end_beam_stack =   2.0)
run_individual(event_no =103, dphase = 'PKiKP',
    start_buff_stack = -10.0, end_buff_stack =  10.0,
    start_beam_stack =  -2.0, end_beam_stack =   1.0)
# run_individual(event_no =104, dphase = 'PKiKP',
#     start_buff_stack = -10.0, end_buff_stack =  10.0,
#     start_beam_stack =  -2.0, end_beam_stack =   2.0)
# run_individual(event_no =105, dphase = 'PKiKP',
#     start_buff_stack = -10.0, end_buff_stack =  10.0,
#     start_beam_stack =  -2.0, end_beam_stack =   2.0)
# run_individual(event_no =106, dphase = 'PKiKP',
#     start_buff_stack = -10.0, end_buff_stack =  10.0,
#     start_beam_stack =  -2.0, end_beam_stack =   2.0)
# run_individual(event_no =107, dphase = 'PKiKP',
#     start_buff_stack = -10.0, end_buff_stack =  10.0,
#     start_beam_stack =  -2.0, end_beam_stack =   2.0)
# run_individual(event_no =108, dphase = 'PKiKP',
#     start_buff_stack = -10.0, end_buff_stack =  10.0,
#     start_beam_stack =  -2.0, end_beam_stack =   2.0)
# run_individual(event_no =109, dphase = 'PKiKP',
#     start_buff_stack = -10.0, end_buff_stack =  10.0,
#     start_beam_stack =  -2.0, end_beam_stack =   2.0)
# run_individual(event_no =110, dphase = 'PKiKP',
#     start_buff_stack = -10.0, end_buff_stack =  10.0,
#     start_beam_stack =  -2.0, end_beam_stack =   2.0)
# run_individual(event_no =111, dphase = 'PKiKP',
#     start_buff_stack = -10.0, end_buff_stack =  10.0,
#     start_beam_stack =  -2.5, end_beam_stack =  -0.5)

#%% PKiKP precursor
# run_get_shift_J(
#     start_buff_align = -10, end_buff_align = 10, start_beam_align = -9, end_beam_align = 9,
#     start_buff_stack = -20, end_buff_stack = 20, start_beam_stack = -10, end_beam_stack = 10, event_no = 102, dphase = 'PKiKP')
# run_get_shift_J(
#     start_buff_align = -50, end_buff_align = 10, start_beam_align = -9, end_beam_align = 9,
#     start_buff_stack = -20, end_buff_stack = 20, start_beam_stack = -10, end_beam_stack = 10, event_no = 103, dphase = 'PKiKP')
# run_get_shift_J(
#     start_buff_align = -10, end_buff_align = 10, start_beam_align = -9, end_beam_align = 9,
#     start_buff_stack = -20, end_buff_stack = 20, start_beam_stack = -10, end_beam_stack = 10, event_no = 104, dphase = 'PKiKP')
# run_get_shift_J(
#     start_buff_align = -10, end_buff_align = 10, start_beam_align = -9, end_beam_align = 9,
#     start_buff_stack = -20, end_buff_stack = 20, start_beam_stack = -10, end_beam_stack = 10, event_no = 105, dphase = 'PKiKP')
# run_get_shift_J(
#     start_buff_align = -50, end_buff_align = 10, start_beam_align = 0, end_beam_align = 8,
#     start_buff_stack = -20, end_buff_stack = 20, start_beam_stack = -10, end_beam_stack = 10, event_no = 106, dphase = 'PKiKP')
# run_get_shift_J(
#     start_buff_align = -10, end_buff_align = 10, start_beam_align = -9, end_beam_align = 9,
#     start_buff_stack = -20, end_buff_stack = 20, start_beam_stack = -10, end_beam_stack = 10, event_no = 107, dphase = 'PKiKP')
# run_get_shift_J(
#     start_buff_align = -10, end_buff_align = 10, start_beam_align = -9, end_beam_align = 9,
#     start_buff_stack = -20, end_buff_stack = 20, start_beam_stack = -10, end_beam_stack = 10, event_no = 108, dphase = 'PKiKP')
# run_get_shift_J(
#     start_buff_align = -10, end_buff_align = 10, start_beam_align = -9, end_beam_align = 9,
#     start_buff_stack = -20, end_buff_stack = 20, start_beam_stack = -10, end_beam_stack = 10, event_no = 109, dphase = 'PKiKP')
# run_get_shift_J(
#     start_buff_align = -10, end_buff_align = 10, start_beam_align = -9, end_beam_align = 9,
#     start_buff_stack = -20, end_buff_stack = 20, start_beam_stack = -10, end_beam_stack = 10, event_no = 110, dphase = 'PKiKP')
# run_get_shift_J(
#     start_buff_align = -10, end_buff_align = 10, start_beam_align = -9, end_beam_align = 9,
#     start_buff_stack = -20, end_buff_stack = 20, start_beam_stack = -10, end_beam_stack = 10, event_no = 111, dphase = 'PKiKP')

#%% PcP precursor
# run_get_shift_J(
#     start_buff_align = -10, end_buff_align = 10, start_beam_align = -9, end_beam_align = 9,
#     start_buff_stack = -20, end_buff_stack = 20, start_beam_stack = -10, end_beam_stack = 10, event_no = 102, dphase = 'PcP')
# run_get_shift_J(
#     start_buff_align = -40, end_buff_align = 10, start_beam_align = -3, end_beam_align = 6,
#     start_buff_stack = -20, end_buff_stack = -4, start_beam_stack = -20, end_beam_stack = -4, event_no = 103, dphase = 'PcP')
# run_get_shift_J(
#     start_buff_align = -10, end_buff_align = 10, start_beam_align = -9, end_beam_align = 9,
#     start_buff_stack = -20, end_buff_stack = 20, start_beam_stack = -10, end_beam_stack = 10, event_no = 104, dphase = 'PcP')
# run_get_shift_J(
#     start_buff_align = -10, end_buff_align = 10, start_beam_align = -9, end_beam_align = 9,
#     start_buff_stack = -20, end_buff_stack = 20, start_beam_stack = -10, end_beam_stack = 10, event_no = 105, dphase = 'PcP')
# run_get_shift_J(
#     start_buff_align = -10, end_buff_align = 10, start_beam_align = -9, end_beam_align = 9,
#     start_buff_stack = -20, end_buff_stack = 20, start_beam_stack = -10, end_beam_stack = 10, event_no = 106, dphase = 'PcP')
# run_get_shift_J(
#     start_buff_align = -10, end_buff_align = 10, start_beam_align = -9, end_beam_align = 9,
#     start_buff_stack = -20, end_buff_stack = 20, start_beam_stack = -10, end_beam_stack = 10, event_no = 107, dphase = 'PcP')
# run_get_shift_J(
#     start_buff_align = -10, end_buff_align = 10, start_beam_align = -9, end_beam_align = 9,
#     start_buff_stack = -20, end_buff_stack = 20, start_beam_stack = -10, end_beam_stack = 10, event_no = 108, dphase = 'PcP')
# run_get_shift_J(
#     start_buff_align = -10, end_buff_align = 10, start_beam_align = -9, end_beam_align = 9,
#     start_buff_stack = -20, end_buff_stack = 20, start_beam_stack = -10, end_beam_stack = 10, event_no = 109, dphase = 'PcP')
# run_get_shift_J(
#     start_buff_align = -10, end_buff_align = 10, start_beam_align = -9, end_beam_align = 9,
#     start_buff_stack = -20, end_buff_stack = 20, start_beam_stack = -10, end_beam_stack = 10, event_no = 110, dphase = 'PcP')
# run_get_shift_J(
#     start_buff_align = -10, end_buff_align = 10, start_beam_align = -9, end_beam_align = 9,
#     start_buff_stack = -20, end_buff_stack = 20, start_beam_stack = -10, end_beam_stack = 10, event_no = 111, dphase = 'PcP')

#%% P precursor
# run_get_shift_J(
#     start_buff_align = -10, end_buff_align = 10, start_beam_align = -9, end_beam_align = 9,
#     start_buff_stack = -20, end_buff_stack = 20, start_beam_stack = -10, end_beam_stack = 10, event_no = 102, dphase = 'P')
# run_get_shift_J(
#     start_buff_align = -10, end_buff_align = 10, start_beam_align = -9, end_beam_align = 9,
#     start_buff_stack = -20, end_buff_stack = 20, start_beam_stack = -10, end_beam_stack = 10, event_no = 103, dphase = 'P')
# run_get_shift_J(
#     start_buff_align = -10, end_buff_align = 10, start_beam_align = -9, end_beam_align = 9,
#     start_buff_stack = -20, end_buff_stack = 20, start_beam_stack = -10, end_beam_stack = 10, event_no = 104, dphase = 'P')
# run_get_shift_J(
#     start_buff_align = -10, end_buff_align = 10, start_beam_align = -9, end_beam_align = 9,
#     start_buff_stack = -20, end_buff_stack = 20, start_beam_stack = -10, end_beam_stack = 10, event_no = 105, dphase = 'P')
# run_get_shift_J(
#     start_buff_align = -10, end_buff_align = 10, start_beam_align = -9, end_beam_align = 9,
#     start_buff_stack = -20, end_buff_stack = 20, start_beam_stack = -10, end_beam_stack = 10, event_no = 106, dphase = 'P')
# run_get_shift_J(
#     start_buff_align = -10, end_buff_align = 10, start_beam_align = -9, end_beam_align = 9,
#     start_buff_stack = -20, end_buff_stack = 20, start_beam_stack = -10, end_beam_stack = 10, event_no = 107, dphase = 'P')
# run_get_shift_J(
#     start_buff_align = -10, end_buff_align = 10, start_beam_align = -9, end_beam_align = 9,
#     start_buff_stack = -20, end_buff_stack = 20, start_beam_stack = -10, end_beam_stack = 10, event_no = 108, dphase = 'P')
# run_get_shift_J(
#     start_buff_align = -10, end_buff_align = 10, start_beam_align = -9, end_beam_align = 9,
#     start_buff_stack = -20, end_buff_stack = 20, start_beam_stack = -10, end_beam_stack = 10, event_no = 109, dphase = 'P')
# run_get_shift_J(
#     start_buff_align = -10, end_buff_align = 10, start_beam_align = -9, end_beam_align = 9,
#     start_buff_stack = -20, end_buff_stack = 20, start_beam_stack = -10, end_beam_stack = 10, event_no = 110, dphase = 'P')
# run_get_shift_J(
#     start_buff_align = -10, end_buff_align = 10, start_beam_align = -9, end_beam_align = 9,
#     start_buff_stack = -20, end_buff_stack = 20, start_beam_stack = -10, end_beam_stack = 10, event_no = 111, dphase = 'P')

