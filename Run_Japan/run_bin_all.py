#!/usr/bin/env python3
# John Vidale 4/2020

import os

#%% Import functions
from run_bin_one      import run_bin_one

#%% ICS
# run_bin_one(event_no = 101, dphase = 'PKiKP', start_buff_stack = -25.0, end_buff_stack = 25.0, JST = True)
# run_bin_one(event_no = 111, dphase = 'PKiKP', start_buff_stack = -25.0, end_buff_stack = 25.0, JST = True)
# run_bin_one(event_no = 101, dphase = 'PcP', start_buff_stack = -25.0, end_buff_stack = 25.0, JST = True)
# run_bin_one(event_no = 111, dphase = 'PcP', start_buff_stack = -25.0, end_buff_stack = 25.0, JST = True)
# run_bin_one(event_no = 152, dphase = 'PKIKP', start_buff_stack = -25.0, end_buff_stack = 25.0, JST = False)
# run_bin_one(event_no = 151, dphase = 'PKIKP', start_buff_stack = -25.0, end_buff_stack = 25.0, JST = False) # too emergent to see much
run_bin_one(event_no = 150, dphase = 'PKIKP', start_buff_stack = -25.0, end_buff_stack = 25.0, JST = False) # too emergent to see much
