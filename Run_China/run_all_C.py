#!/usr/bin/env python3
# John Vidale 4/2020

import os
os.environ['PATH'] += os.pathsep + '/usr/local/bin'
os.chdir('/Users/vidale/Documents/GitHub/Array_codes/Run_ICS')

#%% Import functions
from run_align_C      import run_get_shift_C

#run_each_C_PKiKP(start_beam =  0.0, end_beam = 10.0, start_buff = -20,  end_buff = 20,   event_no = 200, dphase = 'PKiKP')
#run_each_C_PKiKP(start_beam =  0.0, end_beam = 10.0, start_buff = -20,  end_buff = 20,   event_no = 201, dphase = 'PKiKP')
#run_each_C_PKiKP(start_beam =  0.0, end_beam = 10.0, start_buff = -20,  end_buff = 20,   event_no = 202, dphase = 'PKiKP') # a little PKiKP
#run_each_C_PKiKP(start_beam =  4.0, end_beam =  7.0, start_buff = -20,  end_buff = 20,   event_no = 203, dphase = 'PKiKP') # a little PKiKP
#run_each_C_PKiKP(start_beam =  1.5, end_beam =  4.0, start_buff = -20,  end_buff = 20,   event_no = 204, dphase = 'PKiKP')
#run_each_C_PKiKP(start_beam =  2.0, end_beam =  4.3, start_buff = -20,  end_buff = 20,   event_no = 205, dphase = 'PKiKP') # good PKiKP
run_get_shift_C(
    start_buff_align = -10, end_buff_align = 30, start_beam_align = 1, end_beam_align = 5,
    start_buff_stack =   0, end_buff_stack = 80, start_beam_stack = 20, end_beam_stack = 60, event_no = 205, dphase = 'PKiKP')
#run_each_C_PKiKP(start_beam =  0.0, end_beam =  3.0, start_buff = -20,  end_buff = 20,   event_no = 206, dphase = 'PKiKP') # moderate PKiKP
#run_each_C_PKiKP(start_beam =  3.0, end_beam =  7.0, start_buff = -20,  end_buff = 20,   event_no = 207, dphase = 'PKiKP') # moderate PKiKP
#run_each_C_PKiKP(start_beam =  0.0, end_beam = 10.0, start_buff = -20,  end_buff = 20,   event_no = 208, dphase = 'PKiKP')
#run_each_C_PKiKP(start_beam =  4.0, end_beam =  8.0, start_buff = -20,  end_buff = 20,   event_no = 209, dphase = 'PKiKP') # moderate PKiKP
#run_each_C_PKiKP(start_beam =  8.0, end_beam = 13.0, start_buff = -20,  end_buff = 20,   event_no = 210, dphase = 'PKiKP')
#run_each_C_PKiKP(start_beam =  0.0, end_beam = 10,   start_buff = -20,  end_buff = 20,   event_no = 211, dphase = 'PKiKP')
#run_each_C_PKiKP(start_beam =  0.0, end_beam = 10,   start_buff = -20,  end_buff = 20,   event_no = 212, dphase = 'PKiKP')
#run_each_C_PKiKP(start_beam =  0.0, end_beam = 10,   start_buff = -20,  end_buff = 20,   event_no = 213, dphase = 'PKiKP')
#run_each_C_PKiKP(start_beam =  0.0, end_beam = 10,   start_buff = -20,  end_buff = 20,   event_no = 214, dphase = 'PKiKP')
#run_each_C_PKiKP(start_beam =  7.5, end_beam =  9.5, start_buff = -20,  end_buff = 20,   event_no = 215, dphase = 'PKiKP') # a little PKiKP
#run_each_C_PKiKP(start_beam =  1.5, end_beam =  4.5, start_buff = -20,  end_buff = 20,   event_no = 216, dphase = 'PKiKP')
#run_each_C_PKiKP(start_beam =  0.0, end_beam =  5,   start_buff = -20,  end_buff = 20,   event_no = 217, dphase = 'PKiKP') # moderate PKiKP
#run_each_C_PKiKP(start_beam =  0.0, end_beam = 10,   start_buff = -20,  end_buff = 20,   event_no = 218, dphase = 'PKiKP')
#run_each_C_PKiKP(start_beam =  0.0, end_beam = 10,   start_buff = -20,  end_buff = 20,   event_no = 219, dphase = 'PKiKP')
#run_each_C_PKiKP(start_beam =  0.0, end_beam = 10,   start_buff = -20,  end_buff = 20,   event_no = 220, dphase = 'PKiKP')
#run_each_C_PKiKP(start_beam =  0.0, end_beam =  5,   start_buff = -20,  end_buff = 20,   event_no = 221, dphase = 'PKiKP')
#run_each_C_PKiKP(start_beam =  0.0, end_beam = 10,   start_buff = -20,  end_buff = 20,   event_no = 222, dphase = 'PKiKP') # decent PKiKP
#run_each_C_PKiKP(start_beam =  0.0, end_beam = 10,   start_buff = -20,  end_buff = 20,   event_no = 223, dphase = 'PKiKP')
#run_each_C_PKiKP(start_beam =  0.0, end_beam = 10,   start_buff = -20,  end_buff = 20,   event_no = 224, dphase = 'PKiKP')
#run_each_C_PKiKP(start_beam =  0.0, end_beam = 10,   start_buff = -20,  end_buff = 20,   event_no = 225, dphase = 'PKiKP')
#run_each_C_PKiKP(start_beam =  0.0, end_beam = 10,   start_buff = -20,  end_buff = 20,   event_no = 226, dphase = 'PKiKP')
#run_each_C_PKiKP(start_beam =  0.0, end_beam = 10,   start_buff = -20,  end_buff = 20,   event_no = 227, dphase = 'PKiKP')

#run_each_C_PcP( start_beam =  0.0, end_beam = 10.0, start_buff = -20,  end_buff = 20,   event_no = 200, dphase = 'PcP')
#run_each_C_PcP( start_beam =  0.0, end_beam = 10.0, start_buff = -20,  end_buff = 20,   event_no = 201, dphase = 'PcP')
#run_each_C_PcP( start_beam =  0.0, end_beam = 10.0, start_buff = -20,  end_buff = 20,   event_no = 202, dphase = 'PcP')
#run_each_C_PcP( start_beam =  5.0, end_beam =  8.0, start_buff = -20,  end_buff = 20,   event_no = 203, dphase = 'PcP')
#run_each_C_PcP( start_beam =  0.0, end_beam = 10.0, start_buff = -20,  end_buff = 20,   event_no = 204, dphase = 'PcP')
#run_each_C_PcP( start_beam =  1.0, end_beam =  4.0, start_buff = -20,  end_buff = 20,   event_no = 205, dphase = 'PcP')
#run_each_C_PcP( start_beam =  0.0, end_beam =  2.0, start_buff = -20,  end_buff = 20,   event_no = 206, dphase = 'PcP')
#run_each_C_PcP( start_beam =  5.0, end_beam =  9.0, start_buff = -20,  end_buff = 20,   event_no = 207, dphase = 'PcP')
#run_each_C_PcP( start_beam =  3.0, end_beam =  6.0, start_buff = -20,  end_buff = 20,   event_no = 208, dphase = 'PcP')
#run_each_C_PcP( start_beam =  0.0, end_beam =  5.0, start_buff = -20,  end_buff = 20,   event_no = 209, dphase = 'PcP')
#run_each_C_PcP( start_beam =  0.0, end_beam = 10.0, start_buff = -20,  end_buff = 20,   event_no = 210, dphase = 'PcP')
#run_each_C_PcP( start_beam =  0.0, end_beam = 10.0, start_buff = -20,  end_buff = 20,   event_no = 211, dphase = 'PcP')
#run_each_C_PcP( start_beam =  0.0, end_beam =  3.0, start_buff = -20,  end_buff = 20,   event_no = 212, dphase = 'PcP')
#run_each_C_PcP( start_beam =  0.0, end_beam =  3.0, start_buff = -20,  end_buff = 20,   event_no = 213, dphase = 'PcP')
#run_each_C_PcP( start_beam =  0.0, end_beam =  5.0, start_buff = -20,  end_buff = 20,   event_no = 214, dphase = 'PcP')
#run_each_C_PcP( start_beam =  5.0, end_beam =  7.5, start_buff = -20,  end_buff = 20,   event_no = 215, dphase = 'PcP')
#run_each_C_PcP( start_beam =  3.0, end_beam =  6.0, start_buff = -20,  end_buff = 20,   event_no = 216, dphase = 'PcP')
#run_each_C_PcP( start_beam =  4.5, end_beam =  7.0, start_buff = -20,  end_buff = 20,   event_no = 217, dphase = 'PcP')
#run_each_C_PcP( start_beam =  0.0, end_beam =  2.5, start_buff = -20,  end_buff = 20,   event_no = 218, dphase = 'PcP')
#run_each_C_PcP( start_beam =  2.5, end_beam =  8.5, start_buff = -20,  end_buff = 20,   event_no = 219, dphase = 'PcP')
#run_each_C_PcP( start_beam =  0.0, end_beam = 10.0, start_buff = -20,  end_buff = 20,   event_no = 220, dphase = 'PcP')
#run_each_C_PcP( start_beam =  2.0, end_beam =  5.0, start_buff = -20,  end_buff = 20,   event_no = 221, dphase = 'PcP')
#run_each_C_PcP( start_beam =  8.0, end_beam = 11.0, start_buff = -20,  end_buff = 20,   event_no = 222, dphase = 'PcP')
#run_each_C_PcP( start_beam =  3.0, end_beam = 5.0, start_buff = -20,  end_buff = 20,   event_no = 223, dphase = 'PcP')
#run_each_C_PcP( start_beam =  0.0, end_beam = 10.0, start_buff = -20,  end_buff = 20,   event_no = 224, dphase = 'PcP')
#run_each_C_PcP( start_beam =  2.0, end_beam =  5.5, start_buff = -20,  end_buff = 20,   event_no = 225, dphase = 'PcP')
#run_each_C_PcP( start_beam =  0.0, end_beam = 10.0, start_buff = -20,  end_buff = 20,   event_no = 226, dphase = 'PcP')
#run_each_C_PcP( start_beam =  0.0, end_beam = 10.0, start_buff = -20,  end_buff = 20,   event_no = 227, dphase = 'PcP')

# run_get_shift_C(start_buff = -10, end_buff = 70, start_beam = 1, end_beam = 5, start_beam = 20, end_beam = 60, event_no = 205, dphase = 'PKiKP')
