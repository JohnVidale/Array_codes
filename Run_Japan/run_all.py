#!/usr/bin/env python3
# John Vidale 4/2020

import os
import time

#%% Import functions
from run_individual      import run_individual
start_time_wc = time.time()

#%% ICS
# run_individual(event_no = 102, dphase = 'PKiKP', start_buff_stack = 0.0, end_buff_stack =  100.0, start_beam_stack =  20.0, end_beam_stack =  60.0)
# run_individual(event_no = 103, dphase = 'PKiKP', start_buff_stack = 0.0, end_buff_stack =  100.0, start_beam_stack =  20.0, end_beam_stack =  60.0)
# run_individual(event_no = 104, dphase = 'PKiKP', start_buff_stack = 0.0, end_buff_stack =  100.0, start_beam_stack =  20.0, end_beam_stack =  60.0)
# run_individual(event_no = 105, dphase = 'PKiKP', start_buff_stack = 0.0, end_buff_stack =  100.0, start_beam_stack =  20.0, end_beam_stack =  60.0)
# run_individual(event_no = 106, dphase = 'PKiKP', start_buff_stack = 0.0, end_buff_stack =  100.0, start_beam_stack =  20.0, end_beam_stack =  60.0)
# run_individual(event_no = 107, dphase = 'PKiKP', start_buff_stack = 0.0, end_buff_stack =  100.0, start_beam_stack =  20.0, end_beam_stack =  60.0)
# run_individual(event_no = 108, dphase = 'PKiKP', start_buff_stack = 0.0, end_buff_stack =  100.0, start_beam_stack =  20.0, end_beam_stack =  60.0)
# run_individual(event_no = 109, dphase = 'PKiKP', start_buff_stack = 0.0, end_buff_stack =  100.0, start_beam_stack =  20.0, end_beam_stack =  60.0)
# run_individual(event_no = 110, dphase = 'PKiKP', start_buff_stack = 0.0, end_buff_stack =  100.0, start_beam_stack =  20.0, end_beam_stack =  60.0)
# run_individual(event_no = 111, dphase = 'PKiKP', start_buff_stack = 0.0, end_buff_stack =  100.0, start_beam_stack =  20.0, end_beam_stack =  60.0)

# run_individual(event_no = 112, dphase = 'PKiKP', start_buff_stack = 0.0, end_buff_stack =  100.0, start_beam_stack =  20.0, end_beam_stack =  60.0)
# run_individual(event_no = 113, dphase = 'PKiKP', start_buff_stack = 0.0, end_buff_stack =  100.0, start_beam_stack =  20.0, end_beam_stack =  60.0)
# run_individual(event_no = 114, dphase = 'PKiKP', start_buff_stack = 0.0, end_buff_stack =  100.0, start_beam_stack =  20.0, end_beam_stack =  60.0)
# run_individual(event_no = 115, dphase = 'PKiKP', start_buff_stack = 0.0, end_buff_stack =  100.0, start_beam_stack =  20.0, end_beam_stack =  60.0)
# run_individual(event_no = 116, dphase = 'PKiKP', start_buff_stack = 0.0, end_buff_stack =  100.0, start_beam_stack =  20.0, end_beam_stack =  60.0)
# run_individual(event_no = 117, dphase = 'PKiKP', start_buff_stack = 0.0, end_buff_stack =  100.0, start_beam_stack =  20.0, end_beam_stack =  60.0)
# run_individual(event_no = 118, dphase = 'PKiKP', start_buff_stack = 0.0, end_buff_stack =  100.0, start_beam_stack =  20.0, end_beam_stack =  60.0)
# run_individual(event_no = 119, dphase = 'PKiKP', start_buff_stack = 0.0, end_buff_stack =  100.0, start_beam_stack =  20.0, end_beam_stack =  60.0)
# run_individual(event_no = 120, dphase = 'PKiKP', start_buff_stack = 0.0, end_buff_stack =  100.0, start_beam_stack =  20.0, end_beam_stack =  60.0)

#%% PKiKP
# run_individual(event_no = 101, dphase = 'PKiKP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =   0.0, end_beam_stack =   5.0)
# run_individual(event_no = 102, dphase = 'PKiKP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =  -1.0, end_beam_stack =   4.0)
# run_individual(event_no = 103, dphase = 'PKiKP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =  -4.0, end_beam_stack =   1.0)
# run_individual(event_no = 104, dphase = 'PKiKP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =  -2.0, end_beam_stack =   4.0)
# run_individual(event_no = 105, dphase = 'PKiKP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =   0.0, end_beam_stack =   4.0)
# run_individual(event_no = 106, dphase = 'PKiKP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =   2.0, end_beam_stack =   5.0)
# run_individual(event_no = 107, dphase = 'PKiKP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =   0.0, end_beam_stack =   3.0)
# run_individual(event_no = 108, dphase = 'PKiKP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =  -2.0, end_beam_stack =   3.0)
# run_individual(event_no = 109, dphase = 'PKiKP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =  -2.0, end_beam_stack =   4.0)
# run_individual(event_no = 110, dphase = 'PKiKP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =  -2.0, end_beam_stack =   4.0)
# run_individual(event_no = 111, dphase = 'PKiKP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =  -4.0, end_beam_stack =   1.0)
# run_individual(event_no = 150, dphase = 'PKIKP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =   2.0, end_beam_stack =   6.0, JST = False)
run_individual(event_no = 152, dphase = 'PKIKP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =   0.0, end_beam_stack =   5.0, JST = False)

# run_individual(event_no = 112, dphase = 'PKiKP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =   6.0, end_beam_stack =   8.0)
# run_individual(event_no = 113, dphase = 'PKiKP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =  11.0, end_beam_stack =  15.0)
# run_individual(event_no = 114, dphase = 'PKiKP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =   2.0, end_beam_stack =   8.0)
# run_individual(event_no = 115, dphase = 'PKiKP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =  -5.0, end_beam_stack =   0.0)
# run_individual(event_no = 116, dphase = 'PKiKP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =   1.0, end_beam_stack =   9.0)
# run_individual(event_no = 117, dphase = 'PKiKP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =   4.0, end_beam_stack =   9.0)
# run_individual(event_no = 118, dphase = 'PKiKP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =   2.0, end_beam_stack =  11.0)
# run_individual(event_no = 119, dphase = 'PKiKP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =   3.0, end_beam_stack =  11.0)
# run_individual(event_no = 120, dphase = 'PKiKP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =   3.0, end_beam_stack =   7.0)

#%% PcP
# run_individual(event_no = 101, dphase = 'PcP', start_buff_stack =  -5.0, end_buff_stack =  15.0, start_beam_stack =  2.0, end_beam_stack =  8.0)
# run_individual(event_no = 102, dphase = 'PcP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack = -15.0, end_beam_stack =  25.0)
# run_individual(event_no = 102, dphase = 'PcP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =   5.0, end_beam_stack =   8.0)
# run_individual(event_no = 103, dphase = 'PcP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =  -3.0, end_beam_stack =   0.0)
# run_individual(event_no = 104, dphase = 'PcP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =   2.0, end_beam_stack =   6.0)
# run_individual(event_no = 105, dphase = 'PcP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =   2.0, end_beam_stack =   7.0)
# run_individual(event_no = 106, dphase = 'PcP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =   4.0, end_beam_stack =   9.0)
# run_individual(event_no = 107, dphase = 'PcP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =   3.0, end_beam_stack =   7.0)
# run_individual(event_no = 108, dphase = 'PcP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =   3.0, end_beam_stack =   7.0)
# run_individual(event_no = 109, dphase = 'PcP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =   3.0, end_beam_stack =   7.0)
# run_individual(event_no = 110, dphase = 'PcP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =   2.0, end_beam_stack =   5.0)
# run_individual(event_no = 111, dphase = 'PcP', start_buff_stack   = -20.0, end_buff_stack = 30.0, start_beam_stack =  -6.0, end_beam_stack =  -2.0)

# run_individual(event_no = 112, dphase =   'PcP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =  12.0, end_beam_stack =  15.0)
# run_individual(event_no = 113, dphase =   'PcP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =   7.0, end_beam_stack =  11.0)
# run_individual(event_no = 114, dphase =   'PcP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =   5.0, end_beam_stack =   8.0)
# run_individual(event_no = 115, dphase =   'PcP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =  -6.0, end_beam_stack =  -1.0)
# run_individual(event_no = 116, dphase =   'PcP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =   1.0, end_beam_stack =   8.0)
# run_individual(event_no = 117, dphase =   'PcP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =   2.0, end_beam_stack =   9.0)
# run_individual(event_no = 118, dphase =   'PcP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =   6.0, end_beam_stack =  15.0)
# run_individual(event_no = 119, dphase =   'PcP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =   9.0, end_beam_stack =  15.0)
# run_individual(event_no = 120, dphase =   'PcP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =   6.0, end_beam_stack =  18.0)

elapsed_time_wc = time.time() - start_time_wc
print(f'This job took {elapsed_time_wc:.1f} seconds')
os.system('say "All done"')