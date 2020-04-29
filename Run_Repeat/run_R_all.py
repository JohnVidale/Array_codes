#!/usr/bin/env python3
# John Vidale 4/2020

import os
os.environ['PATH'] += os.pathsep + '/usr/local/bin'
os.chdir('/Users/vidale/Documents/GitHub/Array_codes/Run_Repeat')

#%% Import functions
from run_R_each_duet import runR_each_duet
from run_R_each_solo import runR_each_solo

#runR_each_solo(start_buff = -40,   end_buff = 45, start_beam =  -2,   end_beam = 5,   event_no = 153, do_decimate = 0)

#runR_each_solo(start_buff = -10,   end_buff = 80, start_beam =  20,   end_beam = 60,   event_no = 8)

runR_each_duet(start_buff = -20,   end_buff = 40, start_beam = -2,   end_beam = 1,   event1_no = 7,   event2_no = 8)
#runR_each_duet(start_buff = -10,   end_buff = 80, start_beam =  20,   end_beam = 60,   event1_no = 7,   event2_no = 8)
