#!/usr/bin/env python3
# John Vidale 4/2020

import os
os.environ['PATH'] += os.pathsep + '/usr/local/bin'
os.chdir('/Users/vidale/Documents/GitHub/Array_codes')

#%% Import functions
from run_each_J_PKiKP import run_eachPKiKP
from run_each_J_PcP import run_eachPcP
from run_each_J_ICS import run_eachICS

#run_eachPKiKP(start_buff =  1,   end_buff = 3,   event_no = 102)
#run_eachPKiKP(start_buff = -2,   end_buff = 0,   event_no = 103)
#run_eachPKiKP(start_buff = -1,   end_buff = 1,   event_no = 104)
#run_eachPKiKP(start_buff =  0,   end_buff = 2,   event_no = 105)
#run_eachPKiKP(start_buff =  2,   end_buff = 4,   event_no = 106)
#run_eachPKiKP(start_buff = -1,   end_buff = 1,   event_no = 107)
#run_eachPKiKP(start_buff =  0,   end_buff = 2,   event_no = 108)
#run_eachPKiKP(start_buff = -1,   end_buff = 1,   event_no = 109)
#run_eachPKiKP(start_buff = -0.5, end_buff = 1.5, event_no = 110)
#run_eachPKiKP(start_buff = -2,   end_buff = 0,   event_no = 111)

#run_eachICS(start_buff = 40,   end_buff = 100,   event_no = 102)
#run_eachICS(start_buff = 40,   end_buff = 100,   event_no = 103)
#run_eachICS(start_buff = 40,   end_buff = 100,   event_no = 104)
#run_eachICS(start_buff = 40,   end_buff = 100,   event_no = 105)
#run_eachICS(start_buff = 40,   end_buff = 100,   event_no = 106)
#run_eachICS(start_buff = 40,   end_buff = 100,   event_no = 107)
#run_eachICS(start_buff = 40,   end_buff = 100,   event_no = 108)
#run_eachICS(start_buff = 40,   end_buff = 100,   event_no = 109)
#run_eachICS(start_buff = 40,   end_buff = 100,   event_no = 110)
#run_eachICS(start_buff = -2,   end_buff = 0,   event_no = 111)

#run_eachPcP(start_buff =  1,    end_buff = 5,   event_no = 102, rel_time = 1)
#run_eachPcP(start_buff = -0.5,    end_buff = 3,   event_no = 103, rel_time = 1)
#run_eachPcP(start_buff = 0.5,    end_buff = 3.5,   event_no = 104, rel_time = 1)
#run_eachPcP(start_buff =  0,    end_buff = 3,   event_no = 105, rel_time = 1)
#run_eachPcP(start_buff =  4,  end_buff = 7, event_no = 106, rel_time = 1)
#run_eachPcP(start_buff =  0.5,   end_buff = 3.5,   event_no = 107, rel_time = 1)
#run_eachPcP(start_buff = -1,   end_buff = 3,   event_no = 108, rel_time = 1)
#run_eachPcP(start_buff = 1,   end_buff = 4,   event_no = 109, rel_time = 1)
#run_eachPcP(start_buff = 1,   end_buff = 4,   event_no = 110, rel_time = 1)
#run_eachPcP(start_buff = -1,   end_buff = 2,   event_no = 111, rel_time = 1)
