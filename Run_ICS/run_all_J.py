#!/usr/bin/env python3
# John Vidale 4/2020

import os
os.environ['PATH'] += os.pathsep + '/usr/local/bin'
os.chdir('/Users/vidale/Documents/GitHub/Array_codes')

#%% Import functions
from run_each_J_PKiKP import run_eachPKiKP
from run_each_J_PcP import run_eachPcP
from run_each_J_ICS import run_eachICS

#run_eachPKiKP(start_buff = -20,  end_buff = 25, start_beam =  1,   end_beam = 3,   event_no = 102)
#run_eachPKiKP(start_buff = -20,  end_buff = 25, start_beam = -2,   end_beam = 0,   event_no = 103)
#run_eachPKiKP(start_buff = -20,  end_buff = 25, start_beam = -1,   end_beam = 1,   event_no = 104)
#run_eachPKiKP(start_buff = -20,  end_buff = 25, start_beam =  0,   end_beam = 2,   event_no = 105)
#run_eachPKiKP(start_buff = -20,  end_buff = 25, start_beam =  2,   end_beam = 4,   event_no = 106)
#run_eachPKiKP(start_buff = -20,  end_buff = 25, start_beam = -1,   end_beam = 1,   event_no = 107)
#run_eachPKiKP(start_buff = -20,  end_buff = 25, start_beam =  0,   end_beam = 2,   event_no = 108)
#run_eachPKiKP(start_buff = -20,  end_buff = 25, start_beam = -1,   end_beam = 1,   event_no = 109)
#run_eachPKiKP(start_buff = -20,  end_buff = 25, start_beam = -0.5, end_beam = 1.5, event_no = 110)
#run_eachPKiKP(start_buff = -20,  end_buff = 25, start_beam = -2,   end_beam = 0,   event_no = 111)

run_eachICS(start_buff = 10,   end_buff = 70, start_beam = 20,   end_beam = 60,   event_no = 102)
run_eachICS(start_buff = 10,   end_buff = 70, start_beam = 20,   end_beam = 60,   event_no = 103)
run_eachICS(start_buff = 10,   end_buff = 70, start_beam = 20,   end_beam = 60,   event_no = 104)
run_eachICS(start_buff = 10,   end_buff = 70, start_beam = 20,   end_beam = 60,   event_no = 105)
run_eachICS(start_buff = 10,   end_buff = 70, start_beam = 20,   end_beam = 60,   event_no = 106)
run_eachICS(start_buff = 10,   end_buff = 70, start_beam = 20,   end_beam = 60,   event_no = 107)
run_eachICS(start_buff = 10,   end_buff = 70, start_beam = 20,   end_beam = 60,   event_no = 108)
run_eachICS(start_buff = 10,   end_buff = 70, start_beam = 20,   end_beam = 60,   event_no = 109)
run_eachICS(start_buff = 10,   end_buff = 70, start_beam = 20,   end_beam = 60,   event_no = 110)
#run_eachICS(start_buff = 10,   end_buff = 70, start_beam = 20,   end_beam = 60,   event_no = 111)

#run_eachPcP(start_buff = -20,   end_buff = 25, start_beam =  1,   end_beam = 5,   event_no = 102)
#run_eachPcP(start_buff = -20,   end_buff = 25, start_beam = -0.5, end_beam = 3,   event_no = 103)
#run_eachPcP(start_buff = -20,   end_buff = 25, start_beam =  0.5, end_beam = 3.5, event_no = 104)
#run_eachPcP(start_buff = -20,   end_buff = 25, start_beam =  0,   end_beam = 3,   event_no = 105)
#run_eachPcP(start_buff = -20,   end_buff = 25, start_beam =  4,   end_beam = 7,   event_no = 106)
#run_eachPcP(start_buff = -20,   end_buff = 25, start_beam =  0.5, end_beam = 3.5, event_no = 107)
#run_eachPcP(start_buff = -20,   end_buff = 25, start_beam = -1,   end_beam = 3,   event_no = 108)
#run_eachPcP(start_buff = -20,   end_buff = 25, start_beam =  1,   end_beam = 4,   event_no = 109)
#run_eachPcP(start_buff = -20,   end_buff = 25, start_beam =  1,   end_beam = 4,   event_no = 110)
#run_eachPcP(start_buff = -20,   end_buff = 25, start_beam = -1,   end_beam = 2,   event_no = 111)
