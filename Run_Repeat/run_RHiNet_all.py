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

runR_each_duet(start_buff = -50,   end_buff =  200, start_beam = -25,   end_beam =  175,   event1_no = 1,   event2_no = 2, dphase = 'PKiKP')
runR_each_duet(start_buff = -50,   end_buff =  200, start_beam = -25,   end_beam =  175,   event1_no = 4,   event2_no = 5, dphase = 'PKiKP')
runR_each_duet(start_buff = -50,   end_buff =  200, start_beam = -25,   end_beam =  175,   event1_no = 7,   event2_no = 8, dphase = 'PKiKP')

#runR_each_duet(start_buff = -20,   end_buff =  40, start_beam =   0,   end_beam =  20,   event1_no = 4,   event2_no = 5, dphase = 'PKiKP')
#runR_each_duet(start_buff =   0,   end_buff =  60, start_beam =  20,   end_beam =  40,   event1_no = 4,   event2_no = 5, dphase = 'PKiKP')
#runR_each_duet(start_buff =  20,   end_buff =  80, start_beam =  40,   end_beam =  60,   event1_no = 4,   event2_no = 5, dphase = 'PKiKP')
#runR_each_duet(start_buff =  40,   end_buff = 100, start_beam =  60,   end_beam =  80,   event1_no = 4,   event2_no = 5, dphase = 'PKiKP')
#runR_each_duet(start_buff =  60,   end_buff = 120, start_beam =  80,   end_beam = 100,   event1_no = 4,   event2_no = 5, dphase = 'PKiKP')
#runR_each_duet(start_buff =  80,   end_buff = 140, start_beam = 100,   end_beam = 120,   event1_no = 4,   event2_no = 5, dphase = 'PKiKP')
#runR_each_duet(start_buff = 100,   end_buff = 160, start_beam = 120,   end_beam = 140,   event1_no = 4,   event2_no = 5, dphase = 'PKiKP')
#runR_each_duet(start_buff = 120,   end_buff = 180, start_beam = 140,   end_beam = 160,   event1_no = 4,   event2_no = 5, dphase = 'PKiKP')
#runR_each_duet(start_buff = 140,   end_buff = 200, start_beam = 160,   end_beam = 180,   event1_no = 4,   event2_no = 5, dphase = 'PKiKP')
#runR_each_duet(start_buff = 160,   end_buff = 220, start_beam = 180,   end_beam = 200,   event1_no = 4,   event2_no = 5, dphase = 'PKiKP')
#runR_each_duet(start_buff = 180,   end_buff = 240, start_beam = 200,   end_beam = 220,   event1_no = 4,   event2_no = 5, dphase = 'PKiKP')
#runR_each_duet(start_buff = 200,   end_buff = 260, start_beam = 220,   end_beam = 240,   event1_no = 4,   event2_no = 5, dphase = 'PKiKP')
#
#runR_each_duet(start_buff = -20,   end_buff =  40, start_beam =   0,   end_beam =  20,   event1_no = 7,   event2_no = 8, dphase = 'PKiKP')
#runR_each_duet(start_buff =   0,   end_buff =  60, start_beam =  20,   end_beam =  40,   event1_no = 7,   event2_no = 8, dphase = 'PKiKP')
#runR_each_duet(start_buff =  20,   end_buff =  80, start_beam =  40,   end_beam =  60,   event1_no = 7,   event2_no = 8, dphase = 'PKiKP')
#runR_each_duet(start_buff =  40,   end_buff = 100, start_beam =  60,   end_beam =  80,   event1_no = 7,   event2_no = 8, dphase = 'PKiKP')
#runR_each_duet(start_buff =  60,   end_buff = 120, start_beam =  80,   end_beam = 100,   event1_no = 7,   event2_no = 8, dphase = 'PKiKP')
#runR_each_duet(start_buff =  80,   end_buff = 140, start_beam = 100,   end_beam = 120,   event1_no = 7,   event2_no = 8, dphase = 'PKiKP')
#runR_each_duet(start_buff = 100,   end_buff = 160, start_beam = 120,   end_beam = 140,   event1_no = 7,   event2_no = 8, dphase = 'PKiKP')
#runR_each_duet(start_buff = 120,   end_buff = 180, start_beam = 140,   end_beam = 160,   event1_no = 7,   event2_no = 8, dphase = 'PKiKP')
#runR_each_duet(start_buff = 140,   end_buff = 200, start_beam = 160,   end_beam = 180,   event1_no = 7,   event2_no = 8, dphase = 'PKiKP')
#runR_each_duet(start_buff = 160,   end_buff = 220, start_beam = 180,   end_beam = 200,   event1_no = 7,   event2_no = 8, dphase = 'PKiKP')
#runR_each_duet(start_buff = 180,   end_buff = 240, start_beam = 200,   end_beam = 220,   event1_no = 7,   event2_no = 8, dphase = 'PKiKP')
#runR_each_duet(start_buff = 200,   end_buff = 260, start_beam = 220,   end_beam = 240,   event1_no = 7,   event2_no = 8, dphase = 'PKiKP')

#runR_each_duet(start_buff = -20,   end_buff =  30, start_beam =   0,   end_beam =  10,   event1_no = 4,   event2_no = 5, dphase = 'PKiKP')
#runR_each_duet(start_buff = -10,   end_buff =  40, start_beam =  10,   end_beam =  20,   event1_no = 4,   event2_no = 5, dphase = 'PKiKP')
#runR_each_duet(start_buff =   0,   end_buff =  40, start_beam =  20,   end_beam =  30,   event1_no = 4,   event2_no = 5, dphase = 'PKiKP')
#runR_each_duet(start_buff =  10,   end_buff =  60, start_beam =  30,   end_beam =  40,   event1_no = 4,   event2_no = 5, dphase = 'PKiKP')
#runR_each_duet(start_buff =  20,   end_buff =  70, start_beam =  40,   end_beam =  50,   event1_no = 4,   event2_no = 5, dphase = 'PKiKP')
#runR_each_duet(start_buff =  30,   end_buff =  80, start_beam =  50,   end_beam =  60,   event1_no = 4,   event2_no = 5, dphase = 'PKiKP')
#runR_each_duet(start_buff =  40,   end_buff =  90, start_beam =  60,   end_beam =  70,   event1_no = 4,   event2_no = 5, dphase = 'PKiKP')
#runR_each_duet(start_buff =  50,   end_buff = 100, start_beam =  70,   end_beam =  80,   event1_no = 4,   event2_no = 5, dphase = 'PKiKP')
#runR_each_duet(start_buff =  60,   end_buff = 110, start_beam =  80,   end_beam =  90,   event1_no = 4,   event2_no = 5, dphase = 'PKiKP')
#runR_each_duet(start_buff =  70,   end_buff = 120, start_beam =  90,   end_beam = 100,   event1_no = 4,   event2_no = 5, dphase = 'PKiKP')
#runR_each_duet(start_buff =  80,   end_buff = 130, start_beam = 100,   end_beam = 110,   event1_no = 4,   event2_no = 5, dphase = 'PKiKP')
#runR_each_duet(start_buff =  90,   end_buff = 140, start_beam = 110,   end_beam = 120,   event1_no = 4,   event2_no = 5, dphase = 'PKiKP')
#runR_each_duet(start_buff = 100,   end_buff = 150, start_beam = 120,   end_beam = 130,   event1_no = 4,   event2_no = 5, dphase = 'PKiKP')
#runR_each_duet(start_buff = 110,   end_buff = 160, start_beam = 130,   end_beam = 140,   event1_no = 4,   event2_no = 5, dphase = 'PKiKP')
#runR_each_duet(start_buff = 120,   end_buff = 170, start_beam = 140,   end_beam = 150,   event1_no = 4,   event2_no = 5, dphase = 'PKiKP')
#runR_each_duet(start_buff = 130,   end_buff = 180, start_beam = 150,   end_beam = 160,   event1_no = 4,   event2_no = 5, dphase = 'PKiKP')
#runR_each_duet(start_buff = 140,   end_buff = 190, start_beam = 160,   end_beam = 170,   event1_no = 4,   event2_no = 5, dphase = 'PKiKP')
#runR_each_duet(start_buff = 150,   end_buff = 200, start_beam = 170,   end_beam = 180,   event1_no = 4,   event2_no = 5, dphase = 'PKiKP')
#runR_each_duet(start_buff = 160,   end_buff = 210, start_beam = 180,   end_beam = 190,   event1_no = 4,   event2_no = 5, dphase = 'PKiKP')
#runR_each_duet(start_buff = 170,   end_buff = 220, start_beam = 190,   end_beam = 200,   event1_no = 4,   event2_no = 5, dphase = 'PKiKP')
#runR_each_duet(start_buff = 180,   end_buff = 230, start_beam = 200,   end_beam = 210,   event1_no = 4,   event2_no = 5, dphase = 'PKiKP')
#runR_each_duet(start_buff = 190,   end_buff = 240, start_beam = 210,   end_beam = 220,   event1_no = 4,   event2_no = 5, dphase = 'PKiKP')
#runR_each_duet(start_buff = 200,   end_buff = 250, start_beam = 220,   end_beam = 230,   event1_no = 4,   event2_no = 5, dphase = 'PKiKP')
#runR_each_duet(start_buff = 210,   end_buff = 260, start_beam = 230,   end_beam = 240,   event1_no = 4,   event2_no = 5, dphase = 'PKiKP')
