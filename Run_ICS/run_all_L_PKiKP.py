#!/usr/bin/env python3
# John Vidale 4/2020

import os
os.environ['PATH'] += os.pathsep + '/usr/local/bin'
os.chdir('/Users/vidale/Documents/GitHub/Array_codes/Run_ICS')

#%% Import functions
from run_each_L_PKiKP import run_each_L_PKiKP

#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam = -2,   end_beam =  1,   event_no = 12)
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam =  2,   end_beam = 11,   event_no = 13)
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam = -2,   end_beam =  7,   event_no = 15)
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam =  1,   end_beam =  6,   event_no = 16)
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam =  0,   end_beam = 20,   event_no = 17)
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam =  1,   end_beam =  6,   event_no = 19)
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam = -2,   end_beam =  1,   event_no = 18)
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam = -2,   end_beam =  4,   event_no = 14)
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam =  0,   end_beam = 20,   event_no = 20)
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam = -1.5, end_beam =  0.5, event_no = 21)
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam = -0.5, end_beam =  2.5, event_no = 23)
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam =  5,   end_beam = 12,   event_no = 22)
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam = 17,   end_beam = 21,   event_no = 24) # later window
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam =  1,   end_beam =  2.5, event_no = 25)
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam =  1,   end_beam =  5,   event_no = 26)
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam =  0,   end_beam = 20,   event_no = 27)
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam = -2,   end_beam =  0.5, event_no = 29)
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam = -3,   end_beam =  2,   event_no =  7)
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam = -2,   end_beam =  1,   event_no =  8)
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam = -2,   end_beam =  5,   event_no = 30)  # very weak PKiKP
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam = -0.5, end_beam =  1.5, event_no = 31)
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam = -5,   end_beam = 20,   event_no = 33)  # very weak PKiKP
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam = -2,   end_beam =  1.5, event_no = 34)
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam =  7,   end_beam = 12,   event_no = 32)
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam = -2.5, end_beam =  2,   event_no = 35)
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam =  0.0, end_beam = 20,   event_no = 36)  # very weak PKiKP
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam =  0.0, end_beam = 20,   event_no = 37)  # very weak PKiKP
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam =  9,   end_beam = 12,   event_no = 38)
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam = 11,   end_beam = 14,   event_no = 39)
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam = -1,   end_beam =  2.5, event_no = 40)
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam =  4,   end_beam =  7,   event_no = 41)
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam =  1.5, end_beam =  3.5, event_no = 42)
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam =  0,   end_beam = 20,   event_no = 1)  # no apparent PKiKP
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam =  0,   end_beam = 20,   event_no = 2)  # no apparent PKiKP
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam =  8,   end_beam = 12,   event_no = 43)
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam = -2,   end_beam =  4,   event_no = 44)  # no apparent PKiKP
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam = -0.5, end_beam =  2.5, event_no = 45)
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam =  0,   end_beam = 20,   event_no = 46)  # no apparent PKiKP
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam =  0,   end_beam = 10,   event_no = 4)
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam =  4,   end_beam =  8,   event_no = 5)
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam = -1,   end_beam = 10,   event_no = 47)
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam =  1,   end_beam = 10,   event_no = 48)  # no apparent PKiKP
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam =  2,   end_beam =  6,   event_no = 49)
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam =  1,   end_beam =  3.5, event_no = 50)
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam =  0,   end_beam =  6,   event_no = 51)
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam =  0,   end_beam = 10,   event_no = 52)
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam =  1.5, end_beam =  3.5, event_no = 53)
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam = 11,   end_beam = 14,   event_no = 54)
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam = -1,   end_beam =  2,   event_no = 55)
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam =  6,   end_beam =  8,   event_no = 57)
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam =  5,   end_beam = 10,   event_no = 58)
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam =  1,   end_beam =  4,   event_no = 59)
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam =  5,   end_beam = 10,   event_no = 60)
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam =  2,   end_beam =  9,   event_no = 61)
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam =  0,   end_beam = 10,   event_no = 63)
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam =  1,   end_beam =  3.5, event_no = 64)
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam = -0.5, end_beam =  2,   event_no = 65)
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam = 12,   end_beam = 18,   event_no = 66)
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam =  0,   end_beam =  3,   event_no = 67)
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam =  1,   end_beam =  3,   event_no = 68)
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam =  6,   end_beam =  9,   event_no = 69)
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam =  0,   end_beam =  5,   event_no = 70)  #  seems like other signals
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam =  0,   end_beam = 20,   event_no = 71)
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam =  4,   end_beam =  7,   event_no = 72)
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam = -0.5, end_beam =  4,   event_no = 73)
run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam =  1.5, end_beam =  7,   event_no = 74)
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam =  2,   end_beam =  4,   event_no = 75)
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam = 12,   end_beam = 21,   event_no = 76)
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam =  0,   end_beam = 20,   event_no = 77)
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam =  2,   end_beam =  5,   event_no = 79)
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam =  1,   end_beam =  8,   event_no = 80)
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam =  1,   end_beam =  5,   event_no = 81)
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam =  0,   end_beam = 10,   event_no = 83)

#  DUPES NOT YET PROCESSED AND NAMED
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam = XXXX, end_beam = XXXX, event_no = 56)
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam = XXXX, end_beam = XXXX, event_no = 62)
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam = XXXX, end_beam = XXXX, event_no = 78)
#run_each_L_PKiKP(start_buff = -20,  end_buff = 25, start_beam = XXXX, end_beam = XXXX, event_no = 82)