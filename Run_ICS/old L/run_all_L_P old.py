#!/usr/bin/env python
# input is set of LASA traces from NTS event
# This programs deals with a single event.
# this program tapers, filters, selects range and SNR
# plots against traveltime curves, either raw or reduced against traveltimes
# John Vidale 2/2019

import os
os.environ['PATH'] += os.pathsep + '/usr/local/bin'
os.chdir('/Users/vidale/Documents/GitHub/Array_codes')

#%% Import functions
from run_each_P import run_each

#run_each(event_no = 12) # no P wave for this one
#run_each(event_no = 15)
#run_each(event_no = 16)
#run_each(event_no = 17)
#run_each(event_no = 19)
#run_each(event_no = 18)
#run_each(event_no = 14)
#run_each(event_no = 20)
#run_each(event_no = 21)
#run_each(event_no = 23)
#run_each(event_no = 22)
#run_each(event_no = 24) # later window
#run_each(event_no = 25)
#run_each(event_no = 26)
#run_each(event_no = 27)
#run_each(event_no = 29)
#run_each(event_no =  7)
#run_each(event_no =  8)
#run_each(event_no = 30)  # very weak PKiKP
#run_each(event_no = 31)
#run_each(event_no = 33)  # very weak PKiKP
#run_each(event_no = 34)
#run_each(event_no = 32)
#run_each(event_no = 35)
#run_each(event_no = 36)  # very weak PKiKP
#run_each(event_no = 37)  # very weak PKiKP
#run_each(event_no = 38)
#run_each(event_no = 39)
#run_each(event_no = 40)
#run_each(event_no = 41)
#run_each(event_no = 42)
#run_each(event_no = 1)  # no apparent PKiKP
#run_each(event_no = 2)  # no apparent PKiKP
#run_each(event_no = 43)
#run_each(event_no = 44)  # no apparent PKiKP
#run_each(event_no = 45)
#run_each(event_no = 46)  # no apparent PKiKP
#run_each(event_no = 4)
#run_each(event_no = 52)
#run_each(event_no = 47)
#run_each(event_no = 48)  # no apparent PKiKP
#run_each(event_no = 49)
#run_each(event_no = 50)
#run_each(event_no = 51)
#run_each(event_no = 52)
#run_each(event_no = 53)
#run_each(event_no = 54)
#run_each(event_no = 55)
#run_each(event_no = 57)
#run_each(event_no = 58)
#run_each(event_no = 59)
#run_each(event_no = 60)
#run_each(event_no = 61)
#run_each(event_no = 63)
#run_each(event_no = 13)
#run_each(event_no = 64)
#run_each(event_no = 65)
#run_each(event_no = 66)
#run_each(event_no = 67)
#run_each(event_no = 68)
#run_each(event_no = 69)
#run_each(event_no = 70)  #  seems like other signals
#run_each(event_no = 71)
#run_each(event_no = 72)
#run_each(event_no = 73)
#run_each(event_no = 74)
#run_each(event_no = 75)
#run_each(event_no = 76)
#run_each(event_no = 77)
#run_each(event_no = 79)
#run_each(event_no = 80)
#run_each(event_no = 81)
run_each(event_no = 83)  # too far for p-tau

#  DUPES NOT YET PROCESSED AND NAMED
#run_each(start_buff = XXXX, end_buff = XXXX, event_no = 56)
#run_each(start_buff = XXXX, end_buff = XXXX, event_no = 62)
#run_each(start_buff = XXXX, end_buff = XXXX, event_no = 78)
#run_each(start_buff = XXXX, end_buff = XXXX, event_no = 82)