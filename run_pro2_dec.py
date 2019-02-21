#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
os.chdir('/Users/vidale/Documents/GitHub/Hinet-codes')
from pro2_dec import pro2_decimate

os.environ['PATH'] += os.pathsep + '/usr/local/bin'

#os.chdir('/Users/vidale/Documents/PyCode/Hinet/Reps/Aleutians/Pair1')
#pro2_decimate('event1.txt')
#pro2_decimate('event2.txt')
#
#os.chdir('/Users/vidale/Documents/PyCode/Hinet/Reps/Aleutians/Pair2')
#pro2_decimate('event1.txt')
#pro2_decimate('event2.txt')
#
#os.chdir('/Users/vidale/Documents/PyCode/Hinet/Reps/Peru_Chile/Pair1')
#pro2_decimate('event1.txt')
#pro2_decimate('event2.txt')
#
#os.chdir('/Users/vidale/Documents/PyCode/Hinet/Reps/Peru_Chile/Pair2')
#pro2_decimate('event1.txt')
#pro2_decimate('event2.txt')
#
#os.chdir('/Users/vidale/Documents/PyCode/Hinet/Reps/SSI')
#pro2_decimate('event1.txt')
#pro2_decimate('event2.txt')
#
#os.chdir('/Users/vidale/Documents/PyCode/Hinet/Reps/Sumatra/Pair1')
#pro2_decimate('event1.txt')
#pro2_decimate('event2.txt')
#
#os.chdir('/Users/vidale/Documents/PyCode/Hinet/Reps/Sumatra/Pair2')
#pro2_decimate('event1.txt')
#pro2_decimate('event2.txt')
#
#os.chdir('/Users/vidale/Documents/PyCode/Hinet/Reps/Sumatra/Pair3')
#pro2_decimate('event1.txt')
#pro2_decimate('event2.txt')
#
#os.chdir('/Users/vidale/Documents/PyCode/Hinet/Reps/Tonga')
#pro2_decimate('event1.txt')
#pro2_decimate('event2.txt')

#os.chdir('/Users/vidale/Documents/PyCode/Hinet/Biggies/SSI_2016-05-28')
#pro2_decimate('event.txt')

#os.chdir('/Users/vidale/Documents/PyCode/Hinet/Biggies/SA_2018-04-02')
#pro2_decimate('event.txt')

os.chdir('/Users/vidale/Documents/PyCode/Hinet/Biggies/SSI_2011-12-11')
pro2_decimate('event.txt')
