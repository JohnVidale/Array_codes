#!/usr/bin/env python
# input is set of hinet traces
# this program tapers, filters, selects range and SNR
# plots against traveltime curves, either raw or reduced against traveltimes
# This programs deals with a single event.
# John Vidale 2/2019

import os
os.chdir('/Users/vidale/Documents/PyCode/Codes')
from pro3a_sort_plot_pair import pro3pair
from pro3b_sort_plot_singlet import pro3singlet

os.environ['PATH'] += os.pathsep + '/usr/local/bin'

#os.chdir('/Users/vidale/Documents/PyCode/Hinet/Reps/Aleutians/Pair1')
#pro3singlet(eq_file1 = 'event1.txt', eq_file2 = 'event2.txt', start_buff = 100,
#			end_buff = 300, plot_scale_fac = 0.1,qual_threshold = 3,
#			dphase = 'P', dphase2 = 'pP', dphase3 = 'PcP', dphase4 = 'sP',
#			corr_threshold = 0.7, freq_min = 1, freq_max = 3, min_dist = 25, max_dist = 44)

os.chdir('/Users/vidale/Documents/PyCode/Hinet/Reps/Aleutians/Pair1')
pro3pair(eq_file1 = 'event1.txt', eq_file2 = 'event2.txt', start_buff = 100,
			end_buff = 300, plot_scale_fac = 0.1,qual_threshold = 3,
			dphase = 'P', dphase2 = 'pP', dphase3 = 'PcP', dphase4 = 'sP',
			corr_threshold = 0.7, freq_min = 1, freq_max = 3, min_dist = 25, max_dist = 44)