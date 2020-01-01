#!/usr/bin/env python
# input is set of hinet traces
# this program tapers, filters, selects range and SNR
# plots against traveltime curves, either raw or reduced against traveltimes
# This programs deals with a single event.
# John Vidale 2/2019

import os
os.chdir('/Users/vidale/Documents/GitHub/Hinet-codes')
from pro4_get_shifts import pro4statics

os.environ['PATH'] += os.pathsep + '/usr/local/bin'

#os.chdir('/Users/vidale/Documents/PyCode/Hinet/Biggies/SSI_2016-05-28')
#pro4statics(eq_file = 'event.txt', stat_corr = 1, ref_trace = 'N.SZW',
#			dphase = 'PKIKP', dphase2 = 'PKiKP', dphase3 = 'PKIKP', dphase4 = 'PKiKP',
#			start_corr_win = -1, end_corr_win = 3, plot_scale_fac = 0.05,
#			qual_threshold = 0, corr_threshold = 0, max_time_shift = 2)

os.chdir('/Users/vidale/Documents/PyCode/LASA/')
pro4statics(eq_file = 'event10.txt', ref_trace = 'A041z',
			dphase = 'PKIKP', dphase2 = 'PKiKP', dphase3 = 'PKIKP', dphase4 = 'PKiKP',
			start_corr_win = 0, end_corr_win = 23, plot_scale_fac = 0.05,
			qual_threshold = 0, corr_threshold = 0, max_time_shift = 2, ARRAY = 1,
			min_dist = 163.9, max_dist = 165.5)