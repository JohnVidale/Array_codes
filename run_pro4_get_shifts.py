#!/usr/bin/env python
# input is set of hinet traces
# this program tapers, filters, selects range and SNR
# plots against traveltime curves, either raw or reduced against traveltimes
# This programs deals with a single event.
# John Vidale 2/2019

import os
os.chdir('/Users/vidale/Documents/PyCode/Codes')
from pro4_get_shifts import pro4statics

os.environ['PATH'] += os.pathsep + '/usr/local/bin'

os.chdir('/Users/vidale/Documents/PyCode/Hinet/Biggies/SSI_2016-05-28_vgood')
pro4statics(eq_file = 'event.txt', stat_corr = 1, ref_trace = 'N.SZW',
			dphase = 'PKIKP', dphase2 = 'PKiKP', dphase3 = 'PKIKP', dphase4 = 'PKiKP',
			start_corr_win = -1, end_corr_win = 3, plot_scale_fac = 0.05,
			qual_threshold = 0, corr_threshold = 0, max_time_shift = 2)