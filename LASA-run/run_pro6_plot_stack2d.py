#!/usr/bin/env python
# input is set of hinet traces
# this program tapers, filters, selects range and SNR
# plots against traveltime curves, either raw or reduced against traveltimes
# This programs deals with a single event.
# John Vidale 2/2019

import os
os.chdir('/Users/vidale/Documents/GitHub/Hinet-codes')
from pro6_read_plot_stack import pro6plotstack

os.environ['PATH'] += os.pathsep + '/usr/local/bin'

os.chdir('/Users/vidale/Documents/PyCode/LASA')
pro6plotstack(eq_file = 'event1.txt', plot_scale_fac = 0.05,
			  slowR_lo = -0.05, slowR_hi = 0.05, slowT_lo = -0.05, slowT_hi = 0.05,
			  slow_delta = 0.005, start_buff = 10, end_buff = 70,
			  fig_index = 301, plot_dyn_range = 1000, snaptime = 50, snaps=15)