#!/usr/bin/env python
# input is set of hinet traces
# this program tapers, filters, selects range and SNR
# plots against traveltime curves, either raw or reduced against traveltimes
# This programs deals with a single event.
# John Vidale 2/2019

import os
os.chdir('/Users/vidale/Documents/GitHub/Hinet-codes')
from pro5b_stack2d import pro5stack2d

os.environ['PATH'] += os.pathsep + '/usr/local/bin'

os.chdir('/Users/vidale/Documents/PyCode/LASA')
pro5stack2d(eq_file = 'event1.txt', plot_scale_fac = 0.05,
			  slowR_lo = -0.05, slowR_hi = 0.05, slowT_lo = -0.05, slowT_hi = 0.05,
			  slow_delta = 0.005, start_buff = 10, end_buff = 250,
			  norm = 1, global_norm_plot = 1,
			  ARRAY=1, decimate_fac = 5)

#pro5stack2d(eq_file = 'event1.txt', plot_scale_fac = 0.05,
#			  slowR_lo = -0.05, slowR_hi = 0.05, slowT_lo = -0.05, slowT_hi = 0.05,
#			  slow_delta = 0.01, start_buff = 10, end_buff = 250,
#			  envelope = 1, norm = 1, global_norm_plot = 1,
#			  ARRAY=1, fig_index = 301, decimate_fac = 5)