#!/usr/bin/env python
# input is set of hinet traces
# this program tapers, filters, selects range and SNR
# plots against traveltime curves, either raw or reduced against traveltimes
# This programs deals with a single event.
# John Vidale 2/2019

import os
os.chdir('/Users/vidale/Documents/GitHub/Hinet-codes')
from pro5a_stack import pro5stack

os.environ['PATH'] += os.pathsep + '/usr/local/bin'

os.chdir('/Users/vidale/Documents/PyCode/LASA')
pro5stack(ARRAY = 1, eq_file = 'event1.txt', plot_scale_fac = 0.05, slow_lo = -0.05, slow_hi = 0.05,
			  slow_delta = 0.002, start_buff = 10, end_buff = 250,
			  envelope = 1, plot_dyn_range = 1000,
			  norm = 1, global_norm_plot = 1, color_plot = 1, fig_index = 301)