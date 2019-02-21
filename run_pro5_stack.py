#!/usr/bin/env python
# input is set of hinet traces
# this program tapers, filters, selects range and SNR
# plots against traveltime curves, either raw or reduced against traveltimes
# This programs deals with a single event.
# John Vidale 2/2019

import os
os.chdir('/Users/vidale/Documents/PyCode/Codes')
from pro5_stack import pro5stack

os.environ['PATH'] += os.pathsep + '/usr/local/bin'

os.chdir('/Users/vidale/Documents/PyCode/Hinet/Biggies/SA_2018-04-02')
pro5stack(eq_file = 'event.txt', plot_scale_fac = 0.05, slow_lo = -0.05, slow_hi = 0.05,
			  slow_delta = 0.001, start_time = -50, end_time = 500,
			  ref_lat = 36.3, ref_lon = 138.5, envelope = 1,
			  norm = 1, global_norm_plot = 1, color_plot = 1, fig_index = 301)