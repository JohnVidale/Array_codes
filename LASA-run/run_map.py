#!/usr/bin/env python
# input is set of hinet traces
# this program tapers, filters, selects range and SNR
# plots against traveltime curves, either raw or reduced against traveltimes
# This programs deals with a single event.
# John Vidale 2/2019

import os

os.environ['PATH'] += os.pathsep + '/usr/local/bin'
os.chdir('/Users/vidale/Documents/GitHub/Array_codes')

#%% Import functions
from pro12_map_events    import map_plot
from pro13_PKiKP_comp    import map_slo_plot
os.chdir('/Users/vidale/Documents/PyCode/LASA')

#%% Common parameters
min_dist = 0
max_dist = 180

#%% Do it
map_plot(min_dist = min_dist, max_dist = max_dist)
#map_slo_plot(min_dist = min_dist, max_dist = max_dist)