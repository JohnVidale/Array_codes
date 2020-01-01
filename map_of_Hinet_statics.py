##!/usr/bin/env python
##  code to align traces by shifting, and record time shifts for station statics
##  John Vidale, 2/2019
#
#import sys # don't show any warnings
#import warnings
#
#if not sys.warnoptions:
#    warnings.simplefilter("ignore")
#
import os
import numpy as np
import matplotlib.pyplot as plt
from pylab import *

#%% Get Hinet station location file
#sta_file = '/Users/vidale/Documents/PyCode/Codes/Hinet_station/hinet_master_list.txt'
sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/hinet_master_list.txt'
with open(sta_file, 'r') as file:
	lines = file.readlines()
print(str(len(lines)) + ' lines read from hinet_master_list')

# Load station coords into arrays
station_index = range(len(lines))
st_names = []
st_dist  = []
st_lats  = []
st_lons  = []
st_shift = []
st_corr  = []
for ii in station_index:
	line = lines[ii]
	split_line = line.split()
	st_names.append(split_line[0])
	st_dist.append(float(split_line[1]))
	st_lats.append(float(split_line[2]))
	st_lons.append(float(split_line[3]))
	st_shift.append(float(split_line[4]))
	st_corr.append(float(split_line[5]))

#%% Make plot

colors = st_dist[:]
#colors = st_corr[:]
sizes = 10 * (st_dist[:])

figure(4)
plt.scatter(st_lons, st_lats, c=colors, s=100, alpha = 0.3, cmap = 'gist_rainbow')
#plot(st_lons, st_lats, '.')
#plt.scatter(st_lons[:], st_lats[:], c=1, s=100, alpha=0.3, cmap='gist_rainbow')
#plt.scatter(shifts[:,0], shifts[:,1], c=colors, s=100, alpha=0.3, cmap='gist_rainbow')
#plt.scatter(shifts[:,0], shifts[:,1], c=colors, s=100, alpha=0.3, cmap='gist_rainbow')
#plot(shifts[:,0], shifts[:,1], '.')
xlabel('longitude (°N)')
ylabel('latitude (°E)')
title('Map of statics')
plt.colorbar();

os.system('say "Done"')
