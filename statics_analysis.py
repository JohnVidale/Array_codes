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
import numpy as np
from pylab import *

def clean_station_name(name):
	return name.lower() + 'h'
#%% Get Hinet station location file

file_label = '2s.txt'
sta_file = '/Users/vidale/Documents/PyCode/Codes/Hinet_station/hinet_list'
with open(sta_file, 'r') as file:
	lines = file.readlines()
# Load station coords into arrays
station_index = range(781)
st_lats  = []
st_lons  = []
st_deps  = []
st_names = []
sta_master = {}
statics = {}

for ii in station_index:
	line = lines[ii]
	split_line = line.split()
	sta_master[split_line[0]] = {'lat': float(split_line[1]), 'lon': float(split_line[2]), 'elev': float(split_line[3]),
							     'dist1': np.nan, 'shift1': np.nan, 'corr1': np.nan,
								 'dist2': np.nan, 'shift2': np.nan, 'corr2': np.nan,
								 'dist3': np.nan, 'shift3': np.nan, 'corr3': np.nan}

#stations['n.agwh']
#Out[2]: {'lat': '43.0842', 'lon': '140.8199', 'elev': '-77'}
#
#stations['n.agwh']['lat']
#Out[3]: '43.0842'

shift_file1 = '/Users/vidale/Documents/PyCode/Hinet/Statics/HA2016-05-28pro4_PKIKP.txt'
with open(shift_file1, 'r') as file:
	lines = file.readlines()

for line in lines:
	split_line = line.split()
	name = clean_station_name(split_line[0])
	statics[name] = {'dist1': float(split_line[1]), 'shift1': float(split_line[4]), 'corr1': float(split_line[5])}
	sta_master[name].update(statics[name])

shift_file2 = '/Users/vidale/Documents/PyCode/Hinet/Statics/HA2011-12-11pro4_PKIKP.txt'
with open(shift_file2, 'r') as file:
	lines = file.readlines()

for line in lines:
	split_line = line.split()
	name = clean_station_name(split_line[0])
	statics[name] = {'dist2': float(split_line[1]), 'shift2': float(split_line[4]), 'corr2': float(split_line[5])}
	sta_master[name].update(statics[name])

shift_file3 = '/Users/vidale/Documents/PyCode/Hinet/Statics/HA2018-04-02pro4_PKIKP.txt'
with open(shift_file3, 'r') as file:
	lines = file.readlines()

for line in lines:
	split_line = line.split()
	name = clean_station_name(split_line[0])
	statics[name] = {'dist3': float(split_line[1]), 'shift3': float(split_line[4]), 'corr3': float(split_line[5])}
	sta_master[name].update(statics[name])


#sta_master['n.agwh']
#Out[13]:
#{'lat': 43.0842,
# 'lon': 140.8199,
# 'elev': -77.0,
# 'dist1': 150.29,
# 'shift1': 1.233,
# 'corr1': 0.472,
# 'dist2': nan,
# 'shift2': nan,
# 'corr2': nan}

shifts = np.array([(v['shift3'], v['shift2']) for k,v in sta_master.items()])
corrs  = np.array([(v['corr1' ], v['corr2' ]) for k,v in sta_master.items()])
dists  = np.array([(v['dist1' ], v['dist2' ]) for k,v in sta_master.items()])
#lonlat  = np.array([(v['lon' ], v['lat' ]) for k,v in sta_master.items()])

colors = corrs[:,0]
#colors = shifts[:,1]
sizes = 10 * (dists[:,0]-145)
#sizes = 10 * dists[:,0]

figure(3)
plt.scatter(shifts[:,0], shifts[:,1], c=colors, s=100, alpha=0.3, cmap='gist_rainbow')
#plt.scatter(shifts[:,0], shifts[:,1], c=colors, s=100, alpha=0.3, cmap='gist_rainbow')
#plt.scatter(shifts[:,0], shifts[:,1], c=colors, s=100, alpha=0.3, cmap='gist_rainbow')
#plot(shifts[:,0], shifts[:,1], '.')
xlabel('SA_2s.txt')
ylabel('SSI')
title('Shift statics')
plt.colorbar();


