#!/usr/bin/env python
# reads in, removes baseline, and differences statics, makes binned plot
# John Vidale 9/2021

#%% Import
#%% -- Functions
from obspy import UTCDateTime
from obspy import Stream
from obspy import read
from obspy.geodetics import gps2dist_azimuth
import numpy as np
import os
import sys
from obspy.taup import TauPyModel
import matplotlib.pyplot as plt
import time
import statistics
import copy
from termcolor import colored
model = TauPyModel(model='iasp91')

print(colored('Running compile statics', 'cyan'))
start_time_wc = time.time()


#%% -- Station locations and statics

full_sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/sta_hinet.txt'
with open(full_sta_file, 'r') as file:
    linesf = file.readlines()
# Load station coords into arrays
station_index = range(len(linesf))
station_num = len(linesf)
full_st_names = []
for ii in station_index:
    line = linesf[ii]
    split_line = line.split()
    full_st_names.append(split_line[0])
# make lower case and remove trailing "h""
for ii in station_index:
    this_name = full_st_names[ii]
    this_name_truc = this_name[0:5]
    full_st_names[ii]  = this_name_truc.upper()

wei_sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/sta_statics_hinet.txt'
with open(wei_sta_file, 'r') as file:
    linesw = file.readlines()
# Load station coords into arrays
station_index = range(len(linesw))
w_st_names = []
w_st_shift = []
w_st_corr  = []
for ii in station_index:
    line = linesw[ii]
    split_line = line.split()
    w_st_names.append(split_line[0])
    w_st_shift.append(split_line[4])
    w_st_corr.append(split_line[5])

peter2_sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/Tuned_statics/20020210.txt'
with open(peter2_sta_file, 'r') as file:
    lines2 = file.readlines()
# Load station coords into arrays
station_index = range(len(lines2))
print('event 2 has ' + str(len(lines2)-1))
peter2_st_names = []
peter2_st_static = []
for ii in station_index:
    line = lines2[ii]
    split_line = line.split()
    peter2_st_names.append(split_line[1])
    peter2_st_static.append(split_line[7])
# make lower case and remove trailing "h""
for ii in station_index:
    this_name = peter2_st_names[ii]
    peter2_st_names[ii]  = this_name[0:5]

peter4_sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/Tuned_statics/20040911.txt'
with open(peter4_sta_file, 'r') as file:
    lines4 = file.readlines()
# Load station coords into arrays
station_index = range(len(lines4))
print('event 4 has ' + str(len(lines4)-1))
peter4_st_names = []
peter4_st_static = []
for ii in station_index:
    line = lines4[ii]
    split_line = line.split()
    peter4_st_names.append(split_line[1])
    peter4_st_static.append(split_line[7])
# make lower case and remove trailing "h""
for ii in station_index:
    this_name = peter4_st_names[ii]
    peter4_st_names[ii]  = this_name[0:5]

peter9_sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/Tuned_statics/20150610.txt'
with open(peter9_sta_file, 'r') as file:
    lines9 = file.readlines()
# Load station coords into arrays
station_index = range(len(lines9))
peter9_st_names = []
peter9_st_static = []
for ii in station_index:
    line = lines9[ii]
    split_line = line.split()
    peter9_st_names.append(split_line[1])
    peter9_st_static.append(split_line[7])
# make lower case and remove trailing "h""
for ii in station_index:
    this_name = peter9_st_names[ii]
    peter9_st_names[ii]  = this_name[0:5]

peter2_sta_newfile = '/Users/vidale/Documents/GitHub/Array_codes/Files/Tuned_statics/20020210_new.txt'
with open(peter2_sta_newfile, 'r') as file:
    lines2new = file.readlines()
station_index = range(len(lines2new)-1)  # skip first line with labels
peter2_st_newnames  = []
peter2_st_newstatic = []
peter2_st_dist      = []
for ii in station_index:
    line = lines2new[ii+1]
    split_line = line.split()
    peter2_st_newnames.append(split_line[2])
    peter2_st_newstatic.append(split_line[3])
    peter2_st_dist.append(split_line[8])

# make lower case and remove trailing "h""
for ii in station_index:
    this_name = peter2_st_newnames[ii]
    peter2_st_newnames[ii]  = this_name[0:5]

peter4_sta_newfile = '/Users/vidale/Documents/GitHub/Array_codes/Files/Tuned_statics/20040911_new.txt'
with open(peter4_sta_newfile, 'r') as file:
    lines4new = file.readlines()
station_index = range(len(lines4new)-1)  # skip first line with labels
peter4_st_newnames  = []
peter4_st_newstatic = []
peter4_st_dist      = []
for ii in station_index:
    line = lines4new[ii+1]
    split_line = line.split()
    peter4_st_newnames.append(split_line[2])
    peter4_st_newstatic.append(split_line[3])
    peter4_st_dist.append(split_line[8])

# make lower case and remove trailing "h""
for ii in station_index:
    this_name = peter4_st_newnames[ii]
    peter4_st_newnames[ii]  = this_name[0:5]

peter9_sta_newfile = '/Users/vidale/Documents/GitHub/Array_codes/Files/Tuned_statics/20150610_new.txt'
with open(peter9_sta_newfile, 'r') as file:
    lines9new = file.readlines()
station_index = range(len(lines9new)-1)  # skip first line with labels
peter9_st_newnames  = []
peter9_st_newstatic = []
peter9_st_dist      = []
for ii in station_index:
    line = lines9new[ii+1]
    split_line = line.split()
    peter9_st_newnames.append(split_line[2])
    peter9_st_newstatic.append(split_line[3])
    peter9_st_dist.append(split_line[8])

# make lower case and remove trailing "h""
for ii in station_index:
    this_name = peter9_st_newnames[ii]
    peter9_st_newnames[ii]  = this_name[0:5]

w_sta_found  = 0
p2_sta_found = 0
p4_sta_found = 0
p9_sta_found = 0
p2_sta_newfound = 0
p4_sta_newfound = 0
p9_sta_newfound = 0

w_corr     = []
w_coef     = []
p2_corr    = []
p4_corr    = []
p9_corr    = []
p2_newcorr = []
p4_newcorr = []
p9_newcorr = []
p2_dist    = []
p4_dist    = []
p9_dist    = []

for this_name in full_st_names: # traces one by one, find lat-lon
    if this_name in w_st_names:  # find station in station list
        ii = w_st_names.index(this_name)
        w_sta_found += 1

        corr = float(w_st_shift[ii])
        coef = float(w_st_corr[ii])
    else:
        corr = np.nan
        coef = np.nan
    w_corr.append(corr)
    w_coef.append(coef)

    if this_name in peter2_st_names:  # find station in station list
        ii = peter2_st_names.index(this_name)
        p2_sta_found += 1
        corr = float(peter2_st_static[ii])
    else:
        corr = np.nan
    p2_corr.append(corr)

    if this_name in peter4_st_names:  # find station in station list
        ii = peter4_st_names.index(this_name)
        p4_sta_found += 1
        corr = float(peter4_st_static[ii])
    else:
        corr = np.nan
    p4_corr.append(corr)

    if this_name in peter9_st_names:  # find station in station list
        ii = peter9_st_names.index(this_name)
        p9_sta_found += 1
        corr = float(peter9_st_static[ii])
    else:
        corr = np.nan
    p9_corr.append(corr)

    if this_name in peter2_st_newnames:  # find station in station list
        ii = peter2_st_newnames.index(this_name)
        p2_sta_newfound += 1
        corr = float(peter2_st_newstatic[ii])
        dist = float(peter2_st_dist[ii])
    else:
        corr = np.nan
        dist = np.nan
    p2_newcorr.append(corr)
    p2_dist.append(dist)

    if this_name in peter4_st_newnames:  # find station in station list
        ii = peter4_st_newnames.index(this_name)
        p4_sta_newfound += 1
        corr = float(peter4_st_newstatic[ii])
        dist = float(peter4_st_dist[ii])
    else:
        corr = np.nan
        dist = np.nan
    p4_newcorr.append(corr)
    p4_dist.append(dist)

    if this_name in peter9_st_newnames:  # find station in station list
        ii = peter9_st_newnames.index(this_name)
        p9_sta_newfound += 1
        corr = float(peter9_st_newstatic[ii])
        dist = float(peter9_st_dist[ii])
    else:
        corr = np.nan
        dist = np.nan
    p9_newcorr.append(corr)
    p9_dist.append(dist)

# for i in range(station_num):
    # p9_newcorr[i] += -0.1  # ignoring initial numbers
    # p9_newcorr[i] += +1.72 # minus initial numbers
    # p9_newcorr[i] += -2 # plus initial numbers
    # p2_newcorr[i] += +1.4 # minus initial numbers

print('    ' + str(len(linesf)) + ' stations read from ' + full_sta_file)
print('    ' + str(len(linesw)) + ' coarse station statics read from ' + wei_sta_file)
print('    ' + str(len(lines9)) + ' station statics read from ' + peter9_sta_file)
print('    ' + str(len(lines2new)) + ' station statics read from ' + peter2_sta_newfile)
print('    ' + str(len(lines4new)) + ' station statics read from ' + peter4_sta_newfile)
print('    ' + str(len(lines9new)) + ' station statics read from ' + peter9_sta_newfile)

# print(str(w_sta_found) + ' W stations in common')
print(str(p2_sta_found) + ' P2 stations in common')
print(str(p4_sta_found) + ' P4 stations in common')
print(str(p9_sta_found) + ' P9 stations in common')
print(str(p2_sta_newfound) + ' P2new stations in common')
print(str(p4_sta_newfound) + ' P4new stations in common')
print(str(p9_sta_newfound) + ' P9new stations in common')

p9_sta_found = 0
p9_sta_newfound = 0

print(str(p9_sta_found) + ' P9new stations after threshold')
print(str(p2_sta_newfound) + ' P2 stations after threshold')
print(str(p4_sta_newfound) + ' P4 stations after threshold')
print(str(p9_sta_newfound) + ' P9 stations after threshold')

print(str(station_num))

"""
for i in range(station_num):
    print(full_st_names[i] + ' ' + str(i) + ' ' +
          str(p1_corr[i]) + ' ' + str(p2_corr[i]) + ' ' + str(p3_corr[i]) + ' ' +
          str(p4_corr[i]) + ' ' + str(p5_corr[i]))
"""

# corr_median = []
# for i in range(station_num):
#     x = np.nanmedian([p1_corr[i], p2_corr[i], p3_corr[i], p4_corr[i], p5_corr[i]])
#     corr_median.append(x)

# for i in range(station_num):
#     print(full_st_names[i] + ' ' + str(corr_median[i]) + ' ' +
#           str(p1_corr[i]) + ' ' + str(p2_corr[i]) + ' ' + str(p3_corr[i]) + ' ' +
#           str(p4_corr[i]) + ' ' + str(p5_corr[i]))

corr2_2nd = copy.copy(p2_newcorr)
corr4_2nd = copy.copy(p4_newcorr)
corr2_1st = copy.copy(p2_corr)
corr4_1st = copy.copy(p4_corr)

diff = []
for i in range(station_num):
    # diff.append((corr2_2nd[i] + corr2_1st[i]) - (corr4_2nd[i] + corr4_1st[i]))
    diff.append((corr4_2nd[i] - corr4_1st[i]) - (w_corr[i]))

print('P2 max = ' + str(max(p2_corr)) + ' min = ' + str(min(p2_corr)))
print('P4 max = ' + str(max(p4_corr)) + ' min = ' + str(min(p4_corr)))
print('P2new max = ' + str(max(p2_newcorr)) + ' min = ' + str(min(p2_newcorr)))
print('P4new max = ' + str(max(p4_newcorr)) + ' min = ' + str(min(p4_newcorr)))

dist_min = 150  # limits on distance
dist_max = 160
for i in range(station_num):
    if p2_dist[i] < dist_min or p2_dist[i] > dist_max:
        diff[i] = np.nan

stat_thresh = 10  # limit for rejecting outliers
for i in range(station_num):
    if abs(diff[i]) > stat_thresh:
        diff[i] = np.nan

plt.close('all')
plt.figure(50,figsize=(5,1.5))
diff_mean = np.nanmean(diff)
print(f'Mean, stddev {diff_mean:.2f},{np.nanstd(diff - diff_mean):.2f}')
a = np.hstack(diff - diff_mean)
_ = plt.hist(a, bins='auto')  # arguments are passed to np.histogram
# plt.title("2002 - 2004 new + old files, 150-160°")
plt.title("2004 vs Wei, 150-160°")
plt.xlim(-0.6,0.6)
plt.show()

elapsed_time_wc = time.time() - start_time_wc
print(f'    This job took   {elapsed_time_wc:.1f}   seconds')
os.system('say "comparing statics"')