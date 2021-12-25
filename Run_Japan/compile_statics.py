#!/usr/bin/env python
# Reads Peter Nelson (Grand student) Hinet event static files
# 9/30/2021

#%% Import
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

#    import sys # don't show any warnings
#    import warnings
#
#    if not sys.warnoptions:
#        warnings.simplefilter("ignore")

print(colored('Running compile statics', 'cyan'))
start_time_wc = time.time()

#%% -- Read station file
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

#%% -- Read Wei statics file
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

#%% -- Read fine statics files
peter1_sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/Tuned_statics/20010629.txt'
with open(peter1_sta_file, 'r') as file:
    lines1 = file.readlines()
# Load station coords into arrays
station_index = range(len(lines1))
peter1_st_names = []
peter1_st_static = []
for ii in station_index:
    line = lines1[ii]
    split_line = line.split()
    peter1_st_names.append(split_line[1])
    peter1_st_static.append(split_line[7])

# make lower case and remove trailing "h""
for ii in station_index:
    this_name = peter1_st_names[ii]
    peter1_st_names[ii]  = this_name[0:5]

peter2_sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/Tuned_statics/20020210.txt'
with open(peter2_sta_file, 'r') as file:
    lines2 = file.readlines()
# Load station coords into arrays
station_index = range(len(lines2))
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

peter3_sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/Tuned_statics/20030727.txt'
with open(peter3_sta_file, 'r') as file:
    lines3 = file.readlines()
# Load station coords into arrays
station_index = range(len(lines3))
peter3_st_names = []
peter3_st_static = []
for ii in station_index:
    line = lines3[ii]
    split_line = line.split()
    peter3_st_names.append(split_line[1])
    peter3_st_static.append(split_line[7])

# make lower case and remove trailing "h""
for ii in station_index:
    this_name = peter3_st_names[ii]
    peter3_st_names[ii]  = this_name[0:5]

peter4_sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/Tuned_statics/20040911.txt'
with open(peter4_sta_file, 'r') as file:
    lines4 = file.readlines()
# Load station coords into arrays
station_index = range(len(lines4))
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

peter5_sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/Tuned_statics/20050612.txt'
with open(peter5_sta_file, 'r') as file:
    lines5 = file.readlines()
# Load station coords into arrays
station_index = range(len(lines5))
peter5_st_names = []
peter5_st_static = []
for ii in station_index:
    line = lines5[ii]
    split_line = line.split()
    peter5_st_names.append(split_line[1])
    peter5_st_static.append(split_line[7])

# make lower case and remove trailing "h""
for ii in station_index:
    this_name = peter5_st_names[ii]
    peter5_st_names[ii]  = this_name[0:5]

peter6_sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/Tuned_statics/20111122.txt'
with open(peter6_sta_file, 'r') as file:
    lines6 = file.readlines()
# Load station coords into arrays
station_index = range(len(lines6))
peter6_st_names = []
peter6_st_static = []
for ii in station_index:
    line = lines6[ii]
    split_line = line.split()
    peter6_st_names.append(split_line[1])
    peter6_st_static.append(split_line[7])

# make lower case and remove trailing "h""
for ii in station_index:
    this_name = peter6_st_names[ii]
    peter6_st_names[ii]  = this_name[0:5]

peter7_sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/Tuned_statics/20111211.txt'
with open(peter7_sta_file, 'r') as file:
    lines7 = file.readlines()
# Load station coords into arrays
station_index = range(len(lines7))
peter7_st_names = []
peter7_st_static = []
for ii in station_index:
    line = lines7[ii]
    split_line = line.split()
    peter7_st_names.append(split_line[1])
    peter7_st_static.append(split_line[7])

# make lower case and remove trailing "h""
for ii in station_index:
    this_name = peter7_st_names[ii]
    peter7_st_names[ii]  = this_name[0:5]

peter8_sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/Tuned_statics/20140924.txt'
with open(peter8_sta_file, 'r') as file:
    lines8 = file.readlines()
# Load station coords into arrays
station_index = range(len(lines8))
peter8_st_names = []
peter8_st_static = []
for ii in station_index:
    line = lines8[ii]
    split_line = line.split()
    peter8_st_names.append(split_line[1])
    peter8_st_static.append(split_line[7])

# make lower case and remove trailing "h""
for ii in station_index:
    this_name = peter8_st_names[ii]
    peter8_st_names[ii]  = this_name[0:5]

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

#%% -- Read initial statics files
peter1_sta_newfile = '/Users/vidale/Documents/GitHub/Array_codes/Files/Tuned_statics/20010629_new.txt'
with open(peter1_sta_newfile, 'r') as file:
    lines1new = file.readlines()
station_index = range(len(lines1new)-1)  # skip first line with labels
peter1_st_newnames  = []
peter1_st_newstatic = []
peter1_st_dist      = []
for ii in station_index:
    line = lines1new[ii+1]
    split_line = line.split()
    peter1_st_newnames.append(split_line[2])
    peter1_st_newstatic.append(split_line[3])
    peter1_st_dist.append(split_line[8])

# make lower case and remove trailing "h""
for ii in station_index:
    this_name = peter1_st_newnames[ii]
    peter1_st_newnames[ii]  = this_name[0:5]

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

peter3_sta_newfile = '/Users/vidale/Documents/GitHub/Array_codes/Files/Tuned_statics/20030727_new.txt'
with open(peter3_sta_newfile, 'r') as file:
    lines3new = file.readlines()
station_index = range(len(lines3new)-1)  # skip first line with labels
peter3_st_newnames  = []
peter3_st_newstatic = []
peter3_st_dist      = []
for ii in station_index:
    line = lines3new[ii+1]
    split_line = line.split()
    peter3_st_newnames.append(split_line[2])
    peter3_st_newstatic.append(split_line[3])
    peter3_st_dist.append(split_line[8])

# make lower case and remove trailing "h""
for ii in station_index:
    this_name = peter3_st_newnames[ii]
    peter3_st_newnames[ii]  = this_name[0:5]

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

peter5_sta_newfile = '/Users/vidale/Documents/GitHub/Array_codes/Files/Tuned_statics/20050612_new.txt'
with open(peter5_sta_newfile, 'r') as file:
    lines5new = file.readlines()
station_index = range(len(lines5new)-1)  # skip first line with labels
peter5_st_newnames  = []
peter5_st_newstatic = []
peter5_st_dist      = []
for ii in station_index:
    line = lines5new[ii+1]
    split_line = line.split()
    peter5_st_newnames.append(split_line[2])
    peter5_st_newstatic.append(split_line[3])
    peter5_st_dist.append(split_line[8])

# make lower case and remove trailing "h""
for ii in station_index:
    this_name = peter5_st_newnames[ii]
    peter5_st_newnames[ii]  = this_name[0:5]

peter6_sta_newfile = '/Users/vidale/Documents/GitHub/Array_codes/Files/Tuned_statics/20111122_new.txt'
with open(peter6_sta_newfile, 'r') as file:
    lines6new = file.readlines()
station_index = range(len(lines6new)-1)  # skip first line with labels
peter6_st_newnames  = []
peter6_st_newstatic = []
peter6_st_dist      = []
for ii in station_index:
    line = lines6new[ii+1]
    split_line = line.split()
    peter6_st_newnames.append(split_line[2])
    peter6_st_newstatic.append(split_line[3])
    peter6_st_dist.append(split_line[8])

# make lower case and remove trailing "h""
for ii in station_index:
    this_name = peter6_st_newnames[ii]
    peter6_st_newnames[ii]  = this_name[0:5]

peter7_sta_newfile = '/Users/vidale/Documents/GitHub/Array_codes/Files/Tuned_statics/20111211_new.txt'
with open(peter7_sta_newfile, 'r') as file:
    lines7new = file.readlines()
station_index = range(len(lines7new)-1)  # skip first line with labels
peter7_st_newnames = []
peter7_st_newstatic = []
peter7_st_dist      = []
for ii in station_index:
    line = lines7new[ii+1]
    split_line = line.split()
    peter7_st_newnames.append(split_line[2])
    peter7_st_newstatic.append(split_line[3])
    peter7_st_dist.append(split_line[8])

# make lower case and remove trailing "h""
for ii in station_index:
    this_name = peter7_st_newnames[ii]
    peter7_st_newnames[ii]  = this_name[0:5]

peter8_sta_newfile = '/Users/vidale/Documents/GitHub/Array_codes/Files/Tuned_statics/20140924_new.txt'
with open(peter8_sta_newfile, 'r') as file:
    lines8new = file.readlines()
station_index = range(len(lines8new)-1)  # skip first line with labels
peter8_st_newnames  = []
peter8_st_newstatic = []
peter8_st_dist      = []
for ii in station_index:
    line = lines8new[ii+1]
    split_line = line.split()
    peter8_st_newnames.append(split_line[2])
    peter8_st_newstatic.append(split_line[3])
    peter8_st_dist.append(split_line[8])

# make lower case and remove trailing "h""
for ii in station_index:
    this_name = peter8_st_newnames[ii]
    peter8_st_newnames[ii]  = this_name[0:5]

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

#%% Match statics to station list
w_sta_found  = 0
p1_sta_found = 0
p2_sta_found = 0
p3_sta_found = 0
p4_sta_found = 0
p5_sta_found = 0
p6_sta_found = 0
p7_sta_found = 0
p8_sta_found = 0
p9_sta_found = 0
p1_sta_newfound = 0
p2_sta_newfound = 0
p3_sta_newfound = 0
p4_sta_newfound = 0
p5_sta_newfound = 0
p6_sta_newfound = 0
p7_sta_newfound = 0
p8_sta_newfound = 0
p9_sta_newfound = 0

w_corr = []
w_coef = []
p1_corr = []
p2_corr = []
p3_corr = []
p4_corr = []
p5_corr = []
p6_corr = []
p7_corr = []
p8_corr = []
p9_corr = []
p1_newcorr = []
p2_newcorr = []
p3_newcorr = []
p4_newcorr = []
p5_newcorr = []
p6_newcorr = []
p7_newcorr = []
p8_newcorr = []
p9_newcorr = []
p1_dist    = []
p2_dist    = []
p3_dist    = []
p4_dist    = []
p5_dist    = []
p6_dist    = []
p7_dist    = []
p8_dist    = []
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

    if this_name in peter1_st_names:  # find station in station list
        ii = peter1_st_names.index(this_name)
        p1_sta_found += 1
        corr = float(peter1_st_static[ii])
    else:
        corr = np.nan
    p1_corr.append(corr)

    if this_name in peter2_st_names:  # find station in station list
        ii = peter2_st_names.index(this_name)
        p2_sta_found += 1
        corr = float(peter2_st_static[ii])
    else:
        corr = np.nan
    p2_corr.append(corr)

    if this_name in peter3_st_names:  # find station in station list
        ii = peter3_st_names.index(this_name)
        p3_sta_found += 1
        corr = float(peter3_st_static[ii])
    else:
        corr = np.nan
    p3_corr.append(corr)

    if this_name in peter4_st_names:  # find station in station list
        ii = peter4_st_names.index(this_name)
        p4_sta_found += 1
        corr = float(peter4_st_static[ii])
    else:
        corr = np.nan
    p4_corr.append(corr)

    if this_name in peter5_st_names:  # find station in station list
        ii = peter5_st_names.index(this_name)
        p5_sta_found += 1
        corr = float(peter5_st_static[ii])
    else:
        corr = np.nan
    p5_corr.append(corr)

    if this_name in peter6_st_names:  # find station in station list
        ii = peter6_st_names.index(this_name)
        p6_sta_found += 1
        corr = float(peter6_st_static[ii])
    else:
        corr = np.nan
    p6_corr.append(corr)

    if this_name in peter7_st_names:  # find station in station list
        ii = peter7_st_names.index(this_name)
        p7_sta_found += 1
        corr = float(peter7_st_static[ii])
    else:
        corr = np.nan
    p7_corr.append(corr)

    if this_name in peter8_st_names:  # find station in station list
        ii = peter8_st_names.index(this_name)
        p8_sta_found += 1
        corr = float(peter8_st_static[ii])
    else:
        corr = np.nan
    p8_corr.append(corr)

    if this_name in peter9_st_names:  # find station in station list
        ii = peter9_st_names.index(this_name)
        p9_sta_found += 1
        corr = float(peter9_st_static[ii])
    else:
        corr = np.nan
    p9_corr.append(corr)

    if this_name in peter1_st_newnames:  # find station in station list
        ii = peter1_st_newnames.index(this_name)
        p1_sta_newfound += 1
        corr = float(peter1_st_newstatic[ii])
        dist = float(peter1_st_dist[ii])
    else:
        corr = np.nan
        dist = np.nan
    p1_newcorr.append(corr)
    p1_dist.append(dist)

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

    if this_name in peter3_st_newnames:  # find station in station list
        ii = peter3_st_newnames.index(this_name)
        p3_sta_newfound += 1
        corr = float(peter3_st_newstatic[ii])
        dist = float(peter3_st_dist[ii])
    else:
        corr = np.nan
        dist = np.nan
    p3_newcorr.append(corr)
    p3_dist.append(dist)

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

    if this_name in peter5_st_newnames:  # find station in station list
        ii = peter5_st_newnames.index(this_name)
        p5_sta_newfound += 1
        corr = float(peter5_st_newstatic[ii])
        dist = float(peter5_st_dist[ii])
    else:
        corr = np.nan
        dist = np.nan
    p5_newcorr.append(corr)
    p5_dist.append(dist)

    if this_name in peter6_st_newnames:  # find station in station list
        ii = peter6_st_newnames.index(this_name)
        p6_sta_newfound += 1
        corr = float(peter6_st_newstatic[ii])
        dist = float(peter6_st_dist[ii])
    else:
        corr = np.nan
        dist = np.nan
    p6_newcorr.append(corr)
    p6_dist.append(dist)

    if this_name in peter7_st_newnames:  # find station in station list
        ii = peter7_st_newnames.index(this_name)
        p7_sta_newfound += 1
        corr = float(peter7_st_newstatic[ii])
        dist = float(peter7_st_dist[ii])
    else:
        corr = np.nan
        dist = np.nan
    p7_newcorr.append(corr)
    p7_dist.append(dist)

    if this_name in peter8_st_newnames:  # find station in station list
        ii = peter8_st_newnames.index(this_name)
        p8_sta_newfound += 1
        corr = float(peter8_st_newstatic[ii])
        dist = float(peter8_st_dist[ii])
    else:
        corr = np.nan
        dist = np.nan
    p8_newcorr.append(corr)
    p8_dist.append(dist)

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

#%% Print out stats from associating statics
print(str(len(linesf)) + ' stations read from ' + full_sta_file)
print('W stations in common ' + str(w_sta_found))
print('P1: old stations '+ str(len(lines1)) + ' found ' + str(p1_sta_found) + ' new stations '+ str(len(lines1new)-1) + ' found ' + str(p1_sta_newfound))
print('P2: old stations '+ str(len(lines2)) + ' found ' + str(p2_sta_found) + ' new stations '+ str(len(lines2new)-1) + ' found ' + str(p2_sta_newfound))
print('P3: old stations '+ str(len(lines3)) + ' found ' + str(p3_sta_found) + ' new stations '+ str(len(lines3new)-1) + ' found ' + str(p3_sta_newfound))
print('P4: old stations '+ str(len(lines4)) + ' found ' + str(p4_sta_found) + ' new stations '+ str(len(lines4new)-1) + ' found ' + str(p4_sta_newfound))
print('P5: old stations '+ str(len(lines5)) + ' found ' + str(p5_sta_found) + ' new stations '+ str(len(lines5new)-1) + ' found ' + str(p5_sta_newfound))
print('P6: old stations '+ str(len(lines6)) + ' found ' + str(p6_sta_found) + ' new stations '+ str(len(lines6new)-1) + ' found ' + str(p6_sta_newfound))
print('P7: old stations '+ str(len(lines7)) + ' found ' + str(p7_sta_found) + ' new stations '+ str(len(lines7new)-1) + ' found ' + str(p7_sta_newfound))
print('P8: old stations '+ str(len(lines8)) + ' found ' + str(p8_sta_found) + ' new stations '+ str(len(lines8new)-1) + ' found ' + str(p8_sta_newfound))
print('P9: old stations '+ str(len(lines9)) + ' found ' + str(p9_sta_found) + ' new stations '+ str(len(lines9new)-1) + ' found ' + str(p9_sta_newfound))

# print(str(len(w_corr)) + ' number of w_cor values')
# print(str(len(w_coef)) + ' number of w_coef values')
# print(str(len(p1_corr)) + ' number of p1_cor values')
# print(str(len(p2_corr)) + ' number of p2_cor values')
# print(str(len(p3_corr)) + ' number of p3_cor values')
# print(str(len(p4_corr)) + ' number of p4_cor values')
# print(str(len(p5_corr)) + ' number of p5_cor values')
# print(str(len(p6_corr)) + ' number of p6_cor values')
# print(str(len(p7_corr)) + ' number of p7_cor values')
# print(str(len(p8_corr)) + ' number of p8_cor values')
# print(str(len(p9_corr)) + ' number of p9_cor values')
# print(str(len(p1_newcorr)) + ' number of p1_cor values')

# for i in range(station_num):
#     print(full_st_names[i] + ' ' + str(i) + ' ' + str(w_corr[i]) + ' ' + str(w_coef[i]) + ' ' +
#           str(p1_corr[i]) + ' ' + str(p2_corr[i]) + ' ' + str(p3_corr[i]) + ' ' +
#           str(p4_corr[i]) + ' ' + str(p5_corr[i]) + ' ' + str(p6_corr[i]) + ' ' +
#           str(p7_corr[i]) + ' ' + str(p8_corr[i]) + ' ' + str(p9_corr[i]))

plt.close('all')

"""
#%% Scatter - Compare w/ & w/o refinement against Wei
year = '2015'
plot_lim = 2
corr_1st = copy.copy(p9_newcorr) # make event easy to change with four numbers
corr_2nd = copy.copy(p9_corr)
corr_Wei = copy.copy(w_corr)
dist     = copy.copy(p9_dist)

mean_1st = np.nanmean(corr_1st) # find mean
mean_2nd = np.nanmean(corr_2nd)
mean_Wei = np.nanmean(corr_Wei)
corr_1st = np.hstack(corr_1st - mean_1st) # subtract mean
corr_2nd = np.hstack(corr_2nd - mean_2nd)
corr_Wei = np.hstack(corr_Wei - mean_Wei)

p_sum = np.subtract(corr_1st,corr_2nd)
print('old ' + str(len(corr_1st)) + ' new ' + str(len(corr_2nd)) + ' sum ' + str(len(p_sum)))
print('old ' + str(mean_1st) + ' new ' + str(mean_2nd) + ' sum ' + str(mean_Wei))

dist_min = 150  # limits on distance
dist_max = 160
in_range = 0
for i in range(station_num):
    if dist[i] < dist_min or dist[i] > dist_max:
        corr_Wei[i] = np.nan
        corr_1st[i] = np.nan
        corr_2nd[i] = np.nan
    else:
        in_range += 1
print(str(in_range) + ' out of ' + str(station_num) + ' are in range')

plt.figure(100,figsize=(4,4))
plt.xlim(-plot_lim,plot_lim)
plt.ylim(-plot_lim,plot_lim)
plt.xlabel('My statics (s)')
plt.ylabel('Peter statics (s)')
plt.title(year + ' old file')
# for i in range(station_num):
#     p1_corr[i] *= -1
plt.scatter(corr_Wei,corr_2nd)
plt.savefig('/Users/vidale/Desktop/' + year + '_old')

plt.figure(101,figsize=(4,4))
plt.xlim(-plot_lim,plot_lim)
plt.ylim(-plot_lim,plot_lim)
plt.xlabel('My statics (s)')
plt.ylabel('Peter statics (s)')
plt.title(year + ' new file')
plt.scatter(corr_Wei,corr_1st)
plt.savefig('/Users/vidale/Desktop/' + year + '_new')

plt.figure(102,figsize=(4,4))
plt.xlim(-plot_lim,plot_lim)
plt.ylim(-plot_lim,plot_lim)
plt.xlabel('My statics (s)')
plt.ylabel('Peter statics (s)')
plt.title(year + ' diff')
plt.scatter(corr_Wei,p_sum)
plt.savefig('/Users/vidale/Desktop/' + year + '_diff')
"""

#%% Scatter - Compare two of Peter's statics
year = '2001 & 2015'
plot_lim = 2
corr_1_1st = copy.copy(p1_newcorr) # make event easy to change with four numbers
corr_1_2nd = copy.copy(p1_corr)
corr_2_1st = copy.copy(p9_newcorr)
corr_2_2nd = copy.copy(p9_corr)
dist     = copy.copy(p1_dist)

dist_min = 150  # limits on distance
dist_max = 160
in_range = 0
for i in range(station_num):
    if dist[i] < dist_min or dist[i] > dist_max:
        corr_1_1st[i] = np.nan
        corr_1_2nd[i] = np.nan
        corr_2_1st[i] = np.nan
        corr_2_2nd[i] = np.nan
    else:
        in_range += 1
print(str(in_range) + ' out of ' + str(station_num) + ' are in range')

thresh_out = 0
stat_thresh = 3  # limit for rejecting outliers
for i in range(station_num):
    if abs(corr_1_1st[i]) > stat_thresh:
        corr_1_1st[i] = np.nan
        thresh_out += 1
    if abs(corr_1_2nd[i]) > stat_thresh:
        corr_1_2nd[i] = np.nan
        thresh_out += 1
    if abs(corr_2_1st[i]) > stat_thresh:
        corr_2_1st[i] = np.nan
        thresh_out += 1
    if abs(corr_2_2nd[i]) > stat_thresh:
        corr_2_2nd[i] = np.nan
        thresh_out += 1
print(str(thresh_out) + ' out of ' + str(station_num) + ' are rejected from threshold')

mean_1_1st = np.nanmean(corr_1_1st) # find mean
mean_1_2nd = np.nanmean(corr_1_2nd)
mean_2_1st = np.nanmean(corr_2_1st) # find mean
mean_2_2nd = np.nanmean(corr_2_2nd)
corr_1_1st = np.hstack(corr_1_1st - mean_1_1st) # subtract mean
corr_1_2nd = np.hstack(corr_1_2nd - mean_1_2nd)
corr_2_1st = np.hstack(corr_2_1st - mean_2_1st) # subtract mean
corr_2_2nd = np.hstack(corr_2_2nd - mean_2_2nd)

p_sub_1 = np.subtract(corr_1_1st,corr_1_2nd)
p_sub_2 = np.subtract(corr_2_1st,corr_2_2nd)
p_add_1 = np.add(corr_1_1st,corr_1_2nd)
p_add_2 = np.add(corr_2_1st,corr_2_2nd)
# print('old ' + str(len(corr_1st)) + ' new ' + str(len(corr_2nd)) + ' sum ' + str(len(p_sum)))
# print('old ' + str(mean_1st) + ' new ' + str(mean_2nd) + ' sum ' + str(mean_Wei))

plt.figure(100,figsize=(4,4))
plt.xlim(-plot_lim,plot_lim)
plt.ylim(-plot_lim,plot_lim)
plt.xlabel('My statics (s)')
plt.ylabel('Peter statics (s)')
plt.title(year + ' simple')
# for i in range(station_num):
#     p1_corr[i] *= -1
plt.scatter(corr_1_1st,corr_2_1st)
plt.savefig('/Users/vidale/Desktop/' + year + ' simple')

plt.figure(101,figsize=(4,4))
plt.xlim(-plot_lim,plot_lim)
plt.ylim(-plot_lim,plot_lim)
plt.xlabel('My statics (s)')
plt.ylabel('Peter statics (s)')
plt.title(year + ' sub')
plt.scatter(p_sub_1,p_sub_2)
plt.savefig('/Users/vidale/Desktop/' + year + ' sub')

plt.figure(102,figsize=(4,4))
plt.xlim(-plot_lim,plot_lim)
plt.ylim(-plot_lim,plot_lim)
plt.xlabel('My statics (s)')
plt.ylabel('Peter statics (s)')
plt.title(year + ' add')
plt.scatter(p_add_1,p_add_2)
plt.savefig('/Users/vidale/Desktop/' + year + ' add')
"""
#%% Histogram - process numbers
year1 = '2002'
year2 = '2004'
tit_label = ' minus'
event1_1st = copy.copy(p2_newcorr) # make a copy to mess with
event1_2nd = copy.copy(p2_corr)
event2_1st = copy.copy(p4_newcorr)
event2_2nd = copy.copy(p4_corr)
Wei        = copy.copy(w_corr)
dist       = copy.copy(p2_dist)

mean_1_1st = np.nanmean(event1_1st) # find mean
mean_1_2nd = np.nanmean(event1_1st)
mean_2_1st = np.nanmean(event2_1st)
mean_2_2nd = np.nanmean(event2_1st)
mean_Wei   = np.nanmean(Wei)
event1_1st = np.hstack(event1_1st - mean_1_1st) # subtract mean
event1_2nd = np.hstack(event1_2nd - mean_1_2nd)
event2_1st = np.hstack(event2_1st - mean_2_1st)
event2_2nd = np.hstack(event2_2nd - mean_2_2nd)
Wei        = np.nanmean(Wei - mean_Wei)

diff = []
for i in range(station_num):
    diff.append((event1_1st[i] - event1_2nd[i]) - (event2_1st[i] - event2_2nd[i]))
    # diff.append(event1_1st[i] - (w_corr[i]))
    # diff.append((event1_1st[i] - event1_2nd[i]) - (w_corr[i]))
    # diff.append(event1_1st[i] - event2_1st[i])

print('event1_1st max = ' + str(max(event1_1st)) + ' min = ' + str(min(event1_1st)))
print('event1_2nd max = ' + str(max(event1_2nd)) + ' min = ' + str(min(event1_2nd)))
print('event2_1st max = ' + str(max(event2_1st)) + ' min = ' + str(min(event2_1st)))
print('event2_2nd max = ' + str(max(event2_2nd)) + ' min = ' + str(min(event2_2nd)))

dist_min = 150  # limits on distance
dist_max = 160
for i in range(station_num):
    if dist[i] < dist_min or dist[i] > dist_max:
        diff[i] = np.nan

stat_thresh = 10  # limit for rejecting outliers
for i in range(station_num):
    if abs(diff[i]) > stat_thresh:
        diff[i] = np.nan

#%% Histogram plots
plt.figure(50,figsize=(4,4))
diff_mean = np.nanmean(diff)
print(f'Mean, stddev {diff_mean:.2f},{np.nanstd(diff - diff_mean):.2f}')
a = np.hstack(diff - diff_mean)
_ = plt.hist(a, bins='auto')  # arguments are passed to np.histogram
# plt.title("2002 - 2004 new + old files, 150-160°")
plt.title(year1 + ' vs ' + year2 + ' 150-160°' + tit_label)
plt.xlim(-0.6,0.6)
plt.xlabel('time difference (s)', fontsize = 12)
plt.ylabel('count', fontsize = 12)
plt.show()
plt.savefig('/Users/vidale/Desktop/' + year1 + '_' + year2 + tit_label)
"""

elapsed_time_wc = time.time() - start_time_wc
print(f'    This job took   {elapsed_time_wc:.1f}   seconds')
os.system('say "statics done"')