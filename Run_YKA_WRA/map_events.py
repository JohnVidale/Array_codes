#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 29 10:44:55 2020
Makes a map of the LASA events
@author: vidale
"""
import numpy as np
import matplotlib.pyplot as plt

sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Run_YKA_WRA/map.txt'
with open(sta_file, 'r') as file:
    lines = file.readlines()
print('    ' + str(len(lines)) + ' events read from ' + sta_file)
# Load station coords into arrays
event_count = len(lines)
station_index = range(event_count)
event_index = np.zeros(event_count)
comp_name    = []
comp_change  = []
comp_ev1     = np.zeros(event_count)
comp_ev2     = np.zeros(event_count)
comp_lat     = np.zeros(event_count)
comp_lon     = np.zeros(event_count)
comp_depth   = np.zeros(event_count)
for ii in station_index:
    line = lines[ii]
    split_line = line.split()
    # print('input station ' + split_line[0])
    comp_name.append(   split_line[0])
    comp_change.append( split_line[1])
    comp_ev1[ii]   = float(split_line[2])
    comp_ev2[ii]   = float(split_line[3])
    comp_lat[ii]   = float(split_line[4])
    comp_lon[ii]   = float(split_line[5])
    comp_depth[ii] = float(split_line[6])
    # print(line)
    # print(lines[ii])
    # print('lat is ' + str(comp_lat[ii]) + ' lon is  ' + str(comp_lon[ii]))

#%% plot data vs prediction
fig_index = 1
plt.figure(1, figsize=(5,10))
# fig, ax = plt.subplots(1)

#ax = fig.gca()
# ax.set_xticks(np.arange(0, 360, 30))
# ax.set_yticks(np.arange(0, 100, 30))

#    FILL NUMBERS
print('lat0 is ' + str(comp_lat[0]) + ' lon0 is  ' + str(comp_lon[0]))
min_lat = np.min(comp_lat)
max_lat = np.max(comp_lat)
min_lon = np.min(comp_lon)
max_lon = np.max(comp_lon)
# print('lat is ' + str(min_lat) + ' to ' + str(max_lat))
# print('lon is ' + str(min_lon) + ' to ' + str(max_lon))

plt.xlim(min_lon - 0.05 * (max_lon - min_lon),max_lon + 0.05 * (max_lon - min_lon))
plt.ylim(min_lat - 0.05 * (max_lat - min_lat),max_lat + 0.05 * (max_lat - min_lat))
for ii in station_index:
    plt.scatter(comp_lon[ii], comp_lat[ii], c='k', s=100, alpha=1, marker='.')
    plt.text(comp_lon[ii] + 0.01 * (max_lon - min_lon), comp_lat[ii] + 0.01 * (max_lat - min_lat), comp_name[ii])
plt.grid()
plt.rc('grid', linestyle="-", color='black')
plt.xlabel('Longitude (°)')
plt.ylabel('Latitude')
plt.title('Repeater events')
# plt.legend([])
# plt.colorbar()
plt.show()
