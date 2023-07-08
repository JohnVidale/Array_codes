#!/usr/bin/env python
# plot predicted vs observed PcP, ScP, PKiKP slownesses
# John Vidale 7/2022
# -*- coding: utf-8 -*-
"""
Created on June 7, 2022
Compares predicted and observed PKiKP slownesses
@author: vidale
"""
import numpy as np
import matplotlib.pyplot as plt

print('Starting')

sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/slowness_Hinet.txt'
with open(sta_file, 'r') as file:
    lines = file.readlines()
event_count = len(lines)

print(str(event_count) + ' lines read from ' + sta_file)
# Load station coords into arrays
station_index   = range(event_count)
event_names     = []
event_index     = np.zeros(event_count)
event_PE        = np.zeros(event_count)
event_PN        = np.zeros(event_count)
event_OE        = np.zeros(event_count)
event_ON        = np.zeros(event_count)

for ii in station_index:   # read file
    line = lines[ii]
    split_line = line.split()

    event_names.append(split_line[0])
    event_index[ii]  = float(split_line[1])
    event_PE[ii]     = float(split_line[2])
    event_PN[ii]     = float(split_line[3])
    event_OE[ii]     = float(split_line[4])
    event_ON[ii]     = float(split_line[5])

print(event_names[0])

#    fig_index = 4
fig = plt.figure()
ax = fig.add_subplot(111)

c = ax.scatter(event_PE, event_PN, color='blue', s=100, alpha=0.75, label = 'Predicted')
c = ax.scatter(event_OE, event_ON,  color='red', s=100, alpha=0.75, label = 'Observed')
# c = ax.scatter(0.02, 0.01,  color='purple', s=100, alpha=0.75)
for ii in range(event_count):   # read file
    if ii == 0 or ii == 20 or ii == 40:
        if event_names[ii] == 'PcP':
            c = ax.plot([event_PE[ii], event_OE[ii]], [event_PN[ii], event_ON[ii]], color='black', label = 'PcP', linewidth = 4)
        elif event_names[ii] == 'ScP':
            c = ax.plot([event_PE[ii], event_OE[ii]], [event_PN[ii], event_ON[ii]], color='orange', label = 'ScP', linewidth = 4)
        elif event_names[ii] == 'PKiKP':
            c = ax.plot([event_PE[ii], event_OE[ii]], [event_PN[ii], event_ON[ii]], color='purple', label = 'PKiKP', linewidth = 4)
    else:
        if event_names[ii] == 'PcP':
            c = ax.plot([event_PE[ii], event_OE[ii]], [event_PN[ii], event_ON[ii]], color='black', linewidth = 4)
        elif event_names[ii] == 'ScP':
            c = ax.plot([event_PE[ii], event_OE[ii]], [event_PN[ii], event_ON[ii]], color='orange', linewidth = 4)
        elif event_names[ii] == 'PKiKP':
            c = ax.plot([event_PE[ii], event_OE[ii]], [event_PN[ii], event_ON[ii]], color='purple', linewidth = 4)
circle1 = plt.Circle((0, 0), 0.019, color='black', fill=False)
ax.add_artist(circle1)  #inner core limit
circle1 = plt.Circle((0, 0), 0.040, color='black', fill=False)
ax.add_artist(circle1)  #inner core limit
plt.xlabel('East Slowness (s/km)')
plt.ylabel('North Slowness (s/km)')
plt.legend(loc="upper left")

# ax.set_xmax(0.025)
# ax.set_xmin(-0.025)
# ax.set_ymax(0.025)
# ax.set_ymin(-0.025)
ax.grid(True)
plt.title('Predicted vs observed slowness of PcP, ScP, and PKiKP')
plt.show()
