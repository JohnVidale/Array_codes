#!/usr/bin/env python
# plot predicted vs observed PcP, ScP, PKiKP slownesses
# John Vidale 7/2022
# -*- coding: utf-8 -*-
"""
Created on June 7, 2022
Compares predicted and observed PKiKP slownesses
@author: vidale
"""
import os
import numpy as np
import matplotlib.pyplot as plt

print('Starting')

array_name = 'WRA'    # array - WRA or YKA
phase_name = 'ICS'  # phase - PKiKP or ICS

sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/slowness_' + array_name + '_' + phase_name + '.txt'
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

done_first_fail = False
done_first_bad  = False
done_first_seen = False
# c = ax.scatter(0.02, 0.01,  color='purple', s=100, alpha=0.75)
for ii in range(event_count):   # read file
    if event_names[ii] == 'PKiKP_no':
        if done_first_fail:
            c = ax.scatter(event_PE[ii], event_PN[ii], color='purple', s=100, alpha=0.75)
        else:
            done_first_fail = True
            c = ax.scatter(event_PE[ii], event_PN[ii], color='purple', s=100, alpha=0.75, label = 'Not seen')
    elif event_names[ii] == 'PKiKP_bad':
        if done_first_bad:
            c = ax.scatter(event_PE[ii], event_PN[ii], color='yellow', s=100, alpha=0.75)
        else:
            done_first_bad = True
            c = ax.scatter(event_PE[ii], event_PN[ii], color='yellow', s=100, alpha=0.75, label = 'Corrupt')
    else:
        if done_first_seen:
            c = ax.scatter(event_PE[ii], event_PN[ii], color='blue', s=100, alpha=0.75)
            c = ax.scatter(event_OE[ii], event_ON[ii],  color='red', s=100, alpha=0.75)
        else:
            done_first_seen = True
            c = ax.scatter(event_PE[ii], event_PN[ii], color='blue', s=100, alpha=0.75, label = 'Predicted')
            c = ax.scatter(event_OE[ii], event_ON[ii],  color='red', s=100, alpha=0.75, label = 'Observed')
        c = ax.plot([event_PE[ii], event_OE[ii]], [event_PN[ii], event_ON[ii]], color='black', linewidth = 2)
circle1 = plt.Circle((0, 0), 0.019, color='black', fill=False)
ax.add_artist(circle1)  #inner core limit
circle1 = plt.Circle((0, 0), 0.040, color='black', fill=False)
ax.add_artist(circle1)  #inner core limit
ax.set_aspect('equal')
plt.xlabel('East Slowness (s/km)')
plt.ylabel('North Slowness (s/km)')
plt.legend(loc="upper right")
plt.xlim(-0.025,0.025)
plt.ylim(-0.025,0.025)
plt.title(array_name + ' Predicted vs observed slowness of ' + phase_name)
plt.show()
os.chdir('/Users/vidale/Documents/Research/IC/Plots_hold')
plt.savefig('beam_distortion_' + array_name + '_' + phase_name)
