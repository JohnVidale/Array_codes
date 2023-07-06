#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on May 10th, 2023
Makes a timeline of DF waveform and ddt changes for ILAR & YKA
@author: vidale
"""
import numpy as np
import random
import matplotlib.pyplot as plt
import matplotlib.patches as patches
# import pandas as pd
from obspy import UTCDateTime

sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Run_YKA_WRA/event_pairs_data.txt'
with open(sta_file, 'r') as file:
    lines = file.readlines()
print('    ' + str(len(lines)) + ' events read from ' + sta_file)
# Load station coords into arrays
pair_count = len(lines)
pair_index = range(pair_count)
pair_name    = []
similarity   = np.zeros(pair_count)
YKA_offset1  = np.zeros(pair_count)
YKA_offset2  = np.zeros(pair_count)
change_Y     = np.zeros(pair_count)
change_I     = np.zeros(pair_count)
shift_Y      = np.zeros(pair_count)
shift_I      = np.zeros(pair_count)
shift_I_YS   = np.zeros(pair_count)
event1       = np.zeros(pair_count)
event2       = np.zeros(pair_count)
date1        = np.zeros(pair_count)
date2        = np.zeros(pair_count)
lat          = np.zeros(pair_count)
lon          = np.zeros(pair_count)

#%% Parameters
do_YKA_change     = True
do_ILAR_change    = False
do_YKA_shift      = False
do_ILAR_shift     = False
do_ILAR_YS_shift  = False
which_plots = (do_YKA_change, do_ILAR_change, do_YKA_shift, do_ILAR_shift, do_ILAR_YS_shift)
do_both_sets = True   # include first initial and later data sets
do_only_sim  = False  # include only most similar repetitions

plot_Y_vs_I    = False
plot_Y_vs_I_YS = False
plot_I_vs_I_YS = False

#%%  read in table
for i in pair_index:
    line = lines[i]
    split_line = line.split()
    pair_name.append(   split_line[0])
    similarity[i] = int(split_line[1])   # did YKA  waveform change?
    YKA_offset1[i] = float(split_line[2])   # did YKA  waveform change?
    YKA_offset2[i] = float(split_line[3])   # did YKA  waveform change?
    change_Y[i]   = int(split_line[4])   # did YKA  waveform change?
    change_I[i]   = int(split_line[5])   # did ILAR waveform change?
    shift_Y[i]    = float(split_line[6]) # John's     YKA  time shift estimation
    shift_I[i]    = float(split_line[7]) # John's     ILAR time shift estimation
    shift_I_YS[i] = float(split_line[8]) # Y&S's 2020 YKA  time shift estimation
    event1[i]     = int(split_line[9])
    event2[i]     = int(split_line[10])
    t1            = UTCDateTime(split_line[11])
    t2            = UTCDateTime(split_line[12])
    lat[i]        = float(split_line[13])
    lon[i]        = float(split_line[14])

    date1[i]      = t1.year + t1.month/12.
    date2[i]      = t2.year + t2.month/12.
    # print('Pair is ' + pair_name[i] + ' changes are ' + str(change_Y[i]) + ' ' + str(change_I[i]) +
    #       ' events are ' + str(event1[i]) + ' ' + str(event2[i]) + ' t1 is ' + str(date1[i]) +
    #       ' t2 is ' + str(date2[i]) + ' lat & lon are'  + str(lat[i]) + ' ' + str(lon[i]))

min_index = 0
max_index = pair_count + 1
min_year = 1990
max_year = 2021

#%% plot waveform changes
if do_YKA_change or do_ILAR_change or do_YKA_shift or do_ILAR_shift or do_ILAR_YS_shift:
    cnt = 0
    while cnt < 5:
        cnt += 1
        print(str(cnt) + '  ' + str(which_plots[cnt-1]))
        if which_plots[cnt-1] == True:
            fig_index = cnt
            plt.figure(cnt, figsize=(10,10))
            plt.xlim( min_year,  max_year)
            plt.ylim( min_index, max_index)

            # mark multiplets
            colorword = 'lightyellow'
            ax = plt.gca()
            full = 5; count = 3 # purple-red
            rect = patches.Rectangle([date1[full],pair_count - full],date2[full]-date1[full],count-1, linewidth=1, edgecolor='w', facecolor= colorword)
            ax.add_patch(rect)

            full = 14; count = 3 # blue-magenta
            rect = patches.Rectangle([date1[full],pair_count - full],date2[full]-date1[full],count-1, linewidth=1, edgecolor='w', facecolor= colorword)
            ax.add_patch(rect)

            full = 21; count = 6 # blue-red
            rect = patches.Rectangle([date1[full],pair_count - full],date2[full]-date1[full],count-1, linewidth=1, edgecolor='w', facecolor= colorword)
            ax.add_patch(rect)

            full = 26; count = 3 # blue-magenta
            rect = patches.Rectangle([date1[full],pair_count - full],date2[full]-date1[full],count-1, linewidth=1, edgecolor='w', facecolor= colorword)
            ax.add_patch(rect)

            full = 37; count = 10 # light green
            rect = patches.Rectangle([date1[full],pair_count - full],date2[full]-date1[full],count-1, linewidth=1, edgecolor='w', facecolor= colorword)
            ax.add_patch(rect)

            full = 40; count = 3 # red-purple
            rect = patches.Rectangle([date1[full],pair_count - full],date2[full]-date1[full],count-1, linewidth=1, edgecolor='w', facecolor= colorword)
            ax.add_patch(rect)

            full = 47; count = 3 # magenta-blue
            rect = patches.Rectangle([date1[full],pair_count - full],date2[full]-date1[full],count-1, linewidth=1, edgecolor='w', facecolor='lightgray')
            ax.add_patch(rect)

            full = 53; count = 6 # orange-blue
            rect = patches.Rectangle([date1[full],pair_count - full],date2[full]-date1[full],count-1, linewidth=1, edgecolor='w', facecolor= colorword)
            ax.add_patch(rect)

            full = 59; count = 6 # red=orange
            rect = patches.Rectangle([date1[full],pair_count - full],date2[full]-date1[full],count-1, linewidth=1, edgecolor='w', facecolor= colorword)
            ax.add_patch(rect)

            full = 64; count = 3 # pink
            rect = patches.Rectangle([date1[full],pair_count - full],date2[full]-date1[full],count-1, linewidth=1, edgecolor='w', facecolor= colorword)
            ax.add_patch(rect)

            full = 71; count = 3 # blue
            rect = patches.Rectangle([date1[full],pair_count - full],date2[full]-date1[full],count-1, linewidth=1, edgecolor='w', facecolor= colorword)
            ax.add_patch(rect)

            full = 76; count = 3 # magenta-red
            rect = patches.Rectangle([date1[full],pair_count - full],date2[full]-date1[full],count-1, linewidth=1, edgecolor='w', facecolor='lightgray')
            ax.add_patch(rect)

            full = 86; count = 10 # teal
            rect = patches.Rectangle([date1[full],pair_count - full],date2[full]-date1[full],count-1, linewidth=1, edgecolor='w', facecolor= colorword)
            ax.add_patch(rect)

            # ax = plt.gca()
            full = 96; count = 10 # purple
            rect = patches.Rectangle([date1[full],pair_count - full],date2[full]-date1[full],count-1, linewidth=1, edgecolor='w', facecolor='lightgray')
            ax.add_patch(rect)

            # divide clusters
            plt.plot([1990,2021], [97.5,97.5], c='0.2', alpha=1, marker='.',linewidth=1.0, linestyle='dotted')
            plt.plot([1990,2021], [61.5,61.5], c='0.2', alpha=1, marker='.',linewidth=1.0, linestyle='dotted')
            plt.plot([1990,2021], [51.5,51.5], c='0.2', alpha=1, marker='.',linewidth=1.0, linestyle='dotted')
            plt.plot([1990,2021], [31.5,31.5], c='0.2', alpha=1, marker='.',linewidth=1.0, linestyle='dotted')
            plt.plot([1990,2021], [ 2.5, 2.5], c='0.2', alpha=1, marker='.',linewidth=1.0, linestyle='dotted')

    #%% plot requested timelines
            for i in pair_index:
                if pair_name[i][0] == 'P':
                    pair_str = pair_name[i][1:3]
                else:
                    pair_str = pair_name[i][0:3]
                index_num = int(pair_str)
                if cnt == 1 or cnt == 2:
                    if ((similarity[i] == 1) or (do_only_sim == False)):
                        if cnt == 1:
                            change = change_Y[i]
                        elif cnt == 2:
                            change = change_I[i]
                        if (index_num < 50) or do_both_sets:
                            if cnt == 1 or cnt == 2:
                                ii = pair_count - i
                                x = [date1[i], date2[i]]
                                line_index = [ii,ii]
                                if (index_num < 50) or do_both_sets:
                                    if change == 2:
                                        plt.plot(x, line_index, c='r', alpha=1, marker='.',linewidth=3.0)
                                    elif change == 1:
                                        plt.plot(x, line_index, c='y', alpha=1, marker='.',linewidth=3.0)
                                    elif change == 0:
                                        plt.plot(x, line_index, c='g', alpha=1, marker='.',linewidth=3.0)
                                    elif change == -1:
                                        plt.plot(x, line_index, c='0.8', alpha=1, marker='.',linewidth=2.0)
                                    elif change == -2:
                                        plt.plot(x, line_index, c='0.8', alpha=1, marker='.',linewidth=2.0)
                                    elif change == -4:
                                        plt.plot(x, line_index, c='c', alpha=1, marker='.',linewidth=3.0)
                                    plt.text(date1[i] - 1, ii, pair_name[i])
                if cnt == 3 or cnt == 4 or cnt == 5:
                    if ((similarity[i] == 1) or (do_only_sim == False)):
                        if cnt == 3:
                            shift = shift_Y[i]
                        elif cnt == 4:
                            shift = shift_I[i]
                        elif cnt == 5:
                            shift = shift_I_YS[i]
                        if shift > (-1):
                            if (index_num < 50) or do_both_sets:
                                ii = pair_count - i
                                x = [date1[i], date2[i]]
                                line_index = [ii,ii]
                                if   shift < -0.025 and shift > -0.5:
                                    plt.plot(x, line_index, c='r', alpha=1, marker='.',linewidth=3.0)
                                elif shift > -0.025 and shift <  0.025:
                                    plt.plot(x, line_index, c='c', alpha=1, marker='.',linewidth=3.0)
                                elif shift >  0.025 and shift <  0.5:
                                    plt.plot(x, line_index, c='b', alpha=1, marker='.',linewidth=3.0)
                                plt.text(date1[i] - 3, ii, pair_name[i] + ' ' + str(shift))

            plt.xlabel('Year', fontsize=20)
            plt.ylabel('South to North', fontsize=20)
            if cnt ==1:
                plt.title('YKA waveform variation over the years', fontsize=25)
            elif cnt ==2:
                plt.title('ILAR waveform variation over the years', fontsize=25)
            elif cnt ==3:
                plt.title('YKA ddt over the years', fontsize=25)
            plt.legend([])
plt.show()

# %% plot one ddt estimate against another

if plot_Y_vs_I      == True:
    fig_index = 6
    plt.figure(6, figsize=(10,10))
    plt.xlim(-0.155,0.155)
    plt.ylim(-0.155,0.155)
    for i in pair_index:
        # print(str(i) + ' Y ' + str(shift_Y[i]) + ' I ' + str(shift_I[i]) + ' I_YS ' + str(shift_I_YS[i]))
        r1 = (random.random() - 0.5)/200  # jitter to enable recognition of identical points
        r2 = (random.random() - 0.5)/200
        if (shift_Y[i] > -1 and shift_I[i] > -1) and (similarity[i] == 1 or do_only_sim == False):
            plt.scatter( shift_Y[i] + r1, shift_I[i] + r2, c='b', s=50, alpha=1, marker='.')
            print(str(i) + ' Y ' + str(shift_Y[i]) + ' I ' + str(shift_I[i]) + ' plot it. ')
        elif (shift_Y[i] == -1 and shift_I[i] == -1) and (similarity[i] == 1 or do_only_sim == False):
            print(str(i) + ' Y ' + str(shift_Y[i]) + ' I ' + str(shift_I[i]) + ' nothing to do. ')
        elif (shift_Y[i] == -1) and (similarity[i] == 1 or do_only_sim == False):
            print(str(i) + ' Y ' + str(shift_Y[i]) + ' I ' + str(shift_I[i]) + ' plot it on x axis. ')
            plt.scatter(              r1, shift_I[i] + r2, c='r', s=50, alpha=1, marker='.')
        elif (shift_I[i] == -1) and (similarity[i] == 1 or do_only_sim == False):
            print(str(i) + ' Y ' + str(shift_Y[i]) + ' I ' + str(shift_I[i]) + ' plot it on y axis. ')
            plt.scatter( shift_Y[i] + r1,              r2, c='g', s=50, alpha=1, marker='.')
    plt.grid()
    plt.rc('grid', linestyle="-", color='black')
    plt.xlabel('YKA shift estimate (s)')
    plt.ylabel('ILAR shift estimate (s)')
    plt.title('YKA vs ILAR ddt comparison')

if plot_Y_vs_I_YS   == True:
    fig_index = 7
    plt.figure(7, figsize=(10,10))
    plt.xlim(-0.155,0.155)
    plt.ylim(-0.155,0.155)
    for i in pair_index:
        r1 = (random.random() - 0.5)/200
        r2 = (random.random() - 0.5)/200
        if (shift_Y[i] > -1 and shift_I_YS[i] > -1) and (similarity[i] == 1 or do_only_sim == False):
            plt.scatter( shift_Y[i] + r1, shift_I_YS[i] + r2, c='b', s=50, alpha=1, marker='.')
            print(str(i) + ' Y ' + str(shift_Y[i]) + ' I ' + str(shift_I_YS[i]) + ' plot it. ')
        elif (shift_Y[i] == -1 and shift_I_YS[i] == -1) and (similarity[i] == 1 or do_only_sim == False):
            print(str(i) + ' Y ' + str(shift_Y[i]) + ' I ' + str(shift_I_YS[i]) + ' nothing to do. ')
        elif (shift_Y[i] == -1) and (similarity[i] == 1 or do_only_sim == False):
            print(str(i) + ' Y ' + str(shift_Y[i]) + ' I ' + str(shift_I_YS[i]) + ' plot it on x axis. ')
            plt.scatter(              r1, shift_I_YS[i] + r2, c='r', s=50, alpha=1, marker='.')
        elif (shift_I_YS[i] == -1) and (similarity[i] == 1 or do_only_sim == False):
            print(str(i) + ' Y ' + str(shift_Y[i]) + ' I ' + str(shift_I_YS[i]) + ' plot it on y axis. ')
            plt.scatter( shift_Y[i] + r1,              r2, c='g', s=50, alpha=1, marker='.')
    plt.grid()
    plt.rc('grid', linestyle="-", color='black')
    plt.xlabel('YKA shift estimate (s)')
    plt.ylabel('ILAR-YS shift estimate (s)')
    plt.title('YKA vs ILAR-YS ddt comparison')

if plot_I_vs_I_YS   == True:
    fig_index = 8
    plt.figure(8, figsize=(10,10))
    plt.xlim(-0.155,0.155)
    plt.ylim(-0.155,0.155)
    for i in pair_index:
        r1 = (random.random() - 0.5)/200
        r2 = (random.random() - 0.5)/200
        if (shift_I[i] > -1 and shift_I_YS[i] > -1) and (similarity[i] == 1 or do_only_sim == False):
            plt.scatter( shift_I[i] + r1, shift_I_YS[i] + r2, c='b', s=50, alpha=1, marker='.')
            print(str(i) + ' Y ' + str(shift_I[i]) + ' I ' + str(shift_I_YS[i]) + ' plot it. ')
        elif (shift_I[i] == -1 and shift_I_YS[i] == -1) and (similarity[i] == 1 or do_only_sim == False):
            print(str(i) + ' Y ' + str(shift_I[i]) + ' I ' + str(shift_I_YS[i]) + ' nothing to do. ')
        elif (shift_I[i] == -1) and (similarity[i] == 1 or do_only_sim == False):
            print(str(i) + ' Y ' + str(shift_I[i]) + ' I ' + str(shift_I_YS[i]) + ' plot it on x axis. ')
            plt.scatter(              r1, shift_I_YS[i] + r2, c='r', s=50, alpha=1, marker='.')
        elif (shift_I_YS[i] == -1) and (similarity[i] == 1 or do_only_sim == False):
            print(str(i) + ' Y ' + str(shift_I[i]) + ' I ' + str(shift_I_YS[i]) + ' plot it on y axis. ')
            plt.scatter( shift_I[i] + r1,              r2, c='g', s=50, alpha=1, marker='.')
        # if (shift_I[i] > -1 and shift_I[i] > -1) and (similarity[i] == 1) or (do_only_sim == False):
        #     plt.scatter( shift_I[i], shift_I_YS[i], c='b', s=100, alpha=1, marker='.')
        # elif (shift_I[i] == -1 and shift_I[i]  > -1) and (similarity[i] == 1) or (do_only_sim == False):
        #     plt.scatter(          0, shift_I[i], c='r', s=100, alpha=1, marker='.')
        # elif (shift_I[i] >  -1 and shift_I[i] == -1) and (similarity[i] == 1) or (do_only_sim == False):
        #     plt.scatter( shift_Y[i],          0, c='g', s=100, alpha=1, marker='.')
    plt.grid()
    plt.rc('grid', linestyle="-", color='black')
    plt.xlabel('ILAR shift estimate (s)')
    plt.ylabel('ILAR-YS shift estimate (s)')
    plt.title('ILAR vs ILAR-YS ddt comparison')

plt.show()
