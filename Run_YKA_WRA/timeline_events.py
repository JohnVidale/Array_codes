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
do_YKA_change     = False
do_ILAR_change    = False

do_YKA_shift      = False
do_ILAR_shift     = False
do_ILAR_YS_shift  = False
which_plots = (do_YKA_change, do_ILAR_change, do_YKA_shift, do_ILAR_shift, do_ILAR_YS_shift)
do_both_sets = True   # include first initial and later data sets
do_only_sim  = False  # include only most similar repetitions
do_label = False

use_N             = True
use_S             = True
NSsplit           = 60

plot_Y_vs_I    = False
plot_Y_vs_I_YS = False
plot_I_vs_I_YS = False

do_YKA_change2     = False # 0 to 1 1st to 2nd event
do_ILAR_change2    = False
which_plots2 = (do_YKA_change2, do_ILAR_change2)

do_YKA_change3     = False # matches and mismatches in 1st date
do_ILAR_change3    = False
blowup = True
which_plots3 = (do_YKA_change3, do_ILAR_change3)

do_YKA_change4     = True # matches and mismatchs by interval
do_ILAR_change4 = True
combine = True
if do_only_sim:
    time_blank = 5
else:
    time_blank = 0
which_plots4 = (do_YKA_change4, do_ILAR_change4)

do_YKA_change5     = False # start and stops compared to prediction
do_ILAR_change5    = False # Miaki plot
which_plots5 = (do_YKA_change5, do_ILAR_change5)

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
max_year = 2024

#%% plot waveform changes
if do_YKA_change or do_ILAR_change or do_YKA_shift or do_ILAR_shift or do_ILAR_YS_shift:
    cnt = 0
    while cnt < 5:
        cnt += 1
        print(str(cnt) + '  ' + str(which_plots[cnt-1]))
        if which_plots[cnt-1] == True:
            fig_index = cnt
            plt.figure(cnt, figsize=(10,15))
            plt.xlim( min_year,  max_year)
            plt.ylim( min_index, max_index)

            #%% --mark multiplets
            color_non_rev = 'whitesmoke'
            color_rev = 'whitesmoke'
            ax = plt.gca()
            ax.tick_params(right=True, labelright=True, top=True, labeltop=True)
            full = 5; count = 3 # purple-red
            rect = patches.Rectangle([date1[full],pair_count - full],date2[full]-date1[full],count-1, linewidth=1, edgecolor='w', facecolor= color_non_rev)
            ax.add_patch(rect)

            full = 14; count = 3 # blue-magenta
            rect = patches.Rectangle([date1[full],pair_count - full],date2[full]-date1[full],count-1, linewidth=1, edgecolor='w', facecolor= color_non_rev)
            ax.add_patch(rect)

            full = 25; count = 10 # blue-red
            rect = patches.Rectangle([date1[full],pair_count - full],date2[full]-date1[full],count-1, linewidth=1, edgecolor='w', facecolor= color_non_rev)
            ax.add_patch(rect)

            full = 30; count = 3 # dark green
            rect = patches.Rectangle([date1[full],pair_count - full],date2[full]-date1[full],count-1, linewidth=1, edgecolor='w', facecolor= color_non_rev)
            ax.add_patch(rect)

            full = 46; count = 15 # light green
            rect = patches.Rectangle([date1[full],pair_count - full],date2[full]-date1[full],count-1, linewidth=1, edgecolor='w', facecolor=color_rev)
            ax.add_patch(rect)

            full = 49; count = 3 # red-purple
            rect = patches.Rectangle([date1[full],pair_count - full],date2[full]-date1[full],count-1, linewidth=1, edgecolor='w', facecolor= color_non_rev)
            ax.add_patch(rect)

            full = 59; count = 6 # magenta-blue
            rect = patches.Rectangle([date1[full],pair_count - full],date2[full]-date1[full],count-1, linewidth=1, edgecolor='w', facecolor=color_rev)
            ax.add_patch(rect)

            full = 65; count = 6 # orange-blue
            rect = patches.Rectangle([date1[full],pair_count - full],date2[full]-date1[full],count-1, linewidth=1, edgecolor='w', facecolor= color_non_rev)
            ax.add_patch(rect)

            full = 75; count = 10 # red=orange
            rect = patches.Rectangle([date1[full],pair_count - full],date2[full]-date1[full],count-1, linewidth=1, edgecolor='w', facecolor= color_non_rev)
            ax.add_patch(rect)

            full = 80; count = 3 # pink
            rect = patches.Rectangle([date1[full],pair_count - full],date2[full]-date1[full],count-1, linewidth=1, edgecolor='w', facecolor= color_non_rev)
            ax.add_patch(rect)

            full = 89; count = 6 # dark blue-blue
            rect = patches.Rectangle([date1[full],pair_count - full],date2[full]-date1[full],count-1, linewidth=1, edgecolor='w', facecolor= color_non_rev)
            ax.add_patch(rect)

            full = 99; count = 10 # blue
            rect = patches.Rectangle([date1[full],pair_count - full],date2[full]-date1[full],count-1, linewidth=1, edgecolor='w', facecolor= color_non_rev)
            ax.add_patch(rect)

            full = 105; count = 6 # dark blue
            rect = patches.Rectangle([date1[full],pair_count - full],date2[full]-date1[full],count-1, linewidth=1, edgecolor='w', facecolor=color_rev)
            ax.add_patch(rect)

            full = 109; count = 3 # magenta-red
            rect = patches.Rectangle([date1[full],pair_count - full],date2[full]-date1[full],count-1, linewidth=1, edgecolor='w', facecolor=color_rev)
            ax.add_patch(rect)

            full = 119; count = 10 # teal
            rect = patches.Rectangle([date1[full],pair_count - full],date2[full]-date1[full],count-1, linewidth=1, edgecolor='w', facecolor=color_rev)
            ax.add_patch(rect)

            # ax = plt.gca()
            full = 140; count = 21 # purple
            rect = patches.Rectangle([date1[full],pair_count - full],date2[full]-date1[full],count-1, linewidth=1, edgecolor='w', facecolor=color_rev)
            ax.add_patch(rect)

            # divide clusters
            plt.plot([min_year,max_year], [141.5,141.5], c='0.2', alpha=1, marker='.',linewidth=1.0, linestyle='dotted')
            plt.plot([min_year,max_year], [ 96.5, 96.5], c='0.2', alpha=1, marker='.',linewidth=1.0, linestyle='dotted')
            plt.plot([min_year,max_year], [ 83.5, 83.5], c='0.2', alpha=1, marker='.',linewidth=1.0, linestyle='dotted')
            plt.plot([min_year,max_year], [ 59.5, 59.5], c='0.2', alpha=1, marker='.',linewidth=1.0, linestyle='dotted')
            plt.plot([min_year,max_year], [  2.5,  2.5], c='0.2', alpha=1, marker='.',linewidth=1.0, linestyle='dotted')

            # mark timelines
            plt.plot([2001  ,2001  ], [0,pair_count], c='0.2', alpha=1, marker='.',linewidth=1.0, linestyle='dotted')
            plt.plot([2004  ,2004  ], [0,pair_count], c='0.2', alpha=1, marker='.',linewidth=1.0, linestyle='dotted')
            plt.plot([2014  ,2014  ], [0,pair_count], c='0.2', alpha=1, marker='.',linewidth=1.0, linestyle='dotted')
            plt.plot([2020.3,2020.3], [0,pair_count], c='0.2', alpha=1, marker='.',linewidth=1.0, linestyle='dotted')

    #%% --plot requested timelines
            for i in pair_index:
                if pair_name[i][0] == 'P':
                    pair_str = pair_name[i][1:3]
                else:
                    pair_str = pair_name[i][0:3]
                index_num = int(pair_str)
                if cnt == 1 or cnt == 2:
                    if ((similarity[i] == 1) or (do_only_sim == False)):
                        if (use_N and lat[i] > -57) or (use_S and lat[i] <= -57):
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
                                            plt.plot(x, line_index, c='lightgreen', alpha=1, markersize = 10, marker='.',linewidth=3.0)
                                        elif change == 1:
                                            plt.plot(x, line_index, c='deepskyblue', alpha=1, markersize = 10, marker='.',linewidth=3.0)
                                        elif change == 0:
                                            plt.plot(x, line_index, c='crimson', alpha=1, markersize = 10, marker='.',linewidth=3.0)
                                        elif change == -1:
                                            plt.plot(x, line_index, c='0.8', alpha=1, markersize = 5, marker='.',linewidth=2.0)
                                        elif change == -2:
                                            plt.plot(x, line_index, c='0.8', alpha=1, markersize = 5, marker='.',linewidth=2.0)
                                        elif change == -4:
                                            plt.plot(x, line_index, c='c', alpha=1, markersize = 5, marker='.',linewidth=3.0)
                                        plt.text(date1[i] - 1, ii, pair_name[i])
                if cnt == 3 or cnt == 4 or cnt == 5:
                    if ((similarity[i] == 1) or (do_only_sim == False)):
                        if (use_N and lat[i] > -57) or (use_S and lat[i] <= -57):
                            if cnt == 3:
                                shift = shift_Y[i]
                            elif cnt == 4:
                                shift = shift_I[i]
                            elif cnt == 5:
                                shift = shift_I_YS[i]
                            else:
                                shift = 0 # to avoid eval complaint
                                print('cnt not 3, 4, nor 5')
                                exit(-1)
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

            plt.xlabel('Date', fontsize=20)
            plt.ylabel('Repeating pairs arranged South to North', fontsize=20)
            # if cnt ==1:
            #     plt.title('YKA waveform variation over the years', fontsize=25)
            # elif cnt ==2:
            #     plt.title('ILAR waveform variation over the years', fontsize=25)
            # elif cnt ==3:
            #     plt.title('YKA ddt over the years', fontsize=25)
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
    plt.figure(8, figsize=(5,10))
    plt.xlim(-0.06,0.06)
    plt.ylim(-0.06,0.16)
    # plt.xlim(-0.05,0.1)
    # plt.ylim(-0.05,0.1)
    for i in pair_index:
        r1 = (random.random() - 0.5)/200
        r2 = (random.random() - 0.5)/200
        if (shift_I[i] > -1 and shift_I_YS[i] > -1) and (similarity[i] == 1 or do_only_sim == False):
            plt.scatter( shift_I[i] + r1, shift_I_YS[i] + r2, c='b', s=200, alpha=1, marker='.')
            print(str(i) + ' Y ' + str(shift_I[i]) + ' I ' + str(shift_I_YS[i]) + ' plot it. ')
        elif (shift_I[i] == -1 and shift_I_YS[i] == -1) and (similarity[i] == 1 or do_only_sim == False):
            print(str(i) + ' Y ' + str(shift_I[i]) + ' I ' + str(shift_I_YS[i]) + ' nothing to do. ')
        elif (shift_I[i] == -1) and (similarity[i] == 1 or do_only_sim == False):
            print(str(i) + ' Y ' + str(shift_I[i]) + ' I ' + str(shift_I_YS[i]) + ' plot it on x axis. ')
            plt.scatter(              r1, shift_I_YS[i] + r2, c='r', s=150, alpha=1, marker='.')
        elif (shift_I_YS[i] == -1) and (similarity[i] == 1 or do_only_sim == False):
            print(str(i) + ' Y ' + str(shift_I[i]) + ' I ' + str(shift_I_YS[i]) + ' plot it on y axis. ')
            plt.scatter( shift_I[i] + r1,              r2, c='g', s=150, alpha=1, marker='.')

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

#%% plot temporal connections
if do_YKA_change2 or do_ILAR_change2:
    cnt = 0
    while cnt < 2:
        cnt += 1
        print(str(cnt) + '  ' + str(which_plots2[cnt-1]))
        if which_plots2[cnt-1] == True:
            fig_index = cnt + 8
            plt.figure(fig_index, figsize=(15,10))
            plt.xlim( min_year,  max_year)
            plt.ylim( -0.1, 1)

            # divide clusters
            plt.plot([2001,2001], [-0.1, 1.1], c='0.2', alpha=1, marker='.',linewidth=1.0, linestyle='dotted')
            plt.plot([2004,2004], [-0.1, 1.1], c='0.2', alpha=1, marker='.',linewidth=1.0, linestyle='dotted')
            plt.plot([2010,2010], [-0.1, 1.1], c='0.2', alpha=1, marker='.',linewidth=1.0, linestyle='dotted')

            for i in pair_index:
                if pair_name[i][0] == 'P':
                    pair_str = pair_name[i][1:3]
                else:
                    pair_str = pair_name[i][0:3]
                index_num = int(pair_str)
                if cnt == 1:
                    change = change_Y[i]
                elif cnt == 2:
                    change = change_I[i]
                if (date2[i] - date1[i] < 5) and change == 2 :
                    print('quick change ' + pair_name[i])
                if ((index_num < 50) or do_both_sets) and ((date2[i] - date1[i] > 5) or change == 2):
                    if ((similarity[i] == 1) or (do_only_sim == False)):
                        if (use_N and lat[i] > -57) or (use_S and lat[i] <= -57):
                            xy1 = [date1[i],date2[i]]
                            xy2 = [0,1]
                            if change == 2:
                                plt.plot(xy1, xy2, c='r', alpha=1, marker='.',linewidth=3.0)
                            elif change == 1:
                                plt.plot(xy1, xy2, c='y', alpha=1, marker='.',linewidth=3.0)
                            elif change == 0:
                                plt.plot(xy1, xy2, c='g', alpha=1, marker='.',linewidth=3.0)
                            elif change == -1:
                                plt.plot(xy1, xy2, c='0.8', alpha=1, marker='.',linewidth=2.0)
                            elif change == -2:
                                plt.plot(xy1, xy2, c='0.8', alpha=1, marker='.',linewidth=2.0)
                            plt.text(date1[i], -0.05, pair_name[i], rotation = 'vertical')
                            plt.text(date2[i],  1.01, pair_name[i], rotation = 'vertical')

            plt.xlabel('Year', fontsize=20)
            plt.ylabel('1st -> 2nd event in pair', fontsize=20)
            ax = plt.gca()
            ax.tick_params(left=False, labelleft=False, top=True, labeltop=True)
            if cnt ==1:
                plt.title('YKA temporal connections over the years', fontsize=25)
            elif cnt ==2:
                plt.title('ILAR temporal connections over the years', fontsize=25)
            plt.legend([])
plt.show()

#%% plot temporal connections by year of first event
if do_YKA_change3 or do_ILAR_change3:
    cnt = 0
    while cnt < 2:
        cnt += 1
        print(str(cnt) + '  ' + str(which_plots3[cnt-1]))
        if which_plots3[cnt-1] == True:
            fig_index = cnt + 10
            plt.figure(fig_index, figsize=(12,12))
            if cnt == 1:
                plt.xlim(1990, 2024)
                plt.ylim(2024, 1990)
            if cnt == 2:
                plt.xlim(1995, 2024)
                if blowup == True:
                    plt.ylim(2008, 2002)
                else:
                    plt.ylim(2019, 1995)

            # divide clusters
            plt.plot([2001,2001], [min_year,  max_year], c='0.2', alpha=1, marker='.',linewidth=1.0, linestyle='dotted')
            plt.plot([2004,2004], [min_year,  max_year], c='0.2', alpha=1, marker='.',linewidth=1.0, linestyle='dotted')
            plt.plot([2010,2010], [min_year,  max_year], c='0.2', alpha=1, marker='.',linewidth=1.0, linestyle='dotted')

            for i in pair_index:
                if pair_name[i][0] == 'P':
                    pair_str = pair_name[i][1:3]
                else:
                    pair_str = pair_name[i][0:3]
                index_num = int(pair_str)
                if cnt == 1:
                    change = change_Y[i]
                elif cnt == 2:
                    change = change_I[i]
                if ((similarity[i] == 1) or (do_only_sim == False)) and ((index_num < 50) or do_both_sets):
                    if (use_N and lat[i] > -57) or (use_S and lat[i] <= -57):
                        xy1 = [date1[i],date2[i]]
                        xy2 = [date1[i],date1[i]]
                        if change == 0:
                            plt.plot(xy1, xy2, c='g', alpha=1, marker='.',linewidth=3.0, markersize = '20')
                        elif change == 1:
                            plt.plot(xy1, xy2, c='y', alpha=1, marker='.',linewidth=1.0, markersize = '10')
                        elif change == 2:
                            plt.plot(xy1, xy2, c='r', alpha=1, marker='.',linewidth=3.0, markersize = '20')
                        # r1 = (random.random() - 0.5)/200  # jitter to enable recognition of identical points
                        # r2 = (random.random() - 0.5)/200
                        if change >= 0:
                            plt.text(date2[i],date1[i], pair_name[i])

            plt.xlabel('Year', fontsize=20)
            plt.ylabel('Date of 1st event in pair', fontsize=20)
            ax = plt.gca()
            ax.tick_params(right=True, labelright=True, top=True, labeltop=True)
            if cnt ==1:
                plt.title('YKA matches and mismatches', fontsize=25)
            elif cnt ==2:
                plt.title('ILAR matches and mismatches', fontsize=25)
            plt.legend([])
plt.show()

#%% color and plot temporal connections against years of separation
if do_YKA_change4 or do_ILAR_change4:
    cnt = 0
    while cnt < 2:
        cnt += 1
        if cnt == 1:
            print('YKA  ' + str(cnt) + '  ' + str(which_plots4[cnt-1]))
        elif cnt == 2:
            print('ILAR ' + str(cnt) + '  ' + str(which_plots4[cnt-1]))
        else:
            print('cnt should have been 1 or 2 ' + str(cnt) + '  ' + str(which_plots4[cnt-1]))
        if which_plots4[cnt-1] == True:
            fig_index = cnt + 12
            if cnt == 1 or combine == False:
                plt.figure(fig_index, figsize=(12,12))
                # stupidity to get legend right
                plt.plot(2020, 2021, c = 'crimson',  alpha=1, marker='.',linewidth=2.0, markersize = '18', label = "no change")
                plt.plot(2020, 2021, c = 'lightgreen',  alpha=1, marker='.',linewidth=1.5, markersize = '12', label = "moderate change")
                plt.plot(2020, 2021, c = 'deepskyblue', alpha=1, marker='.',linewidth=1.5, markersize = '12', label = "strong change")

            if cnt == 1 or combine:
                if do_only_sim:
                    max_sep = 28
                else:
                    max_sep = 31
                plt.xlim(1990, 2024)
                plt.ylim(0, max_sep)
            else:
                max_sep = 20
                plt.xlim(1995, 2024)
                plt.ylim(0, max_sep)

            for i in pair_index:
                date_diff = date2[i] - date1[i]
                if pair_name[i][0] == 'P':
                    pair_str = pair_name[i][1:3]
                else:
                    pair_str = pair_name[i][0:3]
                index_num = int(pair_str)
                if cnt == 1:
                    change = change_Y[i]
                elif cnt == 2:
                    change = change_I[i]
                if ((similarity[i] == 1) or (do_only_sim == False)) and ((index_num < 50) or do_both_sets) and date_diff > time_blank:
                # if ((similarity[i] == 1) or (do_only_sim == False)) and ((index_num < 50) or do_both_sets):
                    if (use_N and lat[i] > -57) or (use_S and lat[i] <= -57):
                        if combine and cnt == 2:
                            date_diff += 0.4
                        xy1 = [date1[i],date2[i]]
                        xy2 = [date_diff,date_diff]
                        if change == 0:
                            plt.plot(xy1, xy2, c='crimson', alpha=1, marker='.',linewidth=2.0, markersize = '18', label = "no change")
                        if change == 1:
                            plt.plot(xy1, xy2, c='lightgreen', alpha=1, marker='.',linewidth=1.5, markersize = '12', label = "moderate change")
                        if change == 2:
                            plt.plot(xy1, xy2, c='deepskyblue', alpha=1, marker='.',
                                     linewidth=1.5, markersize='12', label="strong change")
                        # r1 = (random.random() - 0.5)/200  # jitter to enable recognition of identical points
                        # r2 = (random.random() - 0.5)/200
                        if do_label:
                            if change >= 0:
                                plt.text(date1[i],date_diff, pair_name[i])
                                plt.text(date2[i],date_diff, pair_name[i])

            if cnt == 1 and combine:
                plt.xlabel('Date (years))', fontsize=20)
                plt.ylabel('Interval between repeaters (years)', fontsize=20)
                ax = plt.gca()
                # ax.tick_params(right=True, labelright=True, top=True, labeltop=True)
                # plt.title('Both arrays - same vs changing waveforms', fontsize=25)
                rect = patches.Rectangle([1990,0], 34, time_blank, linewidth=time_blank/2.0, edgecolor='w', facecolor= 'lightgray')
                ax.add_patch(rect)
            elif combine == False:
                plt.xlabel('Year', fontsize=20)
                plt.ylabel('1st -> 2nd event in pair', fontsize=20)
                ax = plt.gca()
                ax.tick_params(right=True, labelright=True, top=True, labeltop=True)
                if cnt ==1:
                    plt.title('YKA dates vs time spans', fontsize=25)
                elif cnt ==2:
                    plt.title('ILAR dates vs time spans', fontsize=25)
            # fit lines
            xy1 = [2000.5,2008.5]
            xy2 = [28,0]
            plt.plot(xy1, xy2, c='m', alpha=1, marker='.',linewidth=2.0, markersize = '12')
            xy1 = [2008.5, 2024]
            xy2 = [0,22]
            plt.plot(xy1, xy2, c='m', alpha=1, marker='.',linewidth=2.0, markersize = '12')
            ax.tick_params(axis='both', labelsize=18)
            plt.grid()
            plt.rc('grid', linestyle="-", color='black')
            if do_only_sim:
                plt.legend(['similar', 'somewhat similar', 'different'],loc = 'lower right', fontsize = 20)
            else:
                plt.legend(['similar', 'somewhat similar', 'different'],loc = 'upper right', fontsize = 20)
# plt.show()
#%% plot color temporal misfit vs years of separation
if do_YKA_change5 or do_ILAR_change5:
    print('Combine ' + str(combine))
    cnt = 0  # cnt = 1 for YKA, 2 for ILAR
    while cnt < 2:
        cnt += 1
        print(str(cnt) + '  ' + str(which_plots5[cnt-1]))
        if which_plots5[cnt-1] == True:
            fig_index = cnt + 14
            if cnt == 1 or combine == False:
                plt.figure(fig_index, figsize=(12,12))
                # stupidity to get legend right
                plt.scatter(2020, 2021, c='crimson'   , s=400, alpha=1, marker='.', label = "no change")
                plt.scatter(2020, 2021, c='lightgreen', s=400, alpha=1, marker='.', label = "moderate change")
                plt.scatter(2020, 2021, c='deepskyblue', s=400,
                            alpha=1, marker='.', label="strong change")

            plt.xlim(0, 15)
            plt.ylim(time_blank, 31)

            for i in pair_index:
                date_diff = date2[i] - date1[i]
                if pair_name[i][0] == 'P':
                    pair_str = pair_name[i][1:3]
                else:
                    pair_str = pair_name[i][0:3]
                index_num = int(pair_str)
                if cnt == 1:
                    change = change_Y[i]
                elif cnt == 2:
                    change = change_I[i]
                if ((similarity[i] == 1) or (do_only_sim == False)) and ((index_num < 50) or do_both_sets) and date_diff > time_blank:
                # if ((similarity[i] == 1) or (do_only_sim == False)) and ((index_num < 50) or do_both_sets):
                    if (use_N and lat[i] > -57) or (use_S and lat[i] <= -57):
                        t1_pred = 2008.5 - date_diff*(8/28)
                        t2_pred = 2008.5 + date_diff*(15.5/22)
                        t1_diff = abs(date1[i]-t1_pred)
                        t2_diff = abs(date2[i]-t2_pred)
                        if combine and cnt == 2:
                            date_diff += 0.3
                        if change == 0:
                            plt.scatter((t1_diff + t2_diff)/2, date_diff, c='crimson'   , s=300, alpha=1, marker='.')
                        if change == 1:
                            plt.scatter((t1_diff + t2_diff)/2, date_diff, c='lightgreen', s=300, alpha=1, marker='.')
                        if change == 2:
                            plt.scatter((t1_diff + t2_diff)/2, date_diff, c='deepskyblue'   , s=300, alpha=1, marker='.')
                        if do_label:
                            if (change == 0 or change == 1 or change == 2) and cnt == 1:
                                plt.text(t1_diff, date_diff, pair_name[i])

            if cnt == 1 and combine:
                plt.xlabel('Average deviation from model (years)', fontsize=20)
                plt.ylabel('Interval between repeaters (years)', fontsize=20)
                ax = plt.gca()
                # ax.tick_params(right=True, labelright=True, top=True, labeltop=True)
                # plt.title('Model misfit vs time separation', fontsize=25)
                rect = patches.Rectangle([1990,0], 34, time_blank, linewidth=time_blank/2.0, edgecolor='w', facecolor= 'lightgray')
                ax.add_patch(rect)
            elif combine == False:
                plt.xlabel('Year', fontsize=20)
                plt.ylabel('1st -> 2nd event in pair', fontsize=20)
                ax = plt.gca()
                ax.tick_params(right=True, labelright=True, top=True, labeltop=True)
                if cnt ==1:
                    plt.title('YKA dates vs time spans', fontsize=25)
                elif cnt ==2:
                    plt.title('ILAR dates vs time spans', fontsize=25)
            ax.tick_params(axis='both', labelsize=18)
            plt.grid()
            plt.rc('grid', linestyle="-", color='black')
            plt.legend(['similar waveform', 'somewhat similar', 'different'],loc = 'upper right', fontsize = 20)
plt.show()
