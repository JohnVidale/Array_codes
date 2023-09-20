#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on May 10th, 2023
Makes a timeline of DF waveform and ddt changes for ILAR & YKA
@author: vidale
"""
#%%
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
int_of_simY  = np.zeros(pair_count)
int_of_simI  = np.zeros(pair_count)
consensus    = np.zeros(pair_count)
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
do_val = False
zoomer = True

use_N             = True
use_S             = True
NSsplit           = -57

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

do_YKA_change4     = False # matches and mismatchs by interval
do_ILAR_change4    = False
combine4 = True
which_plots4 = (do_YKA_change4, do_ILAR_change4)

do_YKA_change5     = False # start and stops compared to prediction
do_ILAR_change5    = False # Miaki plot
combine5 = True
which_plots5 = (do_YKA_change5, do_ILAR_change5)

do_YKA_change6  = True  # matches and mismatchs by interval
do_ILAR_change6 = True
combine6 = True
which_plots6 = (do_YKA_change6, do_ILAR_change6)

do_ILAR_change7 = False
do_YKA_change7  = False

do_ILAR_change8 = False
do_YKA_change8  = False

do_YKA_change9  = True  # matches and mismatchs by interval
do_ILAR_change9 = True
combine9 = True
which_plots9 = (do_YKA_change9, do_ILAR_change9)

do_cons_change10 = True

#%%  read in table
for i in pair_index:
    line = lines[i]
    split_line = line.split()
    pair_name.append(   split_line[0])
    similarity[i]  = int(split_line[1])   # one of most similar events?
    int_of_simY[i]  = float(split_line[2])   # interval of similarity for ILAR
    int_of_simI[i] = float(split_line[3])   # interval of similarity for YKA
    consensus[i]   = int(split_line[4])   # believe ILAR, then YKA
    YKA_offset2[i] = float(split_line[5])   # PKIKP time shift on YKA
    change_Y[i]   = int(split_line[6])   # did YKA  waveform change?
    change_I[i]   = int(split_line[7])   # did ILAR waveform change?
    shift_Y[i]    = float(split_line[8]) # John's     YKA  time shift estimation
    shift_I[i]    = float(split_line[9]) # John's     ILAR time shift estimation
    shift_I_YS[i] = float(split_line[10]) # Y&S's 2020 YKA  time shift estimation
    event1[i]     = int(split_line[11])
    event2[i]     = int(split_line[12])
    t1            = UTCDateTime(split_line[13])
    t2            = UTCDateTime(split_line[14])
    lat[i]        = float(split_line[15])
    lon[i]        = float(split_line[16])

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
                        if (use_N and lat[i] > -57) or (use_S and lat[i] <= NSsplit):
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
# plt.show()

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

# plt.show()

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
                        if (use_N and lat[i] > NSsplit) or (use_S and lat[i] <= NSsplit):
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
# plt.show()

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
                    if (use_N and lat[i] > NSsplit) or (use_S and lat[i] <= NSsplit):
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
# plt.show()

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
            if cnt == 1 or combine4 == False:
                plt.figure(fig_index, figsize=(12,12))
                # stupidity to get legend right
                plt.plot(2020, 2021, c = 'crimson',  alpha=1, marker='.',linewidth=2.0, markersize = '18', label = "no change")
                plt.plot(2020, 2021, c = 'lightgreen',  alpha=1, marker='.',linewidth=1.5, markersize = '12', label = "moderate change")
                plt.plot(2020, 2021, c = 'deepskyblue', alpha=1, marker='.',linewidth=1.5, markersize = '12', label = "strong change")

            if cnt == 1 or combine4:
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
                if ((similarity[i] == 1) or (do_only_sim == False)) and ((index_num < 50) or do_both_sets):
                # if ((similarity[i] == 1) or (do_only_sim == False)) and ((index_num < 50) or do_both_sets):
                    if (use_N and lat[i] > NSsplit) or (use_S and lat[i] <= NSsplit):
                        marker_type = '.'
                        if combine4 and cnt == 2:
                            date_diff += 0.4
                            marker_type = '*'
                        xy1 = [date1[i],date2[i]]
                        xy2 = [date_diff,date_diff]
                        if change == 0:
                            plt.plot(xy1, xy2, c='crimson', alpha=1, marker=marker_type,linewidth=2.0, markersize = '18', label = "no change")
                            plt.plot(xy1[1], xy2[1], c='crimson', alpha=1, marker=marker_type,markersize='18', markeredgecolor='yellow', markeredgewidth=2)
                        if change == 1:
                            plt.plot(xy1, xy2, c='lightgreen', alpha=1, marker=marker_type,linewidth=1.5, markersize = '12', label = "moderate change")
                            plt.plot(xy1[1], xy2[1], c='lightgreen', alpha=1, marker=marker_type,markersize='18', markeredgecolor='yellow', markeredgewidth=2)
                        if change == 2:
                            plt.plot(xy1, xy2, c='deepskyblue', alpha=1, marker=marker_type,linewidth=1.5, markersize='12', label="strong change")
                            plt.plot(xy1[1], xy2[1], c='deepskyblue', alpha=1, marker=marker_type,
                                     markersize='18', markeredgecolor='yellow', markeredgewidth=2)
                        # r1 = (random.random() - 0.5)/200  # jitter to enable recognition of identical points
                        # r2 = (random.random() - 0.5)/200
                        if do_label:
                            if change >= 0:
                                plt.text(date1[i],date_diff, pair_name[i])
                                plt.text(date2[i],date_diff, pair_name[i])

            if cnt == 1 and combine4:
                plt.xlabel('Date (years)', fontsize=20)
                plt.ylabel('Interval between repeaters (years)', fontsize=20)
                ax = plt.gca()
            elif combine4 == False:
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
    print('Combine ' + str(combine5))
    cnt = 0  # cnt = 1 for YKA, 2 for ILAR
    while cnt < 2:
        cnt += 1
        print(str(cnt) + '  ' + str(which_plots5[cnt-1]))
        if which_plots5[cnt-1] == True:
            fig_index = cnt + 14
            if cnt == 1 or combine5 == False:
                plt.figure(fig_index, figsize=(12,12))
                # stupidity to get legend right
                plt.scatter(2020, 2021, c='crimson'   , s=400, alpha=1, marker='.', label = "no change")
                plt.scatter(2020, 2021, c='lightgreen', s=400, alpha=1, marker='.', label = "moderate change")
                plt.scatter(2020, 2021, c='deepskyblue', s=400,
                            alpha=1, marker='.', label="strong change")

            plt.xlim(0, 15)
            plt.ylim(0, 31)

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
                if ((similarity[i] == 1) or (do_only_sim == False)) and ((index_num < 50) or do_both_sets):
                # if ((similarity[i] == 1) or (do_only_sim == False)) and ((index_num < 50) or do_both_sets):
                    if (use_N and lat[i] > NSsplit) or (use_S and lat[i] <= NSsplit):
                        t1_pred = 2008.5 - date_diff*(8/28)
                        t2_pred = 2008.5 + date_diff*(15.5/22)
                        t1_diff = abs(date1[i]-t1_pred)
                        t2_diff = abs(date2[i]-t2_pred)
                        if combine5 and cnt == 2:
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

            if cnt == 1 and combine5:
                plt.xlabel('Average deviation from model (years)', fontsize=20)
                plt.ylabel('Interval between repeaters (years)', fontsize=20)
                ax = plt.gca()
            elif combine5 == False:
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

#%% color and plot temporal connections against years of separation
if do_YKA_change6 or do_ILAR_change6:
    cnt = 0
    while cnt < 2:
        cnt += 1
        if cnt == 1:
            print('YKA  ' + str(cnt) + '  ' + str(which_plots6[cnt-1]))
        elif cnt == 2:
            print('ILAR ' + str(cnt) + '  ' + str(which_plots6[cnt-1]))
        else:
            print('cnt should have been 1 or 2 ' +
                  str(cnt) + '  ' + str(which_plots6[cnt-1]))
        if which_plots6[cnt-1] == True:
            fig_index = cnt + 13
            if cnt == 1 or combine6 == False:
                if zoomer:
                    plt.figure(fig_index, figsize=(19*0.6, 12*0.6))
                else:
                    plt.figure(fig_index, figsize=(12, 12))
                # stupidity to get legend right
                if zoomer:
                    plt.plot(2010, c='crimson', marker='.',markersize='12', label="no change")
                    plt.plot(2010, c='lightgreen', marker='.',markersize='12', label="moderate change")
                    plt.plot(2010, c='deepskyblue', marker='.',markersize='12', label="strong change")
                else:
                    plt.plot(1995, c='crimson', marker='.',markersize='12', label="no change")
                    plt.plot(1995, c='lightgreen', marker='.',markersize='12', label="moderate change")
                    plt.plot(1995, c='deepskyblue', marker='.',markersize='12', label="strong change")

            if zoomer:
                minx = 2005
                maxx = 2024
                miny = 2000
                maxy = 2012
            else:
                minx = 1990
                maxx = 2024
                miny = 1990
                maxy = 2024
            plt.xlim(minx, maxx)
            plt.ylim(miny, maxy)

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
                    if ((similarity[i] == 1) or (do_only_sim == False)) and ((index_num < 50) or do_both_sets):
                        if (use_N and lat[i] > NSsplit) or (use_S and lat[i] <= NSsplit) or cnt == 2:
                            marker_type = '.'
                            d1 = date1[i]
                            d2 = date2[i]
                            if combine6 and cnt == 2:
                                if zoomer:
                                    d1 = date1[i] + 0.2
                                    d2 = date2[i] - 0.2
                                else:
                                    d1 = date1[i] + 0.4
                                    d2 = date2[i] - 0.4
                                marker_type = '*'
                            if change == 0:
                                plt.plot(d2, d1, c='crimson', alpha=1, marker=marker_type,linewidth=2.0, markersize='12', label="no change")
                            if change == 1:
                                plt.plot(d2, d1, c='lightgreen', alpha=1, marker=marker_type,linewidth=1.5, markersize='12', label="moderate change")
                            if change == 2:
                                plt.plot(d2, d1, c='deepskyblue', alpha=1, marker=marker_type,linewidth=1.5, markersize='12', label="strong change")
                            if do_label:
                                if change >= 0:
                                    if date1[i] > miny and date1[i] < maxy and date2[i] > minx and date2[i] < maxx:
                                        plt.text(
                                            d2, d1, pair_name[i])

            if cnt == 1 and combine6:
                plt.xlabel('Date of 2nd event (years)', fontsize=20)
                plt.ylabel('Date of 1st event (years)', fontsize=20)
                ax = plt.gca()
                plt.title('Original: Both, dates vs time spans', fontsize=25)
            elif combine6 == False:
                ax = plt.gca()
                ax.tick_params(right=True, labelright=True,top=True, labeltop=True)
                if cnt == 1:
                    plt.title('Original: YKA dates vs time spans', fontsize=25)
                elif cnt == 2:
                    plt.title('Original: ILAR dates vs time spans', fontsize=25)
            # diagonal line
            if zoomer:
                xy1 = [2005, 2012]
                xy2 = [2005, 2012]
                plt.plot(xy1, xy2, c='purple', alpha=1, marker='.',linewidth=2.0, markersize='12')
                xy1 = [2010, 2022]
                xy2 = [2000, 2012]
                plt.plot(xy1, xy2, c='purple', alpha=1, marker='.',linewidth=2.0, markersize='12')
                rect = patches.Rectangle([2004, 2004], 15, 10, linewidth=2.0, edgecolor='orange', facecolor='none')
                ax.add_patch(rect)
                plt.figtext(0.20, 0.52, 'little change?',c='orange', fontsize=20)
                plt.figtext(0.15, 0.64, 'dt = 0 years',c='purple', fontsize=18)
                plt.figtext(0.15, 0.12, 'dt = 10 years',c='purple', fontsize=18)
                plt.figtext(0.56, 0.12, 'slope = 1', fontsize=20)
                plt.figtext(0.74, 0.22, 'slope = 2.5', fontsize=20)
                xy1 = [2008.5, 2024]
                xy2 = [2008.5, 2002.5]
                plt.plot(xy1, xy2, 'k--', alpha=1, marker='.',linewidth=2.0, markersize='12')
                xy1 = [2010, 2020]
                xy2 = [2010, 2000]
                plt.plot(xy1, xy2, 'k--', alpha=1, marker='.',linewidth=1.0, markersize='12')
                plt.xticks(range(2005, 2024, 5))
            else:
                rect = patches.Rectangle([2004, 2004], 15, 15, linewidth=2.0, edgecolor='orange', facecolor='none')
                ax.add_patch(rect)
                # rect = patches.Rectangle([2005, 2000], 19, 12, linewidth=2.0, edgecolor='gray', facecolor='none')
                # ax.add_patch(rect)
                xy1 = [1990, 2024]
                xy2 = [1990, 2024]
                plt.plot(xy1, xy2, c='purple', alpha=1, marker='.',linewidth=2.0, markersize='12')
                xy1 = [1994, 2024]
                xy2 = [1990, 2020]
                plt.plot(xy1, xy2, c='purple', alpha=1, marker='.',linewidth=2.0, markersize='12')
                plt.figtext(0.18, 0.35, 'dt = 0 years',c='purple', fontsize=20)
                plt.figtext(0.18, 0.12, 'dt = 4 years',c='purple', fontsize=20)
                plt.figtext(0.50, 0.70, 'little change?',c='orange', fontsize=20)
                # plt.figtext(0.50, 0.55, 'inset', fontsize=20)
                # plt.figtext(1995, 2000,'dt = 0 years', fontsize=20)  # inscrutable "transform" option for plot co-ords
                # plt.figtext(1995, 1992,'dt = 10 years', fontsize=20)
            ax = plt.gca()
            ax.tick_params(axis='both', labelsize=18)
            plt.grid()
            plt.rc('grid', linestyle="-", color='black')
            plt.legend(['similar', 'somewhat similar', 'different'],loc='upper left', fontsize=20)
#%% ILAR - 1st event date vs 2nd event date
if do_ILAR_change7:
    fig_index = 25
    plt.figure(fig_index, figsize=(22*0.4, 28*0.4))
    # stupidity to get legend right
    plt.plot(2005, c='whitesmoke', marker='.',markersize='12')
    plt.plot(2005, c='silver', marker='.', markersize='12')
    plt.plot(2005, c='limegreen', marker='.',markersize='12')
    plt.plot(2005, c='gold', marker='.',markersize='12')
    plt.plot(2005, c='orange', marker='.',markersize='12')
    plt.plot(2005, c='red', marker='.',markersize='12')

    minx = 2002
    maxx = 2024
    miny = 1996
    maxy = 2024
    plt.xlim(minx, maxx)
    plt.ylim(miny, maxy)

    for i in pair_index:
        if pair_name[i][0] == 'P':
            pair_str = pair_name[i][1:3]
        else:
            pair_str = pair_name[i][0:3]
        same_dur = int_of_simI[i]
        marker_type = '.'
        if same_dur == -2:
            plt.plot(date2[i], date1[i], c='whitesmoke', alpha=1, marker=marker_type,linewidth=1.5, markersize='18', label="moderate change")
        elif same_dur == -1:
            plt.plot(date2[i], date1[i], c='silver', alpha=1, marker=marker_type,linewidth=1.5, markersize='18', label="moderate change")
        elif same_dur == 0:
            plt.plot(date2[i], date1[i], c='limegreen', alpha=1, marker=marker_type,linewidth=2.0, markersize='18', label="no change")
        elif same_dur >= 0 and same_dur <= 3:
            plt.plot(date2[i], date1[i], c='gold', alpha=1, marker=marker_type,linewidth=1.5, markersize='18', label="moderate change")
        elif same_dur > 3 and same_dur < 20:
            plt.plot(date2[i], date1[i], c='orange', alpha=1, marker=marker_type,linewidth=1.5, markersize='18', label="strong change")
        elif same_dur >= 20:
            plt.plot(date2[i], date1[i], c='red', alpha=1, marker=marker_type,linewidth=1.5, markersize='18', label="strong change")

        if do_label:
            if same_dur >= 0:
                plt.text(date2[i], date1[i], pair_str)
        if do_val:
            if same_dur >= 0:
                if same_dur == 99:
                    plt.text(date2[i], date1[i], 'same')
                else:
                    plt.text(date2[i], date1[i], str(int(same_dur)))
                # plt.text(date2[i], date1[i], 'test', verticalalignment='bottom')

    plt.xlabel('Date of 2nd event (years)', fontsize=20)
    plt.ylabel('Date of 1st event (years)', fontsize=20)
    ax = plt.gca()
    ax.tick_params(right=True, labelright=False,top=True, labeltop=False)
    plt.title('7: new ILAR waveform similarity', fontsize=25)

    # time boxes    
    rect = patches.Rectangle([2004, 2004], 15, 15, linewidth=2.0, edgecolor='gray', facecolor='none')
    ax.add_patch(rect)
    
    # diagonal line
    xy1 = [1990, 2024]
    xy2 = [1990, 2024]
    plt.plot(xy1, xy2, c='purple', alpha=1, marker='.',linewidth=2.0, markersize='12')
    xy1 = [2000, 2024]
    xy2 = [1996, 2020]
    plt.plot(xy1, xy2, c='purple', alpha=1, marker='.',linewidth=2.0, markersize='12')
    
    plt.figtext(0.33, 0.53, 'dt = 0 years',c='purple', fontsize=16)
    plt.figtext(0.19, 0.12, 'dt = 4 years',c='purple', fontsize=16)
    plt.figtext(0.30, 0.60, 'little change?',c='gray', fontsize=16)
    ax = plt.gca()
    ax.tick_params(axis='both', labelsize=16)
    plt.xticks(range(2000, 2024, 5))
    plt.yticks(range(2000, 2024, 5))
    plt.grid()
    plt.rc('grid', linestyle="-", color='black')
    plt.legend(["no data", "noisy", "different", "same for up to 3s", "same for 3 to 20s", "same throughout"],loc='upper left', fontsize=16)

#%% YKA - 1st event date vs 2nd event date
if do_YKA_change7:
    fig_index = 26
    plt.figure(fig_index, figsize=(22*0.4, 28*0.4))
    # stupidity to get legend right
    plt.plot(1995, c='whitesmoke', marker='.',markersize='12')
    plt.plot(1995, c='silver', marker='.', markersize='12')
    plt.plot(1995, c='limegreen', marker='.',markersize='12')
    plt.plot(1995, c='red', marker='.',markersize='12')

    minx = 1990
    maxx = 2024
    miny = 1990
    maxy = 2024
    plt.xlim(minx, maxx)
    plt.ylim(miny, maxy)

    for i in pair_index:
        if pair_name[i][0] == 'P':
            pair_str = pair_name[i][1:3]
        else:
            pair_str = pair_name[i][0:3]
        same_dur = int_of_simY[i]
        marker_type = '.'
        if same_dur == -2:
            plt.plot(date2[i], date1[i], c='whitesmoke', alpha=1, marker=marker_type,linewidth=1.5, markersize='18', label="moderate change")
        elif same_dur == -1:
            plt.plot(date2[i], date1[i], c='silver', alpha=1, marker=marker_type,linewidth=1.5, markersize='18', label="moderate change")
        elif same_dur == 0:
            plt.plot(date2[i], date1[i], c='limegreen', alpha=1, marker=marker_type,linewidth=2.0, markersize='18', label="no change")
        elif same_dur >= 20:
            plt.plot(date2[i], date1[i], c='red', alpha=1, marker=marker_type,linewidth=1.5, markersize='18', label="strong change")

        if do_label:
            if same_dur >= 0:
                plt.text(date2[i], date1[i], pair_str)
        if do_val:
            if same_dur >= 0:
                if same_dur == 99:
                    plt.text(date2[i], date1[i], 'same')
                else:
                    plt.text(date2[i], date1[i], str(int(same_dur)))
                # plt.text(date2[i], date1[i], 'test', verticalalignment='bottom')

    plt.xlabel('Date of 2nd event (years)', fontsize=20)
    plt.ylabel('Date of 1st event (years)', fontsize=20)
    ax = plt.gca()
    ax.tick_params(right=True, labelright=False,top=True, labeltop=False)
    plt.title('7: new YKA waveform similarity', fontsize=25)

    # time boxes    
    rect = patches.Rectangle([2004, 2004], 15, 15, linewidth=2.0, edgecolor='gray', facecolor='none')
    ax.add_patch(rect)
    
    # diagonal line
    xy1 = [1990, 2024]
    xy2 = [1990, 2024]
    plt.plot(xy1, xy2, c='purple', alpha=1, marker='.',linewidth=2.0, markersize='12')
    xy1 = [1994, 2024]
    xy2 = [1990, 2020]
    plt.plot(xy1, xy2, c='purple', alpha=1, marker='.',linewidth=2.0, markersize='12')
    
    plt.figtext(0.20, 0.38, 'dt = 0 years',c='purple', fontsize=16)
    plt.figtext(0.25, 0.12, 'dt = 4 years',c='purple', fontsize=16)
    plt.figtext(0.45, 0.60, 'little change?',c='gray', fontsize=16)
    ax = plt.gca()
    ax.tick_params(axis='both', labelsize=16)
    plt.xticks(range(1990, 2024, 5))
    plt.yticks(range(1990, 2024, 5))
    plt.grid()
    plt.rc('grid', linestyle="-", color='black')
    plt.legend(["no data", "noisy", "different", "same"],loc='upper left', fontsize=16)

#%% ILAR - 1st event date vs 2nd event date
if do_ILAR_change8:
    fig_index = 27
    plt.figure(fig_index, figsize=(22*0.4, 28*0.4))
    # stupidity to get legend right
    plt.plot(2005, c='whitesmoke', marker='.',markersize='12')
    plt.plot(2005, c='silver', marker='.', markersize='12')
    plt.plot(2005, c='limegreen', marker='.',markersize='12')
    plt.plot(2005, c='gold', marker='.',markersize='12')
    plt.plot(2005, c='red', marker='.',markersize='12')

    minx = 2002
    maxx = 2024
    miny = 1996
    maxy = 2024
    plt.xlim(minx, maxx)
    plt.ylim(miny, maxy)

    for i in pair_index:
        if pair_name[i][0] == 'P':
            pair_str = pair_name[i][1:3]
        else:
            pair_str = pair_name[i][0:3]
        same_dur = change_I[i]
        marker_type = '.'
        if same_dur == -2:
            plt.plot(date2[i], date1[i], c='whitesmoke', alpha=1, marker=marker_type,linewidth=1.5, markersize='18', label="moderate change")
        elif same_dur == -1:
            plt.plot(date2[i], date1[i], c='silver', alpha=1, marker=marker_type,linewidth=1.5, markersize='18', label="moderate change")
        elif same_dur == 2:
            plt.plot(date2[i], date1[i], c='limegreen', alpha=1, marker=marker_type,linewidth=2.0, markersize='18', label="no change")
        elif same_dur == 1:
            plt.plot(date2[i], date1[i], c='gold', alpha=1, marker=marker_type,linewidth=1.5, markersize='18', label="moderate change")
        elif same_dur == 0:
            plt.plot(date2[i], date1[i], c='red', alpha=1, marker=marker_type,linewidth=1.5, markersize='18', label="strong change")

        if do_label:
            if same_dur >= 0:
                plt.text(date2[i], date1[i], pair_str)
        if do_val:
            if same_dur >= 0:
                if same_dur == 99:
                    plt.text(date2[i], date1[i], 'same')
                else:
                    plt.text(date2[i], date1[i], str(int(same_dur)))
                # plt.text(date2[i], date1[i], 'test', verticalalignment='bottom')

    plt.xlabel('Date of 2nd event (years)', fontsize=20)
    plt.ylabel('Date of 1st event (years)', fontsize=20)
    ax = plt.gca()
    ax.tick_params(right=True, labelright=False,top=True, labeltop=False)
    plt.title('8: original ILAR waveform similarity', fontsize=25)

    # time boxes    
    rect = patches.Rectangle([2004, 2004], 15, 15, linewidth=2.0, edgecolor='gray', facecolor='none')
    ax.add_patch(rect)
    
    # diagonal line
    xy1 = [1990, 2024]
    xy2 = [1990, 2024]
    plt.plot(xy1, xy2, c='purple', alpha=1, marker='.',linewidth=2.0, markersize='12')
    xy1 = [2000, 2024]
    xy2 = [1996, 2020]
    plt.plot(xy1, xy2, c='purple', alpha=1, marker='.',linewidth=2.0, markersize='12')
    
    plt.figtext(0.33, 0.53, 'dt = 0 years',c='purple', fontsize=16)
    plt.figtext(0.19, 0.12, 'dt = 4 years',c='purple', fontsize=16)
    plt.figtext(0.30, 0.60, 'little change?',c='gray', fontsize=16)
    ax = plt.gca()
    ax.tick_params(axis='both', labelsize=16)
    plt.xticks(range(2000, 2024, 5))
    plt.yticks(range(2000, 2024, 5))
    plt.grid()
    plt.rc('grid', linestyle="-", color='black')
    plt.legend(["no data", "noisy", "different", "somewhat similar", "similar"],loc='upper left', fontsize=16)

#%% ILAR - 1st event date vs 2nd event date
if do_YKA_change8:
    fig_index = 28
    plt.figure(fig_index, figsize=(22*0.4, 28*0.4))
    # stupidity to get legend right
    plt.plot(1995, c='whitesmoke', marker='.',markersize='12')
    plt.plot(1995, c='silver', marker='.', markersize='12')
    plt.plot(1995, c='limegreen', marker='.',markersize='12')
    plt.plot(1995, c='gold', marker='.',markersize='12')
    plt.plot(1995, c='red', marker='.',markersize='12')

    minx = 1990
    maxx = 2024
    miny = 1990
    maxy = 2024
    plt.xlim(minx, maxx)
    plt.ylim(miny, maxy)

    for i in pair_index:
        if pair_name[i][0] == 'P':
            pair_str = pair_name[i][1:3]
        else:
            pair_str = pair_name[i][0:3]
        same_dur = change_Y[i]
        marker_type = '.'
        if same_dur == -2:
            plt.plot(date2[i], date1[i], c='whitesmoke', alpha=1, marker=marker_type,linewidth=1.5, markersize='18', label="no data")
        elif same_dur == -1:
            plt.plot(date2[i], date1[i], c='silver', alpha=1, marker=marker_type,linewidth=1.5, markersize='18', label="noisy")
        elif same_dur == 2:
            plt.plot(date2[i], date1[i], c='limegreen', alpha=1, marker=marker_type,linewidth=1.5, markersize='18', label="different")
        elif same_dur == 1:
            plt.plot(date2[i], date1[i], c='gold', alpha=1, marker=marker_type,linewidth=1.5, markersize='18', label="same for up to 3s")
        elif same_dur == 0:
            plt.plot(date2[i], date1[i], c='red', alpha=1, marker=marker_type,linewidth=2.0, markersize='18', label="same throughout")

        if do_label:
            if same_dur >= 0:
                plt.text(date2[i], date1[i], pair_str)
        if do_val:
            if same_dur >= 0:
                if same_dur == 99:
                    plt.text(date2[i], date1[i], 'same')
                else:
                    plt.text(date2[i], date1[i], str(int(same_dur)))
                # plt.text(date2[i], date1[i], 'test', verticalalignment='bottom')

    plt.xlabel('Date of 2nd event (years)', fontsize=20)
    plt.ylabel('Date of 1st event (years)', fontsize=20)
    ax = plt.gca()
    ax.tick_params(right=True, labelright=False,top=True, labeltop=False)
    plt.title('8: original YKA waveform similarity', fontsize=25)

    # time boxes    
    rect = patches.Rectangle([2004, 2004], 15, 15, linewidth=2.0, edgecolor='gray', facecolor='none')
    ax.add_patch(rect)
    
    # diagonal line
    xy1 = [1990, 2024]
    xy2 = [1990, 2024]
    plt.plot(xy1, xy2, c='purple', alpha=1, marker='.',linewidth=2.0, markersize='12')
    xy1 = [1994, 2024]
    xy2 = [1990, 2020]
    plt.plot(xy1, xy2, c='purple', alpha=1, marker='.',linewidth=2.0, markersize='12')
    
    plt.figtext(0.20, 0.38, 'dt = 0 years',c='purple', fontsize=16)
    plt.figtext(0.25, 0.12, 'dt = 4 years',c='purple', fontsize=16)
    plt.figtext(0.45, 0.60, 'little change?',c='gray', fontsize=16)
    ax = plt.gca()
    ax.tick_params(axis='both', labelsize=16)
    plt.xticks(range(1990, 2024, 5))
    plt.yticks(range(1990, 2024, 5))
    plt.grid()
    plt.rc('grid', linestyle="-", color='black')
    plt.legend(["no data", "noisy", "different", "somewhat similar", "similar"],loc='upper left', fontsize=16)

#%% color and plot temporal connections against years of separation
if do_YKA_change9 or do_ILAR_change9:
    cnt = 0
    while cnt < 2:
        cnt += 1
        if cnt == 1:
            print('YKA  ' + str(cnt) + '  ' + str(which_plots9[cnt-1]))
        elif cnt == 2:
            print('ILAR ' + str(cnt) + '  ' + str(which_plots9[cnt-1]))
        else:
            print('cnt should have been 1 or 2 ' +
                  str(cnt) + '  ' + str(which_plots9[cnt-1]))
        if which_plots9[cnt-1] == True:
            fig_index = 51
            if cnt == 1 or combine9 == False:
                if zoomer:
                    plt.figure(fig_index, figsize=(12, 6))
                else:
                    plt.figure(fig_index, figsize=(12, 12))
                # stupidity to get legend right
                plt.plot(2005, c='gainsboro', marker='.',markersize='12')
                plt.plot(2005, c='silver', marker='.', markersize='12')
                plt.plot(2005, c='limegreen', marker='.',markersize='12')
                # plt.plot(2005, c='gold', marker='.',markersize='12')
                # plt.plot(2005, c='orange', marker='.',markersize='12')
                plt.plot(2005, c='red', marker='.',markersize='12')

            if zoomer:
                minx = 2003
                maxx = 2024
                miny = 2002
                maxy = 2012
            else:
                minx = 1990
                maxx = 2024
                miny = 1990
                maxy = 2024
            plt.xlim(minx, maxx)
            plt.ylim(miny, maxy)

            for i in pair_index:
                if pair_name[i][0] == 'P':
                    pair_str = pair_name[i][1:3]
                else:
                    pair_str = pair_name[i][0:3]
                index_num = int(pair_str)
                if cnt == 1:
                    change = int_of_simY[i]
                elif cnt == 2:
                    change = int_of_simI[i]
                if ((similarity[i] == 1) or (do_only_sim == False)) and ((index_num < 50) or do_both_sets):
                    if ((similarity[i] == 1) or (do_only_sim == False)) and ((index_num < 50) or do_both_sets):
                        if (use_N and lat[i] > NSsplit) or (use_S and lat[i] <= NSsplit) or cnt == 2:
                            marker_type = '.'
                            marker_size = 18
                            d1 = date1[i]
                            d2 = date2[i]
                            if combine6 and cnt == 2:
                                if zoomer:
                                    d1 = date1[i] + 0.25
                                    d2 = date2[i] - 0.25
                                else:
                                    d1 = date1[i] + 0.4
                                    d2 = date2[i] - 0.4
                                marker_type = '*'
                                marker_size = 12
                            if change == -2:
                                plt.plot(d2, d1, c='gainsboro', alpha=1, marker=marker_type,linewidth=1.5, markersize=marker_size, label="moderate change")
                            elif change == -1:
                                plt.plot(d2, d1, c='silver', alpha=1, marker=marker_type,linewidth=1.5, markersize=marker_size, label="moderate change")
                            elif change == 0:
                                plt.plot(d2, d1, c='limegreen', alpha=1, marker=marker_type,linewidth=2.0, markersize=marker_size, label="no change")
                            # elif change >= 0 and change <= 3:
                            #     plt.plot(d2, d1, c='gold', alpha=1, marker=marker_type,linewidth=1.5, markersize=marker_size, label="moderate change")
                            # elif change > 3 and change < 20:
                            #     plt.plot(d2, d1, c='orange', alpha=1, marker=marker_type,linewidth=1.5, markersize=marker_size, label="strong change")
                            elif change >= 0 and change <= 3:
                                plt.plot(d2, d1, c='limegreen', alpha=1, marker=marker_type,linewidth=1.5, markersize=marker_size, label="moderate change")
                            elif change > 3 and change < 20:
                                plt.plot(d2, d1, c='limegreen', alpha=1, marker=marker_type,linewidth=1.5, markersize=marker_size, label="strong change")
                            elif change >= 20:
                                plt.plot(d2, d1, c='red', alpha=1, marker=marker_type,linewidth=1.5, markersize=marker_size, label="strong change")
                            if do_label:
                                if change >= 0:
                                    if date1[i] > miny and date1[i] < maxy and date2[i] > minx and date2[i] < maxx:
                                        plt.text(
                                            d2, d1, pair_name[i])

            if cnt == 1 and combine6:
                plt.xlabel('2nd event (date)', fontsize=20)
                plt.ylabel('1st event (date)', fontsize=20)
                ax = plt.gca()
                if not zoomer:
                    plt.title('Both arrays', fontsize=25)
            elif combine6 == False:
                ax = plt.gca()
                ax.tick_params(right=True, labelright=True,top=True, labeltop=True)
                if cnt == 1:
                    plt.title('YKA', fontsize=25)
                elif cnt == 2:
                    plt.title('ILAR', fontsize=25)
            # diagonal lines
            rect = patches.Rectangle([2004, 2004], 15.5, 15.5, linewidth=2.0, edgecolor='orange', facecolor='none')
            ax.add_patch(rect)
            rect = patches.Rectangle([2009, 2009], 10.5, 10.5, linewidth=2.0, edgecolor='darkorange', facecolor='none')
            ax.add_patch(rect)
            xy1 = [1990, 2024]
            xy2 = [1990, 2024]
            plt.plot(xy1, xy2, c='purple', alpha=1, marker='.',linewidth=2.0, markersize='12')
            xy1 = [1994, 2024]
            xy2 = [1990, 2020]
            plt.plot(xy1, xy2, c='purple', alpha=1, marker='.',linewidth=2.0, markersize='12')
            if zoomer:
                ax.xaxis.set_ticks(np.arange(2005, 2021, 5))
                ax.yaxis.set_ticks(np.arange(2005, 2011, 5))
                plt.figtext(0.22, 0.35, 'dt = 0 years',c='purple', fontsize=20)
                plt.figtext(0.26, 0.12, 'dt = 4 years',c='purple', fontsize=20)
                plt.figtext(0.18, 0.45, 'no ILAR change',c='orange', fontsize=20)
                plt.figtext(0.35, 0.72, 'no YKA & ILAR change',c='darkorange', fontsize=20)
            else:
                plt.figtext(0.18, 0.35, 'dt = 0 years',c='purple', fontsize=20)
                plt.figtext(0.26, 0.12, 'dt = 4 years',c='purple', fontsize=20)
                plt.figtext(0.45, 0.65, 'no ILAR change',c='orange', fontsize=20)
                plt.figtext(0.57, 0.75, 'no YKA & ILAR change',c='darkorange', fontsize=20)
            ax = plt.gca()
            ax.tick_params(axis='both', labelsize=18)
            plt.grid()
            plt.rc('grid', linestyle="-", color='black')
            plt.legend(["no data", "noisy", "different", "similar"],loc='upper left', fontsize=16)
            # plt.legend(["no data", "noisy", "different", "same for up to 3s", "same for 3 to 20s", "same throughout"],loc='upper left', fontsize=16)

#%% ILAR - 1st event date vs 2nd event date
if do_cons_change10:
    fig_index = 29
    plt.figure(fig_index, figsize=(22*0.4, 22*0.4))
    # stupidity to get legend right
    plt.plot(1995, c='silver', marker='.', markersize='12')
    plt.plot(1995, c='limegreen', marker='.',markersize='12')
    plt.plot(1995, c='red', marker='.',markersize='12')

    minx = 1990
    maxx = 2024
    miny = 1990
    maxy = 2024
    plt.xlim(minx, maxx)
    plt.ylim(miny, maxy)

    for i in pair_index:
        if pair_name[i][0] == 'P':
            pair_str = pair_name[i][1:3]
        else:
            pair_str = pair_name[i][0:3]
        same_dur = consensus[i]
        marker_type = '.'
        if same_dur == -1:
            plt.plot(date2[i], date1[i], c='silver', alpha=1, marker=marker_type,linewidth=1.5, markersize='18')
        elif same_dur == 0:
            plt.plot(date2[i], date1[i], c='limegreen', alpha=1, marker=marker_type,linewidth=1.5, markersize='18')
        elif same_dur == 1:
            plt.plot(date2[i], date1[i], c='red', alpha=1, marker=marker_type,linewidth=2.0, markersize='18')

        if do_label:
            if same_dur >= 0:
                plt.text(date2[i], date1[i], pair_str)
        if do_val:
            if same_dur >= 0:
                if same_dur == 99:
                    plt.text(date2[i], date1[i], 'same')
                else:
                    plt.text(date2[i], date1[i], str(int(same_dur)))
                # plt.text(date2[i], date1[i], 'test', verticalalignment='bottom')

    plt.xlabel('2nd event (date)', fontsize=20)
    plt.ylabel('1st event (date)', fontsize=20)
    ax = plt.gca()
    ax.tick_params(right=True, labelright=False,top=True, labeltop=False)
    plt.title('Best guess', fontsize=25)

    # time boxes    
    rect = patches.Rectangle([2004, 2004], 15, 15, linewidth=2.0, edgecolor='gray', facecolor='none')
    ax.add_patch(rect)
    
    # diagonal line
    xy1 = [1990, 2024]
    xy2 = [1990, 2024]
    plt.plot(xy1, xy2, c='purple', alpha=1, marker='.',linewidth=2.0, markersize='12')
    xy1 = [1994, 2024]
    xy2 = [1990, 2020]
    plt.plot(xy1, xy2, c='purple', alpha=1, marker='.',linewidth=2.0, markersize='12')
    
    plt.figtext(0.20, 0.38, 'dt = 0 years',c='purple', fontsize=16)
    plt.figtext(0.25, 0.12, 'dt = 4 years',c='purple', fontsize=16)
    plt.figtext(0.45, 0.60, 'little change',c='gray', fontsize=16)
    ax = plt.gca()
    ax.tick_params(axis='both', labelsize=16)
    plt.xticks(range(1990, 2024, 5))
    plt.yticks(range(1990, 2024, 5))
    plt.grid()
    plt.rc('grid', linestyle="-", color='black')
    plt.legend(["no reliable estimate", "different", "similar"], fontsize=16)

# show all plots
plt.show()
