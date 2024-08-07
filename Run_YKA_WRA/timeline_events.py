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
import pandas as pd
import os

# look up pair of earthquakes and time shifts in pairs
def search_df(df, column, value, partial_match=True):
    df = df.astype({column:'string'})
    if partial_match:
        return df.loc[df[column].str.contains(value, na=False)]
    else:
        return df.loc[df[column] == value]

df = pd.read_excel('/Users/vidale/Documents/GitHub/Array_codes/Files/ICevents_full.xlsx', sheet_name='pairs')

pair_count = len(df.index)
pair_range = range(pair_count)
pair_index  = np.zeros(pair_count)
pair_name    = []
multiplet    = []
event1       = np.zeros(pair_count)
event2       = np.zeros(pair_count)
date1        = np.zeros(pair_count)
date2        = np.zeros(pair_count)
lat          = np.zeros(pair_count)
lon          = np.zeros(pair_count)
depth        = np.zeros(pair_count)

# shift_Y      = np.zeros(pair_count)
# shift_I      = np.zeros(pair_count)
# shift_I_YS   = np.zeros(pair_count)
similarity   = np.zeros(pair_count)
PKPpre_sim   = np.zeros(pair_count)
# int_of_simY  = np.zeros(pair_count)
int_of_simI  = np.zeros(pair_count)
# change_Y     = np.zeros(pair_count)
# change_I     = np.zeros(pair_count)
consensus    = np.zeros(pair_count)
Ymatch       = np.zeros(pair_count)
Imatch       = np.zeros(pair_count)

#%%  read in table
for i in range(len(df.index)):
    # print('Doing ', str(i))
    # Load station coords into arrays
    pair_index[i]  = int(        df['index'].iloc[i])        # index of pair
    pair_name.append(            df['label'].iloc[i])        # name of pair
    multiplet.append(            df['multiplet'].iloc[i])    # multiplet letter or "no"
    event1[i]      = int(        df['index1'].iloc[i])
    event2[i]      = int(        df['index2'].iloc[i])
    t1             = UTCDateTime(df['date1'].iloc[i])
    t2             = UTCDateTime(df['date2'].iloc[i])
    date1[i]       = t1.year + t1.month/12.
    date2[i]       = t2.year + t2.month/12.
    lat[i]         = float(      df['lat'].iloc[i])
    lon[i]         = float(      df['lon'].iloc[i])
    depth[i]       = float(      df['depth'].iloc[i])
    
    # shift_Y[i]     = float(      df['Y_DF'].iloc[i])         # John's YKA     time shift due to rotation estimation
    # shift_I[i]     = float(      df['I_DF'].iloc[i])         # John's ILAR    time shift due to rotation estimation
    # shift_I_YS[i ] = float(      df['I_Y_S'].iloc[i])        # Y&S's 2020 YKA time shift due to rotation estimation
    similarity[i]  = int(        df['repeat_quality'].iloc[i])   # one of most similar event pairs?
    PKPpre_sim[i] = float(       df['PKPpre_sim'].iloc[i])   #       Whether PKPprecursor is same or different
    # int_of_simY[i] = float(      df['Y_sim_time'].iloc[i])   #       duration of similarity interval for YKA
    int_of_simI[i] = float(      df['I_sim_time'].iloc[i])   #       duration of similarity interval for ILAR
    # change_Y[i]    = int(        df['paper1_change_Y'].iloc[i])     # 3-level did YKA waveform change?
    # change_I[i]    = int(        df['paper1_change_I'].iloc[i])     # 3-level did ILARwaveform change?
    consensus[i]   = int(        df['paper2_Consensus'].iloc[i])    # according to rule - believe ILAR, then YKA
    Ymatch[i]      = float(      df['paper2_Y_match'].iloc[i])      # fixed 2-level interval of similarity for ILAR
    Imatch[i]      = float(      df['paper2_I_match'].iloc[i])      # fixed 2-level interval of similarity for ILAR

#%% Parameters
# Variation or not on N-S vs date
do_YKA_change     = False
do_ILAR_change    = False
do_YKA_shift      = False
do_ILAR_shift     = False

# Time shift on N-S vs date
do_ILAR_YS_shift  = False
which_plots = (do_YKA_change, do_ILAR_change, do_YKA_shift, do_ILAR_shift, do_ILAR_YS_shift)
do_both_sets = True   # include first initial and later data sets
do_only_sim  = False  # include only most similar repetitions

do_label = False
do_val = False
zoomer = False

use_N             = True
use_S             = True
NSsplit           = -57

# Time shifts against each other
plot_Y_vs_I    = False
plot_Y_vs_I_YS = False
plot_I_vs_I_YS = False
do_YKA_change2     = False # 0 to 1 1st to 2nd event
do_ILAR_change2    = False
which_plots2 = (do_YKA_change2, do_ILAR_change2)

# V diagram, matches and mismatchs by interval
do_YKA_change4     = False
do_ILAR_change4    = False
combine4 = True
which_plots4 = (do_YKA_change4, do_ILAR_change4)

# Miaki plot, misfit of start and stops compared to prediction
do_YKA_change5     = False
do_ILAR_change5    = False
combine5 = True
which_plots5 = (do_YKA_change5, do_ILAR_change5)

# YKA similarity, 3 levels of matches, year1 vs year2
do_ILAR_change8 = False  # individual arrays
do_YKA_change8  = False

# Zoomed in plot, both arrays, 2 levels of matches, year1 vs year2, Fig S1 in paper
do_YKA_change9  = False
do_ILAR_change9 = False
combine9 = False  # produces only-YKA  plot if false and YKA True  and ILAR false
                  # produces only-ILAR plot if false and YKA False and ILAR True
                  # shows both if all three are True
which_plots9 = (do_YKA_change9, do_ILAR_change9)

# Best guess plot, year1 vs year2, Fig 2 in paper
do_cons_change10 = False

# Zoomed in plot, best, 2 levels of matches, year1 vs year2, Fig 5 in paper
do_YKA_change11  = False
do_ILAR_change11 = True
label_int = False
# combine11 = False  # produces only-YKA  plot if false and YKA True  and ILAR false
#                   # produces only-ILAR plot if false and YKA False and ILAR True
#                   # shows both if all three are True
# which_plots11 = (do_YKA_change11, do_ILAR_change11)

do_change12 = False # YKA PKPprecursor similarity on 2nd date vs 1st date

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
            for i in pair_range:
                pair_str = pair_index[i]

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
            if cnt ==1:
                plt.title('YKA waveform variation over the years', fontsize=25)
            elif cnt ==2:
                plt.title('ILAR waveform variation over the years', fontsize=25)
            elif cnt ==3:
                plt.title('YKA ddt over the years', fontsize=25)
            plt.legend([])
# plt.show()

# %% plot one ddt estimate against another

if plot_Y_vs_I      == True:
    fig_index = 6
    plt.figure(6, figsize=(10,10))
    plt.xlim(-0.155,0.155)
    plt.ylim(-0.155,0.155)
    for i in pair_range:
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
    for i in pair_range:
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
    for i in pair_range:
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

            for i in pair_range:
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

            for i in pair_range:
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

            for i in pair_range:
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

    for i in pair_range:
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

    for i in pair_range:
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
                plt.plot(2005, c='red', marker='.',markersize='12')

            if zoomer:
                minx = 2003
                maxx = 2024
                miny = 2002
                maxy = 2012
            else:
                minx = 1990
                maxx = 2025
                miny = 1990
                maxy = 2025
            plt.xlim(minx, maxx)
            plt.ylim(miny, maxy)

            for i in pair_range:
                if pair_name[i][0] == 'P':
                    pair_str = pair_name[i][1:3]
                else:
                    pair_str = pair_name[i][0:3]
                index_num = int(pair_str)
                # if cnt == 1:
                #     change = int_of_simY[i]
                # elif cnt == 2:
                    # change = int_of_simI[i]
                if cnt == 1:
                    change = Ymatch[i]
                elif cnt == 2:
                    change = Imatch[i]
                if ((similarity[i] == 1) or (do_only_sim == False)) and ((index_num < 50) or do_both_sets):
                    if ((similarity[i] == 1) or (do_only_sim == False)) and ((index_num < 50) or do_both_sets):
                        if (use_N and lat[i] > NSsplit) or (use_S and lat[i] <= NSsplit) or cnt == 2:
                            marker_type = '.'
                            marker_size = 18
                            d1 = date1[i]
                            d2 = date2[i]
                            if combine9 and cnt == 2:
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
                            # elif change == 0:
                            #     plt.plot(d2, d1, c='limegreen', alpha=1, marker=marker_type,linewidth=2.0, markersize=marker_size, label="no change")
                            # elif change >= 0 and change <= 3:
                            #     plt.plot(d2, d1, c='limegreen', alpha=1, marker=marker_type,linewidth=1.5, markersize=marker_size, label="moderate change")
                            # elif change > 3 and change < 20:
                            #     plt.plot(d2, d1, c='limegreen', alpha=1, marker=marker_type,linewidth=1.5, markersize=marker_size, label="strong change")
                            # elif change >= 20:
                            #     plt.plot(d2, d1, c='red', alpha=1, marker=marker_type,linewidth=1.5, markersize=marker_size, label="strong change")
                            elif change == 0:
                                plt.plot(d2, d1, c='limegreen', alpha=1, marker=marker_type,linewidth=2.0, markersize=marker_size, label="no change")
                            elif change == 1:
                                plt.plot(d2, d1, c='red', alpha=1, marker=marker_type,linewidth=1.5, markersize=marker_size, label="change")
                            if do_label:
                                if date1[i] > miny and date1[i] < maxy and date2[i] > minx and date2[i] < maxx:
                                    plt.text(date2[i], date1[i], int(pair_index[i]))

            if cnt == 1 and combine9:
                plt.xlabel('2nd event (date)', fontsize=20)
                plt.ylabel('1st event (date)', fontsize=20)
                ax = plt.gca()
                if not zoomer:
                    # plt.title('9: Both arrays', fontsize=25)
                    junk = 0
            elif combine9 == False:
                ax = plt.gca()
                ax.tick_params(right=True, labelright=True,top=True, labeltop=True)
                if cnt == 1:
                    plt.title('YKA', fontsize=25)
                elif cnt == 2:
                    plt.title('ILAR', fontsize=25)
            # diagonal lines
            rect = patches.Rectangle([2004, 2004], 14.5, 14.5, linewidth=2.0, edgecolor='orange', facecolor='none')
            ax.add_patch(rect)
            rect = patches.Rectangle([2009, 2009],  9.5,  9.5, linewidth=2.0, edgecolor='darkorange', facecolor='none')
            ax.add_patch(rect)

            # quirky array legend
            if zoomer:
                rect = patches.Rectangle([2004.35, 2006.5], 2.5, 1.9, linewidth=2.0, edgecolor='lightgray', facecolor='none')
                marker_type = '*'
                marker_size = 12
                plt.plot(2005, 2015, c='black', alpha=1, marker=marker_type,linewidth=1.5, markersize=marker_size, label="moderate change")
                marker_type = '.'
                plt.figtext(0.212, 0.541, 'ILAR',c='black', fontsize=16)
                marker_size = 18
                plt.plot(2005, 2007, c='black', alpha=1, marker=marker_type,linewidth=1.5, markersize=marker_size, label="moderate change")
                plt.figtext(0.212, 0.48, 'YKA',c='black', fontsize=16)
            else:
                rect = patches.Rectangle([1990.6, 2014.2], 4, 3, linewidth=1.5, edgecolor='lightgray', facecolor='none')
                marker_type = '*'
                marker_size = 12
                plt.plot(1991.5, 2016.3, c='black', alpha=1, marker=marker_type,linewidth=1.5, markersize=marker_size, label="moderate change")
                marker_type = '.'
                plt.figtext(0.18, 0.685, 'ILAR',c='black', fontsize=16)
                marker_size = 18
                plt.plot(1991.5, 2015, c='black', alpha=1, marker=marker_type,linewidth=1.5, markersize=marker_size, label="moderate change")
                plt.figtext(0.18, 0.655, 'YKA',c='black', fontsize=16)

            ax.add_patch(rect)
            xy1 = [1990, 2025]
            xy2 = [1990, 2025]
            plt.plot(xy1, xy2, c='purple', alpha=1, marker='.',linewidth=2.0, markersize='12')
            xy1 = [1994, 2025]
            xy2 = [1990, 2020]
            plt.plot(xy1, xy2, c='purple', alpha=1, marker='.',linewidth=2.0, markersize='12')
            if zoomer:
                ax.xaxis.set_ticks(np.arange(2005, 2021, 5))
                ax.yaxis.set_ticks(np.arange(2005, 2011, 5))
            #     plt.figtext(0.22, 0.35, 'dt = 0 years',c='purple', fontsize=20)
            #     plt.figtext(0.26, 0.12, 'dt = 4 years',c='purple', fontsize=20)
            #     plt.figtext(0.18, 0.45, 'no ILAR change',c='orange', fontsize=20)
            #     plt.figtext(0.35, 0.72, 'no YKA & ILAR change',c='darkorange', fontsize=20)
            else:
                plt.figtext(0.18, 0.35, 'dt = 0 years',c='purple', fontsize=20)
                plt.figtext(0.26, 0.12, 'dt = 4 years',c='purple', fontsize=20)
                # plt.figtext(0.45, 0.65, 'no ILAR change',c='orange', fontsize=20)
                # plt.figtext(0.57, 0.75, 'no YKA & ILAR change',c='darkorange', fontsize=20)
            ax = plt.gca()
            ax.tick_params(axis='both', labelsize=18)
            plt.grid()
            plt.rc('grid', linestyle="-", color='black')
            plt.legend(["no data", "noisy", "different", "similar"],loc='upper left', fontsize=16)
            # plt.legend(["no data", "noisy", "different", "same for up to 3s", "same for 3 to 20s", "same throughout"],loc='upper left', fontsize=16)
    os.chdir('/Users/vidale/Documents/Research/IC/Plots_hold')
    if zoomer:
        plt.savefig('all_plot_zoom' + '.png')
    else:
        plt.savefig('all_plot' + '.png')

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

    for i in pair_range:
        same_dur = consensus[i]
        marker_type = '.'
        if same_dur == -1:
            plt.plot(date2[i], date1[i], c='silver', alpha=1, marker=marker_type,linewidth=1.5, markersize='18')
        elif same_dur == 0:
            plt.plot(date2[i], date1[i], c='limegreen', alpha=1, marker=marker_type,linewidth=1.5, markersize='18')
        elif same_dur == 1:
            plt.plot(date2[i], date1[i], c='red', alpha=1, marker=marker_type,linewidth=2.0, markersize='18')

        if do_label:
            # if same_dur >= 0:
                plt.text(date2[i], date1[i], int(pair_index[i]))
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
    plt.title('10 Best guess', fontsize=25)

    # time boxes    
    rect = patches.Rectangle([2004, 2004], 14.5, 14.5, linewidth=2.0, edgecolor='gray', facecolor='none')
    ax.add_patch(rect)
    
    # diagonal line
    xy1 = [1990, 2024]
    xy2 = [1990, 2024]
    plt.plot(xy1, xy2, c='purple', alpha=1, marker='.',linewidth=2.0, markersize='12')
    xy1 = [1996, 2024]
    xy2 = [1990, 2018]
    plt.plot(xy1, xy2, c='purple', alpha=1, marker='.',linewidth=2.0, markersize='12')
    
    plt.figtext(0.20, 0.38, 'dt = 0 years',c='purple', fontsize=16)
    plt.figtext(0.25, 0.12, 'dt = 6 years',c='purple', fontsize=16)
    plt.figtext(0.45, 0.60, 'little change',c='gray', fontsize=16)
    ax = plt.gca()
    ax.tick_params(axis='both', labelsize=16)
    plt.xticks(range(1990, 2024, 5))
    plt.yticks(range(1990, 2024, 5))
    plt.grid()
    plt.rc('grid', linestyle="-", color='black')
    plt.legend(["no reliable estimate", "different", "similar"], fontsize=16)
    os.chdir('/Users/vidale/Documents/Research/IC/Plots_hold')
    plt.savefig('best_plot' + '.png')

#%% color and plot temporal connections against years of separation
if do_ILAR_change11 or do_YKA_change11:
    print('made it into do_11')
    fig_index = 52
    plt.figure(fig_index, figsize=(12, 12))
    # stupidity to get legend right
    plt.plot(2005, c='silver', marker='.', markersize='12')
    plt.plot(2005, c='limegreen', marker='.',markersize='12')
    plt.plot(2005, c='red', marker='.',markersize='12')

    minx = 1997
    maxx = 2025
    miny = 1997
    maxy = 2025
    plt.xlim(minx, maxx)
    plt.ylim(miny, maxy)

    for i in pair_range:
        best = False
        best_list = [30,31,34,36,37,38,39,41,43,67,79,109,114,123,319,321,334,339,341,355]
        best = pair_index[i] in best_list

        if pair_name[i][0] == 'P':
            pair_str = pair_name[i][1:3]
        else:
            pair_str = pair_name[i][0:3]
        index_num = int(pair_str)
        if do_ILAR_change11 and (do_YKA_change11 == False):
            change = Imatch[i]
        if (do_ILAR_change11 == False) and do_YKA_change11:
            change = Ymatch[i]
        if ((similarity[i] == 1) or (do_only_sim == False)):
            if ((similarity[i] == 1) or (do_only_sim == False)):
                if (use_N and lat[i] > NSsplit) or (use_S and lat[i] <= NSsplit) or cnt == 2:
                    marker_type = '.'
                    marker_size = 25
                    d1 = date1[i]
                    d2 = date2[i]
                    # if change == -2:
                    #     plt.plot(d2, d1, c='gainsboro', alpha=1, marker=marker_type,linewidth=1.5, markersize=marker_size, label="moderate change")
                    if change == -1:
                        plt.plot(d2, d1, c='silver', alpha=1, marker=marker_type,linewidth=1.5, markersize=marker_size)
                    elif change == 0:
                        plt.plot(d2, d1, c='limegreen', alpha=1, marker=marker_type,linewidth=2.0, markersize=marker_size)
                    elif change == 1:
                        plt.plot(d2, d1, c='red', alpha=1, marker=marker_type,linewidth=1.5, markersize=marker_size)
                    # elif change == 1 and best == False:
                    #     plt.plot(d2, d1, c='red', alpha=1, marker=marker_type,linewidth=1.5, markersize=marker_size)
                    # elif change == 1 and best == True:
                    #     plt.plot(d2, d1, c='black', alpha=1, marker=marker_type,linewidth=1.5, markersize=marker_size)
                    if (date1[i] > miny) and (date1[i] < maxy) and (date2[i] > minx) and (date2[i] < maxx):
                        if do_label and change >= 0:
                            plt.text(date2[i], date1[i], int(pair_index[i]), verticalalignment='bottom')
                        if label_int and change >= 0:
                            plt.text(date2[i], date1[i], int(int_of_simI[i]), verticalalignment='top')

        plt.xlabel('2nd event (date)', fontsize=20)
        plt.ylabel('1st event (date)', fontsize=20)
    ax = plt.gca()
    ax.tick_params(right=True, labelright=True,top=True, labeltop=False)
    if do_ILAR_change11 and (do_YKA_change11 == False):
        plt.title('11 ILAR', fontsize=25)
    if (do_ILAR_change11 == False) and do_YKA_change11:
        plt.title('11 YKA', fontsize=25)
    # diagonal lines
    rect = patches.Rectangle([2004, 2004], 16.3, 16.3, linewidth=2.0, edgecolor='orange', facecolor='none')
    ax.add_patch(rect)

    # quirky array legend
    if zoomer:
        rect = patches.Rectangle([2004.35, 2006.5], 2.5, 1.9, linewidth=2.0, edgecolor='lightgray', facecolor='none')
    else:
        rect = patches.Rectangle([1990.6, 2014.2], 4, 3, linewidth=1.5, edgecolor='lightgray', facecolor='none')

    ax.add_patch(rect)
    xy1 = [1997, 2025]
    xy2 = [1997, 2025]
    plt.plot(xy1, xy2, c='purple', alpha=1, marker='.',linewidth=2.0, markersize='12')
    xy1 = [2003, 2025]
    xy2 = [1997, 2018]
    plt.plot(xy1, xy2, c='purple', alpha=1, marker='.',linewidth=2.0, markersize='12')
    if zoomer:
        ax.xaxis.set_ticks(np.arange(2005, 2021, 5))
        ax.yaxis.set_ticks(np.arange(2005, 2011, 5))
    else:
        plt.figtext(0.16, 0.3, 'dt = 0 years',c='purple', fontsize=20)
        plt.figtext(0.33, 0.125, 'dt = 6 years',c='purple', fontsize=20)
    ax = plt.gca()
    ax.tick_params(axis='both', labelsize=18)
    plt.grid()
    plt.rc('grid', linestyle="-", color='black')
    plt.legend(["noisy", "different", "similar"],loc='upper left', fontsize=16)
    os.chdir('/Users/vidale/Documents/Research/IC/Plots_hold')
    if zoomer:
        plt.savefig('Iall_plot_zoom' + '.png')
    else:
        plt.savefig('Iall_plot' + '.png')

#%% YKA PKPprecursor similarity on 2nd date vs 1st date
if do_change12:
    fig_index = 53
    plt.figure(fig_index, figsize=(22*0.4, 22*0.4))
    # stupidity to get legend right
    plt.plot(1995, c='silver', marker='.', markersize='12')
    plt.plot(1995, c='red', marker='.',markersize='12')
    plt.plot(1995, c='limegreen', marker='.',markersize='12')
    plt.plot(1995, c='black', marker='.',markersize='12')

    minx = 1990
    maxx = 2024
    miny = 1990
    maxy = 2024
    plt.xlim(minx, maxx)
    plt.ylim(miny, maxy)

    for i in pair_range:
        pre_sim = PKPpre_sim[i]
        mult    = multiplet[i]
        marker_type = '.'
        if pre_sim == -1:
            plt.plot(date2[i], date1[i], c='silver',    alpha=1, marker=marker_type,linewidth=1.5, markersize='18')
        elif pre_sim == 0:
            plt.plot(date2[i], date1[i], c='red', alpha=1, marker=marker_type,linewidth=1.5, markersize='18')
        elif pre_sim == 1:
            plt.plot(date2[i], date1[i], c='limegreen',       alpha=1, marker=marker_type,linewidth=2.0, markersize='18')
        elif pre_sim == 2:
            plt.plot(date2[i], date1[i], c='black',     alpha=1, marker=marker_type,linewidth=2.0, markersize='18')

        if do_val:
            if pre_sim >= 0:
                plt.text(date2[i], date1[i], int(pair_index[i]), verticalalignment='top')
                plt.text(date2[i], date1[i], mult, verticalalignment='bottom')

    plt.xlabel('2nd event (date)', fontsize=20)
    plt.ylabel('1st event (date)', fontsize=20)
    ax = plt.gca()
    ax.tick_params(right=True, labelright=False,top=True, labeltop=False)
    plt.title('12 PKPprecursor similarity', fontsize=25)
    
    # diagonal line
    xy1 = [1990, 2024]
    xy2 = [1990, 2024]
    plt.plot(xy1, xy2, c='purple', alpha=1, marker='.',linewidth=2.0, markersize='12')
    xy1 = [1996, 2024]
    xy2 = [1990, 2018]
    plt.plot(xy1, xy2, c='purple', alpha=1, marker='.',linewidth=2.0, markersize='12')
    
    plt.figtext(0.20, 0.38, 'dt = 0 years',c='purple', fontsize=16)
    plt.figtext(0.25, 0.12, 'dt = 6 years',c='purple', fontsize=16)
    ax = plt.gca()
    ax.tick_params(axis='both', labelsize=16)
    plt.xticks(range(1990, 2024, 5))
    plt.yticks(range(1990, 2024, 5))
    plt.grid()
    plt.rc('grid', linestyle="-", color='black')
    plt.legend(["noisy", "same", "similar", "different"], fontsize=16)
    os.chdir('/Users/vidale/Documents/Research/IC/Plots_hold')
    plt.savefig('PKPpre_plot' + '.png')

# show all plots
plt.show()
