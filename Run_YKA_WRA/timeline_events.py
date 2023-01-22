#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Dec 27th, 2022
heck global/ILcMakes a timeline PKIKP waveform changes for ILAR & YKA
@author: vidale
"""
import numpy as np
import matplotlib.pyplot as plt
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
change_Y     = np.zeros(pair_count)
change_I     = np.zeros(pair_count)
shift_Y     = np.zeros(pair_count)
shift_I     = np.zeros(pair_count)
event1       = np.zeros(pair_count)
event2       = np.zeros(pair_count)
date1        = np.zeros(pair_count)
date2        = np.zeros(pair_count)

do_YKA_change  = True
do_ILAR_change = True
do_YKA_shift   = False
do_ILAR_shift  = False

for i in pair_index:
    line = lines[i]
    split_line = line.split()
    # print('input station ' + split_line[0])
    pair_name.append(   split_line[0])
    change_Y[i] = int(split_line[1])
    change_I[i] = int(split_line[2])
    shift_Y[i] = float(split_line[3])
    shift_I[i] = float(split_line[4])
    event1[i]   = int(split_line[5])
    event2[i]   = int(split_line[6])
    t1          = UTCDateTime(split_line[7])
    t2          = UTCDateTime(split_line[8])
    date1[i] = t1.year + t1.month/12.
    date2[i] = t2.year + t2.month/12.
    # print('Pair is ' + pair_name[i] + ' changes are ' + str(change_Y[i]) + ' ' + str(change_I[i]) + ' events are ' + str(event1[i]) + ' ' + str(event2[i]) + ' t1 is ' + str(date1[i]) + ' t2 is ' + str(date2[i]))

#%% plot data vs prediction

# #    FILL NUMBERS
min_index = 0
max_index = pair_count +1
min_year = 1990
max_year = 2021

if do_YKA_change:
    fig_index = 1
    plt.figure(1, figsize=(10,10))
    plt.xlim( min_year,  max_year)
    plt.ylim( min_index, max_index)

    for i in pair_index:
        ii = pair_count - i
        x = [date1[i], date2[i]]
        line_index = [ii,ii]
        if change_Y[i] == 2:
            plt.plot(x, line_index, c='r', alpha=1, marker='.',linewidth=3.0)
        elif change_Y[i] == 1:
            plt.plot(x, line_index, c='y', alpha=1, marker='.',linewidth=3.0)
        elif change_Y[i] == 0:
            plt.plot(x, line_index, c='g', alpha=1, marker='.',linewidth=3.0)
        elif change_Y[i] == -1:
            plt.plot(x, line_index, c='0.8', alpha=1, marker='.',linewidth=2.0)
        elif change_Y[i] == -2:
            plt.plot(x, line_index, c='0.8', alpha=1, marker='.',linewidth=2.0)
        elif change_Y[i] == -4:
            plt.plot(x, line_index, c='c', alpha=1, marker='.',linewidth=3.0)
        # plt.text(comp_lon[ii] + 0.01 * (max_lon - min_lon), comp_lat[ii] + 0.01 * (max_lat - min_lat), comp_name[ii])
        plt.text(date1[i] - 1, ii, pair_name[i])
    plt.xlabel('Year', fontsize=20)
    plt.ylabel('South to North', fontsize=20)
    plt.title('YKA waveform variation over the years', fontsize=25)
    plt.legend([])

if do_ILAR_change:
    fig_index = 2
    plt.figure(2, figsize=(10,10))
    plt.xlim( min_year,  max_year)
    plt.ylim( min_index, max_index)

    for i in pair_index:
        ii = pair_count - i
        x = [date1[i], date2[i]]
        line_index = [ii,ii]
        if change_I[i] == 2:
            plt.plot(x, line_index, c='r', alpha=1, marker='.',linewidth=3.0)
        elif change_I[i] == 1:
            plt.plot(x, line_index, c='y', alpha=1, marker='.',linewidth=3.0)
        elif change_I[i] == 0:
            plt.plot(x, line_index, c='g', alpha=1, marker='.',linewidth=3.0)
        elif change_I[i] == -1:
            plt.plot(x, line_index, c='0.8', alpha=1, marker='.',linewidth=2.0)
        elif change_I[i] == -2:
            plt.plot(x, line_index, c='0.8', alpha=1, marker='.',linewidth=2.0)
        elif change_I[i] == -4:
            plt.plot(x, line_index, c='c', alpha=1, marker='.',linewidth=3.0)
        plt.text(date1[i] - 1, ii, pair_name[i])
        # plt.text(comp_lon[ii] + 0.01 * (max_lon - min_lon), comp_lat[ii] + 0.01 * (max_lat - min_lat), comp_name[ii])

    # plt.grid()
    # plt.rc('grid', linestyle="-", color='black')
    plt.xlabel('Year', fontsize=20)
    plt.ylabel('South to North', fontsize=20)
    plt.title('ILAR waveform variation over the years', fontsize=25)
    plt.legend([])

if do_ILAR_shift:
    fig_index = 4
    plt.figure(4, figsize=(10,10))
    plt.xlim( min_year,  max_year)
    plt.ylim( min_index, max_index)

    for i in pair_index:
        ii = pair_count - i
        x = [date1[i], date2[i]]
        line_index = [ii,ii]
        if shift_I[i] < 0 and shift_I[i] > -0.5:
            plt.plot(x, line_index, c='r', alpha=1, marker='.',linewidth=3.0)
        elif shift_I[i] == 0:
            plt.plot(x, line_index, c='y', alpha=1, marker='.',linewidth=3.0)
        elif shift_I[i] > 0:
            plt.plot(x, line_index, c='g', alpha=1, marker='.',linewidth=3.0)
        elif shift_I[i] == -2:
            plt.plot(x, line_index, c='0.6', alpha=1, marker='.',linewidth=2.0)
        elif shift_I[i] == -3:
            plt.plot(x, line_index, c='0.9', alpha=1, marker='.',linewidth=2.0)
        elif shift_I[i] == -4:
            plt.plot(x, line_index, c='c', alpha=1, marker='.',linewidth=3.0)
        plt.text(date1[i] - 1, ii, pair_name[i])
        # plt.text(comp_lon[ii] + 0.01 * (max_lon - min_lon), comp_lat[ii] + 0.01 * (max_lat - min_lat), comp_name[ii])

    # plt.grid()
    # plt.rc('grid', linestyle="-", color='black')
    plt.xlabel('Year', fontsize=20)
    plt.ylabel('South to North', fontsize=20)
    plt.title('ILAR ddt over the years', fontsize=25)
    plt.legend([])
# # plt.colorbar()
plt.show()
