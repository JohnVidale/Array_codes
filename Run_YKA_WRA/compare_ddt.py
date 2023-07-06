#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 29 10:44:55 2020
Makes a map of the LASA events
@author: vidale
"""
import numpy as np
import matplotlib.pyplot as plt

plot_John   = False
plot_Ruoyan = False
plot_Wei    = False
plot_final  = False
plot_Y_vs_I = True

#%% plot data vs prediction
fig_index = 1
plt.figure(1, figsize=(10,5))
# fig, ax = plt.subplots(1)

#ax = fig.gca()
# ax.set_xticks(np.arange(0, 360, 30))
# ax.set_yticks(np.arange(0, 100, 30))

#    FILL NUMBERS
# print('lat0 is ' + str(comp_lat[0]) + ' lon0 is  ' + str(comp_lon[0]))

plt.xlim(-0.05,0.2)
plt.ylim(-0.05,0.075)

# black maybe reliable
# red impossible to pick
# blue(?) unreliable estimate
if plot_John   == True:
    plt.scatter( 0.153, 0.010, c='b', s=100, alpha=1, marker='.') # P16
    plt.scatter( 0.076, 0.040, c='b', s=100, alpha=1, marker='.') # P17
    plt.scatter( 0.072, 0.020, c='b', s=100, alpha=1, marker='.') # P18
    plt.scatter( 0.018, 0.020, c='b', s=100, alpha=1, marker='.') # P30
    plt.scatter( 0.012, 0.020, c='k', s=100, alpha=1, marker='.') # P31
    plt.scatter(-0.012,-0.010, c='k', s=100, alpha=1, marker='.') # P43
    plt.scatter(-0.002, 0.000, c='k', s=100, alpha=1, marker='.') # P34
    plt.scatter(-0.039, 0.000, c='r', s=100, alpha=1, marker='.') # P75
    plt.scatter( 0.040, 0.060, c='b', s=100, alpha=1, marker='.') # P78

    plt.scatter( 0.079, 0.000, c='r', s=100, alpha=1, marker='.') # P86
    plt.scatter( 0.000, 0.000, c='r', s=100, alpha=1, marker='.') # 100

    plt.scatter( 0.054, 0.000, c='r', s=100, alpha=1, marker='.') # P56
    plt.scatter( 0.058, 0.040, c='k', s=100, alpha=1, marker='.') # P59
    plt.scatter( 0.000,-0.010, c='b', s=100, alpha=1, marker='.') # 106
    plt.scatter( 0.051, 0.040, c='b', s=100, alpha=1, marker='.') # P97
    plt.scatter(-0.017, 0.000, c='r', s=100, alpha=1, marker='.') # P28
    plt.scatter( 0.000,-0.010, c='k', s=100, alpha=1, marker='.') # P79
    plt.scatter( 0.021, 0.015, c='k', s=100, alpha=1, marker='.') # P80
    plt.scatter(-0.020,-0.010, c='k', s=100, alpha=1, marker='.') # P40
    plt.scatter( 0.011, 0.020, c='k', s=100, alpha=1, marker='.') # P51
    plt.scatter(-0.013,-0.015, c='k', s=100, alpha=1, marker='.') # P52
    plt.scatter( 0.006, 0.010, c='k', s=100, alpha=1, marker='.') # P53
    plt.scatter( 0.030, 0.020, c='k', s=100, alpha=1, marker='.') # 104
    plt.scatter( 0.009, 0.000, c='k', s=100, alpha=1, marker='.') # 105

if plot_Ruoyan == True:
    plt.scatter( 0.153, 0.030, c='m', s=100, alpha=1, marker='.') # P16
    plt.scatter( 0.076, 0.040, c='m', s=100, alpha=1, marker='.') # P17
    plt.scatter( 0.072, 0.020, c='m', s=100, alpha=1, marker='.') # P18
    plt.scatter( 0.018, 0.040, c='m', s=100, alpha=1, marker='.') # P30
    plt.scatter( 0.012, 0.040, c='m', s=100, alpha=1, marker='.') # P31
    plt.scatter(-0.012,-0.010, c='m', s=100, alpha=1, marker='.') # P43
    plt.scatter(-0.002,-0.010, c='m', s=100, alpha=1, marker='.') # P34
    plt.scatter(-0.039,-0.030, c='m', s=100, alpha=1, marker='.') # P75
    plt.scatter( 0.040, 0.030, c='m', s=100, alpha=1, marker='.') # P78

    plt.scatter( 0.079, 0.060, c='m', s=100, alpha=1, marker='.') # P86
    plt.scatter( 0.000, 0.020, c='m', s=100, alpha=1, marker='.') # 100

    plt.scatter( 0.054, 0.000, c='r', s=100, alpha=1, marker='.') # P56
    plt.scatter( 0.058, 0.060, c='m', s=100, alpha=1, marker='.') # P59
    plt.scatter( 0.000, 0.000, c='r', s=100, alpha=1, marker='.') # 106
    plt.scatter( 0.051, 0.050, c='m', s=100, alpha=1, marker='.') # P97
    plt.scatter(-0.017, 0.000, c='r', s=100, alpha=1, marker='.') # P28
    plt.scatter( 0.000, 0.000, c='r', s=100, alpha=1, marker='.') # P79
    plt.scatter( 0.021, 0.000, c='m', s=100, alpha=1, marker='.') # P80
    plt.scatter(-0.020, 0.020, c='m', s=100, alpha=1, marker='.') # P40
    plt.scatter( 0.011,-0.020, c='m', s=100, alpha=1, marker='.') # P51
    plt.scatter(-0.013, 0.000, c='m', s=100, alpha=1, marker='.') # P52
    plt.scatter( 0.006, 0.000, c='m', s=100, alpha=1, marker='.') # P53
    plt.scatter( 0.030, 0.000, c='r', s=100, alpha=1, marker='.') # 104
    plt.scatter( 0.009, 0.000, c='m', s=100, alpha=1, marker='.') # 105

if plot_Wei    == True:
    plt.scatter( 0.153, 0.010, c='g', s=100, alpha=1, marker='.') # P16
    plt.scatter( 0.076, 0.040, c='g', s=100, alpha=1, marker='.') # P17
    plt.scatter( 0.072, 0.020, c='g', s=100, alpha=1, marker='.') # P18
    plt.scatter( 0.018, 0.020, c='g', s=100, alpha=1, marker='.') # P30
    plt.scatter( 0.012, 0.040, c='g', s=100, alpha=1, marker='.') # P31
    plt.scatter(-0.012, 0.020, c='g', s=100, alpha=1, marker='.') # P43
    plt.scatter(-0.002, 0.000, c='g', s=100, alpha=1, marker='.') # P34
    plt.scatter(-0.039, 0.000, c='r', s=100, alpha=1, marker='.') # P75
    plt.scatter( 0.040, 0.060, c='g', s=100, alpha=1, marker='.') # P78

    plt.scatter( 0.079, 0.000, c='r', s=100, alpha=1, marker='.') # P86
    plt.scatter( 0.000, 0.000, c='r', s=100, alpha=1, marker='.') # 100

    plt.scatter( 0.054, 0.000, c='r', s=100, alpha=1, marker='.') # P56
    plt.scatter( 0.058, 0.040, c='g', s=100, alpha=1, marker='.') # P59
    plt.scatter( 0.000,-0.010, c='g', s=100, alpha=1, marker='.') # 106
    plt.scatter( 0.051, 0.040, c='g', s=100, alpha=1, marker='.') # P97
    plt.scatter(-0.017, 0.000, c='r', s=100, alpha=1, marker='.') # P28
    plt.scatter( 0.000,-0.010, c='g', s=100, alpha=1, marker='.') # P79
    plt.scatter( 0.021, 0.010, c='g', s=100, alpha=1, marker='.') # P80
    plt.scatter(-0.020,-0.010, c='g', s=100, alpha=1, marker='.') # P40
    plt.scatter( 0.011,-0.020, c='g', s=100, alpha=1, marker='.') # P51
    plt.scatter(-0.013,-0.010, c='g', s=100, alpha=1, marker='.') # P52
    plt.scatter( 0.006, 0.010, c='g', s=100, alpha=1, marker='.') # P53
    plt.scatter( 0.030, 0.020, c='g', s=100, alpha=1, marker='.') # 104
    plt.scatter( 0.009, 0.000, c='g', s=100, alpha=1, marker='.') # 105

if plot_final    == True:
    # most similar
    plt.scatter(-0.002, 0.00, c='b', s=100, alpha=1, marker='.') # P34
    plt.scatter( 0.018, 0.00, c='r', s=100, alpha=1, marker='.') # P30
    plt.scatter(-0.012,-0.01, c='b', s=100, alpha=1, marker='.') # P43
    plt.scatter( 0.076, 0.04, c='c', s=100, alpha=1, marker='.') # P17
    plt.scatter( 0.072, 0.00, c='r', s=100, alpha=1, marker='.') # P18
    plt.scatter( 0.040, 0.04, c='b', s=100, alpha=1, marker='.') # P78
    plt.scatter( 0.079, 0.00, c='r', s=100, alpha=1, marker='.') # P86
    plt.scatter( 0.021, 0.01, c='b', s=100, alpha=1, marker='.') # P80
    plt.scatter( 0.000, 0.00, c='b', s=100, alpha=1, marker='.') # P79
    plt.scatter( 0.000, 0.00, c='b', s=100, alpha=1, marker='.') # 106
    plt.scatter( 0.058, 0.04, c='b', s=100, alpha=1, marker='.') # P59
    plt.scatter( 0.030, 0.00, c='r', s=100, alpha=1, marker='.') # 104
    plt.scatter( 0.006, 0.00, c='b', s=100, alpha=1, marker='.') # P53
    plt.scatter( 0.009, 0.02, c='b', s=100, alpha=1, marker='.') # 105
    plt.scatter(-0.013,-0.01, c='b', s=100, alpha=1, marker='.') # P52
    plt.scatter( 0.011, 0.02, c='c', s=100, alpha=1, marker='.') # P51

    # less similar
    plt.scatter( 0.153, 0.00, c='m', s=100, alpha=1, marker='.') # P16
    plt.scatter( 0.012, 0.00, c='m', s=100, alpha=1, marker='.') # P31
    plt.scatter(-0.039, 0.00, c='m', s=100, alpha=1, marker='.') # P75
    plt.scatter( 0.000, 0.00, c='m', s=100, alpha=1, marker='.') # 100
    plt.scatter( 0.054, 0.00, c='m', s=100, alpha=1, marker='.') # P56
    plt.scatter( 0.051, 0.00, c='m', s=100, alpha=1, marker='.') # P97
    plt.scatter(-0.017, 0.00, c='m', s=100, alpha=1, marker='.') # P28
    plt.scatter(-0.020, 0.00, c='m', s=100, alpha=1, marker='.') # P40

if plot_Y_vs_I    == True:
    # most similar
    plt.scatter( 0.00,-0.02, c='b', s=100, alpha=1, marker='.') # P34
    plt.scatter( 0.00, 0.00, c='r', s=100, alpha=1, marker='.') # P30
#    plt.scatter(-0.01,-0.01, c='b', s=100, alpha=1, marker='.') # P43
    plt.scatter( 0.04, 0.00, c='c', s=100, alpha=1, marker='.') # P17
    plt.scatter( 0.00,-0.01, c='r', s=100, alpha=1, marker='.') # P18
    plt.scatter( 0.04,-0.02, c='b', s=100, alpha=1, marker='.') # P78
    plt.scatter( 0.00,-0.02, c='r', s=100, alpha=1, marker='.') # P86
    plt.scatter( 0.05, 0.04, c='b', s=100, alpha=1, marker='.') # P25
    plt.scatter( 0.01, 0.03, c='b', s=100, alpha=1, marker='.') # P80
#    plt.scatter( 0.00, 0.00, c='b', s=100, alpha=1, marker='.') # P79
    plt.scatter( 0.00, 0.02, c='b', s=100, alpha=1, marker='.') # 106
    plt.scatter( 0.04, 0.00, c='b', s=100, alpha=1, marker='.') # P59
    plt.scatter( 0.00,-0.02, c='r', s=100, alpha=1, marker='.') # 104
    plt.scatter( 0.00, 0.02, c='b', s=100, alpha=1, marker='.') # P53
    plt.scatter( 0.02,-0.03, c='b', s=100, alpha=1, marker='.') # 105
    plt.scatter(-0.01, 0.02, c='b', s=100, alpha=1, marker='.') # P52
    plt.scatter( 0.02,-0.01, c='c', s=100, alpha=1, marker='.') # P51

xpoints = np.array([-0.05, 0.2])
ypoints = np.array([-0.05, 0.2])
plt.plot(xpoints, ypoints, linestyle = 'dotted')

    # plt.text(comp_lon[ii] + 0.01 * (max_lon - min_lon), comp_lat[ii] + 0.01 * (max_lat - min_lat), comp_name[ii])
plt.grid()
plt.rc('grid', linestyle="-", color='black')
plt.xlabel('Yang & Song estimate (s)')
plt.ylabel('Wang & Vidale estimate (s)')
plt.title('ddt comparison')
# plt.legend([])
# plt.colorbar()
plt.show()
