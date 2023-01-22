#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 29 10:44:55 2020
Makes a map of the LASA events
@author: vidale
"""
import numpy as np
import matplotlib.pyplot as plt

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
# unreliable estimate
plt.scatter( 0.153, 0.010, c='b', s=100, alpha=1, marker='.')
plt.scatter( 0.076, 0.040, c='b', s=100, alpha=1, marker='.')
plt.scatter( 0.072, 0.020, c='b', s=100, alpha=1, marker='.')
plt.scatter( 0.018, 0.020, c='b', s=100, alpha=1, marker='.')
plt.scatter(-0.012, 0.020, c='k', s=100, alpha=1, marker='.')
plt.scatter(-0.002,-0.010, c='k', s=100, alpha=1, marker='.')
plt.scatter(-0.039, 0.000, c='k', s=100, alpha=1, marker='.')
plt.scatter( 0.040, 0.000, c='r', s=100, alpha=1, marker='.')

plt.scatter( 0.079, 0.000, c='r', s=100, alpha=1, marker='.')
plt.scatter( 0.000, 0.000, c='r', s=100, alpha=1, marker='.')

plt.scatter( 0.054, 0.000, c='r', s=100, alpha=1, marker='.')
plt.scatter( 0.058, 0.040, c='k', s=100, alpha=1, marker='.')
plt.scatter( 0.000,-0.010, c='b', s=100, alpha=1, marker='.')
plt.scatter( 0.051, 0.040, c='b', s=100, alpha=1, marker='.')
plt.scatter(-0.017, 0.000, c='r', s=100, alpha=1, marker='.')
plt.scatter( 0.000,-0.010, c='k', s=100, alpha=1, marker='.')
plt.scatter( 0.021, 0.015, c='k', s=100, alpha=1, marker='.')
plt.scatter(-0.020,-0.010, c='k', s=100, alpha=1, marker='.')
plt.scatter( 0.011, 0.020, c='k', s=100, alpha=1, marker='.')
plt.scatter(-0.013,-0.015, c='k', s=100, alpha=1, marker='.')
plt.scatter( 0.006, 0.010, c='k', s=100, alpha=1, marker='.')
plt.scatter( 0.030, 0.020, c='k', s=100, alpha=1, marker='.')
plt.scatter( 0.009, 0.000, c='k', s=100, alpha=1, marker='.')

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
