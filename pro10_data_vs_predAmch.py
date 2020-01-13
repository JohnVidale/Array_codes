#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 17:21:27 2019

@author: vidale
"""
import numpy as np
import matplotlib.pyplot as plt


#%% plot data vs prediction
fig, ax = plt.subplots(1)

#ax = fig.gca()
ax.set_xticks(np.arange(-1, 1, 0.25))
ax.set_yticks(np.arange(-0.3, 0.2, 0.1))

obs = [0.17, 0.10, -0.20, -0.30, -0.20, -0.18, -0.17]
#%% but the signal is averaging 1.2 Hz, hence time delat, which assumed 1 Hz, needs reduction
obs = (1/1) * np.array(obs)
pred = [0.5, 0.7, -0.6, -0.6, -0.4, -0.5, 0.0]
#pred = -1 * pred
plt.scatter(pred,obs, c='k', s=100, alpha=1, marker='.')
plt.plot(np.unique(pred), np.poly1d(np.polyfit(pred, obs, 1))(np.unique(pred)))
plt.grid()
plt.rc('grid', linestyle="-", color='black')
plt.xlabel('Predicted time shift for 1° rotation (s)')
plt.ylabel('Observed time shift over 2 years (s)')
plt.title('Amchitka - Observed vs predicted time shifts')
plt.show()