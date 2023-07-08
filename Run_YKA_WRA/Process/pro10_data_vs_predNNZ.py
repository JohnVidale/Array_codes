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
ax.set_yticks(np.arange(-0.2, 0.2, 0.05))

# new
obs = [0.05,  0.10,  0.03, -0.04, -0.20,  0.10, -0.05, -0.18, -0.10,  0.20, -0.05,  -0.10,  0.05,  0.00]
pred = [0.3, -0.6, -0.1,  0.2,  0.5, -0.5,  0.1,  0.5,  0.5, -0.5,  0.4,  0.7, -0.5, 0.3]
# old
#obs = [-0.06, 0.01, 0.12, -0.04, -0.12, -0.20,  0.10, -0.05, -0.20, -0.15,  0.10, -0.12, -0.05]
#pred = [-0.1,-0.0, 0.6, -0.3, -0.3, -0.7, 0.5, -0.3, -0.8, -0.7, 0.4, -0.35, -0.2]
#pred = -1 * pred
#%% but the signal is averaging 1.2 Hz, hence time delat, which assumed 1 Hz, needs reduction
obs = (1/1.2) * np.array(obs)
plt.scatter(pred,obs, c='k', s=100, alpha=1, marker='.')
plt.plot(np.unique(pred), np.poly1d(np.polyfit(pred, obs, 1))(np.unique(pred)))
plt.grid()
plt.rc('grid', linestyle="-", color='black')
plt.xlabel('Predicted time shift for 1° rotation (s)')
plt.ylabel('Observed time shift over 3 years (s)')
plt.title('Observed vs predicted time shifts')
plt.show()