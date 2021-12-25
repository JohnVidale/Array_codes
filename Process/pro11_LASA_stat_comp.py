#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 17:21:27 2019

@author: vidale
"""
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt


# Engdahl-Felix (1971) numbers from plot
shiftF = [ 0.0,-0.3, -1.0, -0.6, -0.6] # 0 doesn't exist
shiftE = [ 0.0, 0.0, -0.5, -0.6, -0.2]
shiftD = [ 0.0, 0.3, -0.4, -0.2, -0.1]
shiftC = [ 0.0, 0.1, -0.3, -0.3,  0.2]
shiftB = [ 0.0, 0.0, -0.2,  0.0,  0.2]
shiftA = [ 0.0, 0.0,  0.0,  0.0,  0.0] # only A == 1 element exists

# test to see effect on beam, not yet implemented
TshiftF = [ 0.0, 0.0,  0.0, 0.0, 0.0] # 0 doesn't exist in station numbering
TshiftE = [ 0.0, 0.0,  0.0, 0.0, 0.0]
TshiftD = [ 0.0, 0.0,  0.0, 0.0, 0.0]
TshiftC = [ 0.0, 0.0,  0.0, 0.0, 0.0]
TshiftB = [ 0.0, 0.0,  0.0, 0.0, 0.0]
TshiftA = [ 0.0, 0.0,  0.0, 0.0, 0.0] # only A == 1 element exists

sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/L_sta_statics.txt'
with open(sta_file, 'r') as file:
	lines = file.readlines()
print(str(len(lines)) + ' stations read from ' + sta_file)
# Load station coords into arrays
station_index = range(len(lines))
st_names = []
st_dist  = []
st_lats  = []
st_lons  = []
st_shift = []
st_corr  = []
st_Eshift = []
for ii in station_index:
	line = lines[ii]
	split_line = line.split()
	st_names.append(split_line[0])
	st_dist.append( split_line[1])
	st_lats.append( split_line[2])
	st_lons.append( split_line[3])
	st_shift.append(split_line[4])
	st_corr.append( split_line[5])

	array_L = split_line[0][0]
	indexE = split_line[0][1]
	if array_L == 'A':
		Eshift = shiftA[int(indexE)]
	if array_L == 'B':
		Eshift = shiftB[int(indexE)]
	if array_L == 'C':
		Eshift = shiftC[int(indexE)]
	if array_L == 'D':
		Eshift = shiftD[int(indexE)]
	if array_L == 'E':
		Eshift = shiftE[int(indexE)]
	if array_L == 'F':
		Eshift = shiftF[int(indexE)]
#	print(str(array_L) + ' ' + str(split_line[0][1]) + ' ' + str(indexE))
#	Ecorr = array_L[int(indexE)]
#	print(str(array_L) + ' ' + str(split_line[0][1]) + ' ' + str(indexE) + ' ' + str(Ecorr))
	st_Eshift.append(Eshift)

for ii in station_index:
	print(st_names[ii] + ' ' +str(st_Eshift[ii]) + ' ' + str(st_shift[ii]))

st_Eshift = np.array(st_Eshift,dtype=np.float)
st_shift  = np.array(st_shift, dtype=np.float)
#%% plot data vs prediction
fig, ax = plt.subplots(1)

ax = fig.gca()

plt.scatter(st_shift,st_Eshift, c='b', s=50, alpha=1, marker='.')
ax.set_xticks(np.arange(-1, 1, 0.5))
ax.set_yticks(np.arange(-1, 1, 0.5))

print(type(st_shift))
print(type(st_Eshift))
slope, intercept, r_value, p_value, std_err = stats.linregress(st_shift,st_Eshift)

line = slope*st_shift+intercept
plt.plot(st_shift, line, 'r', label='y={:.2f}x+{:.2f}'.format(slope,intercept))

plt.plot(np.unique(st_shift), np.poly1d(np.polyfit(st_shift, st_Eshift, 1))(np.unique(st_shift)))
plt.plot(st_shift, np.poly1d(np.polyfit(st_shift, st_Eshift, 1))(st_shift))

plt.grid()

#plt.rc('grid', linestyle="-", color='black')
#plt.yticks(np.arange(-1,1,step=0.5))
plt.ylabel('Engdahl-Felix correction (s)')
plt.xlabel('Correlation computed correction (s)')
plt.title('CC vs Engdahl time shifts')
#plt.gca().invert_yaxis()
plt.show()