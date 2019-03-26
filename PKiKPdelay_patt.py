#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calculate rough estimate of expected pattern of time change from rotation
Created on Sun Mar 24 09:48:49 2019

@author: vidale
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math

velocity = 10    # km/s
IC_depth = 5000  # km
S_R_dist = 6000  # km
nx = 101
ny = nx
mesh_sp = S_R_dist / (nx-1)     # km
print(f'Mesh spacing {mesh_sp:.2f}')
#print(f'Mesh spacing {mesh_sp:.2f}  {it:4d}')

x = range(nx)       # source index
y = range(ny)       # radial index

PKiKP_data  = np.zeros(len(x) * len(y)) # sum of both legs for scatterer in specified location
PKiKP_diff  = np.zeros(len(x) * len(y)) # difference from first arrival scattering at midpoint
PKiKP_data2 = np.zeros(len(x) * len(y)) # difference from moving 20km tranversely
PKiKP_data2diff = np.zeros(len(x) * len(y)) # difference from moving 20km tranversely
PKiKP_data3 = np.zeros(len(x) * len(y)) # difference from moving 20km radially
PKiKP_data3diff = np.zeros(len(x) * len(y)) # difference from moving 20km radially

xx1         = np.zeros(len(x) * len(y)) # x distance in 1st leg, also x coordinate of scatterer
xx2         = np.zeros(len(x) * len(y)) # x distance in 2nd leg
xx1d        = np.zeros(len(x) * len(y)) # x distance in perturbed 1st leg, also x coordinate of perturbed scatterer
xx2d        = np.zeros(len(x) * len(y)) # x distance in perturbed 2nd leg
yy1         = np.zeros(len(x) * len(y)) # y distance in 1st leg, also y coordinate of scatterer
yy2         = np.zeros(len(x) * len(y)) # y distance in 2nd leg
yy1d        = np.zeros(len(x) * len(y)) # y distance in perturbed 1st leg, also y coordinate of perturbed scatterer
yy2d        = np.zeros(len(x) * len(y)) # y distance in perturbed 2nd leg
dr21         = np.zeros(len(x) * len(y)) # distance squared in receiver leg
ds21         = np.zeros(len(x) * len(y)) # distance squared in source leg
dr22         = np.zeros(len(x) * len(y)) # distance squared in receiver leg
ds22         = np.zeros(len(x) * len(y)) # distance squared in source leg
dr23         = np.zeros(len(x) * len(y)) # distance squared in receiver leg
ds23         = np.zeros(len(x) * len(y)) # distance squared in source leg

ICd2 = IC_depth * IC_depth

# original scatterer location
for ix in x:
	for iy in y:
		ii = iy*nx + ix
		xx1[ii] = ix*mesh_sp                # radial     source-scatterer   distance
		xx2[ii] = S_R_dist - xx1[ii]        # radial     receiver-scatterer distance
		yy1[ii] = (iy - ((nx-1)/2))*mesh_sp # y distance from center line
		yy2[ii] =  yy1[ii]                  # distance from center line
		ds21[ii] = xx1[ii]*xx1[ii] + yy1[ii]*yy1[ii] + ICd2
		dr21[ii] = xx2[ii]*xx2[ii] + yy2[ii]*yy2[ii] + ICd2
		PKiKP_data[ii] = (math.sqrt(ds21[ii]) + math.sqrt(dr21[ii]))/velocity

PKiKP_diff = PKiKP_data - min(PKiKP_data)
for ix in x:
	for iy in y:
		ii = iy*nx + ix
		if PKiKP_diff[ii] > 150:
			PKiKP_diff[ii] = 150

grid1 = PKiKP_diff.reshape((nx, ny))

color_map = plt.imshow(grid1, extent=(xx1[0], xx1[-1], yy1[0], yy1[-1]),
           interpolation='nearest', cmap=cm.gist_rainbow)
color_map.set_cmap = cm.gist_rainbow
plt.colorbar()

plt.figure(1)
plt.xlabel('radial distance (km)')
plt.ylabel('tranverse distance (km)')
plt.title('Relative arrival time (s)')

# move scatterer location transversely 10 km
plt.figure(2)
xx1d = xx1 + 10
xx2d = xx2 - 10

for ix in x:
	for iy in y:
		ii = iy*nx + ix
		ds22[ii] = xx1d[ii]*xx1d[ii] + yy1[ii]*yy1[ii] + ICd2
		dr22[ii] = xx2d[ii]*xx2d[ii] + yy2[ii]*yy2[ii] + ICd2
		PKiKP_data2[ii] = (math.sqrt(ds22[ii]) + math.sqrt(dr22[ii]))/velocity
		PKiKP_data2diff[ii] = PKiKP_data2[ii] - PKiKP_data[ii]

grid2 = PKiKP_data2diff.reshape((nx, ny))

color_map = plt.imshow(grid2, extent=(xx1[0], xx1[-1], yy1[0], yy1[-1]),
           interpolation='nearest', cmap=cm.gist_rainbow)
color_map.set_cmap = cm.gist_rainbow
plt.colorbar()

plt.xlabel('radial distance (km)')
plt.ylabel('tranverse distance (km)')
plt.title('Spatial delay pattern for radial shift (s)')

# move scatterer location transversely 10 km
plt.figure(3)
yy1d = yy1 + 10
yy2d = yy2 + 10

for ix in x:
	for iy in y:
		ii = iy*nx + ix
		ds23[ii] = xx1[ii]*xx1[ii] + yy1d[ii]*yy1d[ii] + ICd2
		dr23[ii] = xx2[ii]*xx2[ii] + yy2d[ii]*yy2d[ii] + ICd2
		PKiKP_data3[ii] = (math.sqrt(ds23[ii]) + math.sqrt(dr23[ii]))/velocity
		PKiKP_data3diff[ii] = PKiKP_data3[ii] - PKiKP_data[ii]

grid3 = PKiKP_data3diff.reshape((nx, ny))

color_map = plt.imshow(grid3, extent=(xx1[0], xx1[-1], yy1[0], yy1[-1]),
           interpolation='nearest', cmap=cm.gist_rainbow)
color_map.set_cmap = cm.gist_rainbow
plt.colorbar()

plt.xlabel('radial distance (km)')
plt.ylabel('tranverse distance (km)')
plt.title('Spatial delay pattern for tranverse shift (s)')

# plot delay vs offset

plt.figure(4)

xxx = np.zeros(nx)
for ix in x:
	iy = int((ny - 1) / 2)
	ii = int(iy*nx + ix)
	xxx[ix] = PKiKP_data2[ii]                # radial     source-scatterer   distance
plt.plot(x, xxx)
plt.xlabel('radial distance (km)')
plt.ylabel('time delay (s)')
plt.title('Scattering time vs radial distance (s)')

plt.figure(5)

yyy = np.zeros(ny)
for iy in y:
	ix = int((nx - 1) / 2)
	ii = int(iy*nx + ix)
	yyy[iy] = PKiKP_data3[ii]                # radial     source-scatterer   distance
plt.plot(y, yyy)
plt.xlabel('tranverse distance (km)')
plt.ylabel('time delay (s)')
plt.title('Scattering time vs tranverse distance (s)')

# plot delay vs offset

plt.figure(6)

xxx = np.zeros(nx)
for ix in x:
	iy = int((ny - 1) / 2)
	ii = int(iy*nx + ix)
	xxx[ix] = PKiKP_data2diff[ii]                # radial     source-scatterer   distance
plt.plot(x, xxx)
plt.xlabel('radial distance (km)')
plt.ylabel('time delay (s)')
plt.title('Variation in scattering time with radial distance (s)')

plt.figure(7)

yyy = np.zeros(ny)
for iy in y:
	ix = int((nx - 1) / 2)
	ii = int(iy*nx + ix)
	yyy[iy] = PKiKP_data3diff[ii]                # radial     source-scatterer   distance
plt.plot(y, yyy)
plt.xlabel('tranverse distance (km)')
plt.ylabel('time delay (s)')
plt.title('Variation in scattering time with tranverse distance (s)')

plt.show()

test = 50*101 + 50
print(f'xx1 {xx1[test]:.2f} yy1 {yy1[test]:.2f}')
print(f'xx2 {xx2[test]:.2f} yy2 {yy2[test]:.2f}')

print(f'xx1d {xx1d[test]:.2f} yy1d {yy1d[test]:.2f}')
print(f'xx2d {xx2d[test]:.2f} yy2d {yy2d[test]:.2f}')

print(f'ds21 {ds21[test]:.2f} dr21 {dr21[test]:.2f}')
print(f'ds22 {ds22[test]:.2f} dr22 {dr22[test]:.2f}')
print(f'ds23 {ds23[test]:.2f} dr23 {dr23[test]:.2f}')

print(f'PKiKP_data2diff {PKiKP_data2diff[test]:.2f} PKiKP_data2 {PKiKP_data2[test]:.2f}')
print(f'PKiKP_data3diff {PKiKP_data3diff[test]:.2f} PKiKP_data3 {PKiKP_data3[test]:.2f}')
