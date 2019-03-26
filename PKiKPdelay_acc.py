#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calculate accurate estimate of expected pattern of time change with lag time and back-azimuth from rotation
Started on Sun Mar 24 09:48:49 2019
1.  grid ICB scattering points in lat & lon
2.  for each point, compute angular distance to both source and receiver
3.  compute travel time to twice that distance for PKiKP (divide by two for segment time)
4.  also compute back-azimuth and slowness receiver to scatterer
5.  do same for scattering points 1° rotated in lat & lon
6.  transform scattering point grid (lon, lat, time) into grid of slowness anomaly relative PKiKP (Sr, St, time)

@author: vidale
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math
from obspy.geodetics import gps2dist_azimuth
from obspy.taup import TauPyModel
model = TauPyModel(model='iasp91')

LASA_lat =   46.7    # °N
LASA_lon = -106.22   # °E
# LASA_rad =    0.4    # ° radius to exclude outer rings (°)
NZ_lat   = 73.393
NZ_lon   = 54.929
NZ_depth = 0

min_lat = 10
dlat = 5
dlon = 5

# test obspy routines
distance = gps2dist_azimuth(LASA_lat,LASA_lon,NZ_lat,NZ_lon) # Get traveltimes again, hard to store
dist= distance[0] # distance in m
az1 = distance[1] # azimuth A->B in degrees
az2 = distance[2] # azimuth B->A in degrees)
dist = distance[0]/(1000*111)
arrival = model.get_travel_times(source_depth_in_km=NZ_depth,
								 distance_in_degree=dist,phase_list=['PKiKP'])
arrival_time = arrival[0].time
arrival_rayp = arrival[0].ray_param
arrival_ic   = arrival[0].incident_angle
print(f'Distance {dist:.0f} azimuth {az1:.0f} back-azimuth {az2:.0f}')
print(f'Arrival time {arrival_time:.0f} rayp {arrival_rayp:.0f} incident angle {arrival_ic:.0f}')

# parameterize
tsteps = int((90 - min_lat)/dlat) + 1
r  = np.linspace(     90-min_lat, 0, int(tsteps))
rr = np.linspace(min_lat,      90, int(tsteps))
theta = np.linspace(0, 2.*np.pi, int((360/dlon) + 1))

nx = len(theta) # longitude, every 10 degrees
ny = len(rr) # latitude, every dlat degrees

x = range(nx)       # longitude index
y = range(ny)       # latitude index
nn = nx * ny

xlat        = np.zeros(nn) # lat at each grid point
xlon        = np.zeros(nn) # lon at each grid point

PKiKP_data  = np.zeros(nn) # time sum of both legs for scatterer in specified location
PKiKP_diff  = np.zeros(nn) # difference from first arrival scattering at midpoint
PKiKP_rot   = np.zeros(nn) # time from rotating 1°

dd1         = np.zeros(nn) # distance from scatterer to NZ
dd2         = np.zeros(nn) # distance from scatterer to LASA
baz         = np.zeros(nn) # back-azimuth from scatterer to LASA
rayp        = np.zeros(nn) # ray parameter at LASA
tt1         = np.zeros(nn) # time in 1st leg
tt2         = np.zeros(nn) # time in 2nd leg
tt1d        = np.zeros(nn) # time in 5° rotated 1st leg
tt2d        = np.zeros(nn) # time in 5° rotated 2nd leg

# time and angle for grid of scatterer locations
for ix in x:
	for iy in y:
#		ii = iy*nx + ix
		ii = iy + ix*ny
		xlat[ii] = rr[iy] # specify scatterer
		xlon[ii] = theta[ix]*180/math.pi

		distance1 = gps2dist_azimuth(xlat[ii],xlon[ii],NZ_lat,NZ_lon) # source to scatterer
		distance2 = gps2dist_azimuth(LASA_lat,LASA_lon,xlat[ii],xlon[ii]) # scatterer to receiver
		dd1[ii] = distance1[0]/(1000*111)
		dd2[ii] = distance2[0]/(1000*111)
		if (dd1[ii] < 77 and dd2[ii] < 77):
			arrival = model.get_travel_times(source_depth_in_km=NZ_depth,
									 distance_in_degree = dd1[ii]*2, phase_list=['PKiKP'])
			tt1[ii] = 0.5 * arrival[0].time # twice the  distance, half the traveltime for one half PKiKP

			baz[ii] = distance2[1] # azimuth A->B in degrees
			arrival = model.get_travel_times(source_depth_in_km=NZ_depth,
									 distance_in_degree = dd2[ii]*2, phase_list=['PKiKP'])
			tt2[ii] = 0.5 * arrival[0].time # twice the  distance, half the traveltime for one half PKiKP
			rayp[ii] = arrival[0].ray_param  # ray parameter
		else:
			tt1[ii]        = float('nan')
			tt2[ii]        = float('nan')
			baz[ii]        = float('nan')
			rayp[ii]       = float('nan')
		PKiKP_data[ii] = tt1[ii] + tt2[ii]

first_arrival = np.nanmin(PKiKP_data)
print(f'first arrival {first_arrival:.0f}')
PKiKP_diff = PKiKP_data - first_arrival

for ix in x:
	for iy in y:
		ii  = iy + ix*ny
		if ix == 0:
#			iii = int(iy + ny(nx-1)) # wrap-around case, because the endpoint repeats the initial point
#			iii = int(iy + (ix-1)*ny + nx) # wrap-around case, because the endpoint repeats the initial point
			iii = iy + (nx-2)*ny
		else:
			iii = ii - ny
		PKiKP_rot[ii] = (PKiKP_diff[iii] - PKiKP_diff[ii])/5 # divide because this was 5° shift, corrected to 1°


# arrival relative to first arrival
plt.figure(1)
grid1 = PKiKP_diff.reshape((nx, ny))

R, Theta = np.meshgrid(r, theta)
X1 = R*np.cos(Theta)
X2 = R*np.sin(Theta)

plt.pcolormesh(X1, X2, grid1, cmap=cm.gist_rainbow)
plt.colorbar()
plt.axis('equal')
plt.title('Relative arrival time (s)')

# time shift from 1° rotation
plt.figure(2)
grid2 = PKiKP_rot.reshape((nx, ny))
plt.pcolormesh(X1, X2, grid2, cmap=cm.gist_rainbow)
plt.colorbar()
plt.axis('equal')
plt.title('Time shift from 1° rotation (s)')

# ray parameter
plt.figure(3)
grid3 = rayp.reshape((nx, ny))
plt.pcolormesh(X1, X2, grid3, cmap=cm.gist_rainbow)
plt.colorbar()
plt.axis('equal')
plt.title('Ray parameter')

# ray parameter
plt.figure(4)
grid4 = baz.reshape((nx, ny))
plt.pcolormesh(X1, X2, grid4, cmap=cm.gist_rainbow)
plt.colorbar()
plt.axis('equal')
plt.title('Back azimuth')
plt.show()

#plt.figure(2)
#color_map = plt.imshow(grid1, extent=(xlon[0], xlon[-1], xlat[0], xlat[-1]),
#           interpolation='nearest', cmap=cm.gist_rainbow, origin = 'lower')
#color_map.set_cmap = cm.gist_rainbow
#plt.colorbar()
#
#plt.xlabel('longitude (°)')
#plt.ylabel('latitude (°)')
#plt.title('scatterer to receiver, Relative arrival time (s)')
#plt.show()

'''
# rotate scatterer location transversely 10 km
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
'''