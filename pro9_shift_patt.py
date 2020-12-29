#!/usr/bin/env python3
# calculate predicted timing and timeshift of waves scattered from the ICB
"""
Calculate accurate estimate of expected pattern of time change with lag time and back-azimuth from rotation
Started on Sun Mar 24 09:48:49 2019
1.  grid ICB scattering points in lat & lon
2.  for each point, compute angular distance to both source and receiver
3.  compute travel time to twice that distance for PKiKP (divide by two for segment time)
4.  also compute back-azimuth and slowness receiver to scatterer
5.  do same for scattering points 5° rotated in lat & lon
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
from scipy.interpolate import griddata

source = 0  # 0 Amchitka, 1 is Northern NZ, 2 is southern NZ
mark_long = 1 # 1 will make longitudes of 0, 90, 180, and 270 plain, lat of 45 a divot
LASA_lat =   46.7    # °N
LASA_lon = -106.22   # °E
Am_lat = 51.416
Am_lon = 179.18
# LASA_rad =    0.4    # ° radius to exclude outer rings (°)
NZN_lat   = 73.393
NZN_lon   = 54.929
NZS_lat   = 70.80
NZS_lon   = 53.96
depth = 0

if source == 0:
	s_lat = Am_lat
	s_lon = Am_lon
elif source == 1:
	s_lat = NZN_lat
	s_lon = NZN_lon
elif source == 2:
	s_lat = NZS_lat
	s_lon = NZS_lon

min_lat = 0
dlat = 5
dlon = 5
min_slow = -0.02
max_slow =  0.02
nslow    =  21

# test obspy routines
distance = gps2dist_azimuth(LASA_lat,LASA_lon,s_lat,s_lon) # Get traveltimes again, hard to store
dist   = distance[0] # distance in m
s_bazi = distance[1] # azimuth A->B in degrees
az2    = distance[2] # azimuth B->A in degrees)
dist = distance[0]/(1000*111)
arrival = model.get_travel_times(source_depth_in_km=depth,
								 distance_in_degree=dist,phase_list=['PKiKP'])
arrival_time = arrival[0].time
s_rayp = arrival[0].ray_param # ray parmeter in s/°
arrival_ic   = arrival[0].incident_angle
print(f'Distance {dist:.0f} back-azimuth {s_bazi:.0f}')
print(f'Arrival time {arrival_time:.0f} rayp {s_rayp:.0f} incident angle {arrival_ic:.0f}')

# parameterize
tsteps = int((90 - min_lat)/dlat) + 1
r  = np.linspace(90-min_lat, 0, int(tsteps))
rr = np.linspace(min_lat,   90, int(tsteps))
theta = np.linspace(0, 2.*np.pi, int((360/dlon) + 1))

nx = len(theta) # longitude, every dlon degrees
ny = len(rr) # latitude, every dlat degrees

x = range(nx)       # longitude index
y = range(ny)       # latitude index
nn = nx * ny

xlat        = np.zeros(nn) # lat at each grid point
xlon        = np.zeros(nn) # lon at each grid point

PKiKP_data  = np.zeros(nn) # total scatterer arrival time
PKiKP_diff  = np.zeros(nn) # scatterer arrival time (relative to PKiKP)
PKiKP_rot   = np.zeros(nn) # change in time from rotating 1°
slowT       = np.zeros(nn) # transverse slowness of scattered arrival
slowR       = np.zeros(nn) # radial     slowness of scattered arrival

dd1         = np.zeros(nn) # distance from scatterer to NZ
dd2         = np.zeros(nn) # distance from scatterer to LASA
baz         = np.zeros(nn) # back-azimuth from scatterer to LASA
rayp        = np.zeros(nn) # ray parameter at LASA, s/°, p = r sin(theta)/v
tt1         = np.zeros(nn) # time in 1st leg
tt2         = np.zeros(nn) # time in 2nd leg
tt1d        = np.zeros(nn) # time in 5° rotated 1st leg
tt2d        = np.zeros(nn) # time in 5° rotated 2nd leg

# time and angle for grid of scatterer locations
print('Calculation has ' + str(nx) + ' rows of ' + str(ny) + ' points each.')
for ix in x:
	if ix%10 == 0:
		print('starting row ' + str(ix))
	for iy in y:
#		ii = iy*nx + ix
		ii = iy + ix*ny
		# specify scatterer lat 0 to 90°, lon
		xlat[ii] = rr[iy] # specify scatterer
		xlon[ii] = theta[ix]*180/math.pi

		distance1 = gps2dist_azimuth(xlat[ii],xlon[ii],s_lat,s_lon) # source to scatterer
		distance2 = gps2dist_azimuth(LASA_lat,LASA_lon,xlat[ii],xlon[ii]) # scatterer to receiver
		dd1[ii] = distance1[0]/(1000*111)
		dd2[ii] = distance2[0]/(1000*111)
		if (dd1[ii] < 77 and dd2[ii] < 77):
			arrival = model.get_travel_times(source_depth_in_km=depth,
									 distance_in_degree = dd1[ii]*2, phase_list=['PKiKP'])
			tt1[ii] = 0.5 * arrival[0].time # twice the  distance, half the traveltime for one half PKiKP

			baz[ii] = distance2[1] # azimuth A->B in degrees
			arrival = model.get_travel_times(source_depth_in_km=depth,
									 distance_in_degree = dd2[ii]*2, phase_list=['PKiKP'])
			tt2[ii] = 0.5 * arrival[0].time # twice the  distance, half the traveltime for one half PKiKP
			if mark_long == 1: # just to mark longitude visible, gives wrong answer
				if xlon[ii] == 0 or xlon[ii] == 90 or xlon[ii] == 180 or xlon[ii] == 270 or xlat[ii] == 45:
					tt2[ii] = tt2[ii] + 10
			rayp[ii] = arrival[0].ray_param  # ray parameter
			 # radial and transverse slowness of scattered arrival
			slowR[ii] = (rayp[ii] / 6370) * math.cos((math.pi/180)*(baz[ii] - s_bazi))
			slowT[ii] = (rayp[ii] / 6370) * math.sin((math.pi/180)*(baz[ii] - s_bazi))
		else:
			tt1[ii]        = float('nan')
			tt2[ii]        = float('nan')
			baz[ii]        = float('nan')
			rayp[ii]       = float('nan')
			slowT[ii]      = float('nan')
			slowR[ii]      = float('nan')
		PKiKP_data[ii] = tt1[ii] + tt2[ii]

first_arrival = np.nanmin(PKiKP_data)
print(f'first arrival {first_arrival:.0f}')
PKiKP_diff = PKiKP_data - first_arrival

for ix in x: # numerically differentiate to find change in scattered arrival time
	for iy in y:
		ii  = iy + ix*ny
		if ix == 0:
			iii = iy + (nx-2)*ny
		else:
			iii = ii - ny
		PKiKP_rot[ii] = (PKiKP_diff[iii] - PKiKP_diff[ii])/5 # divide because this was 5° shift, corrected to 1°

for ix in x: # make all arrays nan pattern match PKiKP_diff, necessary for griddata
	for iy in y: # just mainly affects one grid point on the end of each row
		ii  = iy + ix*ny
		if math.isnan(PKiKP_rot[ii]):
			PKiKP_diff[ii] = float('nan')
			slowT[ii]      = float('nan')
			slowR[ii]      = float('nan')

# remove nan from arrays to allow 'linear' option in griddata
slowR_pure = slowR[~np.isnan(slowR)]
slowT_pure = slowT[~np.isnan(slowT)]
PKiKP_diff_pure = PKiKP_diff[~np.isnan(PKiKP_diff)]
PKiKP_rot_pure  = PKiKP_rot[ ~np.isnan(PKiKP_rot )] * (-1)

# project delay and shift onto slowness coordinates
x = np.linspace(min_slow, max_slow,nslow)
y = np.linspace(min_slow, max_slow,nslow)
X, Y = np.meshgrid(x,y)

# plot delay time vs slownesses
fig, ax = plt.subplots(1)
Ti = griddata((slowR_pure, slowT_pure), PKiKP_diff_pure, (X, Y), method='linear')
#plt.pcolormesh(Y, X, Ti, cmap=cm.gist_rainbow)
plt.pcolormesh(Y, X, Ti, cmap=cm.gist_yarg)
plt.colorbar()
plt.scatter(slowT_pure, slowR_pure, c='k', alpha=0.2, marker='.')
plt.axis('equal')
#ax.axis([x.min(), x.max(), y.min(), y.max()])
circle1 = plt.Circle((0, 0), 0.019, color='black', fill=False)
ax.add_artist(circle1)
plt.xlabel('Transverse slowness (s/km)')
plt.ylabel('Radial slowness (s/km)')
plt.title('Scattered wave arrival time (after PKiKP, s)')

# plot shift vs slownesses
fig, ax = plt.subplots(1)
Ti = griddata((slowR_pure, slowT_pure), PKiKP_rot_pure, (X, Y), method='linear')
#mask = (x*x + y*y > 0.019*0.019)
#Ti[mask] = np.nan
#plt.pcolormesh(Y, X, Ti, cmap=cm.gist_rainbow_r)
plt.pcolormesh(Y, X, Ti, cmap=cm.coolwarm)
plt.colorbar()
plt.scatter(slowT_pure, slowR_pure, c='k', alpha=0.2, marker='.')
plt.axis('equal')
circle1 = plt.Circle((0, 0), 0.019, color='black', fill=False)
ax.add_artist(circle1)
plt.xlabel('Transverse slowness (s/km)')
plt.ylabel('Radial slowness (s/km)')
plt.title('Time shift of scattered arrival (s/° rotation)')

#%% more primitive plot projected onto slowness rectangle rather than sphere
#x = np.linspace(min_slow, max_slow,nslow)
#y = np.linspace(min_slow, max_slow,nslow)
#X, Y = np.meshgrid(x,y)
#mask = (x*x + y*y > 0.019*0.019)
##mask = (math.sqrt((x*x) + (y*y))>0.019)
#Ti = griddata((slowR, slowT), PKiKP_rot, (X, Y), method='nearest')
#Ti[mask] = np.nan
#plt.pcolormesh(Y, X, Ti, cmap=cm.gist_rainbow)
#plt.colorbar()
#plt.scatter(slowT, slowR, c='k', alpha=0.2, marker='.')
#plt.axis('equal')
#plt.xlabel('Transverse slowness (s/km)')
#plt.ylabel('Radial slowness (s/km)')
#plt.title('Time shift from 1° rotation (s)')

# arrival relative to first arrival
fig, ax = plt.subplots(1)
grid1 = PKiKP_diff.reshape((nx, ny))

R, Theta = np.meshgrid(r, theta)
X1 = R*np.cos(Theta)
X2 = R*np.sin(Theta)

plt.pcolormesh(X1, X2, grid1, cmap=cm.gist_rainbow)
plt.colorbar()
plt.axis('equal')
circle1 = plt.Circle((0, 0), 75, color='black', fill=False)
ax.add_artist(circle1)
plt.xlabel('°E of lon = 0,180')
plt.ylabel('°N of lon = 90,270')
plt.title('Relative arrival time (s)')

# time shift from 1° rotation
fig, ax = plt.subplots(1)
grid2 = PKiKP_rot.reshape((nx, ny))
plt.pcolormesh(X1, X2, grid2, cmap=cm.gist_rainbow)
plt.colorbar()
plt.axis('equal')
circle1 = plt.Circle((0, 0), 75, color='black', fill=False)
ax.add_artist(circle1)
plt.xlabel('°E of lon = 0,180')
plt.ylabel('°N of lon = 90,270')
plt.title('Time shift from 1° rotation (s)')

# ray parameter
fig, ax = plt.subplots(1)
grid3 = rayp.reshape((nx, ny))
plt.pcolormesh(X1, X2, grid3, cmap=cm.gist_rainbow)
plt.colorbar()
plt.axis('equal')
circle1 = plt.Circle((0, 0), 75, color='black', fill=False)
ax.add_artist(circle1)
plt.xlabel('°E of lon = 0,180')
plt.ylabel('°N of lon = 90,270')
plt.title('Ray parameter')

# back azimuth
fig, ax = plt.subplots(1)
grid4 = baz.reshape((nx, ny))
plt.pcolormesh(X1, X2, grid4, cmap=cm.gist_rainbow)
plt.colorbar()
plt.axis('equal')
circle1 = plt.Circle((0, 0), 75, color='black', fill=False)
ax.add_artist(circle1)
plt.xlabel('°E of lon = 0,180')
plt.ylabel('°N of lon = 90,270')
plt.title('Back azimuth')

# radial slowness
fig, ax = plt.subplots(1)
grid5 = slowR.reshape((nx, ny))
plt.pcolormesh(X1, X2, grid5, cmap=cm.gist_rainbow)
plt.colorbar()
plt.axis('equal')
circle1 = plt.Circle((0, 0), 75, color='black', fill=False)
ax.add_artist(circle1)
plt.xlabel('°E of lon = 0,180')
plt.ylabel('°N of lon = 90,270')
plt.title('Radial slowness')

# transverse slowness
fig, ax = plt.subplots(1)
grid6 = slowT.reshape((nx, ny))
plt.pcolormesh(X1, X2, grid6, cmap=cm.gist_rainbow)
plt.colorbar()
plt.axis('equal')
#circle1 = plt.Circle((0, 0), 75, color='black', fill=False)
#ax.add_artist(circle1)
plt.xlabel('°E of lon = 0,180')
plt.ylabel('°N of lon = 90,270')
plt.title('Transverse slowness')

plt.show()

# flawed lat-lon plot
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