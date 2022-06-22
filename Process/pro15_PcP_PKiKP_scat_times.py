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

#%% Define parameters
source = 3      # 0 Amchitka, 1 is Northern NZ, 2 is southern NZ, 3 is Kawakatsu event, 4 is test case
array = 1       # 0 is LASA, 1 is Hinet
phase = 'PcP'   # phase can be PcP or PKiKP
NS = True       # Use N-S coordinates, if false, use radial and transverse

mark_long = False # 1 will make longitudes of 0, 90, 180, and 270 plain, lat of 45 a divot

if array == 0:   # LASA
    array_lat =   46.7    # °N
    array_lon = -106.22   # °E
elif array == 1:   # Hinet
    array_lat =  36.3  # °N
    array_lon = 138.5  # °E

if source == 0: # 0 Amchitka
    s_lat = 51.416
    s_lon = 179.18
elif source == 1: # 1 is Northern NZ
    s_lat   = 73.393
    s_lon   = 54.929
elif source == 2: # 2 is southern NZ
    s_lat   = 70.80
    s_lon   = 53.96
elif source == 3: # 3 is Kawakatsu event
    s_lat = 21.635
    s_lon = 142.988

if source == 4:
    array_lat = 45   # °N
    array_lon = 0  # °E
    s_lat     = 0
    s_lon     = 0

depth = 0

dlatlon = 2  # grid spacing in degrees
min_slow = -0.04
max_slow =  0.04
nslow    =  81

#%% Compute distance, azimuth, etc.
distance = gps2dist_azimuth(array_lat,array_lon,s_lat,s_lon)
dist   = distance[0] # distance in m
s_bazi = distance[1] # azimuth A->B in degrees
az2    = distance[2] # azimuth B->A in degrees
dist = distance[0]/(1000*111)
arrival = model.get_travel_times(source_depth_in_km=depth, distance_in_degree=dist,phase_list=[phase])
arrival_time = arrival[0].time
print(f'Distance {dist:.0f} back-azimuth {s_bazi:.0f}')
print(f'Arrival time {arrival_time:.0f}')

#%% Define grids
rr = np.linspace(-90, 90, int((180/dlatlon) + 1))    # steps in lat
theta = np.linspace(0, 360, int((360/dlatlon) + 1))  # steps in lon

nx = len(theta) # longitude, every dlon degrees
ny = len(rr) # latitude, every dlat degrees

x = range(nx)       # array with longitude index
y = range(ny)       # array with latitude index
nn = nx * ny

xlat        = np.zeros(nn) # lat at each grid point
xlon        = np.zeros(nn) # lon at each grid point

phase_time  = np.zeros(nn) # total scatterer arrival time
phase_diff  = np.zeros(nn) # scatterer arrival time (relative to phase)
slowT       = np.zeros(nn) # transverse slowness of scattered arrival
slowR       = np.zeros(nn) # radial     slowness of scattered arrival

dd1         = np.zeros(nn) # distance from scatterer to source
dd2         = np.zeros(nn) # distance from scatterer to array
baz         = np.zeros(nn) # back-azimuth from scatterer to array
rayp        = np.zeros(nn) # ray parameter at station, s/°, p = r sin(theta)/v
tt1         = np.zeros(nn) # time in 1st leg
tt2         = np.zeros(nn) # time in 2nd leg

#%% Fill grids
# time and angle for grid of scatterer locations
print('Calculation has ' + str(nx) + ' rows of ' + str(ny) + ' points each.')

is_nan = 0 # counter, in range of phase
no_nan = 0 # counter, out range of phase
for ix in x:  # specify scatterer lon 0 to 360°
    if ix%10 == 0:
        print('starting row ' + str(ix))
    for iy in y: # specify scatterer lat -90 to 90°
        ii = iy + ix*ny
        xlat[ii] = rr[iy] # specify scatterer
        xlon[ii] = theta[ix]

        distance1 = gps2dist_azimuth(xlat[ii],xlon[ii],s_lat,s_lon) # source to scatterer
        distance2 = gps2dist_azimuth(array_lat,array_lon,xlat[ii],xlon[ii]) # scatterer to receiver
        dd1[ii] = distance1[0]/(1000*111)
        dd2[ii] = distance2[0]/(1000*111)
        # beyond 2*77 = 144° for PKiKP and 2*49 = 98°, reflected phases are not returned
        if (phase == 'PKiKP' and (dd1[ii] < 77 and dd2[ii] < 77)) or (phase == 'PcP' and (dd1[ii] < 49 and dd2[ii] < 49)):
            arrival1 = model.get_travel_times(source_depth_in_km=depth, distance_in_degree = dd1[ii]*2, phase_list=[phase])
            arrival2 = model.get_travel_times(source_depth_in_km=depth, distance_in_degree = dd2[ii]*2, phase_list=[phase])
            tt1[ii] = 0.5 * arrival1[0].time # twice the  distance, half the traveltime for one half reflected phase
            tt2[ii] = 0.5 * arrival2[0].time # twice the  distance, half the traveltime for one half reflected phase
            baz[ii] = distance2[1] # azimuth A->B in degrees

            if mark_long: # just to mark longitude visible, gives wrong answer
                if xlon[ii] == 0 or xlon[ii] == 90 or xlon[ii] == 180 or xlon[ii] == 270 or xlat[ii] == 45:
                    tt2[ii] = tt2[ii] + 100

            rayp[ii] = arrival2[0].ray_param  # ray parameter
             # radial and transverse slowness of scattered arrival
            if NS:
                slowR[ii] = (rayp[ii] / 6370) * math.cos((math.pi/180)*baz[ii])
                slowT[ii] = (rayp[ii] / 6370) * math.sin((math.pi/180)*baz[ii])
            else:
                slowR[ii] = (rayp[ii] / 6370) * math.cos((math.pi/180)*(baz[ii] - s_bazi))
                slowT[ii] = (rayp[ii] / 6370) * math.sin((math.pi/180)*(baz[ii] - s_bazi))
            no_nan += 1
        else:
            tt1[ii]        = float('nan')
            tt2[ii]        = float('nan')
            baz[ii]        = float('nan')
            rayp[ii]       = float('nan')
            slowT[ii]      = float('nan')
            slowR[ii]      = float('nan')
            is_nan += 1
        phase_time[ii] = tt1[ii] + tt2[ii]

first_arrival = np.nanmin(phase_time)      # find first arrival
phase_diff    = phase_time - first_arrival # find all relative times

print(f'first arrival {first_arrival:.0f}')
print('no_nan ' + str(no_nan) + ' is_nan ' + str(is_nan))

# remove nan from arrays to allow 'linear' option in griddata
slowR_pure = slowR[~np.isnan(slowR)]
slowT_pure = slowT[~np.isnan(slowT)]
phase_diff_pure = phase_diff[~np.isnan(phase_diff)]

# project delay and shift onto slowness coordinates
x = np.linspace(min_slow, max_slow, nslow)
y = np.linspace(min_slow, max_slow, nslow)
X, Y = np.meshgrid(x,y)

# plot delay time vs slownesses
fig, ax = plt.subplots(1)
Ti = griddata((slowR_pure, slowT_pure), phase_diff_pure, (X, Y), method='linear')
plt.pcolormesh(Y, X, Ti, cmap=cm.gist_rainbow, vmin = 0, vmax = 2)
plt.colorbar()
plt.scatter(slowT_pure, slowR_pure, c='k', alpha=0.2, marker='.')
plt.axis('equal')
#ax.axis([x.min(), x.max(), y.min(), y.max()])
circle1 = plt.Circle((0, 0), 0.019, color='black', fill=False)
ax.add_artist(circle1)
circle2 = plt.Circle((0, 0), 0.040, color='black', fill=False)
ax.add_artist(circle2)
if NS:
    plt.xlabel('East slowness (s/km)')
    plt.ylabel('North slowness (s/km)')
else:
    plt.xlabel('Transverse slowness (s/km)')
    plt.ylabel('Radial slowness (s/km)')
plt.title('Scattered wave arrival time (after phase, s)')