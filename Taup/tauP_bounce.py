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

#%% Set parameters
source = 3  # 0 Amchitka, 1 is Northern NZ, 2 is southern NZ, 3 is symmetric test case
mark_long = 0 # 1 will make longitudes of 0, 90, 180, and 270 plain, lat of 45 a divot
LASA_lat =   46.7    # °N
LASA_lon = -106.22   # °E
Am_lat = 51.416
Am_lon = 179.18
# LASA_rad =    0.4    # ° radius to exclude outer rings (°)
NZN_lat   = 73.393
NZN_lon   = 54.929
NZS_lat   = 70.80
NZS_lon   = 53.96
phase1 = 'PKiKP'
phase2 = 'PKiKP'  # last leg before coming to array
max_Pdiff_dist = 135
depth = 0

#%% Geometry parameters
if source == 0:
    s_lat = Am_lat   # °N
    s_lon = Am_lon   # °E
elif source == 1:
    s_lat = NZN_lat
    s_lon = NZN_lon
elif source == 2:
    s_lat = NZS_lat
    s_lon = NZS_lon
elif source == 3:
    s_lat =   45
    s_lon =    0
    LASA_lat = 0
    LASA_lon = 0

#%% Raypath midpoints: near or far side?
if   (phase1 == 'PcP' or phase1 == 'Pdiff' or phase1 == 'PKiKP' or phase1 == 'ScS' or phase1 == 'P'):
    nearside1 = True
elif (phase1 == 'PKKP' or phase1 == 'PKIKKIKP' or phase1 == 'PKPPKP' or phase1 == 'PKIKPPKIKP'):
    nearside1 = False
else:
    print('Source side raypath not labeled as either near or far side')
    quit()
if   (phase2 == 'PcP' or phase2 == 'Pdiff' or phase2 == 'PKiKP' or phase2 == 'ScS' or phase1 == 'P'):
    nearside2 = True
elif (phase2 == 'PKKP' or phase2 == 'PKIKKIKP' or phase2 == 'PKPPKP' or phase2 == 'PKIKPPKIKP'):
    nearside2 = False
else:
    print('Source side raypath not labeled as either near or far side')
    quit()
print('nearside1 ' + str(nearside1) + ' nearside2 ' + str(nearside2))

# min_lat = -90
# min_lat = 0
dlat = 2.5
dlon = 2.5
min_slow = -0.04
max_slow =  0.04
nslow    =  81
extra_plots = False

#%% Get source-array azimuth for rotation
distance = gps2dist_azimuth(LASA_lat,LASA_lon,s_lat,s_lon)
s_to_array_dist = distance[0]/(1000*111)
s_bazi = distance[1] # azimuth A->B in degrees, for rotating to R, T in plotting

#%% Define scatterer grid in Earth co-ordinates
lat_array = np.linspace(  -90,  90, int((180/dlat) + 1))  # array for latitude
lon_array = np.linspace( -180, 180, int((360/dlon) + 1))  # array for longitude
# lon_array = np.linspace( -90, 90, int((180/dlon) + 1))  # array for longitude
print('Latitude runs from '  + str(min(lat_array)) + ' to ' + str(max(lat_array)))
print('Longitude runs from ' + str(min(lon_array)) + ' to ' + str(max(lon_array)))

nx = len(lon_array)        # longitude, every dlon degrees
ny = len(lat_array)        # latitude,  every dlat degrees
print('Calculation has ' + str(nx) + ' rows of ' + str(ny) + ' points each.')

x = range(nx)              # array with longitude index
y = range(ny)              # array with latitude index
nn = nx * ny

xlat        = np.zeros(nn) # lat at each grid point
xlon        = np.zeros(nn) # lon at each grid point

all_times   = np.zeros(nn) # total scatterer arrival time
slowT       = np.zeros(nn) # transverse slowness of scattered arrival
slowR       = np.zeros(nn) # radial     slowness of scattered arrival

#%% Fill grids
# time and angle for grid of scatterer locations
full_path     = 0
only_one_path = 0
neither_path  = 0
neg_ray2  = 0
mult_arr1 = 0
mult_arr2 = 0
for ix in x:
    if ix%10 == 0: # report calculation progress
        print('starting row ' + str(ix) + ' of ' + str(len(x)))
    for iy in y:
        # turn 2D array to 1D vector
        ii = iy + ix*ny
        # specify scatterer co-ord in this iteration
        xlat[ii] = lat_array[iy]
        xlon[ii] = lon_array[ix]

        distance1 = gps2dist_azimuth(xlat[ii],xlon[ii],s_lat,s_lon) # source to scatterer
        distance2 = gps2dist_azimuth(LASA_lat,LASA_lon,xlat[ii],xlon[ii]) # scatterer to receiver
        dd1 = (distance1[0]*2.)/(1000*111) # double distance to include down and back, meters to degrees
        dd2 = (distance2[0]*2.)/(1000*111)
        # distance back to surface
        arrival1 = model.get_travel_times(source_depth_in_km=depth,distance_in_degree = dd1, phase_list=[phase1])
        arrival2 = model.get_travel_times(source_depth_in_km=depth,distance_in_degree = dd2, phase_list=[phase2])
        #%%  ---- Maybe viable pair of raypaths
        if (len(arrival1) >= 1) and (len(arrival2) >= 1):
            if len(arrival1) > 1:  # are there multiple arrivals for some scatterers?
                mult_arr1 += 1
            if len(arrival2) > 1:
                mult_arr2 += 1
            rayp1 = arrival1[0].ray_param  # ray parameter
            rayp2 = arrival2[0].ray_param
            baz1 = distance1[1]            # azimuth A->B in degrees
            baz2 = distance2[1]
            az1 = distance1[2]            # azimuth B->A in degrees
            az2 = distance2[2]

             # radial and transverse slowness of scattered arrival
            slowR[ii] = (rayp2 / 6370) * math.cos((math.pi/180)*(baz2 - s_bazi))
            slowT[ii] = (rayp2 / 6370) * math.sin((math.pi/180)*(baz2 - s_bazi))
            baz_scat = math.atan2(slowT[ii], slowR[ii]) * 180./math.pi
            if baz_scat < 0:
                baz_scat += 360
            # print(f'Rayp {rayp1:.1f} diff  {baz_scat - baz1:.1f} lat {xlat[ii]:.1f} lon {xlon[ii]:.1f} baz2 {baz1:.1f} baz_scat {baz_scat:.1f}')

            #%% ---- ---- Viable: both ray segments are on either minor or major arc
            ray1_works = (dd1 < 180 and nearside1) or (dd1 > 180 and not nearside1)
            ray2_works = (dd2 < 180 and nearside2) or (dd2 > 180 and not nearside2)
            ray1_pdiff_ok = (phase1 != 'Pdiff') or (dd1 < max_Pdiff_dist)
            ray2_pdiff_ok = (phase2 != 'Pdiff') or (dd2 < max_Pdiff_dist)
            # if (((dd1 < 180 and nearside1) or (dd1 > 180 and not nearside1)) and ((dd2 < 180 and nearside2) or (dd2 > 180 and not nearside2))):
            if ray1_works and ray2_works and ray1_pdiff_ok and ray2_pdiff_ok:
                # print(f'Yes  Rays {rayp1:.1f} {rayp2:.1f} ddistances {dd1:5.1f} {dd2:5.1f} b_az {baz1:5.1f} {baz2:5.1f} az {az1:5.1f} {az2:5.1f} slows {slowR[ii]:6.3f} {slowT[ii]:6.3f} baz_scat {baz_scat:5.1f}')
                full_path += 1
                tt1 = 0.5 * arrival1[0].time   # distance doubled, so traveltime halved
                tt2 = 0.5 * arrival2[0].time
            else:
                only_one_path += 1
                tt1        = float('nan')
                tt2        = float('nan')
            if mark_long == 1: # just to mark longitude visible, gives wrong answer
                if xlon[ii] == 0 or xlon[ii] == 90 or xlon[ii] == 180 or xlon[ii] == 270 or xlat[ii] == 45:
                    tt2[ii] = tt2[ii] + 10
        #%% ---- ---- Not viable
        elif(len(arrival1) >= 1): # only source-scatterer path viable
            only_one_path += 1
            tt1        = float('nan')
            tt2        = float('nan')
             # radial and transverse slowness of scattered arrival
            slowR[ii]      = float('nan')
            slowT[ii]      = float('nan')
        elif(len(arrival2) >= 1): # only scatterer-receiver path viable
            only_one_path += 1
            tt1        = float('nan')
            tt2        = float('nan')
            baz2 = distance2[1] # azimuth A->B in degrees
            rayp2 = arrival2[0].ray_param  # ray parameter
             # radial and transverse slowness of scattered arrival
            slowR[ii] = (rayp2 / 6370) * math.cos((math.pi/180)*(baz2 - s_bazi))
            slowT[ii] = (rayp2 / 6370) * math.sin((math.pi/180)*(baz2 - s_bazi))
        else: # still want slownesses for empty symbol
            neither_path += 1
            tt1        = float('nan')
            tt2        = float('nan')
            slowR[ii]      = float('nan')
            slowT[ii]      = float('nan')
        all_times[ii] = tt1 + tt2

print('Full path ' + str(full_path) + ' only one path ' + str(only_one_path) + ' neither path connects ' + str(neither_path) + ' neg_ray2 ' + str(neg_ray2))
print('instances of multiple arrivals for 1st path: ' + str(mult_arr1) + ' for 2nd path: ' + str(mult_arr2))

# remove nan from arrays to allow 'linear' option in griddata
slowR_pure = slowR[~np.isnan(all_times)] # record slownesses that contribute
slowT_pure = slowT[~np.isnan(all_times)]
slowR_nope = slowR[~np.isnan(slowR)] # record slownesses that had a second leg
slowT_nope = slowT[~np.isnan(slowR)]
calc_times = all_times[~np.isnan(all_times)] # select the times at those slownesses

# project delay and shift onto slowness coordinates
x = np.linspace(min_slow, max_slow,nslow)
y = np.linspace(min_slow, max_slow,nslow)
X, Y = np.meshgrid(x,y)

# plot delay time vs slownesses

# colored symbols
fig, ax = plt.subplots(1)
if len(calc_times > 0):
    t_min=min(calc_times)
    t_max=max(calc_times)
    print('t_min ' + str(t_min) + ' t_max ' + str(t_max))
print('Tslows ' + str(len(slowT_pure)) +' Rslows ' + str(len(slowR_pure)) + ' calcs ' + str(len(calc_times)))
print('nopesT ' + str(len(slowT_nope)) +' nopesR ' + str(len(slowR_nope)))

#%% Colored dots
plt.scatter(slowT_nope, slowR_nope, c = 'k', s=2)
if len(calc_times > 0):
    plt.scatter(slowT_pure, slowR_pure, c=calc_times, cmap=cm.gist_rainbow, s=50, vmin=t_min, vmax=t_max)
    plt.colorbar(label = 'arrival time, s')

    #%% Colored interpolation
    Ti = griddata((slowR_pure, slowT_pure), calc_times, (X, Y), method='linear')
    # plt.pcolormesh(Y, X, Ti, cmap=cm.gist_rainbow)
    # plt.pcolormesh(Y, X, Ti, cmap=cm.gist_yarg)

    # plt.scatter(slowT_nope, slowR_nope, c='k', marker='o')
    # plt.scatter(slowT_nope, slowR_nope)
plt.axis('equal')
#ax.axis([x.min(), x.max(), y.min(), y.max()])
circle1 = plt.Circle((0, 0), 0.0175, color='black', fill=False)
ax.add_artist(circle1)
circle1 = plt.Circle((0, 0), 0.040, color='black', fill=False)
ax.add_artist(circle1)
plt.xlim([-0.04, 0.04])
plt.ylim([-0.04, 0.04])
plt.xlabel('Transverse slowness (s/km)')
plt.ylabel('Radial slowness (s/km)')

if (phase1 == 'Pdiff') or (phase2 == 'Pdiff'):
    plt.title(f'{phase1} to {phase2} arrival time (s) dist = {s_to_array_dist:.1f}, max Pdiff= {max_Pdiff_dist:.0f}')
else:
    plt.title(f'{phase1} to {phase2} scattered wave arrival time (s), dist= {s_to_array_dist:.1f}')