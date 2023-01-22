#!/usr/bin/env python
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from obspy.geodetics import gps2dist_azimuth
from obspy import Stream
import numpy as np
import matplotlib.pyplot as plt
from obspy.taup import TauPyModel
from obspy.taup import plot_travel_times
from obspy.taup.tau import plot_ray_paths
model = TauPyModel(model="iasp91")
from obspy.geodetics import gps2dist_azimuth

ref1_lat = -56
ref1_lon = -27
ref2_lat = -57
ref2_lon = -26
ref3_lat = -59
ref3_lon = -25

sta_name = 'ALE'; ev_lat = 82.5033; ev_lon = -62.35
ref1_distance = gps2dist_azimuth(ref1_lat,ref1_lon,ev_lat,ev_lon)
ref2_distance = gps2dist_azimuth(ref2_lat,ref2_lon,ev_lat,ev_lon)
ref3_distance = gps2dist_azimuth(ref3_lat,ref3_lon,ev_lat,ev_lon)
distance1 = ref1_distance[0]/(1000*111)
distance2 = ref2_distance[0]/(1000*111)
distance3 = ref3_distance[0]/(1000*111)
print(f'{sta_name}   1 {distance1:.1f}   2 {distance2:.1f}   3 {distance3:.1f}  ')

sta_name = 'DAG'; ev_lat =  76.77; ev_lon = -18.65
ref1_distance = gps2dist_azimuth(ref1_lat,ref1_lon,ev_lat,ev_lon)
ref2_distance = gps2dist_azimuth(ref2_lat,ref2_lon,ev_lat,ev_lon)
ref3_distance = gps2dist_azimuth(ref3_lat,ref3_lon,ev_lat,ev_lon)
distance1 = ref1_distance[0]/(1000*111)
distance2 = ref2_distance[0]/(1000*111)
distance3 = ref3_distance[0]/(1000*111)
print(f'{sta_name}   1 {distance1:.1f}   2 {distance2:.1f}   3 {distance3:.1f}  ')

sta_name = 'LVZ'; ev_lat =  67.8979; ev_lon = 34.6514
ref1_distance = gps2dist_azimuth(ref1_lat,ref1_lon,ev_lat,ev_lon)
ref2_distance = gps2dist_azimuth(ref2_lat,ref2_lon,ev_lat,ev_lon)
ref3_distance = gps2dist_azimuth(ref3_lat,ref3_lon,ev_lat,ev_lon)
distance1 = ref1_distance[0]/(1000*111)
distance2 = ref2_distance[0]/(1000*111)
distance3 = ref3_distance[0]/(1000*111)
print(f'{sta_name}   1 {distance1:.1f}   2 {distance2:.1f}   3 {distance3:.1f}  ')

sta_name = 'KIEV'; ev_lat =  50.6944; ev_lon = 29.2083
ref1_distance = gps2dist_azimuth(ref1_lat,ref1_lon,ev_lat,ev_lon)
ref2_distance = gps2dist_azimuth(ref2_lat,ref2_lon,ev_lat,ev_lon)
ref3_distance = gps2dist_azimuth(ref3_lat,ref3_lon,ev_lat,ev_lon)
distance1 = ref1_distance[0]/(1000*111)
distance2 = ref2_distance[0]/(1000*111)
distance3 = ref3_distance[0]/(1000*111)
print(f'{sta_name}   1 {distance1:.1f}   2 {distance2:.1f}   3 {distance3:.1f}  ')

sta_name = 'OBN'; ev_lat =  55.1146; ev_lon = 36.5674
ref1_distance = gps2dist_azimuth(ref1_lat,ref1_lon,ev_lat,ev_lon)
ref2_distance = gps2dist_azimuth(ref2_lat,ref2_lon,ev_lat,ev_lon)
ref3_distance = gps2dist_azimuth(ref3_lat,ref3_lon,ev_lat,ev_lon)
distance1 = ref1_distance[0]/(1000*111)
distance2 = ref2_distance[0]/(1000*111)
distance3 = ref3_distance[0]/(1000*111)
print(f'{sta_name}   1 {distance1:.1f}   2 {distance2:.1f}   3 {distance3:.1f}  ')

sta_name = 'ARU'; ev_lat =  56.4302; ev_lon = 58.5625
ref1_distance = gps2dist_azimuth(ref1_lat,ref1_lon,ev_lat,ev_lon)
ref2_distance = gps2dist_azimuth(ref2_lat,ref2_lon,ev_lat,ev_lon)
ref3_distance = gps2dist_azimuth(ref3_lat,ref3_lon,ev_lat,ev_lon)
distance1 = ref1_distance[0]/(1000*111)
distance2 = ref2_distance[0]/(1000*111)
distance3 = ref3_distance[0]/(1000*111)
print(f'{sta_name}   1 {distance1:.1f}   2 {distance2:.1f}   3 {distance3:.1f}  ')

sta_name = 'AAK'; ev_lat =  42.6375; ev_lon = 74.4942
ref1_distance = gps2dist_azimuth(ref1_lat,ref1_lon,ev_lat,ev_lon)
ref2_distance = gps2dist_azimuth(ref2_lat,ref2_lon,ev_lat,ev_lon)
ref3_distance = gps2dist_azimuth(ref3_lat,ref3_lon,ev_lat,ev_lon)
distance1 = ref1_distance[0]/(1000*111)
distance2 = ref2_distance[0]/(1000*111)
distance3 = ref3_distance[0]/(1000*111)
print(f'{sta_name}   1 {distance1:.1f}   2 {distance2:.1f}   3 {distance3:.1f}  ')

sta_name = 'NIL'; ev_lat =  33.6506; ev_lon = 73.2686
ref1_distance = gps2dist_azimuth(ref1_lat,ref1_lon,ev_lat,ev_lon)
ref2_distance = gps2dist_azimuth(ref2_lat,ref2_lon,ev_lat,ev_lon)
ref3_distance = gps2dist_azimuth(ref3_lat,ref3_lon,ev_lat,ev_lon)
distance1 = ref1_distance[0]/(1000*111)
distance2 = ref2_distance[0]/(1000*111)
distance3 = ref3_distance[0]/(1000*111)
print(f'{sta_name}   1 {distance1:.1f}   2 {distance2:.1f}   3 {distance3:.1f}  ')

sta_name = 'LSA'; ev_lat =  29.703; ev_lon = 91.127
ref1_distance = gps2dist_azimuth(ref1_lat,ref1_lon,ev_lat,ev_lon)
ref2_distance = gps2dist_azimuth(ref2_lat,ref2_lon,ev_lat,ev_lon)
ref3_distance = gps2dist_azimuth(ref3_lat,ref3_lon,ev_lat,ev_lon)
distance1 = ref1_distance[0]/(1000*111)
distance2 = ref2_distance[0]/(1000*111)
distance3 = ref3_distance[0]/(1000*111)
print(f'{sta_name}   1 {distance1:.1f}   2 {distance2:.1f}   3 {distance3:.1f}  ')

sta_name = 'CHTO'; ev_lat =  18.79; ev_lon = 98.9769
ref1_distance = gps2dist_azimuth(ref1_lat,ref1_lon,ev_lat,ev_lon)
ref2_distance = gps2dist_azimuth(ref2_lat,ref2_lon,ev_lat,ev_lon)
ref3_distance = gps2dist_azimuth(ref3_lat,ref3_lon,ev_lat,ev_lon)
distance1 = ref1_distance[0]/(1000*111)
distance2 = ref2_distance[0]/(1000*111)
distance3 = ref3_distance[0]/(1000*111)
print(f'{sta_name}   1 {distance1:.1f}   2 {distance2:.1f}   3 {distance3:.1f}  ')

sta_name = 'XAN'; ev_lat =  34.0394; ev_lon = 108.9214
ref1_distance = gps2dist_azimuth(ref1_lat,ref1_lon,ev_lat,ev_lon)
ref2_distance = gps2dist_azimuth(ref2_lat,ref2_lon,ev_lat,ev_lon)
ref3_distance = gps2dist_azimuth(ref3_lat,ref3_lon,ev_lat,ev_lon)
distance1 = ref1_distance[0]/(1000*111)
distance2 = ref2_distance[0]/(1000*111)
distance3 = ref3_distance[0]/(1000*111)
print(f'{sta_name}   1 {distance1:.1f}   2 {distance2:.1f}   3 {distance3:.1f}  ')

#%% Plot a few travel time curves
# ax = plot_travel_times(source_depth=0, npoints=200, min_degrees=90, max_degrees=180, verbose=True, phase_list=['PKIKP', 'PKP', 'PKiKP'])

#%%  Plot ray paths
# fig, ax = plt.subplots(subplot_kw=dict(polar=True))
# ax = plot_ray_paths(source_depth=100, ax=ax, legend=True, label_arrivals=True, min_degrees=90, max_degrees=180, fig=fig, phase_list=['PKiKP', 'PKIKP', 'PKP'], npoints=18)
