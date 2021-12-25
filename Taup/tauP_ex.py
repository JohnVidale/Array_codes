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

#%%  Plot all ray paths
# arrivals = model.get_ray_paths(source_depth_in_km=55, distance_in_degree=150)
# ax = arrivals.plot_rays()

#%%  Plot ray paths
#arrivals = model.get_ray_paths(source_depth_in_km=55, distance_in_degree=150)
#arrival = arrivals[0]
# arrivals = model.get_ray_paths(source_depth_in_km=55, distance_in_degree=150, phase_list = ["PKP", "PKIKP", "PKiKP"])
# ax = arrivals.plot_rays()

#%%  Plot one ray path at one distance
# ev_dep = 50.
# dist = 50.
# dphase = 'P'
# arrivals = model.get_ray_paths(source_depth_in_km=ev_dep,distance_in_degree=dist,phase_list=[dphase])
# ax = arrivals.plot_rays()

#%% Plot all travel time curves
# fig, ax = plt.subplots(figsize=(9, 9))
# ax = plot_travel_times(source_depth=10, fig=fig)

#%% Plot a few travel time curves
# ax = plot_travel_times(source_depth=0, npoints=200, min_degrees=0, max_degrees=180, verbose=True, phase_list=['PKiKP', 'ScP', 'PKKP', 'Pdiff', 'PKIKKIKP'])
# ax = plot_travel_times(source_depth=0, ax=ax, fig=fig, min_degrees=0, max_degrees=180, npoints=200, verbose=True)

#%%  Plot travel time curves at one distance
#ev_dep = 50.
#dist = 50.
#dphase = 'P'
#arrivals = model.get_travel_times(source_depth_in_km=ev_dep,distance_in_degree=dist,phase_list=[dphase])

#%%  Plot ray paths
fig, ax = plt.subplots(subplot_kw=dict(polar=True))
ax = plot_ray_paths(source_depth=100, ax=ax, legend=True, label_arrivals=True, fig=fig, phase_list=['PKiKP', 'PKKP', 'Pdiff', 'PKIKKIKP'], npoints=36)