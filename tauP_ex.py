#!/usr/bin/env python
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from obspy.geodetics import gps2dist_azimuth
from obspy import Stream
import numpy as np
import matplotlib.pyplot as plt
from obspy.taup import TauPyModel
from obspy.taup import plot_travel_times
model = TauPyModel(model="iasp91")

#%%  Plot ray paths
#arrivals = model.get_ray_paths(source_depth_in_km=55, distance_in_degree=150)
#arrival = arrivals[0]
#arrivals = model.get_ray_paths(source_depth_in_km=55, distance_in_degree=150,
#							   phase_list = ["PKP", "PKIKP", "PKiKP"])
#ax = arrivals.plot_rays()

#%% another way to plot travel time curves
fig, ax = plt.subplots(figsize=(9, 9))
#ax = plot_travel_times(source_depth=10, phase_list=["PKP", "PKIKP", "PKiKP"],
ax = plot_travel_times(source_depth=0,
					   ax=ax, fig=fig, min_degrees=0, max_degrees=180,
					   npoints=200, verbose=True)

#%%  Plot travel time curves at one distance
#arrivals = model.get_travel_times(source_depth_in_km=55,
#								  distance_in_degree=150, phase_list = ["PKP", "PKIKP", "PKiKP"])
#ax = arrivals.plot_times()
#ev_dep = 50.
#dist = 50.
#dphase = 'P'
#arrivals = model.get_travel_times(source_depth_in_km=ev_dep,distance_in_degree=dist,phase_list=[dphase])
