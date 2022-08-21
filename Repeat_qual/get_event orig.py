#!/usr/bin/env python
#from __future__ import print_function

from obspy import UTCDateTime
from obspy.clients.fdsn import Client, RoutingClient
from obspy import Stream, Trace
import os
import matplotlib.pyplot as plt
from obspy import read_events, read
import time

#%% Get catalog

start_time_wc = time.time()
# etime      = '2017-06-20T12:54:35.91' # event has to be unique and follow time by less than 1 minute
etime      = '1998-04-12T21:33:47.42' # event has to be unique and follow time by less than 1 minute
ev_lon     = -27.00  # South Sandwich Islands
ev_lat     = -56.19
ev_depth   =  77.5

chan_type  = 'EHZ,HHZ,HNZ,HLZ,BHZ' # e.g., BHZ
start_buff = 10   # Pre-event buffer
end_buff   = 1800 # Time after OT recorded
min_dist   = 0      # Minimum distance loaded, degrees
max_dist   = 30    # Maximum distance loaded, degrees
verbose    = 1       # 1 prints networks and stations
refine     = 1        # Use refined location

client = Client('IRIS')
# client = RoutingClient("iris-federator")
t = UTCDateTime(etime)

# refined time and hypocenter

#%% Make inventory of all stations recording this channel
inventory = client.get_stations(longitude=ev_lon, latitude=ev_lat, starttime =
                                t, endtime = t + 3600, minradius=min_dist, maxradius=max_dist,
                                channel=chan_type, level='channel', matchtimeseries=True)
# inventory = client.get_stations(longitude=ev_lon, latitude=ev_lat, starttime =
#                                 t, endtime = t + 3600, minradius=min_dist, maxradius=max_dist,
#                                 channel=chan_type, level='channel', matchtimeseries=True, network = 'IU')
if verbose == 1:
    print(inventory)
print('inventory has ' + str(len(inventory)) + ' networks recording data')

#%% Get waveforms, weed out synthetic waveforms
st = Stream()
cnt = 0
for network in inventory:
    if network.code != 'SY':  # Skip synthetic traces
        for station in network:
            cnt = cnt +1
            if cnt%10 == 0:
                print(str(cnt) + ' stations examined, ' + str(len(st)) + ' traces extracted')
            try:
                s_t = t - start_buff
                e_t = t + end_buff
                st += client.get_waveforms(network.code, station.code, location='*',channel=chan_type, starttime=s_t, endtime = e_t, attach_response=True)
            except:
                pass
print('st has data from ' + str(len(st)) + ' traces')
#%%  How many of which channels
tr = Trace()
hnz_chosen = 0
ehz_chosen = 0
hhz_chosen = 0
hlz_chosen = 0
bhz_chosen = 0
for tr in st:
    if tr.stats.channel == 'HNZ':
        hnz_chosen += 1
    if tr.stats.channel == 'EHZ':
        ehz_chosen += 1
    if tr.stats.channel == 'HHZ':
        hhz_chosen += 1
    if tr.stats.channel == 'HLZ':
        hlz_chosen += 1
    if tr.stats.channel == 'BHZ':
        bhz_chosen += 1

print('Total channels ' + str(len(st)) + ' - HNZ, EHZ, HHZ, HLZ, BHZ have '
       + str(hnz_chosen) + ' ' + str(ehz_chosen) + ' '
       + str(hhz_chosen) + ' ' + str(hlz_chosen) + ' ' + str(bhz_chosen))

#%%
st_cull = Stream()
for tr1 in st:# keep only a single location per station
        reject = 0  # keep only a single location per station
        # for tr2 in st_cull:
            # if ((tr2.stats.network  == tr1.stats.network) &
            #     (tr2.stats.station  == tr1.stats.station)):
            #     reject = 1
        tr1.detrend(type='simple')
        tr1.taper(0.05)
        tr1.filter('bandpass', freqmin=1, freqmax=4, corners=4, zerophase=False)
        if reject == 0:
            st_cull += tr1
print('st_cull has data from ' + str(len(st_cull)) + ' traces')
#%% Plot and write files
fname     = 'H' + etime[:10] + '.mseed'
fig = plt.figure(2)
plt.title(fname)
plt.xlabel('seconds')
plt.ylabel('stations')
st_cull.plot(equal_scale=False, fig=fig)

st_cull.write(       fname,     format = 'MSEED')
elapsed_time_wc = time.time() - start_time_wc
print(f'    Whole binning code took   {elapsed_time_wc:.1f}   seconds')
os.system('say "get event is done"')
