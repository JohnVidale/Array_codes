#!/usr/bin/env python
#from __future__ import print_function
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from obspy import Stream
import os

#%% Get catalog
chan_type = 'BHZ' # e.g., BHZ
start_buff = 10   # Pre-event buffer
end_buff   = 3600 # Time after OT recorded
min_mag = 4.5     # Threshold to find event in IRIS catalog
etime = 'XXXX-XX-XXTXX:XX:00' # event has to be unique and follow time by less than 1 minute
#etime = '2010-03-24T03:06:00' # event has to be unique and follow time by less than 1 minute
min_dist = 0      # Minimum distance loaded, degrees
max_dist = 180    # Maximum distance loaded, degrees
verbose = 0       # 1 prints networks and stations
refine = 0        # Use refined location

client = Client('IRIS')
t = UTCDateTime(etime)
catalog = client.get_events(starttime = t, endtime = t + 60, minmagnitude = min_mag)
print('IRIS epicenter:',catalog)

if refine == 1:
	# refined time and hypocenter
	ev_lon   = XXX.XXX  # South Sandwich Islands
	ev_lat   = XXX.XXX
	ev_depth = XX.X
	t        = UTCDateTime('XXXX-XX-XXTXX:XX:XX.XX')
#	t        = UTCDateTime('XXXX-XX-XXTXX:XX:XX.XX')


	#  Overwrite IRIS location in catalog
	catalog[0].origins[0].longitude = ev_lon
	catalog[0].origins[0].latitude = ev_lat
	catalog[0].origins[0].depth = ev_depth
	catalog[0].origins[0].time = t
	print('refined epicenter:',catalog)
else:
	#  Use IRIS location in catalog
	ev_lon   = catalog[0].origins[0].longitude
	ev_lat   = catalog[0].origins[0].latitude
	ev_depth = catalog[0].origins[0].depth
	t        = catalog[0].origins[0].time

#%% Make inventory of all stations recording this channel
inventory = client.get_stations(longitude=ev_lon, latitude=ev_lat, starttime =
								t, endtime = t + 3600, minradius=min_dist, maxradius=max_dist,
								channel=chan_type, level='channel', matchtimeseries=True)
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
			if cnt%100 == 0:
				print(str(cnt) + ' stations examined, ' + str(len(st)) + ' traces extracted')
			try:
				s_t = t - start_buff
				e_t = t + end_buff
				st += client.get_waveforms(network.code, station.code, location='*',channel=chan_type, starttime=s_t, endtime = e_t, attach_response=True)
			except:
				pass
print('st has data from ' + str(len(st)) + ' traces')
#%%
st_cull = Stream()
for tr1 in st:# keep only a single location per station
		reject = 0  # keep only a single location per station
		for tr2 in st_cull:
			if ((tr2.stats.network  == tr1.stats.network) &
			    (tr2.stats.station  == tr1.stats.station)):
				reject = 1
		if reject == 0:
			st_cull += tr1
print('st_cull has data from ' + str(len(st_cull)) + ' traces')
#%% Plot and write files
fname     = 'H' + etime[:10] + 'wvf_' + '.mseed'
fname_st  = 'H' + etime[:10] + 'sta_' + '.xml'
fname_cat = 'H' + etime[:10] + 'cat_' + '.xml'
#fig = plt.figure(2)
#plt.title(fname)
#plt.xlabel('seconds')
#plt.ylabel('stations')
#st_cull.plot(equal_scale=False, fig=fig)

st_cull.write(       fname,     format = 'MSEED')
inventory.write(fname_st,  format = 'STATIONXML')
catalog.write(  fname_cat, format = 'QUAKEML')
os.system('say "done"')