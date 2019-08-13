#!/usr/bin/env python
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from obspy import Stream
from obspy import read_inventory
from obspy import read_events
from obspy import read
import os

verbose = 0  # more output
#%% Set parameters
etime1 = 'XXXX-XX-XXTXX:XX:00' # event has to be unique and follow time by less than 10 minutes
etime2 = 'XXXX-XX-XXTXX:XX:00' # event has to be unique and follow time by less than 10 minutes
min_mag = 4.5

#%% Read in corrected catalog for accurate OT, lat and lon for distance
fname1     = 'H' + etime1[:10] + 'wvf_' + '.mseed'
fname2     = 'H' + etime2[:10] + 'wvf_' + '.mseed'
fname_st1  = 'H' + etime1[:10] + 'sta_' + '.xml'
fname_st2  = 'H' + etime2[:10] + 'sta_' + '.xml'
fname_cat1 = 'H' + etime1[:10] + 'cat_' + '.xml'
fname_cat2 = 'H' + etime2[:10] + 'cat_' + '.xml'

cat1 = read_events(fname_cat1)
cat2 = read_events(fname_cat2)
print('event1:',cat1)
print('event2:',cat2)
# catalog.plot()

#%% Reload station inventory and waveforms
inventory1 = read_inventory(fname_st1)
inventory2 = read_inventory(fname_st2)
st1 = Stream()
st2 = Stream()
st1=read(fname1)
print('st1 has ' + str(len(st1)) + ' traces')
st2=read(fname2)
print('st2 has ' + str(len(st2)) + ' traces')
if verbose == 1:
	print(inventory1)
	print(inventory2)
	print(st1)
	print(st2)

#%%  find stations recording both events, make common dataset
st1cull = Stream()
st2cull = Stream()
for tr1 in st1:
	for tr2 in st2:
		if ((tr1.stats.network  == tr2.stats.network) &
		    (tr1.stats.station  == tr2.stats.station) &
		    (tr1.stats.location == tr2.stats.location)):
			reject = 0  # keep only a single location per station
			for tr3 in st1cull:
				if ((tr3.stats.network  == tr1.stats.network) &
				    (tr3.stats.station  == tr1.stats.station)):
					reject = 1
			if reject == 0:
				st1cull += tr1
				st2cull += tr2
print('st1cull has ' + str(len(st1cull)) + ' traces')
print('st2cull has ' + str(len(st2cull)) + ' traces')

#  Save culled files
fname1     = 'HC' + etime1[:10] + 'wvf_' + '.mseed'
fname2     = 'HC' + etime2[:10] + 'wvf_' + '.mseed'
fname_st1  = 'HC' + etime1[:10] + 'sta_' + '.xml'
fname_st2  = 'HC' + etime2[:10] + 'sta_' + '.xml'
fname_cat1 = 'HC' + etime1[:10] + 'cat_' + '.xml'
fname_cat2 = 'HC' + etime2[:10] + 'cat_' + '.xml'
st1cull.write(   fname1,     format = 'MSEED')
st2cull.write(   fname2,     format = 'MSEED')
inventory1.write(fname_st1,  format = 'STATIONXML')
inventory2.write(fname_st2,  format = 'STATIONXML')
cat1.write(      fname_cat1, format = 'QUAKEML')
cat2.write(      fname_cat2, format = 'QUAKEML')

os.system('say "Done"')