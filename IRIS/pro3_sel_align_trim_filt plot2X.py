#!/usr/bin/env python
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from obspy import Stream
from obspy import read_inventory
from obspy import read_events
from obspy import read
from obspy.geodetics import gps2dist_azimuth
import numpy as np
import os
from obspy.taup import TauPyModel
import matplotlib.pyplot as plt
model = TauPyModel(model='iasp91')
client = Client('IRIS')

verbose = 0           # more output
rel_time = 1          # timing is relative to a chosen phase, otherwise relative to OT
dphase  = 'P'     # phase to be aligned
dphase2 = 'no'     # another phase to have traveltime plotted
start_buff = 20       # plots start Xs before PKiKP
end_buff   = 20       # plots end Xs before PKiKP
chan_type = 'BHZ'     # was BHZ
taper_frac = .05      #Fraction of window tapered on both ends
signal_dur = 5.       # signal length used in SNR calculation
plot_scale_fac = 2    #  Bigger numbers make each trace amplitude bigger on plot
qual_threshold =  2 # minimum SNR
plot_tt = 1           # plot the traveltimes?
min_dist = 0
max_dist = 90
freq_min = 1
freq_max = 3

#%% Is taper too long compared to noise estimation window?
totalt = start_buff + end_buff
noise_time_skipped = taper_frac * totalt
if noise_time_skipped >= 0.5 * start_buff:
	print('Specified taper of ' + str(taper_frac * totalt) + ' is not big enough compared to available noise estimation window ' + str(start_buff - noise_time_skipped) + '. May not work well.')
	old_taper_frac = taper_frac
	taper_frac = 0.5*start_buff/totalt
	print('Taper reset from ' + str(old_taper_frac * totalt) + ' to ' + str(taper_frac * totalt) + ' seconds.')

if rel_time == 0: # SNR requirement not implemented for unaligned traces
	qual_threshold = 0

# Plot with reduced velocity?
red_plot = 0
red_dist = 55
red_time = 300
red_slow = 7.2 # seconds per degree

#%% Name and read input files
etime1 = 'XXXX-XX-XXTXX:XX:00' # event has to be unique and in specified 10 minutes
etime2 = 'XXXX-XX-XXTXX:XX:00' # event has to be unique and in specified 10 minutes
t1 = UTCDateTime(etime1)
t2 = UTCDateTime(etime2)
refine = 1

fname1     = 'HC' + etime1[:10] + 'wvf_' + '.mseed'
fname2     = 'HC' + etime2[:10] + 'wvf_' + '.mseed'
fname_st1  = 'HC' + etime1[:10] + 'sta_' + '.xml'
fname_st2  = 'HC' + etime2[:10] + 'sta_' + '.xml'
fname_cat1 = 'HC' + etime1[:10] + 'cat_' + '.xml'
fname_cat2 = 'HC' + etime2[:10] + 'cat_' + '.xml'

cat1 = read_events(fname_cat1)
cat2 = read_events(fname_cat2)
print('event1:',cat1)
print('event2:',cat2)

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

if refine == 1:  # Insert different hypocentral coordinates
	# relocation from Jiayuan

	ev_lon1   = XX.XXX
	ev_lat1   = XX.XXX
	ev_depth1 = XX.X
	t1        = UTCDateTime('XXXX-XX-XXTXX:XX:XX.XX')

	ev_lon2   = XX.XXX
	ev_lat2   = XX.XXX
	ev_depth2 = XX.X
	t2        = UTCDateTime('XXXX-XX-XXTXX:XX:XX.XX')

	#  Overwrite IRIS location in catalog
	cat1[0].origins[0].longitude = ev_lon1
	cat1[0].origins[0].latitude = ev_lat1
	cat1[0].origins[0].depth = ev_depth1
	cat1[0].origins[0].time = t1
	cat2[0].origins[0].longitude = ev_lon2
	cat2[0].origins[0].latitude = ev_lat2
	cat2[0].origins[0].depth = ev_depth2
	cat2[0].origins[0].time = t2

else:
	#  Use IRIS location in catalog
	ev_lon1   = cat1[0].origins[0].longitude
	ev_lat1   = cat1[0].origins[0].latitude
	ev_depth1 = cat1[0].origins[0].depth
	t1        = cat1[0].origins[0].time
	ev_lon2   = cat2[0].origins[0].longitude
	ev_lat2   = cat2[0].origins[0].latitude
	ev_depth2 = cat2[0].origins[0].depth
	t2        = cat2[0].origins[0].time

if verbose == 1:
	print('refined event:',cat1)
	print('refined event:',cat2)

#%%
# window and adjust start time to align picked times
st_pickalign1 = Stream()
st_pickalign2 = Stream()
for tr in st1: # traces one by one
	for network in inventory1:  # find lat-lon by searching entire inventory.  Inefficient
		for station in network:
			if (network.code == tr.stats.network) & (station.code == tr.stats.station): # find station in inventory
				stalon = station.longitude # look up lat & lon again to find distance
				stalat = station.latitude
				distance = gps2dist_azimuth(stalat,stalon,ev_lat1,ev_lon1) # Get traveltimes again, hard to store
				tr.stats.distance=distance[0] # distance in km
				dist = distance[0]/(1000*111)
				if min_dist < dist and dist < max_dist:
					ev_dep = ev_depth1
					try:
						arrivals = model.get_travel_times(source_depth_in_km=ev_dep,distance_in_degree=dist,phase_list=[dphase])
						atime = arrivals[0].time
						if rel_time == 1:
							s_t = t1 + atime - start_buff
							e_t = t1 + atime + end_buff
						else:
							s_t = t1 - start_buff
							e_t = t1 + end_buff
						tr.trim(starttime=s_t,endtime = e_t)
						# deduct theoretical traveltime and start_buf from starttime
						if rel_time == 1:
							tr.stats.starttime = tr.stats.starttime - atime
						st_pickalign1 += tr
					except:
						pass

for tr in st2: # traces one by one
	for network in inventory2:  # find lat-lon by searching entire inventory.  Inefficient
		for station in network:
			if (network.code == tr.stats.network) & (station.code == tr.stats.station): # find station in inventory
				stalon = station.longitude # look up lat & lon again to find distance
				stalat = station.latitude
				distance = gps2dist_azimuth(stalat,stalon,ev_lat2,ev_lon2) # Get traveltimes again, hard to store
				tr.stats.distance=distance[0] # distance in km
				dist = distance[0]/(1000*111)
				if min_dist < dist and dist < max_dist:
					ev_dep = ev_depth2
					try:
						arrivals = model.get_travel_times(source_depth_in_km=ev_dep,distance_in_degree=dist,phase_list=[dphase])
						atime = arrivals[0].time
						if rel_time == 1:
							s_t = t2 + atime - start_buff
							e_t = t2 + atime + end_buff
						else:
							s_t = t2 - start_buff
							e_t = t2 + end_buff
						tr.trim(starttime=s_t,endtime = e_t)
						# deduct theoretical traveltime and start_buf from starttime
						if rel_time == 1:
							tr.stats.starttime = tr.stats.starttime - atime
						st_pickalign2 += tr
					except:
						pass
print('After alignment and range selection - event 1: ' + str(len(st_pickalign1)) + ' traces')
print('After alignment and range selection - event 2: ' + str(len(st_pickalign2)) + ' traces')

#%%
#print(st) # at length
if verbose:
	print(st1.__str__(extended=True))
	print(st2.__str__(extended=True))
	if rel_time == 1:
		print(st_pickalign1.__str__(extended=True))
		print(st_pickalign2.__str__(extended=True))


#%%  detrend, taper, filter
st_pickalign1.detrend(type='simple')
st_pickalign2.detrend(type='simple')
st_pickalign1.taper(taper_frac)
st_pickalign2.taper(taper_frac)
st_pickalign1.filter('bandpass', freqmin=freq_min, freqmax=freq_max, corners=2, zerophase=True)
st_pickalign2.filter('bandpass', freqmin=freq_min, freqmax=freq_max, corners=2, zerophase=True)
st_pickalign1.taper(taper_frac)
st_pickalign2.taper(taper_frac)

#%%  Cull further by imposing SNR threshold on both traces
st1good = Stream()
st2good = Stream()
for tr1 in st_pickalign1:
	for tr2 in st_pickalign2:
		if ((tr1.stats.network  == tr2.stats.network) &
		    (tr1.stats.station  == tr2.stats.station)):
# estimate mean noise
			t_noise1_start  = int(len(tr1.data) * taper_frac)
			t_noise2_start  = int(len(tr2.data) * taper_frac)
			t_noise1_end    = int(len(tr1.data) * start_buff/(start_buff + end_buff))
			t_noise2_end    = int(len(tr2.data) * start_buff/(start_buff + end_buff))
			noise1          = np.mean(abs(tr1.data[t_noise1_start:t_noise1_end]))
			noise2          = np.mean(abs(tr2.data[t_noise2_start:t_noise2_end]))
# estimate mean signal
			t_signal1_start = int(len(tr1.data) * start_buff/(start_buff + end_buff))
			t_signal2_start = int(len(tr2.data) * start_buff/(start_buff + end_buff))
			t_signal1_end   = t_signal1_start + int(len(tr1.data) * signal_dur/(start_buff + end_buff))
			t_signal2_end   = t_signal2_start + int(len(tr2.data) * signal_dur/(start_buff + end_buff))
			signal1         = np.mean(abs(tr1.data[t_signal1_start:t_signal1_end]))
			signal2         = np.mean(abs(tr2.data[t_signal2_start:t_signal2_end]))
#			test SNR
			SNR1 = signal1/noise1;
			SNR2 = signal2/noise2;
			if (SNR1 > qual_threshold and SNR2 > qual_threshold):
				st1good += tr1
				st2good += tr2
#print(st1good)
#print(st2good)
print('Above SNR threshold - event 1: ' + str(len(st1good)) + ' traces')
print('Above SNR threshold - event 2: ' + str(len(st2good)) + ' traces')

#%%  get station lat-lon, compute distance for plot, maybe can comment out now?
for tr in st1good:
	for network in inventory1:  # find lat-lon by searching entire inventory.  Inefficient
		for station in network:
			if (network.code == tr.stats.network) & (station.code == tr.stats.station):
				stalon = station.longitude
				stalat = station.latitude
				distance = gps2dist_azimuth(stalat,stalon,ev_lat1,ev_lon1)
				tr.stats.distance=distance[0] # distance in km

for tr in st2good:
	for network in inventory2:  # find lat-lon by searching entire inventory.  Inefficient
		for station in network:
			if (network.code == tr.stats.network) & (station.code == tr.stats.station):
				stalon = station.longitude
				stalat = station.latitude
				distance = gps2dist_azimuth(stalat,stalon,ev_lat2,ev_lon2)
				tr.stats.distance=distance[0] # distance in km

#%%
# plot traces
plt.figure(7)
for tr in st1good:
	dist_offset = tr.stats.distance/(1000*111) # trying for approx degrees
	time = np.arange(len(tr.data)) * tr.stats.delta + (tr.stats.starttime - t1)
	if red_plot == 1:
		shift = red_time + (dist_offset - red_dist) * red_slow
		time = time - shift
	plt.plot(time, (tr.data - np.mean(tr.data))*plot_scale_fac /(tr.data.max()
		- tr.data.min()) + dist_offset, color = 'green')
#plt.xlabel('Time (s)')
#plt.ylabel('Epicentral distance from event (°)')
#plt.title(fname1)
plt.show()

for tr in st2good:
	dist_offset = tr.stats.distance/(1000*111) # trying for approx degrees
	time = np.arange(len(tr.data)) * tr.stats.delta + (tr.stats.starttime - t2)
	if red_plot == 1:
		shift = red_time + (dist_offset - red_dist) * red_slow
		time = time - shift
	time = time # arbitrary eyeballed shift
#	time = time + 2.0 # arbitrary eyeballed shift
	plt.plot(time, (tr.data - np.mean(tr.data))*plot_scale_fac /(tr.data.max()
		- tr.data.min()) + dist_offset, color = 'red')

	#%% Plot traveltime curves
if plot_tt:
	# first traveltime curve
	line_pts = 50
	dist_vec  = np.arange(min_dist, max_dist, (max_dist - min_dist)/line_pts) # distance grid
	time_vec1 = np.arange(min_dist, max_dist, (max_dist - min_dist)/line_pts) # empty time grid of same length (filled with -1000)
	for i in range(0,line_pts):
		arrivals = model.get_travel_times(source_depth_in_km=ev_dep,distance_in_degree
									=dist_vec[i],phase_list=[dphase])
		num_arrivals = len(arrivals)
		found_it = 0
		for j in range(0,num_arrivals):
			if arrivals[j].name == dphase:
				time_vec1[i] = arrivals[j].time
				found_it = 1
		if found_it == 0:
			time_vec1[i] = np.nan
	# second traveltime curve
	if dphase2 != 'no':
		time_vec2 = np.arange(min_dist, max_dist, (max_dist - min_dist)/line_pts) # empty time grid of same length (filled with -1000)
		for i in range(0,line_pts):
			arrivals = model.get_travel_times(source_depth_in_km=ev_dep,distance_in_degree
										=dist_vec[i],phase_list=[dphase2])
			num_arrivals = len(arrivals)
			found_it = 0
			for j in range(0,num_arrivals):
				if arrivals[j].name == dphase2:
					time_vec2[i] = arrivals[j].time
					found_it = 1
			if found_it == 0:
				time_vec2[i] = np.nan

	if rel_time == 1:
		if dphase2 != 'no':
			time_vec2 = time_vec2 - time_vec1
		time_vec1 = time_vec1 - time_vec1
	plt.plot(time_vec1,dist_vec, color = 'blue')
	if dphase2 != 'no':
		plt.plot(time_vec2,dist_vec, color = 'purple')

plt.xlabel('Time (s)')
plt.ylabel('Epicentral distance from event (°)')
plt.title(dphase + ' for ' + fname1[2:12] + ' vs ' + fname2[2:12])
plt.show()
os.system('say "Done"')

#  Save processed files
#fname3 = 'A' + etime1[:10] + 'pro_' + dphase + '.mseed'
#fname4 = 'A' + etime2[:10] + 'pro_' + dphase + '.mseed'
#st1good.write(fname3,format = 'MSEED')
#st2good.write(fname4,format = 'MSEED')
#inventory.write(fname_st,format = 'STATIONXML')

##  Save culled files
#fname1     = 'HC' + etime1[:10] + 'wvf_' + '.mseed'
#fname2     = 'HC' + etime2[:10] + 'wvf_' + '.mseed'
#fname_st1  = 'HC' + etime1[:10] + 'sta_' + '.xml'
#fname_st2  = 'HC' + etime2[:10] + 'sta_' + '.xml'
#fname_cat1 = 'HC' + etime1[:10] + 'cat_' + '.xml'
#fname_cat2 = 'HC' + etime2[:10] + 'cat_' + '.xml'
#
#st1cull.write(   fname1,format = 'MSEED')
#st2cull.write(   fname2,format = 'MSEED')
## inventory and sta_list not updated, no point to re-writing it, fix to shrink lists
#inventory1.write(fname_st1,  format = 'STATIONXML')
#inventory2.write(fname_st2,  format = 'STATIONXML')
#cat1.write(      fname_cat1, format = 'QUAKEML')
#cat2.write(      fname_cat2, format = 'QUAKEML')
#
#os.system('say "Done"')