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
dphase  = 'PKiKP'     # phase to be aligned
dphase2 = 'PKIKP'     # another phase to be plotted
start_buff = 20       # plots start Xs before PKiKP
end_buff   = 20       # plots end Xs before PKiKP
chan_type = 'BHZ'     # was BHZ
taper_frac = .05      #Fraction of window tapered on both ends
signal_dur = 5.       # signal length used in SNR calculation
plot_scale_fac = 2    #  Bigger numbers make each trace amplitude bigger on plot
qual_threshold =  1.5 # minimum SNR
plot_tt = 1           # plot the traveltimes?
min_dist = 120
max_dist = 180
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
etime = 'XXXX-XX-XXTXX:XX:00' # event has to be unique and in specified 10 minutes
#etime = 'XXXX-XX-XXTXX:XX:00' # event has to be unique and in specified 10 minutes
t = UTCDateTime(etime)
refine = 1

fname     = 'HC' + etime[:10] + 'wvf_' + '.mseed'
fname_st  = 'HC' + etime[:10] + 'sta_' + '.xml'
fname_cat = 'HC' + etime[:10] + 'cat_' + '.xml'

cat = read_events(fname_cat)
print('event:',cat)

#%% Reload station inventory and waveforms
inventory = read_inventory(fname_st)
st = Stream()
st=read(fname)
print('st has ' + str(len(st)) + ' traces')
if verbose == 1:
	print(inventory)
	print(st)

if refine == 1:  # Insert different hypocentral coordinates

	ev_lon   = XX.XXX
	ev_lat   = XX.XXX
	ev_depth = XX.X
	t        = UTCDateTime('XXXX-XX-XXTXX:XX:XX.XX')

#	ev_lon   = XX.XXX
#	ev_lat   = XX.XXX
#	ev_depth = XX.X
#	t        = UTCDateTime('XXXX-XX-XXTXX:XX:XX.XX')

	#  Overwrite IRIS location in catalog
	cat[0].origins[0].longitude = ev_lon
	cat[0].origins[0].latitude = ev_lat
	cat[0].origins[0].depth = ev_depth
	cat[0].origins[0].time = t

else:
	#  Use IRIS location in catalog
	ev_lon   = cat[0].origins[0].longitude
	ev_lat   = cat[0].origins[0].latitude
	ev_depth = cat[0].origins[0].depth
	t        = cat[0].origins[0].time

if verbose == 1:
	print('refined event:',cat)

#%%
# window and adjust start time to align picked times
st_pickalign = Stream()
for tr in st: # traces one by one
	for network in inventory:  # find lat-lon by searching entire inventory.  Inefficient
		for station in network:
			if (network.code == tr.stats.network) & (station.code == tr.stats.station): # find station in inventory
				stalon = station.longitude # look up lat & lon again to find distance
				stalat = station.latitude
				distance = gps2dist_azimuth(stalat,stalon,ev_lat,ev_lon) # Get traveltimes again, hard to store
				tr.stats.distance=distance[0] # distance in km
				dist = distance[0]/(1000*111)
				if min_dist < dist and dist < max_dist:
					ev_dep = ev_depth
					try:
						arrivals = model.get_travel_times(source_depth_in_km=ev_dep,distance_in_degree=dist,phase_list=[dphase])
						atime = arrivals[0].time
						if rel_time == 1:
							s_t = t + atime - start_buff
							e_t = t + atime + end_buff
						else:
							s_t = t - start_buff
							e_t = t + end_buff
						tr.trim(starttime=s_t,endtime = e_t)
						# deduct theoretical traveltime and start_buf from starttime
						if rel_time == 1:
							tr.stats.starttime = tr.stats.starttime - atime
						st_pickalign += tr
					except:
						pass
print('After alignment and range selection - event: ' + str(len(st_pickalign)) + ' traces')

#%%
#print(st) # at length
if verbose:
	print(st.__str__(extended=True))
	if rel_time == 1:
		print(st_pickalign.__str__(extended=True))


#%%  detrend, taper, filter
st_pickalign.detrend(type='simple')
st_pickalign.taper(taper_frac)
st_pickalign.filter('bandpass', freqmin=freq_min, freqmax=freq_max, corners=2, zerophase=True)
st_pickalign.taper(taper_frac)

#%%  Cull further by imposing SNR threshold on both traces
stgood = Stream()
for tr in st_pickalign:
# estimate mean noise
	t_noise_start  = int(len(tr.data) * taper_frac)
	t_noise_end    = int(len(tr.data) * start_buff/(start_buff + end_buff))
	noise          = np.mean(abs(tr.data[t_noise_start:t_noise_end]))
# estimate mean signal
	t_signal_start = int(len(tr.data) * start_buff/(start_buff + end_buff))
	t_signal_end   = t_signal_start + int(len(tr.data) * signal_dur/(start_buff + end_buff))
	signal         = np.mean(abs(tr.data[t_signal_start:t_signal_end]))
#			test SNR
	SNR = signal/noise;
	if (SNR > qual_threshold):
		stgood += tr
#print(st1good)
#print(st2good)
print('Above SNR threshold - event: ' + str(len(stgood)) + ' traces')

#%%  get station lat-lon, compute distance for plot
for tr in stgood:
	for network in inventory:  # find lat-lon by searching entire inventory.  Inefficient
		for station in network:
			if (network.code == tr.stats.network) & (station.code == tr.stats.station):
				stalon = station.longitude
				stalat = station.latitude
				distance = gps2dist_azimuth(stalat,stalon,ev_lat,ev_lon)
				tr.stats.distance=distance[0] # distance in km

#%%
# plot traces
plt.figure(8)
for tr in stgood:
	dist_offset = tr.stats.distance/(1000*111) # trying for approx degrees
	time = np.arange(len(tr.data)) * tr.stats.delta + (tr.stats.starttime - t)
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
plt.title(dphase + ' for ' + fname[2:12])
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