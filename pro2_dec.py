#!/usr/bin/env python
# this program only detrends, tapers and decimates
# John Vidale 2/2019

def pro2_decimate(eq_file, decimate_fac = 5):

	from obspy import UTCDateTime
	from obspy import Stream
	from obspy import read
	import os
	import time

	import sys # don't show any warnings
	import warnings

	if not sys.warnoptions:
	    warnings.simplefilter("ignore")

	start_time_wc = time.time()

	#%%
	taper_frac = .05      #Fraction of window tapered on both ends
#	decimate_fac = 10         # 0 if no decimation desired

	#%% input event data with 1-line file of format
	#  event 2016-05-28T09:47:00.000 -56.241 -26.935 78
	file = open(eq_file, 'r')
	lines=file.readlines()
	split_line = lines[0].split()
#			ids.append(split_line[0])  ignore label for now
	t           = UTCDateTime(split_line[1])
	date_label  = split_line[1][0:10]
	print('date_label ' + date_label + ' time ' + str(t))

	#%% Load waveforms
	st = Stream()
	fname = 'HiNet' + date_label + '_wvf' + '.mseed'
	fname_out = 'HD' + date_label + '.mseed'
	st=read(fname)
	nt = len(st[0].data)
	dt = st[0].stats.delta
	print('Read in:\n' + str(len(st)) + ' traces' + ' from file ' + fname +
	   ', ' + str(nt) + ' time pts, time sampling of '
		  + str(dt) + ' and thus duration of ' + str((nt-1)*dt))

		#%%  detrend, taper, decimate
	st.detrend(type='simple')
	st.taper(taper_frac)
	st.decimate(decimate_fac)

	nt = len(st[0].data)
	dt = st[0].stats.delta
	print('After decimation:\n' + str(len(st)) + ' traces written to file ' + fname_out +
	   ', ' + str(nt) + ' time pts, time sampling of '
		  + str(dt) + ' and thus duration of ' + str((nt-1)*dt))

	#  Save processed files
	st.write(fname_out,format = 'MSEED')

	elapsed_time_wc = time.time() - start_time_wc
	print('This job took ' + str(elapsed_time_wc) + ' seconds')

	os.system('say "Done"')