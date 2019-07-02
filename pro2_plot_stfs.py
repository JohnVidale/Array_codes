#!/usr/bin/env python
# this program only detrends, tapers and decimates
# John Vidale 2/2019

def pro2_test(eq_file1, conv_file1, eq_file2, conv_file2):

	from obspy import UTCDateTime
	from obspy import Stream, Trace
	from obspy import read
	import os
	import time
	import numpy as np
	import matplotlib.pyplot as plt

	import sys # don't show any warnings
	import warnings

	if not sys.warnoptions:
	    warnings.simplefilter('ignore')

	start_time_wc = time.time()

	#%%
	taper_frac = .05      #Fraction of window tapered on both ends
#	conv_file = 'HD1971-11-06_stf'
	'''
	#%% input event data with 1-line file of format
	#  event 2016-05-28T09:47:00.000 -56.241 -26.935 78
	file1 = open('EvLocs/' + eq_file1, 'r')

	lines1=file1.readlines()
	split_line1 = lines1[0].split()
#			ids.append(split_line[0])  ignore label for now
	t1           = UTCDateTime(split_line1[1])
	date_label1  = split_line1[1][0:10]
	print('date_label ' + date_label1 + ' time ' + str(t1))

	#%% Load waveforms and convolution trace
	st1        = Stream()
	con_trace1 = Stream()
	st_out1    = Stream()
	tr1 = Trace()

	fname_sel1     = 'Pro_Files/HD' + date_label1 + 'sel.mseed'

	st1        = read(fname_sel1)
	fname1     = conv_file1
	con_trace1 = read(fname1)

	nt1 = len(st1[0].data)
	dt1 = st1[0].stats.delta
	print('Read in:\n' + str(len(st1)) + ' traces' + ' from file ' + fname1 +
	   ', ' + str(nt1) + ' time pts, time sampling of '
		  + str(dt1) + ' and thus duration of ' + str((nt1-1)*dt1))

		#%%  detrend, taper, decimate
	st1.detrend(type='simple')
	st1.taper(taper_frac)
	fig_index = 2
	con_trace1.plot()
	'''
	#%% input event data with 1-line file of format
	#  event 2016-05-28T09:47:00.000 -56.241 -26.935 78
	file2 = open('EvLocs/' + eq_file2, 'r')

	lines2=file2.readlines()
	split_line2 = lines2[0].split()
#			ids.append(split_line[0])  ignore label for now
	t2           = UTCDateTime(split_line2[1])
	date_label2  = split_line2[1][0:10]
	print('date_label ' + date_label2 + ' time ' + str(t2))

	#%% Load waveforms and convolution trace
	st2        = Stream()
	con_trace2 = Stream()
	st_out2    = Stream()
	tr2 = Trace()

	fname_sel2     = 'Pro_Files/HD' + date_label2 + 'sel.mseed'

	st2        = read(fname_sel2)
	fname2     = conv_file2
	con_trace2 = read(fname2)

	nt2 = len(st2[0].data)
	dt2 = st2[0].stats.delta
	print('Read in:\n' + str(len(st2)) + ' traces' + ' from file ' + fname2 +
	   ', ' + str(nt2) + ' time pts, time sampling of '
		  + str(dt2) + ' and thus duration of ' + str((nt2-1)*dt2))

		#%%  detrend, taper, decimate
	st2.detrend(type='simple')
	st2.taper(taper_frac)
	fig_index = 3
	con_trace2.plot()

	'''
	done = 0

	for tr in st: # traces one by one, find lat-lon by searching entire inventory.  Inefficient but cheap
		print('con_trace data has length ' + str(len(con_trace[0].data)))
		print('Tr data has length ' + str(len(tr.data)) + 'con_trace data has length ' + str(len(con_trace[0].data)))
		tr.data = np.convolve(tr.data, con_trace[0].data)
		print('Now, Tr data has length ' + str(len(tr.data)))
		st_out += tr
		done += 1
		if done%50 == 0:
			print('Done stacking ' + str(done) + ' out of ' + str(len(st)) + ' stations.')

	nt = len(st_out[0].data)
	dt = st_out[0].stats.delta
	print('After decimation:\n' + str(len(st_out)) + ' traces written to file ' + fname_sel +
	   ', ' + str(nt) + ' time pts, time sampling of '
		  + str(dt) + ' and thus duration of ' + str((nt-1)*dt))

	#  Save processed files
	st_out.write(fname_sel,format = 'MSEED')
	'''

	elapsed_time_wc = time.time() - start_time_wc
	print('This job took ' + str(elapsed_time_wc) + ' seconds')
	os.system('say "Done"')