#!/usr/bin/env python
# this program only detrends, tapers and decimates
# John Vidale 2/2019

def pro2_test(eq_file1, conv_file1, eq_file2, conv_file2):

	from obspy import UTCDateTime
	from obspy import Stream, Trace
	from obspy import read
#	from obspy.signal import correlate_template
	import obspy
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
	taper_frac = .5      #Fraction of window tapered on both ends
#	conv_file = 'HD1971-11-06_stf'

	#%% Load waveforms and convolution trace
#	con_trace1 = Stream()
#	fname1     = conv_file1
#	con_trace1 = read(fname1)
#
#	con_trace2 = Stream()
#	fname2     = conv_file2
#	con_trace2 = read(fname2)
#
#	con_trace2.append(con_trace1[0])
#	con_trace2.plot(typ = 'relative')

	con_trace1 = Stream()
	con_trace2 = Stream()
	tr = Trace()
	con_trace1 = read(conv_file1)
	con_trace2 = read(conv_file2)
	con_trace1.taper(taper_frac)
	con_trace2.taper(taper_frac)
	con_trace1.normalize()
	con_trace2.normalize()
	tr.data = np.convolve(con_trace2[0].data, con_trace1[0].data)
	tr.stats.delta = 0.1
	tr.plot(type = 'relative')

#	print('got to print statement')
#	os.system('say "made it to here"')
#	os.system('say "and here"')

#	st3 = Stream()
#	st3 = correlate_template(st, con_trace, normalize='none')
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