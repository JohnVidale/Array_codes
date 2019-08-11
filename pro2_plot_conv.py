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

#	import sys # don't show any warnings
#	import warnings
#
#	if not sys.warnoptions:
#	    warnings.simplefilter('ignore')
#
	start_time_wc = time.time()

	taper_frac = .5      #Fraction of window tapered on both ends

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
	tr.normalize()

	tr.stats.delta = 0.1
	tr.trim(starttime = tr.stats.starttime + 10, endtime = tr.stats.endtime - 10)
	con_trace1[0].stats.starttime = tr.stats.starttime
	con_trace2[0].stats.starttime = tr.stats.starttime
	print(str(tr.stats.delta) + ' ' + str(tr.stats.starttime) + ' ' + str(tr.stats.endtime))
	print(str(con_trace1[0].stats.delta) + ' ' + str(con_trace1[0].stats.starttime) + ' ' + str(con_trace1[0].stats.endtime))
	print(str(con_trace2[0].stats.delta) + ' ' + str(con_trace2[0].stats.starttime) + ' ' + str(con_trace2[0].stats.endtime))
	con_trace1[0].stats.channel = 'trace1'
	con_trace2[0].stats.channel = 'trace2'
	tr.stats.channel = 'convolved'
	con_trace1.append(con_trace2[0])
	con_trace1.append(tr)
	print('length of con_trace1 is ' + str(len(con_trace1)))

	sgrams = Stream()
	sgrams += con_trace1[0]
	sgrams += con_trace1[1]
	sgrams += con_trace1[2]
#	con_trace1 += tr
#	con_trace1.plot(type = 'relative')
#	con_trace2.plot(type = 'relative')
	sgrams[0].stats.station = ''
	sgrams[1].stats.station = ''
	sgrams[2].stats.station = ''
	sgrams[0].stats.channel = '1971'
	sgrams[1].stats.channel = '1969'
	sgrams[2].stats.channel = 'convolved'
	sgrams.filter('bandpass', freqmin=1, freqmax=2, corners=4, zerophase=True)
	sgrams.normalize()


	sgrams.plot(size = (800,600))

	elapsed_time_wc = time.time() - start_time_wc
	print('This job took ' + str(elapsed_time_wc) + ' seconds')
	os.system('say "Done"')