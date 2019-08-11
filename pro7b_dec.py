#!/usr/bin/env python
# decimates output of pro6 for quickly read into pro7a
# John Vidale 3/2019

def pro7dec(eq_file1, eq_file2, decimate_fac = 5, ARRAY = 0):

	from obspy import Stream
	from obspy import read
	import os
	import time

	start_time_wc = time.time()

	if ARRAY == 0:
		file = open(eq_file1, 'r')
	elif ARRAY == 1:
		file = open('EvLocs/' + eq_file1, 'r')
	lines=file.readlines()
	split_line = lines[0].split()
#			ids.append(split_line[0])  ignore label for now
	date_label1  = split_line[1][0:10]

	if ARRAY == 0:
		file = open(eq_file2, 'r')
	elif ARRAY == 1:
		file = open('EvLocs/' + eq_file2, 'r')
	lines=file.readlines()
	split_line = lines[0].split()
#			ids.append(split_line[0])  ignore label for now
	date_label2  = split_line[1][0:10]

	#%% Input parameters
	# #%% Get saved event info, also used to name files
	# date_label = '2018-04-02' # date for filename
	if ARRAY == 0:
		fname1 = 'HD' + date_label1 + '_' + date_label2 + '_tshift.mseed'
		fname2 = 'HD' + date_label1 + '_' + date_label2 + '_amp_ave.mseed'
		fname3 = 'HD' + date_label1 + '_' + date_label2 + '_amp_ratio.mseed'
	elif ARRAY == 1:
		fname1 = 'Pro_files/HD' + date_label1 + '_' + date_label2 + '_tshift.mseed'
		fname2 = 'Pro_files/HD' + date_label1 + '_' + date_label2 + '_amp_ave.mseed'
		fname3 = 'Pro_files/HD' + date_label1 + '_' + date_label2 + '_amp_ratio.mseed'
	st        = Stream()
	amp_ave   = Stream()
	amp_ratio = Stream()
	st        = read(fname1)
	amp_ave   = read(fname2)
	amp_ratio = read(fname3)
	print('Read in: ' + str(len(st)) + '  ' + str(len(amp_ave)) + '  ' + str(len(amp_ratio)) + ' traces for st, amp_ave, amp_ratio, time sampling is '
	   + str(st[0].stats.delta), ' number of time points is ' + str(len(st[0].data)))

	for i in range(len(st)):  # loop over traces
		st[i].decimate(decimate_fac)
		amp_ave[i].decimate(decimate_fac)
		amp_ratio[i].decimate(decimate_fac)

	elapsed_time_wc = time.time() - start_time_wc
	print('Decimation took ' + str(elapsed_time_wc) + ' seconds')

	#  Save processed files
	if ARRAY == 0:
		fname = 'HD' + date_label1 + '_' + date_label2 + '_tshift_dec.mseed'
	elif ARRAY == 1:
		fname = 'Pro_Files/HD' + date_label1 + '_' + date_label2 + '_tshift_dec.mseed'
	st.write(fname,format = 'MSEED')

	if ARRAY == 0:
		fname = 'HD' + date_label1 + '_' + date_label2 + '_amp_ave_dec.mseed'
	elif ARRAY == 1:
		fname = 'Pro_Files/HD' + date_label1 + '_' + date_label2 + '_amp_ave_dec.mseed'
	amp_ave.write(fname,format = 'MSEED')

	if ARRAY == 0:
		fname = 'HD' + date_label1 + '_' + date_label2 + '_amp_ratio_dec.mseed'
	elif ARRAY == 1:
		fname = 'Pro_Files/HD' + date_label1 + '_' + date_label2 + '_amp_ratio_dec.mseed'
	amp_ratio.write(fname,format = 'MSEED')

	print('Wrote out: time sampling is ' + str(st[0].stats.delta) + ' number of time points is ' + str(len(st[0].data)))

	elapsed_time_wc = time.time() - start_time_wc
	print('This job took ' + str(elapsed_time_wc) + ' seconds')
	os.system('say "Done"')