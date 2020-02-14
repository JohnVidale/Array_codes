#!/usr/bin/env python
# looks at envelope stack differential parameters from pro6
# Reads in tdiff, ave_amp, amp_ratio computed from a pair of events
# window by quality signals
# plots snapshots at a range of lag times
# plus plots sections at 4 radial slownesses (0.0, 0.05, 0.01, 0.015)
# plus plots 1 transverse slowness at 0 slowness
# John Vidale 2/2019

def pro7plotstack3(eq_file1, eq_file2, plot_scale_fac = 0.05, slow_delta = 0.0005,
			  slowR_lo = -0.1, slowR_hi = 0.1, slowT_lo = -0.1, slowT_hi = 0.1,
			  start_buff = -50, end_buff = 50,
			  R_slow_plot = 0.06, T_slow_plot = 0.0,
			  snaptime = 8, snaps = 10, tdiff_clip = -1,
			  plot_dyn_range = 1000, fig_index = 401, skip_T = 0, skip_R = 0, skip_snaps = 0,
			  decimate_fac = 0, in_dec = 0, ref_phase = 'blank', ARRAY = 0,
			  min_rat = 0.6, max_rat = 1.8, min_amp = 0.3):

	from obspy import read
	import numpy as np
	import os
	import matplotlib.pyplot as plt
	import time
	from obspy import UTCDateTime
	from obspy import Stream

	print('Running pro7b_plot_stack')

	start_time_wc = time.time()

	if ARRAY == 0:
		file = open(eq_file1, 'r')
	elif ARRAY == 1:
		file = open('EvLocs/' + eq_file1, 'r')
	lines=file.readlines()
	split_line = lines[0].split()
#			ids.append(split_line[0])  ignore label for now
	t1           = UTCDateTime(split_line[1])
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
		if in_dec == 0:
			fname1 = 'HD' + date_label1 + '_' + date_label2 + '_tshift.mseed'
			fname2 = 'HD' + date_label1 + '_' + date_label2 + '_amp_ave.mseed'
			fname3 = 'HD' + date_label1 + '_' + date_label2 + '_amp_ratio.mseed'
		else:
			fname1 = 'HD' + date_label1 + '_' + date_label2 + '_tshift_dec.mseed'
			fname2 = 'HD' + date_label1 + '_' + date_label2 + '_amp_ave_dec.mseed'
			fname3 = 'HD' + date_label1 + '_' + date_label2 + '_amp_ratio_dec.mseed'
	elif ARRAY == 1:
		if in_dec == 0:
			fname1 = 'Pro_files/HD' + date_label1 + '_' + date_label2 + '_tshift.mseed'
			fname2 = 'Pro_files/HD' + date_label1 + '_' + date_label2 + '_amp_ave.mseed'
			fname3 = 'Pro_files/HD' + date_label1 + '_' + date_label2 + '_amp_ratio.mseed'
		else:
			fname1 = 'Pro_files/HD' + date_label1 + '_' + date_label2 + '_tshift_dec.mseed'
			fname2 = 'Pro_files/HD' + date_label1 + '_' + date_label2 + '_amp_ave_dec.mseed'
			fname3 = 'Pro_files/HD' + date_label1 + '_' + date_label2 + '_amp_ratio_dec.mseed'

	tdiff     = Stream()
	amp_ave   = Stream()
	amp_ratio = Stream()
	tdiff     = read(fname1)
	amp_ave   = read(fname2)
	amp_ratio = read(fname3)
	print('Read in: ' + str(len(tdiff)) + '  ' + str(len(amp_ave)) + '  ' + str(len(amp_ratio)) + ' traces for tdiff, amp_ave, amp_ratio')

	for i in range(len(tdiff)):  # loop over traces, probably already decimated in previous step, pro7dec
			if decimate_fac != 0:
				tdiff[i].decimate(decimate_fac)
				amp_ave[i].decimate(decimate_fac)
				amp_ratio[i].decimate(decimate_fac)

	elapsed_time_wc = time.time() - start_time_wc
	print('Decimation took ' + str(elapsed_time_wc) + ' seconds')

	nt = len(tdiff[0].data)
	dt = tdiff[0].stats.delta
	print('After decimation, first trace has : ' + str(nt) + ' time pts, time sampling of '
		  + str(dt) + ' and thus duration of ' + str((nt-1)*dt))

	#%% Make grid of slownesses for
	slowR_n = int(1 + (slowR_hi - slowR_lo)/slow_delta)  # number of radial slownesses
	slowT_n = int(1 + (slowT_hi - slowT_lo)/slow_delta)  # number of tranverse slownesses
	stack_nt = int(1 + ((end_buff - start_buff)/dt))  # number of time points
	print(str(slowT_n) + ' trans slownesses, hi and lo are ' + str(slowT_hi) + '  ' + str(slowT_lo))
	# In English, stack_slows = range(slow_n) * slow_delta - slow_lo
	a1R = range(slowR_n)
	a1T = range(slowT_n)
	stack_Rslows = [(x * slow_delta + slowR_lo) for x in a1R]
	stack_Tslows = [(x * slow_delta + slowT_lo) for x in a1T]
	print(str(slowR_n) + ' radial slownesses, hi and lo are ' + str(slowR_hi) + '  ' + str(slowR_lo))
	print('Input trace starttime ' + str(tdiff[0].stats.starttime))

	#%% mask out bad points
	global_max = 0
	for slow_i in range(len(amp_ave)): # find global max of ave_amp
		local_max = max(amp_ave[slow_i].data)
		if local_max > global_max:
			global_max = local_max

	for slow_i in range(len(tdiff)): # ignore less robust points
		for it in range(nt):
			if ((amp_ratio[slow_i].data[it] < min_rat) or (amp_ratio[slow_i].data[it] > max_rat) or (amp_ave[slow_i].data[it] < (min_amp * global_max))):
				tdiff[slow_i].data[it] = np.nan

	#%% Slice near transverse slowness T_slow_plot
	if skip_R != 1:
		lowest_Tslow = 1000000
		for slow_i in range(slowT_n):
			if abs(stack_Tslows[slow_i] - T_slow_plot) < lowest_Tslow:
				lowest_Tindex = slow_i
				lowest_Tslow = abs(stack_Tslows[slow_i] - T_slow_plot)

		print(f'{lowest_Tindex:4d} is T slow nearest {T_slow_plot:.3f}, difference is {lowest_Tslow:.3f}')

		# Select only stacks with that slowness for Transverse plot
		centralR_st = Stream()
		for slowR_i in range(slowR_n):
			centralR_st += tdiff[slowR_i*slowT_n + lowest_Tindex]

	#%% Slice near radial slowness R_slow_plot
	if skip_T != 1:
		lowest_Rslow = 1000000
		for slow_i in range(slowR_n):
			if abs(stack_Rslows[slow_i] - R_slow_plot) < lowest_Rslow:
				lowest_Rindex = slow_i
				lowest_Rslow = abs(stack_Rslows[slow_i] - R_slow_plot)

		print(f'{lowest_Rindex:4d} is R slow nearest 0.005, difference is {lowest_Rslow:.3f}')

		# Select only stacks with that slowness for Radial plot
		centralT_st = Stream()
		for slowT_i in range(slowT_n):
			centralT_st += tdiff[lowest_Rindex*slowT_n + slowT_i]

#%%
	#%% compute timing time series
	ttt = (np.arange(len(tdiff[0].data)) * tdiff[0].stats.delta + start_buff) # in units of seconds

#%%  Plotting
# Regular radial-time stack
	if skip_R != 1:
		stack_array = np.zeros((slowR_n,stack_nt))

		for it in range(stack_nt):  # check points one at a time
			for slowR_i in range(slowR_n):  # for this station, loop over slownesses
				num_val = centralR_st[slowR_i].data[it]
				stack_array[slowR_i, it] = num_val

		y, x = np.mgrid[slice(stack_Rslows[0], stack_Rslows[-1] + slow_delta, slow_delta),
					 slice(ttt[0], ttt[-1] + dt, dt)]

		fig, ax = plt.subplots(1, figsize=(10,3))
		c = ax.pcolormesh(x, y, stack_array, cmap=plt.cm.coolwarm, vmin=-tdiff_clip, vmax=tdiff_clip)
		fig.subplots_adjust(bottom=0.2)
		ax.axis([x.min(), x.max(), y.min(), y.max()])
		fig.colorbar(c, ax=ax)
		plt.xlabel('Time (s)')
		plt.ylabel('Radial slowness (s/km)')
		plt.title(ref_phase + ' Time lag at ' + str(T_slow_plot) + ' s/km transverse slowness, ' + fname1[12:22] + ' ' + fname1[23:33])
		plt.show()

		fig_index += 1
#%%  Transverse-time stacks
	if skip_T != 1:
		stack_array = np.zeros((slowT_n,stack_nt))

		for it in range(stack_nt):  # check points one at a time
			for slowT_i in range(slowT_n):  # for this station, loop over slownesses
				num_val = centralT_st[slowT_i].data[it]
				stack_array[slowT_i, it] = num_val

		y, x = np.mgrid[slice(stack_Tslows[0], stack_Tslows[-1] + slow_delta, slow_delta),
					 slice(ttt[0], ttt[-1] + dt, dt)]

		fig, ax = plt.subplots(1, figsize=(10,3))
		fig.subplots_adjust(bottom=0.2)
		c = ax.pcolormesh(x, y, stack_array, cmap=plt.cm.coolwarm, vmin=-tdiff_clip, vmax=tdiff_clip)
		ax.axis([x.min(), x.max(), y.min(), y.max()])
		fig.colorbar(c, ax=ax)
		plt.xlabel('Time (s)')
		plt.ylabel('Transverse slowness (s/km)')
		plt.title(ref_phase + ' Time lag at ' + str(R_slow_plot) + ' s/km radial slowness, ' + fname1[12:22] + ' ' + fname1[23:33])
		plt.show()

		fig_index += 1
#%% R-T stack time difference
	if skip_snaps == 0:
		stack_slice = np.zeros((slowR_n,slowT_n))
		for snap_num in range(snaps):
			fig_index += 1
			it = int((snaptime - start_buff)/dt) + snap_num
			for slowR_i in range(slowR_n):  # loop over radial slownesses
				for slowT_i in range(slowT_n):  # loop over transverse slownesses
					index = slowR_i*slowT_n + slowT_i
					num_val = tdiff[index].data[it]
					stack_slice[slowR_i, slowT_i] = num_val

			y1, x1 = np.mgrid[slice(stack_Rslows[0], stack_Rslows[-1] + slow_delta, slow_delta),
						 slice(stack_Tslows[0], stack_Tslows[-1] + slow_delta, slow_delta)]

			fig, ax = plt.subplots(1)
			c = ax.pcolormesh(x1, y1, stack_slice, cmap=plt.cm.coolwarm, vmin=-tdiff_clip, vmax=tdiff_clip)
			ax.axis([x1.min(), x1.max(), y1.min(), y1.max()])
			fig.colorbar(c, ax=ax)
			circle1 = plt.Circle((0, 0), 0.019, color='black', fill=False)
			ax.add_artist(circle1)
			plt.xlabel('T Slowness (s/km)')
			plt.ylabel('R Slowness (s/km)')
			plt.title(ref_phase + ' T-R plot of time lag at rel time ' + str(snaptime + snap_num*dt) + '  ' + fname1[12:22] + ' ' + fname1[23:33])
			plt.show()

#%% R-T stack amplitude
	if skip_snaps == 0:
		stack_slice = np.zeros((slowR_n,slowT_n))
		for snap_num in range(snaps):
			fig_index += 1
			it = int((snaptime - start_buff)/dt) + snap_num
			for slowR_i in range(slowR_n):  # loop over radial slownesses
				for slowT_i in range(slowT_n):  # loop over transverse slownesses
					index = slowR_i*slowT_n + slowT_i
					num_val = amp_ave[index].data[it]
					stack_slice[slowR_i, slowT_i] = num_val

			y1, x1 = np.mgrid[slice(stack_Rslows[0], stack_Rslows[-1] + slow_delta, slow_delta),
						 slice(stack_Tslows[0], stack_Tslows[-1] + slow_delta, slow_delta)]

			fig, ax = plt.subplots(1)
			c = ax.pcolormesh(x1, y1, stack_slice, cmap=plt.cm.gist_rainbow_r)
			ax.axis([x1.min(), x1.max(), y1.min(), y1.max()])
			fig.colorbar(c, ax=ax)
			circle1 = plt.Circle((0, 0), 0.019, color='black', fill=False)
			ax.add_artist(circle1)
			plt.xlabel('T Slowness (s/km)')
			plt.ylabel('R Slowness (s/km)')
			plt.title(ref_phase + ' T-R plot of amplitude at rel time ' + str(snaptime + snap_num*dt) + '  ' + fname1[12:22] + ' ' + fname1[23:33])
			plt.show()

	#  Save processed files
#	fname = 'HD' + date_label + '_slice.mseed'
#	stack.write(fname,format = 'MSEED')

	elapsed_time_wc = time.time() - start_time_wc
	print('This job took ' + str(elapsed_time_wc) + ' seconds')
	os.system('say "Done"')