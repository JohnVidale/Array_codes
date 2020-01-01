#!/usr/bin/env python
# Program to read event data from Hinet and remove instrument response
# pre-filters, takes UTC request, returns UTC data
# John Vidale 2/2019
def pro1get(eq_file):

	import HinetPy as hinet
	import datetime as dt
	import obspy as obs
	from obspy import Stream
	from obspy import read
	import os
	os.environ['PATH'] += os.pathsep + '/usr/local/bin'
#	print(os.environ)

	#%% Get Hinet station location file
	sta_file = '/Users/vidale/Documents/GitHub/Array_codes/Files/hinet_sta.txt'
	with open(sta_file, 'r') as file:
		lines = file.readlines()
	# Load station coords into arrays
	#station_index = range(781)
	station_index = range(780)
	st_lats  = []
	st_lons  = []
	st_deps  = []
	stations = []
	for ii in station_index:
		line = lines[ii]
		split_line = line.split()
		stations.append(split_line[0])
		stations[ii] = stations[ii].upper()
	#	name_truc_cap  = this_name_truc.upper()
		st_lats.append( split_line[1]) # not used
		st_lons.append( split_line[2])
		st_deps.append( split_line[3])

	dtformat = '%Y-%m-%dT%H:%M:%S.%f'
	offset = 0 # seconds
	duration = 30 # minutes
	pre_filter = (0.1, 0.2, 30, 48)
	components = ['N', 'E', 'Z']

	path = eq_file
	with open(path, 'r') as file:
		lines = file.readlines()

	ids = []
	times = []
	timesUTC = []
	for line in lines:
		split_line = line.split()
		ids.append(split_line[0])
		offset += 9*60*60 # convert UTC -> JST, Japan 9 hours later than Greenwich
		times.append(   dt.datetime.strptime(split_line[1], dtformat) + dt.timedelta(seconds=offset))
		timesUTC.append(dt.datetime.strptime(split_line[1], dtformat)) # keep UTC version for output

	client = hinet.Client('jvidale', 'h2L6VQvLarP3')
	for i in range(len(times)):  # loop over list of requested times
		outdir = 'data/{}'.format(ids[i])

		time    = times[i].strftime(dtformat)[:-3] # remove microseconds
		timeUTC = timesUTC[i].strftime(dtformat)[:-3] # remove microseconds
		data, ctable = client.get_waveform('0101', time, duration,  # Hinet is network 0101
										   data='{}.cnt'.format(ids[i]),
										   ctable='{}.ch'.format(ids[i]),
										   outdir='events/win32',
										   cleanup=True)
		print('i is ' + str(i) + ', ctable is ' + ctable)
		if ctable is None:
			continue

		print('stations has ' + str(len(stations)) + ' entries to retrieve')

		for station in stations:
			hinet.win32.extract_sac(data, ctable, suffix='sac', filter_by_name=station, outdir=outdir, with_pz=False)
	#		with open(ctable, 'r') as file:
			with open('/Users/vidale/Documents/GitHub/Array_codes/Files/hinet_station_list_full.txt', 'r') as file:
				while True:
					line = file.readline()
	#				print('line is ' + line + ' in station ' + station)
					if line.split()[1].lower() == station.lower():
						header = []
						header.append(file.readline().split())
						header.append(file.readline().split())
						header.append(file.readline().split())
	#					print('line is ' + line + ' in station ' + station + ' header is ' + header[0][0])
						break

					if not line or line == '':
						break

			for c in range(3):
				component = components[c]
				filename = '{}_'.format(ids[i]) + station + '.' + component + '.sac'
				filename_ic = '{}_'.format(ids[i]) + station + '.' + component + '_IC.sac'
				filepath = outdir + '/' + filename
				filepath_ic = outdir + '/' + filename_ic

				if component is 'Z':
					old = station + '.U.sac'
				else:
					old = station + '.{}.sac'.format(component)

				try:
					os.rename(outdir + '/' + old, filepath)
				except:
					continue

				if float(header[c][12]) == 1.023e-7:
					default_resp = '/Users/vidale/Documents/GitHub/Array_codes/Files/RESP.TYPE3.txt'
					sensitivity_line = 345
				elif float(header[c][12]) == 1.000e-7:
					default_resp = '/Users/vidale/Documents/GitHub/Array_codes/Files/RESP.TYPE2.txt'
					sensitivity_line = 345
				elif float(header[c][12]) == 1.192e-7:
					default_resp = '/Users/vidale/Documents/GitHub/Array_codes/Files/RESP.TYPE1.txt'
					sensitivity_line = 404

				gain = header[c][7]

				with open(default_resp, 'r') as file:
					resp_file_contents = file.readlines()

				# change station
				old_line = resp_file_contents[3]
				idx = old_line.find('XXXH')
				new_line = old_line[:idx]
				new_line += station + '\n'
				resp_file_contents[3] = new_line

				# change channel
				old_line = resp_file_contents[6]
				idx = old_line.find('EHZ')
				new_line = old_line[:idx]
				new_line += component + '\n'
				resp_file_contents[6] = new_line

				# first gain change at 37th line
				default_gain = '2.000094E+02'
				old_line = resp_file_contents[36]
				idx = old_line.find(default_gain)
				new_line = old_line[:idx]
				new_line += gain + '\n'
				resp_file_contents[36] = new_line

				# compute total counts to m/s constant
				line = resp_file_contents[67]
				sensitivity = line.split()[-1]
				counts2mps = float(gain) * float(sensitivity)

				# second gain change
				old_line = resp_file_contents[sensitivity_line - 1]
				new_line = old_line[:-13]
				new_line += '{:1.6E}'.format(counts2mps) + '\n'
				resp_file_contents[sensitivity_line - 1] = new_line

				resp_file = filename[:-3]
				resp_file += 'resp'

				with open(outdir + '/' + resp_file, 'w') as file:
					file.writelines(resp_file_contents)

				stream = obs.read(filepath)
				if c == 2:
					stream.traces[0].stats['channel'] = 'Z'

				stream.detrend('demean')
				stream.detrend('linear')
				stream.detrend('demean')
				stream_ic = stream.copy()

				stream.traces[0].data /= float(counts2mps)
				stream.write(filepath)

				resp = {'filename': outdir + '/' + resp_file, 'date': obs.UTCDateTime(times[i]), 'units': 'VEL'}
				print(resp)
				stream_ic.simulate(paz_remove=None, pre_filt=pre_filter, seedresp=resp)
				stream_ic.write(filepath_ic)

				try:
					os.remove(outdir + '/' + resp_file)
				except FileNotFoundError:
					pass

		st    = Stream()
		st_all = Stream()
		for station in stations:
			fname_in  = outdir + '/' + ids[i] + '_' + station + '.Z_IC.sac'
			# print(fname_in)
			try:
				st = read(fname_in)  # >>> st = read("my_file.sac")
				st.traces[0].stats.starttime -= 9*60*60 # convert JST back to UTC
				st_all += st
			except:  # skip if the station doesn't exist
				pass
		fname_out = 'HiNet' + time[0:10] + '_wvf' + '.mseed'
		fname_out = 'HiNet' + timeUTC[0:10] + '_wvf' + '.mseed'
		st_all.write(fname_out,format = 'MSEED')
		print("\nEvent #{}/{} done...\n".format(i, len(times)))
	os.system('rm -r events data')
	os.system('say "done"')
