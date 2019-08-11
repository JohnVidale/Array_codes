import HinetPy as hinet
import datetime as dt
import obspy as obs
from obspy import Stream
from obspy import read
import os

sta_file = 'hinet_station_list_full.txt'

with open(sta_file, 'r') as file:
	lines = file.readlines()

nsta = len(lines)/4
print('File has ' + str(nsta) + ' stations')
print(lines[0])
print(lines[1])

line0 = lines[0]
line1 = lines[1]

split_line1 = line0.split()
split_line2 = line1.split()
i = range(781)
for ii in i:
	line1 = lines[4*(ii-1)]
	line2 = lines[4*(ii-1) + 1]
	split_line1 = line1.split()
	split_line2 = line2.split()
	name = split_line1[1]
	lat = split_line2[13]
	lon = split_line2[14]
	dep = split_line2[15]
	line_out = name + '	' + lat + '	' + lon + '	' + dep
	f = open('hinet_sta','a')
	f.write(line_out + '\n')
	# print(name + '	' + lat + '	' + lon + '	' + dep)
os.system('say "done"')