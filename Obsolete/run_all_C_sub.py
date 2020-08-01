#!/usr/bin/env python3
# John Vidale 4/2020

import os
os.environ['PATH'] += os.pathsep + '/usr/local/bin'
os.chdir('/Users/vidale/Documents/GitHub/Array_codes')

#%% Import functions
from run_each_C_PKiKP_sub import run_each_C_PKiKP_sub
from run_each_C_ICS_sub   import run_each_C_ICS_sub

ref_loc  =   1 # 0 select stations by distance from epicenter, 1 select stations by distance from ref location
ref_rad  =   4 # radius of stations around ref_loc chosen
ref_lat  =  37
ref_lon  = 108

#run_each_C_PKiKP_sub(start_beam =  1.5, end_beam =  4.0, start_buff = -20,  end_buff = 30,
#					 event_no = 204, ref_loc = ref_loc, ref_rad = ref_rad, ref_lat = ref_lat, ref_lon = ref_lon, dphase = 'PKiKP')
#run_each_C_ICS_sub( start_beam = 20,   end_beam= 60, start_buff = 0,   end_buff = 80,
#				   event_no = 204, ref_loc = ref_loc, ref_rad = ref_rad, ref_lat = ref_lat, ref_lon = ref_lon, dphase = 'PKiKP')

#run_each_C_PKiKP_sub(start_beam =  2.0, end_beam =  4.3, start_buff = -20,  end_buff = 30,
#					 event_no = 205, ref_loc = ref_loc, ref_rad = ref_rad, ref_lat = ref_lat, ref_lon = ref_lon, dphase = 'PKiKP') # good PKiKP
#run_each_C_ICS_sub( start_beam = 20,   end_beam= 60, start_buff = 0,   end_buff = 80,
#				   event_no = 205, ref_loc = ref_loc, ref_rad = ref_rad, ref_lat = ref_lat, ref_lon = ref_lon, dphase = 'PKiKP')

#run_each_C_PKiKP_sub(start_beam =  0.0, end_beam =  3.0, start_buff = -20,  end_buff = 30,
#					 event_no = 206, ref_loc = ref_loc, ref_rad = ref_rad, ref_lat = ref_lat, ref_lon = ref_lon, dphase = 'PKiKP') # moderate PKiKP
#run_each_C_ICS_sub( start_beam = 20,   end_beam= 60, start_buff = 0,   end_buff = 80,
#				   event_no = 206, ref_loc = ref_loc, ref_rad = ref_rad, ref_lat = ref_lat, ref_lon = ref_lon, dphase = 'PKiKP')

#run_each_C_PKiKP_sub(start_beam =  0.0, end_beam =  5,   start_buff = -20,  end_buff = 30,
#					 event_no = 217, ref_loc = ref_loc, ref_rad = ref_rad, ref_lat = ref_lat, ref_lon = ref_lon, dphase = 'PKiKP') # moderate PKiKP
#run_each_C_ICS_sub( start_beam = 20,   end_beam= 60, start_buff = 0,   end_buff = 80,
#				   event_no = 217, ref_loc = ref_loc, ref_rad = ref_rad, ref_lat = ref_lat, ref_lon = ref_lon, dphase = 'PKiKP')

run_each_C_PKiKP_sub(start_beam = 1.0, end_beam =  4.0, start_buff = -20,  end_buff = 30,
					 event_no = 221, ref_loc = ref_loc, ref_rad = ref_rad, ref_lat = ref_lat, ref_lon = ref_lon, dphase = 'PKiKP')
run_each_C_ICS_sub( start_beam = 5,   end_beam= 30, start_buff = -20,   end_buff = 60,
				   event_no = 221, ref_loc = ref_loc, ref_rad = ref_rad, ref_lat = ref_lat, ref_lon = ref_lon, dphase = 'PKiKP')

#run_each_C_PKiKP_sub(start_beam =  0.0, end_beam = 10,   start_buff = -20,  end_buff = 30,
#					 event_no = 226, ref_loc = ref_loc, ref_rad = ref_rad, ref_lat = ref_lat, ref_lon = ref_lon, dphase = 'PKiKP')
#run_each_C_ICS_sub( start_beam = 20,   end_beam= 60, start_buff = 0,   end_buff = 80,
#				   event_no = 226, ref_loc = ref_loc, ref_rad = ref_rad, ref_lat = ref_lat, ref_lon = ref_lon, dphase = 'PKiKP')