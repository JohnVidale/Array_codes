#!/usr/bin/env python3
# John Vidale 4/2020

import os
os.environ['PATH'] += os.pathsep + '/usr/local/bin'
os.chdir('/Users/vidale/Documents/GitHub/Array_codes')

#%% Import functions
from run_each_J_PKiKP_sub import run_each_J_PKiKP_sub
from run_each_J_ICS_sub   import run_each_J_ICS_sub

ref_loc  =   1 # 0 select stations by distance from epicenter, 1 select stations by distance from ref location
ref_rad  = 0.25 # radius of stations around ref_loc chosen
ref_lat  =  34
ref_lon  = 133

#run_each_J_PKiKP_sub(start_buff = -20,  end_buff = 25, start_beam =  1,   end_beam = 3,   event_no = 102,
#				ref_loc = ref_loc, ref_rad = ref_rad, ref_lat = ref_lat, ref_lon = ref_lon, dphase = 'PKiKP')
#run_each_J_ICS_sub(start_buff = 10,   end_buff = 70, start_beam = 20,   end_beam = 60,   event_no = 102,
#				ref_loc = ref_loc, ref_rad = ref_rad, ref_lat = ref_lat, ref_lon = ref_lon, dphase = 'PKiKP')

run_each_J_PKiKP_sub(start_buff = -20,  end_buff = 25, start_beam = -2,   end_beam = 0,   event_no = 103,
				ref_loc = ref_loc, ref_rad = ref_rad, ref_lat = ref_lat, ref_lon = ref_lon, dphase = 'PKiKP')
#run_each_J_ICS_sub(start_buff = -25,   end_buff = 80, start_beam = -5,   end_beam = 60,   event_no = 103,
#				ref_loc = ref_loc, ref_rad = ref_rad, ref_lat = ref_lat, ref_lon = ref_lon, dphase = 'PKiKP')

#run_each_J_PKiKP_sub(start_buff = -20,  end_buff = 25, start_beam = -1,   end_beam = 1,   event_no = 104,
#				ref_loc = ref_loc, ref_rad = ref_rad, ref_lat = ref_lat, ref_lon = ref_lon, dphase = 'PKiKP')
#run_each_J_ICS_sub(start_buff = 10,   end_buff = 70, start_beam = 20,   end_beam = 60,   event_no = 104,
#				ref_loc = ref_loc, ref_rad = ref_rad, ref_lat = ref_lat, ref_lon = ref_lon, dphase = 'PKiKP')

#run_each_J_PKiKP_sub(start_buff = -20,  end_buff = 25, start_beam =  0,   end_beam = 2,   event_no = 105,
#				ref_loc = ref_loc, ref_rad = ref_rad, ref_lat = ref_lat, ref_lon = ref_lon, dphase = 'PKiKP')
#run_each_J_ICS_sub(start_buff = 10,   end_buff = 70, start_beam = 20,   end_beam = 60,   event_no = 105,
#				ref_loc = ref_loc, ref_rad = ref_rad, ref_lat = ref_lat, ref_lon = ref_lon, dphase = 'PKiKP')

#run_each_J_PKiKP_sub(start_buff = -20,  end_buff = 25, start_beam =  2,   end_beam = 4,   event_no = 106,
#				ref_loc = ref_loc, ref_rad = ref_rad, ref_lat = ref_lat, ref_lon = ref_lon, dphase = 'PKiKP')
#run_each_J_ICS_sub(start_buff = 10,   end_buff = 70, start_beam = 20,   end_beam = 60,   event_no = 106,
#				ref_loc = ref_loc, ref_rad = ref_rad, ref_lat = ref_lat, ref_lon = ref_lon, dphase = 'PKiKP')

#run_each_J_PKiKP_sub(start_buff = -20,  end_buff = 25, start_beam = -1,   end_beam = 1,   event_no = 107,
#				ref_loc = ref_loc, ref_rad = ref_rad, ref_lat = ref_lat, ref_lon = ref_lon, dphase = 'PKiKP')
#run_each_J_ICS_sub(start_buff = 10,   end_buff = 70, start_beam = 20,   end_beam = 60,   event_no = 107,
#				ref_loc = ref_loc, ref_rad = ref_rad, ref_lat = ref_lat, ref_lon = ref_lon, dphase = 'PKiKP')

#run_each_J_PKiKP_sub(start_buff = -20,  end_buff = 25, start_beam =  0,   end_beam = 2,   event_no = 108,
#				ref_loc = ref_loc, ref_rad = ref_rad, ref_lat = ref_lat, ref_lon = ref_lon, dphase = 'PKiKP')
#run_each_J_ICS_sub(start_buff = 10,   end_buff = 70, start_beam = 20,   end_beam = 60,   event_no = 108,
#				ref_loc = ref_loc, ref_rad = ref_rad, ref_lat = ref_lat, ref_lon = ref_lon, dphase = 'PKiKP')

#run_each_J_PKiKP_sub(start_buff = -20,  end_buff = 25, start_beam = -1,   end_beam = 1,   event_no = 109,
#				ref_loc = ref_loc, ref_rad = ref_rad, ref_lat = ref_lat, ref_lon = ref_lon, dphase = 'PKiKP')
#run_each_J_ICS_sub(start_buff = 10,   end_buff = 70, start_beam = 20,   end_beam = 60,   event_no = 109,
#				ref_loc = ref_loc, ref_rad = ref_rad, ref_lat = ref_lat, ref_lon = ref_lon, dphase = 'PKiKP')

#run_each_J_PKiKP_sub(start_buff = -20,  end_buff = 25, start_beam = -0.5, end_beam = 1.5, event_no = 110,
#				ref_loc = ref_loc, ref_rad = ref_rad, ref_lat = ref_lat, ref_lon = ref_lon, dphase = 'PKiKP')
#run_each_J_ICS_sub(  start_buff = -25,  end_buff = 80, start_beam = -5,   end_beam = 60,   event_no = 110,
#				ref_loc = ref_loc, ref_rad = ref_rad, ref_lat = ref_lat, ref_lon = ref_lon, dphase = 'PKiKP')

#run_each_J_PKiKP_sub(start_buff = -20,  end_buff = 25, start_beam = -2,   end_beam =  0,   event_no = 111,
#				ref_loc = ref_loc, ref_rad = ref_rad, ref_lat = ref_lat, ref_lon = ref_lon, dphase = 'PKiKP')
#run_each_J_ICS_sub(  start_buff = -25,  end_buff = 80, start_beam = -5,   end_beam = 60,   event_no = 111,
#				ref_loc = ref_loc, ref_rad = ref_rad, ref_lat = ref_lat, ref_lon = ref_lon, dphase = 'PKiKP')
