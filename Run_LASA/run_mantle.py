#!/usr/bin/env python
# input is set of hinet traces
# this program tapers, filters, selects range and SNR
# plots against traveltime curves, either raw or reduced against traveltimes
# This programs deals with a single event.
# John Vidale 2/2019

import os
import matplotlib.pyplot as plt
#%% close plots
plt.close('all')

os.chdir('/Users/vidale/Documents/GitHub/Array_codes/Run_LASA')
from run_mantle_one   import run_mantle_one

# eq_num = 1
# start_time = 1050
# wind_len = 200
# ref_phase = 'PKiKP'
# slow_limit = 0.04
# slow_delta = 0.002
# for i in range(1,98):
for i in range(1,2):
    # run_mantle_one(fig_index = 100, eq_num = i, start_time = 1250, wind_len = 50)
    # run_mantle_one(fig_index = 200, eq_num = i, start_time = 1350, wind_len = 50)
    # run_mantle_one(fig_index = 300, eq_num = i, start_time = 1450, wind_len = 50)
    run_mantle_one(fig_index = 400, eq_num = i, start_time = 1550, wind_len = 50)
    # run_mantle_one(fig_index = 500, eq_num = i, start_time = 1650, wind_len = 50)
    # run_mantle_one(fig_index = 600, eq_num = i, start_time = 1750, wind_len = 50)
    # except:
    #     print(str(i) + ' did not work!')
# run_mantle_one(eq_num = 7, start_time = 1450, wind_len=300, slow_limit = 0.1,
#                ref_phase = 'PKiKP', slow_delta = 0.005, wind_buff = 300)