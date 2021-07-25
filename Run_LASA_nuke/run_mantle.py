#!/usr/bin/env python
# input is set of hinet traces
# this program tapers, filters, selects range and SNR
# plots against traveltime curves, either raw or reduced against traveltimes
# This programs deals with a single event.
# John Vidale 2/2019

# import os
import matplotlib.pyplot as plt
#%% close plots
plt.close('all')

from run_mantle_one   import run_mantle_one

# eq_num = 1
# start_time = 1050
# wind_len = 200
# ref_phase = 'PKiKP'
# slow_limit = 0.04
# slow_delta = 0.002
# for i in range(1,98):
#     run_mantle_one(eq_num = i, start_time = 1300)
#     run_mantle_one(eq_num = i, start_time = 1400)
#     run_mantle_one(eq_num = i, start_time = 1500)
#     run_mantle_one(eq_num = i, start_time = 1600)
#     run_mantle_one(eq_num = i, start_time = 1700)
#     run_mantle_one(eq_num = i, start_time = 1800)
    # except:
    #     print(str(i) + ' did not work!')
run_mantle_one(eq_num = 7, start_time = 1450, wind_len=300, slow_limit = 0.1,
               ref_phase = 'PKiKP', slow_delta = 0.005, wind_buff = 300)