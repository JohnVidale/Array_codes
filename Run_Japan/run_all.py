#!/usr/bin/env python3
# John Vidale 4/2020

import os
import time
import matplotlib.pyplot as plt

#%% Import functions
from run_individual      import run_individual
from run_dual            import run_dual
from run_bin_one      import run_bin_one
from run_individual_get_shifts      import run_individual_get_shifts

start_time_wc = time.time()

plt.close('all')

#%% Recent runs

# Dual run
# PKiKP fine window
# run_dual(eq_num1 = 161, eq_num2 = 162, dphase = 'PKiKP', Zstart_buff = -3, Zend_buff = 10, start_buff = -20.0, end_buff = 20.0,
#           min_dist = 148, max_dist = 155, beam_offset = 0.018, beam_width  = 0.005, beam_step   = 0.00025)

# PKiKP coda fine window
# run_dual(eq_num1 = 161, eq_num2 = 162, dphase = 'PKiKP', Zstart_buff = 2, Zend_buff = 10, start_buff = -20.0, end_buff = 30.0,
#           min_dist = 148, max_dist = 155, beam_offset = 0.018, beam_width  = 0.006, beam_step   = 0.00025)

# PKiKP coda late fine window
# run_dual(eq_num1 = 161, eq_num2 = 162, dphase = 'PKiKP', Zstart_buff = -10, Zend_buff = 70, start_buff = -20.0, end_buff = 80.0,
#           min_dist = 148, max_dist = 155, beam_offset = 0.018, beam_width  = 0.006, beam_step   = 0.00025)

# PKiKP coda late fine window
# run_dual(eq_num1 = 179, eq_num2 = 180, dphase = 'PKIKP', Zstart_buff = -10, Zend_buff = 70, start_buff = -20.0, end_buff = 80.0,
#           min_dist = 0, max_dist = 180, beam_offset = 0.018, beam_width  = 0.006, beam_step   = 0.00025, JST = True)

# PKiKP focused fine window
# run_dual(eq_num1 = 161, eq_num2 = 162, dphase = 'PKiKP', Zstart_buff = 5.5, Zend_buff = 7.5, start_buff = -20.0, end_buff = 30.0,
#           min_dist = 148, max_dist = 155, beam_offset = 0.018, beam_width  = 0.006, beam_step   = 0.00025)

# PKIKP window
# run_dual(eq_num1 = 161, eq_num2 = 162, dphase = 'PKiKP', Zstart_buff = -12, Zend_buff = -6, start_buff = -40.0, end_buff = 30.0,
#           min_dist = 148, max_dist = 155, beam_offset = 0.01, beam_width = 0.012, beam_step   = 0.0005, freq_min = 0.5, freq_max = 2)

# Individual runs
# run_individual(eq_num = 161, dphase = 'PKiKP', start_buff_stack = -50.0, end_buff_stack = 50.0, start_beam_stack = -5, end_beam_stack =  5, min_dist = 148, max_dist = 160)
# run_individual(eq_num = 162, dphase = 'PKiKP', start_buff_stack = -5.0, end_buff_stack = 5.0, start_beam_stack = -20, end_beam_stack =  20, min_dist = 148, max_dist = 155)

# Get shifts
# run_individual_get_shifts(eq_num = 161, dphase = 'PKiKP', start_buff_align = -50, end_buff_align = 50, precursor_shift = -3, signal_dur = 5, min_dist = 148, max_dist = 155)

# Binning runs
# run_bin_one(event_no = 161, start_buff_stack = -20.0, end_buff_stack = 30.0, rel_time = 0)
# run_bin_one(event_no = 162, start_buff_stack = -20.0, end_buff_stack = 30.0, rel_time = 0)

#%% PKIKP
# run_individual(event_no = 161, dphase = 'PKP', start_buff_stack = -40.0, end_buff_stack = 40.0, start_beam_stack = -20, end_beam_stack =  20, JST = False, min_dist = 0, max_dist = 180)

#  The option I'd immediately like to see work
# run_individual(event_no = 161, dphase = 'PKP', start_buff_stack = -150.0, end_buff_stack = 150.0, start_beam_stack = -50, end_beam_stack =  50, JST = False, min_dist = 145, max_dist = 152.4)

# Repeaters
# run_dual(eq_num1 = 157, eq_num2 = 158, dphase = 'PKiKP', wind_buff = 40, wind_len = 20.0, Zstart_buff = -10, JST = False, min_dist = 115, max_dist = 140)
# run_dual(eq_num1 = 161, eq_num2 = 162, just_do_3 = True, dphase = 'PKIKP', start_buff = -100, end_buff = 100.0, Zstart_buff = 0, Zend_buff = 10, JST = False, min_dist = 148, max_dist = 152)

# PKP & PKiKP
# run_individual(event_no = 162, dphase = 'PKP', start_buff_stack = -200.0, end_buff_stack = 400.0, start_beam_stack =  -40.0, end_beam_stack =  20.0, JST = False, min_dist = 148, max_dist = 150, beam_width = 0.015, beam_step = 0.0005, beam_offset = 0.02)
# run_individual(event_no = 161, dphase = 'PKiKP', start_buff_stack = -20.0, end_buff_stack = 20.0, start_beam_stack =  -4.0, end_beam_stack = -2.0, JST = False, min_dist = 148, max_dist = 155, beam_width = 0.01, beam_step = 0.00025, beam_offset = 0.02)

# Kawakatsu event
run_individual(eq_num = 101, dphase = 'PcP', precursor_shift = 2, signal_dur = 4, start_buff_stack = -20, end_buff_stack = 20, start_beam_stack = 0, end_beam_stack = 6, JST = True,
               min_dist = 13.5, max_dist = 24.5, beam_width = 0.02, beam_step = 0.00025, beam_offset = 0.01, R_slow_plot = 0.013, T_slow_plot = -0.004)

# PcP + sPcP
# run_individual(event_no = 165, dphase = 'PcP', start_buff_stack = -75.0, end_buff_stack = 75.0, start_beam_stack = -10.0, end_beam_stack = 30.0, JST = False, min_dist = 42, max_dist = 60, beam_width = 0.03, beam_step = 0.001, beam_offset = 0.03)

# P wave
# run_individual(event_no = 165, dphase = 'P', start_buff_stack = -50.0, end_buff_stack = 50.0, start_beam_stack = -10.0, end_beam_stack = 10.0, JST = False, min_dist = 42, max_dist = 60, beam_width = 0.05, beam_step = 0.0025, beam_offset = 0.05)

#%% ICS
# run_individual(event_no = 102, dphase = 'PKiKP', start_buff_stack = 0.0, end_buff_stack =  100.0, start_beam_stack =  20.0, end_beam_stack =  60.0)
# run_individual(event_no = 103, dphase = 'PKiKP', start_buff_stack = 0.0, end_buff_stack =  100.0, start_beam_stack =  20.0, end_beam_stack =  60.0)
# run_individual(event_no = 104, dphase = 'PKiKP', start_buff_stack = 0.0, end_buff_stack =  100.0, start_beam_stack =  20.0, end_beam_stack =  60.0)
# run_individual(event_no = 105, dphase = 'PKiKP', start_buff_stack = 0.0, end_buff_stack =  100.0, start_beam_stack =  20.0, end_beam_stack =  60.0)
# run_individual(event_no = 106, dphase = 'PKiKP', start_buff_stack = 0.0, end_buff_stack =  100.0, start_beam_stack =  20.0, end_beam_stack =  60.0)
# run_individual(event_no = 107, dphase = 'PKiKP', start_buff_stack = 0.0, end_buff_stack =  100.0, start_beam_stack =  20.0, end_beam_stack =  60.0)
# run_individual(event_no = 108, dphase = 'PKiKP', start_buff_stack = 0.0, end_buff_stack =  100.0, start_beam_stack =  20.0, end_beam_stack =  60.0)
# run_individual(event_no = 109, dphase = 'PKiKP', start_buff_stack = 0.0, end_buff_stack =  100.0, start_beam_stack =  20.0, end_beam_stack =  60.0)
# run_individual(event_no = 110, dphase = 'PKiKP', start_buff_stack = 0.0, end_buff_stack =  100.0, start_beam_stack =  20.0, end_beam_stack =  60.0)
# run_individual(event_no = 111, dphase = 'PKiKP', start_buff_stack = 0.0, end_buff_stack =  100.0, start_beam_stack =  20.0, end_beam_stack =  60.0)

# run_individual(event_no = 112, dphase = 'PKiKP', start_buff_stack = 0.0, end_buff_stack =  100.0, start_beam_stack =  20.0, end_beam_stack =  60.0)
# run_individual(event_no = 113, dphase = 'PKiKP', start_buff_stack = 0.0, end_buff_stack =  100.0, start_beam_stack =  20.0, end_beam_stack =  60.0)
# run_individual(event_no = 114, dphase = 'PKiKP', start_buff_stack = 0.0, end_buff_stack =  100.0, start_beam_stack =  20.0, end_beam_stack =  60.0)
# run_individual(event_no = 115, dphase = 'PKiKP', start_buff_stack = 0.0, end_buff_stack =  100.0, start_beam_stack =  20.0, end_beam_stack =  60.0)
# run_individual(event_no = 116, dphase = 'PKiKP', start_buff_stack = 0.0, end_buff_stack =  100.0, start_beam_stack =  20.0, end_beam_stack =  60.0)
# run_individual(event_no = 117, dphase = 'PKiKP', start_buff_stack = 0.0, end_buff_stack =  100.0, start_beam_stack =  20.0, end_beam_stack =  60.0)
# run_individual(event_no = 118, dphase = 'PKiKP', start_buff_stack = 0.0, end_buff_stack =  100.0, start_beam_stack =  20.0, end_beam_stack =  60.0)
# run_individual(event_no = 119, dphase = 'PKiKP', start_buff_stack = 0.0, end_buff_stack =  100.0, start_beam_stack =  20.0, end_beam_stack =  60.0)
# run_individual(event_no = 120, dphase = 'PKiKP', start_buff_stack = 0.0, end_buff_stack =  100.0, start_beam_stack =  20.0, end_beam_stack =  60.0)

#%% PKiKP
# run_individual(event_no = 101, dphase = 'PKiKP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =   0.0, end_beam_stack =   5.0)
# run_individual(event_no = 102, dphase = 'PKiKP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =  -1.0, end_beam_stack =   4.0)
# run_individual(event_no = 103, dphase = 'PKiKP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =  -4.0, end_beam_stack =   1.0)
# run_individual(event_no = 104, dphase = 'PKiKP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =  -2.0, end_beam_stack =   4.0)
# run_individual(event_no = 105, dphase = 'PKiKP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =   0.0, end_beam_stack =   4.0)
# run_individual(event_no = 106, dphase = 'PKiKP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =   2.0, end_beam_stack =   5.0)
# run_individual(event_no = 107, dphase = 'PKiKP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =   0.0, end_beam_stack =   3.0)
# run_individual(event_no = 108, dphase = 'PKiKP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =  -2.0, end_beam_stack =   3.0)
# run_individual(event_no = 109, dphase = 'PKiKP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =  -2.0, end_beam_stack =   4.0)
# run_individual(event_no = 110, dphase = 'PKiKP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =  -2.0, end_beam_stack =   4.0)
# run_individual(event_no = 111, dphase = 'PKiKP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =  -4.0, end_beam_stack =   1.0)

# run_individual(event_no = 112, dphase = 'PKiKP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =   6.0, end_beam_stack =   8.0)
# run_individual(event_no = 113, dphase = 'PKiKP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =  11.0, end_beam_stack =  15.0)
# run_individual(event_no = 114, dphase = 'PKiKP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =   2.0, end_beam_stack =   8.0)
# run_individual(event_no = 115, dphase = 'PKiKP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =  -5.0, end_beam_stack =   0.0)
# run_individual(event_no = 116, dphase = 'PKiKP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =   1.0, end_beam_stack =   9.0)
# run_individual(event_no = 117, dphase = 'PKiKP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =   4.0, end_beam_stack =   9.0)
# run_individual(event_no = 118, dphase = 'PKiKP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =   2.0, end_beam_stack =  11.0)
# run_individual(event_no = 119, dphase = 'PKiKP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =   3.0, end_beam_stack =  11.0)
# run_individual(event_no = 120, dphase = 'PKiKP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =   3.0, end_beam_stack =   7.0)

#%% PKIKP
# run_individual(event_no = 150, dphase = 'PKIKP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =   2.0, end_beam_stack =   6.0, JST = False)
# run_individual(event_no = 152, dphase = 'PKIKP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =   0.0, end_beam_stack =   5.0, JST = False)

#%% PcP
# run_individual(event_no = 101, dphase = 'PcP', start_buff_stack =  -5.0, end_buff_stack =  15.0, start_beam_stack =  2.0, end_beam_stack =  8.0)
# run_individual(event_no = 102, dphase = 'PcP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack = -15.0, end_beam_stack =  25.0)
# run_individual(event_no = 102, dphase = 'PcP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =   5.0, end_beam_stack =   8.0)
# run_individual(event_no = 103, dphase = 'PcP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =  -3.0, end_beam_stack =   0.0)
# run_individual(event_no = 104, dphase = 'PcP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =   2.0, end_beam_stack =   6.0)
# run_individual(event_no = 105, dphase = 'PcP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =   2.0, end_beam_stack =   7.0)
# run_individual(event_no = 106, dphase = 'PcP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =   4.0, end_beam_stack =   9.0)
# run_individual(event_no = 107, dphase = 'PcP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =   3.0, end_beam_stack =   7.0)
# run_individual(event_no = 108, dphase = 'PcP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =   3.0, end_beam_stack =   7.0)
# run_individual(event_no = 109, dphase = 'PcP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =   3.0, end_beam_stack =   7.0)
# run_individual(event_no = 110, dphase = 'PcP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =   2.0, end_beam_stack =   5.0)
# run_individual(event_no = 111, dphase = 'PcP', start_buff_stack   = -20.0, end_buff_stack = 30.0, start_beam_stack =  -6.0, end_beam_stack =  -2.0)

# run_individual(event_no = 112, dphase =   'PcP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =  12.0, end_beam_stack =  15.0)
# run_individual(event_no = 113, dphase =   'PcP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =   7.0, end_beam_stack =  11.0)
# run_individual(event_no = 114, dphase =   'PcP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =   5.0, end_beam_stack =   8.0)
# run_individual(event_no = 115, dphase =   'PcP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =  -6.0, end_beam_stack =  -1.0)
# run_individual(event_no = 116, dphase =   'PcP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =   1.0, end_beam_stack =   8.0)
# run_individual(event_no = 117, dphase =   'PcP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =   2.0, end_beam_stack =   9.0)
# run_individual(event_no = 118, dphase =   'PcP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =   6.0, end_beam_stack =  15.0)
# run_individual(event_no = 119, dphase =   'PcP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =   9.0, end_beam_stack =  15.0)
# run_individual(event_no = 120, dphase =   'PcP', start_buff_stack = -20.0, end_buff_stack =  30.0, start_beam_stack =   6.0, end_beam_stack =  18.0)

elapsed_time_wc = time.time() - start_time_wc
print(f'This job took {elapsed_time_wc:.1f} seconds')
os.system('say "All done"')