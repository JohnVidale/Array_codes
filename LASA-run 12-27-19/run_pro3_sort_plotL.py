#!/usr/bin/env python
# input is set of hinet traces
# this program tapers, filters, selects range and SNR
# plots against traveltime curves, either raw or reduced against traveltimes
# This programs deals with a single event.
# John Vidale 2/2019

import os

os.environ['PATH'] += os.pathsep + '/usr/local/bin'

os.chdir('/Users/vidale/Documents/GitHub/Array_codes')
from pro3a_sort_plot_pair import pro3pair
from pro3b_sort_plot_singlet import pro3singlet
os.chdir('/Users/vidale/Documents/PyCode/LASA')

# whole network
pro3singlet(ARRAY = 1, stat_corr = 1, eq_file = 'event1.txt',
			start_buff = 10, end_buff = 250,
			plot_scale_fac = 0.1,
			dphase = 'PKiKP', dphase2 = 'PKIKKIKP', dphase3 = 'PP', dphase4 = 'SP',
			freq_min = 1, freq_max = 3,
			min_dist = 58.5, max_dist = 60.5)

# just inner rings
#pro3singlet(ARRAY = 1, stat_corr = 0, eq_file = 'event1.txt', start_buff = 10,
#			end_buff = 70, plot_scale_fac = 0.1, qual_threshold = 0,
#			dphase = 'P', dphase2 = 'PcP', dphase3 = 'PP', dphase4 = 'SP',
#			corr_threshold = 0.7, freq_min = 0.5, freq_max = 2,ref_loc = 1,
#			min_dist = 59.2, max_dist = 59.75)

#North NZ
#pro3pair(ARRAY = 1, stat_corr = 0, eq_file1 = 'event1.txt', eq_file2 = 'event2.txt', start_buff = 150,
#			end_buff = 100, plot_scale_fac = 0.02, qual_threshold = 0,
#			dphase = 'PcP', dphase2 = 'S', dphase3 = 'PP', dphase4 = 'P',
#			corr_threshold = 0.7, freq_min = 1.5, freq_max = 4, min_dist = 59.2, max_dist = 59.73)

#South NZ,
#pro3pair(LASA = 1, stat_corr = 0, eq_file1 = 'event4.txt', eq_file2 = 'event5.txt', start_buff = 100,
#			end_buff = 2000, plot_scale_fac = 0.4, qual_threshold = 0,
#			dphase = 'P', dphase2 = 'PKiKP', dphase3 = 'PKIKPPKIKP', dphase4 = 'PKIKKIKP',
#			corr_threshold = 0.7, freq_min = 0.5, freq_max = 2, min_dist = 61.6, max_dist = 62.2)


# Cannikin
#pro3singlet(ARRAY = 1, stat_corr = 0, eq_file = 'event8.txt', start_buff = 10,
#			end_buff = 70, plot_scale_fac = 0.1, qual_threshold = 0,
#			dphase = 'P', dphase2 = 'PcP', dphase3 = 'PP', dphase4 = 'SP',
#			corr_threshold = 0.7, freq_min = 1, freq_max = 4,ref_loc = 0,
#			min_dist = 46, max_dist = 48.1)


# Antipodal earthquake #10 - BEST
#pro3singlet(ARRAY = 1, stat_corr = 0, eq_file = 'event10.txt', start_buff = 20,
#			end_buff = 40, plot_scale_fac = 0.01, qual_threshold = 0,
#			dphase = 'PKIKP', dphase2 = 'PKIKKIKP', dphase3 = 'PKIKPPKIKP', dphase4 = 'PP',
#			corr_threshold = 0.7, freq_min = 0.5, freq_max = 2,ref_loc = 0,
#			min_dist = 163.9, max_dist = 165.5)

# Antipodal earthquake #11
#pro3singlet(ARRAY = 1, stat_corr = 1, eq_file = 'event11.txt', start_buff = 10,
#			end_buff = 10, plot_scale_fac = 0.01, qual_threshold = 0,
#			dphase = 'PKiKP', dphase2 = 'PKIKKIKP', dphase3 = 'PKIKPPKIKP', dphase4 = 'PP',
#			corr_threshold = 0.7, freq_min = 1, freq_max = 4,ref_loc = 0,
#			min_dist = 126, max_dist = 128.1)

# Antipodal earthquake #9
#pro3singlet(ARRAY = 1, stat_corr = 1, eq_file = 'event9.txt', start_buff = 20,
#			end_buff = 40, plot_scale_fac = 0.03, qual_threshold = 0,
#			dphase = 'PKIKP', dphase2 = 'PKIKKIKP', dphase3 = 'PKIKPPKIKP', dphase4 = 'PP',
#			corr_threshold = 0.7, freq_min = 1, freq_max = 4,ref_loc = 0,
#			min_dist = 159, max_dist = 159.7)
