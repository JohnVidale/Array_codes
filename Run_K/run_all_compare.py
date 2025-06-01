#!/usr/bin/env python
# John Vidale 7/2022

# input is pair of traces
# this program tapers, filters, selects range and SNR
# plots against traveltime curves, either raw or reduced against traveltimes

import sys
import os
import matplotlib.pyplot as plt
# pro_directory = '/Users/vidale/Documents/GitHub/Array_codes/Run_YKA_WRA'
# os.chdir(pro_directory)
sys.path.append('/Users/vidale/Documents/GitHub/Array_codes/Process')
sys.path.append('/Users/vidale/Documents/GitHub/Array_codes/Run_K')
# Print the current working directory (CWD)
cwd = os.getcwd()
print("Run_all current working directory: ", cwd)

from run_compare_ind      import run_compare_ind

run_compare_ind(repeater = 'P500',do_global = False, do_KK = False, do_KUR = True)
# run_compare_ind(repeater = 'P501',do_global = False, do_KK = True, do_KUR = True)
# run_compare_ind(repeater = 'P502',do_global = False, do_KK = True, do_KUR = True)

plt.show()