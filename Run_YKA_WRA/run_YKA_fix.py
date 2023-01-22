#!/usr/bin/env python
# John Vidale 7/2022
# input is pair of traces
# this program tapers, filters, selects range and SNR
# plots against traveltime curves, either raw or reduced against traveltimes

import os
pro_directory = '/Users/vidale/Documents/GitHub/Array_codes/Run_YKA_WRA'
os.chdir(pro_directory)
from fix_YKA_timing      import fix_YKA_timing

# earlier correction
fix_YKA_timing(eq_num=813)

# Fix 0.1s shift from unknown source
# did
# fix_YKA_timing(eq_num=725)
# fix_YKA_timing(eq_num=729)

# skipped
# fix_YKA_timing(eq_num=730)
# fix_YKA_timing(eq_num=731)
# fix_YKA_timing(eq_num=732)

# did nominal correction
# fix_YKA_timing(eq_num=733)
# fix_YKA_timing(eq_num=734)
# fix_YKA_timing(eq_num=735)
# fix_YKA_timing(eq_num=736)
# fix_YKA_timing(eq_num=737)
# fix_YKA_timing(eq_num=738)
# fix_YKA_timing(eq_num=739)
# fix_YKA_timing(eq_num=740)
# fix_YKA_timing(eq_num=741)
# fix_YKA_timing(eq_num=742)
# fix_YKA_timing(eq_num=743)
# fix_YKA_timing(eq_num=745)
# fix_YKA_timing(eq_num=814)
# fix_YKA_timing(eq_num=815)
# fix_YKA_timing(eq_num=816)
# fix_YKA_timing(eq_num=817)
# fix_YKA_timing(eq_num=819)
# fix_YKA_timing(eq_num=822)
# fix_YKA_timing(eq_num=825)
# fix_YKA_timing(eq_num=826)
# fix_YKA_timing(eq_num=828)
# fix_YKA_timing(eq_num=829)
# fix_YKA_timing(eq_num=830)

# Fix clock errors from data logger changeover, 2013
# fix_YKA_timing(eq_num=744) # 2012-12-03, initial group with no errors
# fix_YKA_timing(eq_num=746) # initial group with no errors
# fix_YKA_timing(eq_num=747) # start of group with logger errors 2014-04-24
# fix_YKA_timing(eq_num=831)
# fix_YKA_timing(eq_num=832)
# fix_YKA_timing(eq_num=748)
# fix_YKA_timing(eq_num=749)
# fix_YKA_timing(eq_num=750)
# fix_YKA_timing(eq_num=751)
# fix_YKA_timing(eq_num=833)
# fix_YKA_timing(eq_num=835)
# fix_YKA_timing(eq_num=836)
# fix_YKA_timing(eq_num=837)
# fix_YKA_timing(eq_num=752) # also no errors
# fix_YKA_timing(eq_num=753)
# fix_YKA_timing(eq_num=754)
# fix_YKA_timing(eq_num=839)
# fix_YKA_timing(eq_num=755)
# fix_YKA_timing(eq_num=756)
# fix_YKA_timing(eq_num=840)
# fix_YKA_timing(eq_num=841)
# fix_YKA_timing(eq_num=757)
# fix_YKA_timing(eq_num=758)
# fix_YKA_timing(eq_num=842)
# fix_YKA_timing(eq_num=843)
# fix_YKA_timing(eq_num=844)
# fix_YKA_timing(eq_num=845) # end of group with logger errors 2020-02-14

# fix_YKA_timing(eq_num=846)
# fix_YKA_timing(eq_num=847)
# fix_YKA_timing(eq_num=848)
# fix_YKA_timing(eq_num=849)
# fix_YKA_timing(eq_num=850)
