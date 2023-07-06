#!/usr/bin/env python
# John Vidale 7/2022
# input is pair of traces
# this program tapers, filters, selects range and SNR
# plots against traveltime curves, either raw or reduced against traveltimes

import os
pro_directory = '/Users/vidale/Documents/GitHub/Array_codes/Run_YKA_WRA'
os.chdir(pro_directory)
from fix_YKA_timing      import fix_YKA_timing

# fix_YKA_timing(eq_num=701)
# fix_YKA_timing(eq_num=702)
# fix_YKA_timing(eq_num=703)
# fix_YKA_timing(eq_num=704)
# fix_YKA_timing(eq_num=800)
# fix_YKA_timing(eq_num=801)
# fix_YKA_timing(eq_num=802)
# fix_YKA_timing(eq_num=803)
# fix_YKA_timing(eq_num=705)
# fix_YKA_timing(eq_num=706)

# fix_YKA_timing(eq_num=707)
# fix_YKA_timing(eq_num=708)
# fix_YKA_timing(eq_num=709)
# fix_YKA_timing(eq_num=710)
# fix_YKA_timing(eq_num=711)
# fix_YKA_timing(eq_num=712)
# fix_YKA_timing(eq_num=804)
# fix_YKA_timing(eq_num=713)
# fix_YKA_timing(eq_num=805)
# fix_YKA_timing(eq_num=806)

# fix_YKA_timing(eq_num=714)
# fix_YKA_timing(eq_num=715)
# fix_YKA_timing(eq_num=716)
# fix_YKA_timing(eq_num=717)
# fix_YKA_timing(eq_num=807)
# fix_YKA_timing(eq_num=718)
# fix_YKA_timing(eq_num=719)
# fix_YKA_timing(eq_num=808)
# fix_YKA_timing(eq_num=720)
# fix_YKA_timing(eq_num=721)

# fix_YKA_timing(eq_num=722)
# fix_YKA_timing(eq_num=723)
# fix_YKA_timing(eq_num=724)
# fix_YKA_timing(eq_num=809)
# fix_YKA_timing(eq_num=810)
# fix_YKA_timing(eq_num=725) # Noted in "fix" but not set - Fix 0.1s shift
# fix_YKA_timing(eq_num=726)
# fix_YKA_timing(eq_num=727)
# fix_YKA_timing(eq_num=811)
# fix_YKA_timing(eq_num=728)

# fix_YKA_timing(eq_num=729) # Noted in "fix" but not set - Fix 0.1s shift
# fix_YKA_timing(eq_num=730)
# fix_YKA_timing(eq_num=812)
# fix_YKA_timing(eq_num=731)
# fix_YKA_timing(eq_num=732) # Fix 0.1s shift
# fix_YKA_timing(eq_num=813)
# fix_YKA_timing(eq_num=733) # Fix 0.1s shift
# fix_YKA_timing(eq_num=734) # disabled - Fix 0.1s shift
# fix_YKA_timing(eq_num=814) # Fix 0.1s shift
# fix_YKA_timing(eq_num=735) # Fix 0.1s shift

# fix_YKA_timing(eq_num=815) # Fix 0.1s shift
# fix_YKA_timing(eq_num=736)
# fix_YKA_timing(eq_num=737) # Fix 0.1s shift
# fix_YKA_timing(eq_num=738) # Fix 0.1s shift
# fix_YKA_timing(eq_num=739) # Fix 0.1s shift
# fix_YKA_timing(eq_num=816) # Fix 0.1s shift
# fix_YKA_timing(eq_num=817)
# fix_YKA_timing(eq_num=740) # Fix 0.1s shift
# fix_YKA_timing(eq_num=818)
# fix_YKA_timing(eq_num=819) # Fix 0.1s shift

# fix_YKA_timing(eq_num=820)
# fix_YKA_timing(eq_num=821)
# fix_YKA_timing(eq_num=822)
# fix_YKA_timing(eq_num=741) # Fix 0.1s shift
# fix_YKA_timing(eq_num=823)
# fix_YKA_timing(eq_num=742)
# fix_YKA_timing(eq_num=824)
# fix_YKA_timing(eq_num=743)
# fix_YKA_timing(eq_num=825)
# fix_YKA_timing(eq_num=826) # Fix 0.1s shift

# fix_YKA_timing(eq_num=827)
# fix_YKA_timing(eq_num=828) # all the rest also have "Fix 0.1s shift"
# fix_YKA_timing(eq_num=744) # 2012-12-03, initial group with no errors
# fix_YKA_timing(eq_num=745)
# fix_YKA_timing(eq_num=829)
# fix_YKA_timing(eq_num=830)
# fix_YKA_timing(eq_num=746) # initial group with no errors # Fix clock errors from data logger changeover, 2013   has doubled traces
# fix_YKA_timing(eq_num=747) # start of group with logger errors 2014-04-24
# fix_YKA_timing(eq_num=831)
# fix_YKA_timing(eq_num=832)

# fix_YKA_timing(eq_num=748)
# fix_YKA_timing(eq_num=749)
# fix_YKA_timing(eq_num=750)
# fix_YKA_timing(eq_num=751)
# fix_YKA_timing(eq_num=833)
# fix_YKA_timing(eq_num=834) # no unshifted?
# fix_YKA_timing(eq_num=752) # also no errors
# fix_YKA_timing(eq_num=835)
# fix_YKA_timing(eq_num=836)
# fix_YKA_timing(eq_num=837)

# fix_YKA_timing(eq_num=753)
# fix_YKA_timing(eq_num=754)
# fix_YKA_timing(eq_num=838)
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
fix_YKA_timing(eq_num=845) # end of group with logger errors 2020-02-14
# fix_YKA_timing(eq_num=846)
# fix_YKA_timing(eq_num=847)
# fix_YKA_timing(eq_num=848)
# fix_YKA_timing(eq_num=849)
# fix_YKA_timing(eq_num=850)
