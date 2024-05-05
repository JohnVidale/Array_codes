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
sys.path.append('/Users/vidale/Documents/GitHub/Array_codes/Run_YKA_WRA')
# Print the current working directory (CWD)
cwd = os.getcwd()
print("Run_all current working directory: ", cwd)

from run_compare_ind      import run_compare_ind

# run_compare_ind(repeater = 'P01')
# run_compare_ind(repeater = 'P02') # long
# run_compare_ind(repeater = 'P03')
# run_compare_ind(repeater = 'P04')
# run_compare_ind(repeater = 'P05')
# run_compare_ind(repeater = 'P06') # long
# run_compare_ind(repeater = 'P07')
# run_compare_ind(repeater = 'P08')
run_compare_ind(repeater = 'P09')
# run_compare_ind(repeater = 'P10')
# run_compare_ind(repeater = 'P11')
# run_compare_ind(repeater = 'P12')
# run_compare_ind(repeater = 'P13')
# run_compare_ind(repeater = 'P14')
# run_compare_ind(repeater = 'P15')
# run_compare_ind(repeater = 'P16')
# run_compare_ind(repeater = 'P17')
# run_compare_ind(repeater = 'P18')
# run_compare_ind(repeater = 'P20')
# run_compare_ind(repeater = 'P21')
# run_compare_ind(repeater = 'P24')
# run_compare_ind(repeater = 'P25')
# run_compare_ind(repeater = 'P26')
# run_compare_ind(repeater = 'P27')
# run_compare_ind(repeater = 'P28')
# run_compare_ind(repeater = 'P29')
# run_compare_ind(repeater = 'P30')
# run_compare_ind(repeater = 'P31') # diffy
# run_compare_ind(repeater = 'P32')
# run_compare_ind(repeater = 'P33') # diffy
# run_compare_ind(repeater = 'P34') # diffy
# run_compare_ind(repeater = 'P35')
# run_compare_ind(repeater = 'P36') # diffy
# run_compare_ind(repeater = 'P37') # diffy
# run_compare_ind(repeater = 'P38') # diffy
# run_compare_ind(repeater = 'P39') # diffy
# run_compare_ind(repeater = 'P40')
# run_compare_ind(repeater = 'P41')
# run_compare_ind(repeater = 'P42')
# run_compare_ind(repeater = 'P43')
# run_compare_ind(repeater = 'P44')
# run_compare_ind(repeater = 'P50')
# run_compare_ind(repeater = 'P51') # nothing
# run_compare_ind(repeater = 'P52')
# run_compare_ind(repeater = 'P53')
# run_compare_ind(repeater = 'P54') # long
# run_compare_ind(repeater = 'P55')
# run_compare_ind(repeater = 'P56')
# run_compare_ind(repeater = 'P57')
# run_compare_ind(repeater = 'P58')
# run_compare_ind(repeater = 'P59')
# run_compare_ind(repeater = 'P60')
# run_compare_ind(repeater = 'P61')
# run_compare_ind(repeater = 'P62')
# run_compare_ind(repeater = 'P63') # diffy
# run_compare_ind(repeater = 'P64')
# run_compare_ind(repeater = 'P65')
# run_compare_ind(repeater = 'P66')
# run_compare_ind(repeater = 'P67')
# run_compare_ind(repeater = 'P68')
# run_compare_ind(repeater = 'P69')
# run_compare_ind(repeater = 'P70')
# run_compare_ind(repeater = 'P71')
# run_compare_ind(repeater = 'P72')
# run_compare_ind(repeater = 'P73')
# run_compare_ind(repeater = 'P74')
# run_compare_ind(repeater = 'P75')
# run_compare_ind(repeater = 'P76')
# run_compare_ind(repeater = 'P77')
# run_compare_ind(repeater = 'P78')
# run_compare_ind(repeater = 'P79') # nothing
# run_compare_ind(repeater = 'P80')
# run_compare_ind(repeater = 'P81')
# run_compare_ind(repeater = 'P82')
# run_compare_ind(repeater = 'P83')
# run_compare_ind(repeater = 'P84')
# run_compare_ind(repeater = 'P85')
# run_compare_ind(repeater = 'P86')
# run_compare_ind(repeater = 'P87')
# run_compare_ind(repeater = 'P88')
# run_compare_ind(repeater = 'P89') # Missing
# run_compare_ind(repeater = 'P90')
# run_compare_ind(repeater = 'P91')
# run_compare_ind(repeater = 'P92')
# run_compare_ind(repeater = 'P93')
# run_compare_ind(repeater = 'P94')
# run_compare_ind(repeater = 'P95')
# run_compare_ind(repeater = 'P96')
# run_compare_ind(repeater = 'P97')
# run_compare_ind(repeater = 'P98')
# run_compare_ind(repeater = 'P99')
# run_compare_ind(repeater = 'P101')
# run_compare_ind(repeater = 'P102')
# run_compare_ind(repeater = 'P103')
# run_compare_ind(repeater = 'P104')
# run_compare_ind(repeater = 'P105')
# run_compare_ind(repeater = 'P106')
# run_compare_ind(repeater = 'P107')
# run_compare_ind(repeater = 'P108')
# run_compare_ind(repeater = 'P109') # nothing
# run_compare_ind(repeater = 'P110')
# run_compare_ind(repeater = 'P111')
# run_compare_ind(repeater = 'P112') # nothing
# run_compare_ind(repeater = 'P113')
# run_compare_ind(repeater = 'P114')
# run_compare_ind(repeater = 'P115')
# run_compare_ind(repeater = 'P116')
# run_compare_ind(repeater = 'P117')
# run_compare_ind(repeater = 'P118')
# run_compare_ind(repeater = 'P119')
# run_compare_ind(repeater = 'P120')
# run_compare_ind(repeater = 'P121')
# run_compare_ind(repeater = 'P122')
# run_compare_ind(repeater = 'P123') # nothing
# run_compare_ind(repeater = 'P124')
# run_compare_ind(repeater = 'P125')
# run_compare_ind(repeater = 'P126')
# run_compare_ind(repeater = 'P127')
# run_compare_ind(repeater = 'P128')
# run_compare_ind(repeater = 'P129')
# run_compare_ind(repeater = 'P130')
# run_compare_ind(repeater = 'P131')
# run_compare_ind(repeater = 'P132')
# run_compare_ind(repeater = 'P133')
# run_compare_ind(repeater = 'P134')
# run_compare_ind(repeater = 'P135')
# run_compare_ind(repeater = 'P136')
# run_compare_ind(repeater = 'P137')
# run_compare_ind(repeater = 'P138')
# run_compare_ind(repeater = 'P139')
# run_compare_ind(repeater = 'P140') # nothing
# run_compare_ind(repeater = 'P141')
# run_compare_ind(repeater = 'P142')
# run_compare_ind(repeater = 'P143')
# run_compare_ind(repeater = 'P144')
# run_compare_ind(repeater = 'P145') # nothing
# run_compare_ind(repeater = 'P146')
# run_compare_ind(repeater = 'P147')
# run_compare_ind(repeater = 'P148')
# run_compare_ind(repeater = 'P149') # long
# run_compare_ind(repeater = 'P150')
# run_compare_ind(repeater = 'P151')
# run_compare_ind(repeater = 'P152')

# run_compare_ind(repeater = 'P200') # 863 similar, not identical
# run_compare_ind(repeater = 'P201') # 863 similar, not identical
# run_compare_ind(repeater = 'P222') # 869
# run_compare_ind(repeater = 'P223') # 870
# run_compare_ind(repeater = 'P231') # 873
# run_compare_ind(repeater = 'P232') # 874
# run_compare_ind(repeater = 'P266') # 888
# run_compare_ind(repeater = 'P302') # 889 good
# run_compare_ind(repeater = 'P303') # 889 good
# run_compare_ind(repeater = 'P316') # 904 good
# run_compare_ind(repeater = 'P318') # 904
# run_compare_ind(repeater = 'P319') # 904
# run_compare_ind(repeater = 'P321') # 905 good
# run_compare_ind(repeater = 'P322') # 905 good
# run_compare_ind(repeater = 'P324') # 906
# run_compare_ind(repeater = 'P325') # 906
# run_compare_ind(repeater = 'P329') # 907
# run_compare_ind(repeater = 'P331') # 908
# run_compare_ind(repeater = 'P332') # 908
# run_compare_ind(repeater = 'P333') # 899
# run_compare_ind(repeater = 'P334') # 899
# run_compare_ind(repeater = 'P335') # 899
# run_compare_ind(repeater = 'P336') # 899 & 908
# run_compare_ind(repeater = 'P337') # 904
# run_compare_ind(repeater = 'P338') # 904
# run_compare_ind(repeater = 'P339') # 905
# run_compare_ind(repeater = 'P340') # 906

# run_compare_ind(repeater = 'P350')
# run_compare_ind(repeater = 'P351')
# run_compare_ind(repeater = 'P352')
# run_compare_ind(repeater = 'P353')
# run_compare_ind(repeater = 'P354')
# run_compare_ind(repeater = 'P356')

# run_compare_ind(repeater = 'COMP')

# Sumatra, IV2 is too close
# run_compare_ind(repeater = 'P45')
plt.show()