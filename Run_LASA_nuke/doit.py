#!/usr/bin/env python
# input is set of hinet traces
# this program tapers, filters, selects range and SNR
# plots against traveltime curves, either raw or reduced against traveltimes
# This programs deals with a single event.
# John Vidale 2/2019

from run_pair_func             import runpair

# runpair(Tstart=600, Tend=700, eq_num1=1, eq_num2=2)
# runpair(Tstart=700, Tend=800, eq_num1=1, eq_num2=2)
# runpair(Tstart=800, Tend=900, eq_num1=1, eq_num2=2)
# runpair(Tstart=900, Tend=1000, eq_num1=1, eq_num2=2)
# runpair(Tstart=1000, Tend=1100, eq_num1=1, eq_num2=2)
# runpair(Tstart=1100, Tend=1200, eq_num1=1, eq_num2=2)
# runpair(Tstart=1200, Tend=1300, eq_num1=1, eq_num2=2)
# runpair(Tstart=1300, Tend=1400, eq_num1=1, eq_num2=2)
# runpair(Tstart=1400, Tend=1500, eq_num1=1, eq_num2=2)
# runpair(Tstart=1500, Tend=1600, eq_num1=1, eq_num2=2)
# runpair(Tstart=1600, Tend=1700, eq_num1=1, eq_num2=2)
# runpair(Tstart=1700, Tend=1800, eq_num1=1, eq_num2=2)
# runpair(Tstart=1800, Tend=1900, eq_num1=1, eq_num2=2)
# runpair(Tstart=1900, Tend=2000, eq_num1=1, eq_num2=2)
# runpair(Tstart=2000, Tend=2100, eq_num1=1, eq_num2=2)
# runpair(Tstart=2100, Tend=2200, eq_num1=1, eq_num2=2)
# runpair(Tstart=2200, Tend=2300, eq_num1=1, eq_num2=2)
# runpair(Tstart=2300, Tend=2400, eq_num1=1, eq_num2=2)

# runpair(Tstart=1000, Tend=1300, eq_num1=1, eq_num2=2)
# runpair(Tstart=1100, Tend=1200, eq_num1=1, eq_num2=2)

# runpair(Tstart=600, Tend=700, eq_num1=4, eq_num2=5)
# runpair(Tstart=700, Tend=800, eq_num1=4, eq_num2=5)
# runpair(Tstart=800, Tend=900, eq_num1=4, eq_num2=5)
# runpair(Tstart=900, Tend=1000, eq_num1=4, eq_num2=5)
# runpair(Tstart=1000, Tend=1100, eq_num1=4, eq_num2=5)
# runpair(Tstart=1100, Tend=1200, eq_num1=4, eq_num2=5)
# runpair(Tstart=1200, Tend=1300, eq_num1=4, eq_num2=5)
# runpair(Tstart=1300, Tend=1400, eq_num1=4, eq_num2=5)
# runpair(Tstart=1400, Tend=1500, eq_num1=4, eq_num2=5)
# runpair(Tstart=1500, Tend=1600, eq_num1=4, eq_num2=5)
# runpair(Tstart=1600, Tend=1700, eq_num1=4, eq_num2=5)
# runpair(Tstart=1700, Tend=1800, eq_num1=4, eq_num2=5)
# runpair(Tstart=1800, Tend=1900, eq_num1=4, eq_num2=5)
# runpair(Tstart=1900, Tend=2000, eq_num1=4, eq_num2=5)
# runpair(Tstart=2000, Tend=2100, eq_num1=4, eq_num2=5)
# runpair(Tstart=2100, Tend=2200, eq_num1=4, eq_num2=5)
# runpair(Tstart=2200, Tend=2300, eq_num1=4, eq_num2=5)
# runpair(Tstart=2300, Tend=2400, eq_num1=4, eq_num2=5)

decon78 = False

# runpair(Tstart=500, Tend=600, eq_num1=7, eq_num2=8, decon78 = decon78)
# runpair(Tstart=600, Tend=700, eq_num1=7, eq_num2=8, decon78 = decon78)
# runpair(Tstart=700, Tend=800, eq_num1=7, eq_num2=8, decon78 = decon78)
# runpair(Tstart=800, Tend=900, eq_num1=7, eq_num2=8, decon78 = decon78)
# runpair(Tstart=900, Tend=1000, eq_num1=7, eq_num2=8, decon78 = decon78)
# runpair(Tstart=1000, Tend=1100, eq_num1=7, eq_num2=8, decon78 = decon78)
runpair(Tstart=1100, Tend=1200, eq_num1=7, eq_num2=8, decon78 = decon78)
# runpair(Tstart=1200, Tend=1300, eq_num1=7, eq_num2=8, decon78 = decon78)
# runpair(Tstart=1300, Tend=1400, eq_num1=7, eq_num2=8, decon78 = decon78)
# runpair(Tstart=1400, Tend=1500, eq_num1=7, eq_num2=8, decon78 = decon78)
# runpair(Tstart=1500, Tend=1600, eq_num1=7, eq_num2=8, decon78 = decon78)
# runpair(Tstart=1600, Tend=1700, eq_num1=7, eq_num2=8, decon78 = decon78)
# runpair(Tstart=1700, Tend=1800, eq_num1=7, eq_num2=8, decon78 = decon78)
# runpair(Tstart=1800, Tend=1900, eq_num1=7, eq_num2=8, decon78 = decon78)
# runpair(Tstart=1900, Tend=2000, eq_num1=7, eq_num2=8, decon78 = decon78)
# runpair(Tstart=2000, Tend=2100, eq_num1=7, eq_num2=8, decon78 = decon78)
# runpair(Tstart=2100, Tend=2200, eq_num1=7, eq_num2=8, decon78 = decon78)
# runpair(Tstart=2200, Tend=2300, eq_num1=7, eq_num2=8, decon78 = decon78)
# runpair(Tstart=2300, Tend=2400, eq_num1=7, eq_num2=8, decon78 = decon78)
