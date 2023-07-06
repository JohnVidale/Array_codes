#!/usr/bin/env python3
import matplotlib.pyplot as plt
#%% close plots
# matplotlib inline
import matplotlib.pyplot as plt
plt.style.use('ggplot')
plt.rcParams['figure.figsize'] = 12, 8
import obspy

# ObsPy automatically detects the file format.
st = obspy.read("data/example.mseed")
print(st)

# Fileformat specific information is stored here.
print(st[0].stats.mseed)
