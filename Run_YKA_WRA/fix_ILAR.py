from obspy import read, UTCDateTime, Stream
import os
os.chdir('/Users/vidale/Documents/Research/IC/Mseed/ILAR')
x = '19990326_1427'
fname_in  = x + '_old.mseed'
fname_out = x + '.mseed'
# 1. Read the mseed file
stream = read(fname_in)
for tr in stream:
    tr.stats.network = "IM"  # change network to "IM"

# 2. Specify your reference time
time1 = UTCDateTime("1999-03-26T14:27:44.00")

# 3. Create a new empty stream for filtered traces
filtered_stream = Stream()

# 4. Loop and apply corrected filtering
for tr in stream:
    station = tr.stats.station
    end_offset = tr.stats.endtime - time1

    if station == "IL31":
        print(f"Dropped {station}: station is IL31")
    elif end_offset <= 20 * 60:
        print(f"Dropped {station}: ends only {end_offset/60:.1f} minutes after time1 (needs >20 min)")
    else:
        print(f"Kept {station}: ends {end_offset/60:.1f} minutes after time1")
        filtered_stream.append(tr)

# 5. Write if there are traces left
if len(filtered_stream) > 0:
    filtered_stream.write(fname_out, format="MSEED")
    print(f"✅ Wrote {len(filtered_stream)} traces to {fname_out}")
else:
    print("🚫 No traces matched the criteria. Nothing written.")