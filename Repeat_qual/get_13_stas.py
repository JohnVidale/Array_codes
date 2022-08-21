#!/usr/bin/env python
# John Vidale 7/2022
# gets a specified set of stations for an individual event, saves file

def get_13_stas(eq_num = 401, fig_index = 1):

    from obspy import UTCDateTime
    from obspy.clients.fdsn import Client
    from obspy import Stream, Trace
    import os
    import matplotlib.pyplot as plt
    import time

    #%% Get catalog

    # plt.close('all')

    start_time_wc = time.time()

    # eq_num = 606
    #%% Get saved event info, also used to name files
    print('Opening location for event ' + str(eq_num))
    fname = '/Users/vidale/Documents/Research/IC/EvLocs/event' + str(eq_num) + '.txt'
    file = open(fname, 'r')
    lines = file.readlines()
    split_line = lines[0].split()
    #            ids.append(split_line[0])  ignore label for now
    t           = UTCDateTime(split_line[1])

    #  new lines to match more specific naming
    date_label    = split_line[1][0:10]
    year_label    = split_line[1][0:4]
    month_label   = split_line[1][5:7]
    day_label     = split_line[1][8:10]
    hour_label    = split_line[1][11:13]
    minute_label  = split_line[1][14:16]

    ev_lat      = float(      split_line[2])
    ev_lon      = float(      split_line[3])
    ev_depth    = float(      split_line[4])
    print('event: date_label ' + date_label + ' time ' + str(t) + ' lat '
       + str(ev_lat) + ' lon ' + str( ev_lon) + ' depth ' + str(ev_depth))

    # etime      = '2017-06-20T12:54:35.91' # event has to be unique and follow time by less than 1 minute
    # etime      = '1998-04-12T21:33:47.42' # event has to be unique and follow time by less than 1 minute
    # ev_lon     = -27.00  # South Sandwich Islands
    # ev_lat     = -56.19
    # ev_depth   =  77.5

    chan_type  = 'HNZ,EHZ,HHZ,HLZ,BHZ' # e.g., BHZ
    network_list = 'GT,II,IU'
    # sta_list = 'COLA,SDV,LPAZ,ESK,AAK,ULN,MAJO,CHTO,CTAO,SNZO,SBA,VNDA'
    # sta_list = 'COLA,LPAZ,ULN,MAJO,SBA,VNDA'
    sta_list = 'AAK,ANMO,COLA,CTAO,DUG,ESK,LPAZ,MAJO,NNA,SBA,SDV,SNZO,TUC,ULN,VNDA'
    start_buff = 100   # Pre-event buffer
    end_buff   = 3600 # Time after OT recorded
    min_dist   = 0      # Minimum distance loaded, degrees
    max_dist   = 180    # Maximum distance loaded, degrees
    verbose    = 1       # 1 prints networks and stations
    s_t = t - start_buff
    e_t = t + end_buff

    client = Client('IRIS')
    # client = RoutingClient("iris-federator")
    # t = UTCDateTime(etime)

    # refined time and hypocenter

    #%% Make inventory of all stations recording this channel
    inventory = client.get_stations(longitude=ev_lon, latitude=ev_lat, starttime =
                                    s_t, endtime = e_t, minradius=min_dist, maxradius=max_dist,
                                    station = '*', channel = '*Z', network = network_list,
                                    level='channel', matchtimeseries=True)
    if verbose == 1:
        print(inventory)
    print('inventory has ' + str(len(inventory)) + ' networks recording data')

    #%% Get waveforms, weed out synthetic waveforms
    st = Stream()
    cnt = 0
    for network in inventory:
        if network.code != 'SY':  # Skip synthetic traces
            for station in network:
                cnt = cnt +1
                if cnt%10 == 0:
                    print(str(cnt) + ' stations examined, ' + str(len(st)) + ' traces extracted')
                try:
                    # st += client.get_waveforms(network.code, station.code, location='00',channel=chan_type, starttime=s_t, endtime = e_t, attach_response=True)
                    st += client.get_waveforms(network.code, station.code, location='*',channel=chan_type, starttime=s_t, endtime = e_t, attach_response=True)
                except:
                    pass
    print('st has data from ' + str(len(st)) + ' traces')
    #%%  How many of which channels
    tr = Trace()
    hnz_chosen   = 0
    ehz_chosen   = 0
    hhz_chosen   = 0
    hlz_chosen   = 0
    bhz_chosen   = 0
    other_chosen = 0
    for tr in st:
        print('Station ' + str(tr.stats.station) + '  location: ' + str(tr.stats.location) + '  channel: ' + str(tr.stats.channel) + '  delta: ' + str(tr.stats.delta))
        if tr.stats.channel == 'HNZ':
            hnz_chosen += 1
        elif tr.stats.channel == 'EHZ':
            ehz_chosen += 1
        elif tr.stats.channel == 'HHZ':
            hhz_chosen += 1
        elif tr.stats.channel == 'HLZ':
            hlz_chosen += 1
        elif tr.stats.channel == 'BHZ':
            bhz_chosen += 1
        else:
            other_chosen += 1

    print('Total channels ' + str(len(st)) + ' - HNZ, EHZ, HHZ, HLZ, BHZ, other have '
           + str(hnz_chosen) + ' ' + str(ehz_chosen) + ' '
           + str(hhz_chosen) + ' ' + str(hlz_chosen) + ' ' + str(bhz_chosen) + ' ' + str(other_chosen))

    #%%
    st_cull = Stream()
    for tr in st:# keep only a single location per station
            reject = False  # keep only a single location per station
            if tr.stats.channel != 'BHZ':
                reject = True
            for tr2 in st_cull:
                if ((tr2.stats.network  == tr.stats.network) &
                    (tr2.stats.station  == tr.stats.station)):
                    reject = True
            tr.detrend(type='simple')
            tr.taper(0.05)
            tr.filter('bandpass', freqmin=1, freqmax=4, corners=2, zerophase=False)
            if reject == False:
                print('Station ' + str(tr.stats.station) + '  delta: ' + str(tr.stats.delta))
                tr.resample(100)
                print('Station ' + str(tr.stats.station) + '  delta: ' + str(tr.stats.delta))
                st_cull += tr
                print(f'trace has {len(tr.data)} time pts and {tr.stats.delta} dt, which is {len(tr.data)*tr.stats.delta:.1f} s, trace starts at {tr.stats.starttime}, event at {t}')
    print('st_cull has data from ' + str(len(st_cull)) + ' traces')
    #%% Plot and write files
    fig = plt.figure(fig_index)
    plt.title(fname)
    plt.xlabel('seconds')
    plt.ylabel('stations')
    st_cull.plot(equal_scale=False, fig=fig)

    mseed_name = year_label + month_label + day_label + '_' + hour_label + minute_label
    write_name = '/Users/vidale/Documents/Research/IC/Mseed/Global/' + mseed_name + '.mseed'

    st_cull.write(write_name,     format = 'MSEED')
    elapsed_time_wc = time.time() - start_time_wc
    print(f'    Whole binning code took   {elapsed_time_wc:.1f}   seconds')
    os.system('say "get event program finished"')
