#!/usr/bin/env python
# Read in 2D stacks for two events
# Compute tdiff, ave_amp, amp_ratio
# Plot radial and transverse cuts through stack, plus beam sum
# Write out tdiff, ave_amp results
# John Vidale 3/2019

def pro6_cc_pair(eq_file1, eq_file2, plot_scale_fac = 0.03, slow_delta = 0.0005,
              slowR_lo = -0.1, slowR_hi = 0.1, slowT_lo = -0.1, slowT_hi = 0.1,
              start_buff = 1040, end_buff = 1180, freq_corr = 1.0,
              get_stf = 0, NS = False,
              ARRAY = 0, max_rat = 1.8, min_amp = 0.2, turn_off_black = 0,
              R_slow_plot = 0, T_slow_plot = 0, tdiff_clip = 0.5,
              ref_loc = False, ref_lat = 36.3, ref_lon = 138.5, dphase = 'PKiKP',
              cc_twin = 2, cc_len = 0.5, cc_interp1d = 5, cc_delta = 0.1, cc_thres = 0.8):

    # import obspy
    from obspy import UTCDateTime
    from obspy import read
    from obspy import Stream
    from obspy.geodetics import gps2dist_azimuth
    from obspy.taup import TauPyModel
    from scipy.signal import hilbert
    import numpy as np
    from termcolor import colored
    import matplotlib.pyplot as plt
    import os
    import math
    import time
    import sys
    from pro_proceed_function import cc_measure_tshift
    model = TauPyModel(model='iasp91')
    # from obspy import Trace
    # import obspy.signal
    # import obspy.signal as sign
    # import statistics

#%% Get info
    # get locations
    print(colored('Running pro6_pair_cc', 'cyan'))

    start_time_wc = time.time()

        #%% Input parameters and computed files
    folder_name = '/Users/vidale/Documents/Research/IC/'
    file1 = open(folder_name + 'EvLocs/' + eq_file1, 'r')
    file2 = open(folder_name + 'EvLocs/' + eq_file2, 'r')
    lines1=file1.readlines()
    lines2=file2.readlines()
    split_line1 = lines1[0].split()
    split_line2 = lines2[0].split()
    t1          = UTCDateTime(split_line1[1])
    # t2 = UTCDateTime(split_line2[1])
    date_label1  = split_line1[1][0:10]
    date_label2  = split_line2[1][0:10]
    ev_lat      = float(      split_line1[2])
    ev_lon      = float(      split_line1[3])
    ev_depth    = float(      split_line1[4])
    # date_label = '2018-04-02' # dates in filename

    #  find predicted slowness
    if ref_loc == False:
        if ARRAY == 0:
            ref_lat = 36.3  # °N, around middle of Japan
            ref_lon = 138.5 # °E
        elif ARRAY == 1:
            ref_lat = 46.7  # °N keep only inner rings A-D
            ref_lon = -106.22   # °E
        elif ARRAY == 2: # China set and center
            ref_lat = 38      # °N
            ref_lon = 104.5   # °E
    ref_dist_az = gps2dist_azimuth(ev_lat,ev_lon,ref_lat,ref_lon)
    ref_back_az = ref_dist_az[2]
    ref_dist    = ref_dist_az[0]/(1000*111)  # distance in °

    arrivals1 = model.get_travel_times(source_depth_in_km=ev_depth,distance_in_degree=ref_dist-0.5,phase_list=[dphase])
    arrivals2 = model.get_travel_times(source_depth_in_km=ev_depth,distance_in_degree=ref_dist+0.5,phase_list=[dphase])
    dtime = arrivals2[0].time - arrivals1[0].time
    event_pred_slo  = dtime/111.  # s/km

    # convert to pred rslo and tslo
    if NS == True:
        sin_baz = np.sin(ref_back_az * np.pi /180)
        cos_baz = np.cos(ref_back_az * np.pi /180)
    #  rotate predicted slowness to N and E
        pred_Nslo = event_pred_slo * cos_baz
        pred_Eslo = event_pred_slo * sin_baz
        print(f'Pred N {pred_Nslo:.4f} Pred E {pred_Eslo:.4f}')
    else:
        pred_Nslo = event_pred_slo
        pred_Eslo = 0
        print(f'Pred R {pred_Nslo:.4f} Pred T {pred_Eslo:.4f}')

    #%% -- read files
    # Get saved event info, also used to name files
    # date_label = '2018-04-02' # date for filename
    goto = '/Users/vidale/Documents/Research/IC/Pro_files'
    os.chdir(goto)
    fname1 = 'HD' + date_label1 + '_2dstack.mseed'
    fname2 = 'HD' + date_label2 + '_2dstack.mseed'
    st1 = Stream()
    st2 = Stream()
    st1 = read(fname1)
    st2 = read(fname2)

    tshift        = st1.copy()  # make array for time shift
    tshift_cc     = st1.copy()  # make array for new time shift
    cc            = st1.copy()  # make array for cc coefficient
    amp_ratio     = st1.copy()  # make array for relative amplitude
    amp_ave       = st1.copy()  # make array for average amplitude

    # decimate several arrays to sampling of time shift
    dec_fac = int(round(cc_delta/st1[0].stats.delta))
    if (dec_fac - cc_delta/st1[0].stats.delta)/cc_delta/st1[0].stats.delta > 0.00001:
        print('dec_fac must be an integer, pick more suitable parameters')
        sys.quit()
    print(f'decimation factor {dec_fac:.3f} original sampling  {st1[0].stats.delta:.3f} correlation sampling  {cc_delta:.3f}.')
    #%% Decimate some files from input dt to output dt
    if dec_fac > 1:
        tshift_cc.decimate(dec_fac, no_filter=True)
        cc.decimate(dec_fac, no_filter=True)

    print('Beams read: event1: ' + str(len(st1)) + ' event2: ' + str(len(st2)))
    nt1 = len(st1[0].data)
    nt2 = len(st2[0].data)
    dt1 = st1[0].stats.delta
    dt2 = st2[0].stats.delta
    nt_new = len(cc[0].data)
    print('Event1: 1st trace has ' + str(nt1) + ' time pts, time sampling of '
          + str(dt1) + ' thus duration ' + str((nt1-1)*dt1))
    print('Event2: 1st trace has ' + str(nt2) + ' time pts, time sampling of '
          + str(dt2) + ' thus duration ' + str((nt2-1)*dt2))
    if nt1 != nt2 or dt1 != dt2:
        print('Trouble, nt or dt not does not match')
        sys.exit(-1)

#%% -- Make grid of slownesses
    # count rows and columns
    slowR_n = int(round(1 + (slowR_hi - slowR_lo)/slow_delta))  # number of slownesses
    slowT_n = int(round(1 + (slowT_hi - slowT_lo)/slow_delta))  # number of slownesses
    # enumerate indices for slowness rows and columns
    # In English, stack_slows = range(slow_n) * slow_delta - slow_lo
    a1R = range(slowR_n)
    a1T = range(slowT_n)
    # define R and T slownesses for each beam
    stack_Rslows = [(x * slow_delta + slowR_lo) for x in a1R]
    stack_Tslows = [(x * slow_delta + slowT_lo) for x in a1T]

    #  Loop over slowness
    total_slows = slowR_n * slowT_n
    print(f'{slowR_n} radial slownesses, low {slowR_lo:.4f} high {slowR_hi:.4f}')
    print(f'{slowT_n} transv slownesses, low {slowT_lo:.4f} high {slowT_hi:.4f}  total slows: {total_slows}')


    global_max = 0
    ttt             = start_buff + (np.arange( st1[0].stats.npts) * st1[0].stats.delta) # time array for inputs
    cc_ttt_full     = start_buff + (np.arange((st1[0].stats.npts-1)/dec_fac + 1) * st1[0].stats.delta * dec_fac) # time array for outputs
    print(f'input npts {len(ttt)} end points {ttt[0]:.3f} {ttt[-1]:.3f} ')
    print(f'output npts {len(cc_ttt_full)}  end points {cc_ttt_full[0]:.3f} {cc_ttt_full[-1]:.3f} ')

    doit_once = False

#%% Find envelope, phase, tshift, and global max
    for slow_i in range(total_slows):

        if slow_i % 50 == 0:
            print('Measuring time shifts, ' + str(slow_i) + ' finished slownesses out of ' + str(total_slows))

        if len(st1[slow_i].data) == 0: # test for zero-length traces, indexing errors
            print('Slowness ' + str(slow_i) + ' trace has zero length, problem!')

        seismogram1 = hilbert(st1[slow_i].data)  # make analytic seismograms
        seismogram2 = hilbert(st2[slow_i].data)

        env1 = np.abs(seismogram1) # amplitude and amp ratio
        env2 = np.abs(seismogram2)
        amp_ave[slow_i].data    = 0.5 * (env1 + env2)
        amp_ratio[slow_i].data  = env1/env2

        # old pointwise method
        angle1 = np.angle(seismogram1) # time shift
        angle2 = np.angle(seismogram2)
        d_phase = (angle1 - angle2)
        for it in range(nt1):
            if   d_phase[it] >      math.pi:
                d_phase[it]  -= 2 * math.pi
            elif d_phase[it] < -1 * math.pi:
                d_phase[it]  += 2 * math.pi
            if   d_phase[it] > math.pi or d_phase[it] < -math.pi:
                print(f'Bad d_phase value {d_phase[it]:.2f}  {it:4d}')

        old_shift = False
#%% -- Old time shift estimate
        # if old_shift == True:  # always do this for comparison, almost free
        tshift[slow_i].data     = d_phase/(2*math.pi*freq_corr)

        # Plot new method for two slowness points (SloR, SloT) == (0.01, +/-0.01):
            # for reference - index = slowR_i*slowT_n + slowT_i  =>  I = (Ri * Tn) + Ti
        # re-compute radial and transverse slowness from slowness array index
        # adjust for dominant frequency of 1.2 Hz, not 1 Hz
        SloT = slowT_lo + (int(round(slow_i %  (slowT_n)) * slow_delta))
        SloR = slowR_lo + (int(round(slow_i /  (slowT_n)) * slow_delta))

        if (abs(SloR - 0.01) < 0.00001) and (abs(SloT + 0.01) < 0.00001):
                pt1 = True
        else:
            pt1 = False
        if (abs(SloR - 0.01) < 0.00001) and (abs(SloT - 0.01) < 0.00001):
                pt2 = True
        else:
            pt2 = False

#%% -- New time shift estimate
        if old_shift == False:
            tr1 = st1[slow_i]
            tr2 = st2[slow_i]
            cc_ttt, cc_coef, tshift_new = cc_measure_tshift( tr1=tr1, tr2=tr2, tarr_beg = start_buff,
                              cc_twin=cc_twin, cc_len=cc_len, cc_delta=cc_delta, cc_interp1d=cc_interp1d)

            # 2 plots at strong ICS slownesses
            if (pt1 == True) or ( pt2 == True):
                print(f'slow_i {slow_i} SloR {SloR:.4f} SloT {SloT:.4f}')
                if pt1 == True:
                    fig4 = plt.figure(4)
                    ax1 = fig4.add_subplot(3,1,1)
                    ax2 = fig4.add_subplot(3,1,2)
                    ax3 = fig4.add_subplot(3,1,3)
                    fig4.suptitle(f'(R,T) slowness ({SloR:.2f}, {SloT:.2f}) s/°')
                    plt.subplots_adjust(hspace=0.8)
                if pt2 == True:
                    fig5 = plt.figure(5)
                    ax1 = fig5.add_subplot(3,1,1)
                    ax2 = fig5.add_subplot(3,1,2)
                    ax3 = fig5.add_subplot(3,1,3)
                    fig5.suptitle(f'(R,T) slowness ({SloR:.2f}, {SloT:.2f}) s/°')
                    plt.subplots_adjust(hspace=0.8)
                ax1.title.set_text(f'Correlation, threshold {cc_thres:.2f}')
                ax2.title.set_text('Time shift (s)')
                ax3.title.set_text('Beams')
                plt.subplot(3,1,1)
                plt.xlim(start_buff,end_buff)
                plt.plot(cc_ttt,cc_coef, 'b')  # cross correlation
                plt.subplot(3,1,2)
                plt.xlim(start_buff,end_buff)
                plt.plot(ttt,tshift[slow_i].data,'k')  # old time shift, replaced with new below
                plt.plot(cc_ttt,tshift_new,'b')  # new time shift
                zero_line = tshift_new * 0
                plt.plot(cc_ttt,zero_line,'c')  # new time shift
                plt.subplot(3,1,3)
                plt.xlim(start_buff,end_buff)
                plt.plot(ttt,st1[slow_i].data, 'y')  # traces for 1st event
                plt.plot(ttt,st2[slow_i].data, 'g')  # traces for 2nd event
                # Mark all the cc>0.8 points
                nn = [index for index,x in enumerate(cc_coef) if x >= cc_thres]
                #cc_tshift[nn]=np.nan
                plt.subplot(3,1,1)
                plt.scatter(cc_ttt[nn],cc_coef[nn],color='r',s=20)
                plt.subplot(3,1,2)
                plt.scatter(cc_ttt[nn],tshift_new[nn],color='r',s=20)

            # extend arrays to full length by adding zeroes to ends
            # tshift_cc.decimate(dec_fac, no_filter=True)
            # cc.decimate(dec_fac, no_filter=True)
            zero_fill = int((len(tshift_cc[0].data) - len(tshift_new)) / 2)  # zeroes or NAN to be added to each end
            cc_coef_full    = np.concatenate([np.zeros(zero_fill),    cc_coef, np.zeros(zero_fill + 1)])
            tshift_new_full = np.concatenate([np.zeros(zero_fill), tshift_new, np.zeros(zero_fill + 1)])
            tshift[slow_i].data = tshift_new_full  # tshift just gets  NaNs for plotting, full is file that is saved
            cc[slow_i].data = cc_coef_full

            if doit_once == False:
                print(f'tshift_new_full length is {len(tshift_new_full)} tshift[0].data length is {len(tshift[0].data)}')
                print(f'length of len(tshift_cc[0].data {len(tshift_cc[0].data)} len(tshift_new) {len(tshift_new)}')
                print(f'length of tester {len(cc_coef_full)} zero-fill {zero_fill}')
                print(f'time window starts at {ttt[0]:.2f} ends at {ttt[-1]:.2f}')
                print(f'beam length is {len(tr1.data)} pts, sampling (s) {tr1.stats.delta:.2f}')
                print(f'timeline length is {len(cc_ttt)} pts, sampling (s) {(cc_ttt[1] - cc_ttt[0]):.2f}')
                print(f'coeff length is {len(cc_coef)} pts')
                print(f'target cc and tshift_cc length is {len(cc[0].data)} pts')
                print(f'tshift length is {len(tshift_new)} pts')
                print(f'corr timeline starts at {cc_ttt[0]:.2f}, ends at {cc_ttt[-1]:.2f}')
                print(f'cc_twin {cc_twin:.2f} cc_delta {cc_delta:.2f}')
                # difference in length is cc_twin/2 + cc_delta
                doit_once = True

        local_max = max(abs(amp_ave[slow_i].data))
        if local_max > global_max:
            global_max = local_max

    #%% Decimate some files from input dt to output dt
    if dec_fac > 1:
        amp_ave_dec = amp_ave.copy()
        amp_ave_dec.decimate(dec_fac, no_filter=True)

    #%% -- Flag less robust points as NAN
    for slow_i in range(total_slows):
        tshift[slow_i].stats.delta = tshift[slow_i].stats.delta * dec_fac
    tshift_full = tshift.copy()  # preserve full array for saved time shift file
    for slow_i in range(total_slows): # ignore less robust points
        if slow_i % 200 == 0:
            print('Applying nan to low correlations, ' + str(slow_i) + ' finished slownesses out of ' + str(total_slows))
        if old_shift == True:
            for it in range(nt1):
                if ((amp_ratio[slow_i].data[it] < (1/max_rat)) or (amp_ratio[slow_i].data[it] > max_rat)
                    or (amp_ave[slow_i].data[it] < (min_amp * global_max))):
                    tshift[slow_i].data[it] = np.nan
        elif old_shift == False:
            print(f'cc         .data length is {len(cc[slow_i].data)}          pts')
            print(f'amp_ave_dec.data length is {len(amp_ave_dec[slow_i].data)} pts')
            print(f'amp_ave.    data length is {len(amp_ave[slow_i].data)}     pts')
            print(f'tshift     .data length is {len(tshift[slow_i].data)}      pts')
            print(f'cc_thres    {cc_thres}')
            print(f'min_amp     {min_amp}')
            print(f'global_max  {global_max}')
            print(f'nt_new  {nt_new}')
            for it in range(nt_new):
                if (cc[slow_i].data[it] < cc_thres) or (amp_ave_dec[slow_i].data[it] < (min_amp * global_max)):
                    tshift[slow_i].data[it] = np.nan

#%% Extract slices to plot
    #%% -- Collect T slowness nearest T_slow
    lowest_Tslow = 1000000
    for slow_i in range(slowT_n):
        if abs(stack_Tslows[slow_i] - T_slow_plot) < lowest_Tslow:
            lowest_Tindex = slow_i
            lowest_Tslow = abs(stack_Tslows[slow_i] - T_slow_plot)

    print(f'{slowT_n} T slownesses, index {lowest_Tindex} is closest to requested plot T slowness {T_slow_plot:.4f}, slowness diff there is {lowest_Tslow:.4f} and slowness is {stack_Tslows[lowest_Tindex]:.4f}')
    # Select only stacks with that slowness for radial plot
    centralR_st1 = Stream()
    centralR_st2 = Stream()
    centralR_amp   = Stream()
    centralR_tdiff = Stream()
    for slowR_i in range(slowR_n):
        ii = slowR_i*slowT_n + lowest_Tindex
        centralR_st1 += st1[ii]
        centralR_st2 += st2[ii]
        centralR_amp   += amp_ave[ii]
        centralR_tdiff += tshift[ii]

    #%% -- Collect R slowness nearest R_slow
    lowest_Rslow = 1000000
    for slow_i in range(slowR_n):
        if abs(stack_Rslows[slow_i] - R_slow_plot) < lowest_Rslow:
            lowest_Rindex = slow_i
            lowest_Rslow = abs(stack_Rslows[slow_i] - R_slow_plot)

    print(f'{slowR_n} R slownesses, index {lowest_Rindex} is closest to requested plot R slowness {R_slow_plot:.4f}, slowness diff there is {lowest_Rslow:.4f} and slowness is {stack_Rslows[lowest_Rindex]:.4f}')

    # Select only stacks with that slowness for transverse plot
    centralT_st1 = Stream()
    centralT_st2 = Stream()
    centralT_amp   = Stream()
    centralT_tdiff = Stream()

    for slowT_i in range(slowT_n):
        ii = lowest_Rindex*slowT_n + slowT_i
        centralT_st1 += st1[ii]
        centralT_st2 += st2[ii]
        centralT_amp   += amp_ave[ii]
        centralT_tdiff += tshift[ii]

#%%  Save processed files
    goto = '/Users/vidale/Documents/Research/IC/Pro_Files'
    os.chdir(goto)

    fname = 'HD' + date_label1 + '_' + date_label2 + '_tshift.mseed'
    tshift_full.write(fname,format = 'MSEED')

    fname = 'HD' + date_label1 + '_' + date_label2 + '_amp_ave.mseed'
    amp_ave_dec.write(    fname,format = 'MSEED')

    fname = 'HD' + date_label1 + '_' + date_label2 + '_cc.mseed'
    cc.write(         fname,format = 'MSEED')

    elapsed_time_wc = time.time() - start_time_wc
    print(f'This job took   {elapsed_time_wc:.1f}   seconds')
    os.system('say "Done"')
