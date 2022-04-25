# -*- coding: utf-8 -*-
"""
Created on Tue Jan 18 17:16:39 2022

@author: felix
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import sys
from pandas import read_json, read_csv
from obspy import read, UTCDateTime
from obspy.geodetics.base import gps2dist_azimuth
from obspy.signal.trigger import classic_sta_lta
import utils

mpl.style.use('./style/jupyter-light.mplstyle')

HOUR = 3600
DAY  = 86400
MIN  = 60

def rolling_rms(x, N):
    xc = np.cumsum(abs(x)**2);
    return np.sqrt((xc[N:] - xc[:-N]) / N)

def main():
    # input parameters
    config      = utils.read_config("infrasection.cfg")
    
    # processing queue and meta data
    csv         = read_csv(config["stationlist"], delimiter=" ",
                        header=None, usecols=[1,2])
    requests    = zip(csv[1].values.tolist(),
                    csv[2].values.tolist())
    n           = len(csv[2].values)
    stations    = read_json(config["metafile"])
    target      = read_json(config["targetfile"], typ='series')
    lat         = target["Latitude"]
    lon         = target["Longitude"]
    
    # load x days of data from 5 h prior to the event
    starttime   = UTCDateTime(target["Date"]) - config["leading_hours"]*HOUR
    starttime  -= starttime.microsecond/1000
    offset      = starttime.hour*HOUR + starttime.minute*MIN + starttime.second
    endtime     = starttime + DAY*config["days"]
    
    # resampling rate
    rsr         = config["resampling"]
    rms_len     = int(np.round(config["rms_length"]*rsr))
    xdays       = config["days"]
    ticksperday = config["ticks_per_day"]
    interval    = 24/ticksperday
    
    vels = np.arange(150, 600, 1)
    semblance = np.zeros((int(xdays*DAY*rsr), len(vels)))
    n_semb = np.zeros(len(vels))
    for j, request in enumerate(requests):
        sys.stdout.write(\
    '\rprocessing file {:04d} from {} | {:04.1f} %'.format(j+1,n,
                                                          100*(j+1)/n))
        sys.stdout.flush()
        try:
            file = '{}/*/{}.{}*.mseed'.format(config["data_directory"],
                                              request[0], request[1])
            st = read(file, starttime=starttime, endtime=endtime)
        except:
            # data not on disk, skip request
            continue
        
        # data exists but is empty, skip request
        if len(st) == 0:
            continue
        
        # data processing
        st.detrend(type='linear')
        st.merge(fill_value='interpolate')
        st.trim(starttime=starttime, pad=True, fill_value=0)
        
        # sometimes merging fails, then only keep the longest trace in stream
        while len(st) > 1:
            st.remove(st[np.argmin([len(t) for t in st])])
            
        # bandpass filter
        st.filter('bandpass',
                  freqmax=1/config["filter_up"],
                  freqmin=1/config["filter_low"])
        st.resample(rsr)
        
        tr    = st[0]
        
        # get station for current data stream read meta data
        stat  = stations.at[tr.meta.station, tr.meta.network]
        lat_s = stat['latitude']
        lon_s = stat['longitude']
        dist  = gps2dist_azimuth(lat, lon, lat_s, lon_s)[0]
        
        # time and data vector
        times = utils.get_times(tr, starttime)
        times = times[int(rms_len/2):int(len(times)-(rms_len/2))]
        data  = rolling_rms(tr.data, rms_len)
        
        # add shifted data (^2) according to semblance velocity
        for i, vel in enumerate(sorted(vels)):
            shift = dist / vel
            try:
                i_0 = np.where(times-shift >= 0)[0][0]
            except:
                break
            if len(data[i_0:]) < semblance.shape[0]:
                semblance[:len(data[i_0:]), i] += np.power(data[i_0:], 2)
            else:
                semblance[:, i] += \
                    np.power(data[i_0:semblance.shape[0]+i_0], 2)
            n_semb[i] += 1
    
    # rms semblance
    semblance = np.sqrt(semblance / n_semb)
        
    # initialize figure
    fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(12,10), sharex=True,
                            gridspec_kw={'height_ratios': [.7, .3]})
    
    # create grids and prepare data
    times             = np.arange(starttime.timestamp,
                                  (starttime+xdays*DAY).timestamp, 1/rsr)
    velgrid, timegrid = np.meshgrid(vels, times)
    semblance_norm    = (semblance / (np.max(semblance, axis=None)))
    
    # semblance plotting
    axs[0].pcolormesh(timegrid, velgrid, semblance_norm, shading='nearest')
    
    # detection of semblance maximum
    i_t, i_v = np.where(semblance == np.max(semblance, axis=None))
    i_t, i_v = int(i_t), int(i_v)
    
    # search for origin time (first arrival) within 3 h around semblance max
    search_win  = np.arange(int(i_t-7200*rsr),int(i_t+7200*rsr))
    semb_search = semblance[search_win,i_v]
    
    # "on"-threshold for sta/lta based on 67th percentile of semblance
    i_t_on  = np.where(semb_search > np.percentile(semb_search, 67))[0][0]
    i_t_on += search_win[0]
    
    # calculate sta/lta on 2:24 ratio
    cft       = classic_sta_lta(semblance[:,i_v], 2, 24)
    
    # find first arrival on a 50 sample window around threshold
    i_trigger = np.where(cft[i_t_on-50:i_t_on]<1)[0][-1]+(i_t_on-50)
    
    # string representation of first arrival (calculated origin time of event)
    timefmt = UTCDateTime(times[i_trigger]).strftime('%H:%M')
    
    # indicate maximum semblance
    axs[0].plot(times[i_t], vels[i_v], marker='o', mfc='none', mec='white',
                mew=4, ms=16)
    axs[0].annotate('{:.0f} m/s'.format(float(vels[i_v])),
                  (times[i_t]*0.8, vels[i_v]*1.1), ha='right', va='center',
                  fontsize=22)
    
    # plot velocity trace with maximum semblance to subplot
    axs[1].plot(times, semblance[:, i_v])
    axs[1].annotate(timefmt, (times[i_trigger], semblance[i_trigger,i_v]),
                    xytext=(times[i_trigger]-7000,
                            semblance[i_t_on, i_v]*1.2),
                    arrowprops={'arrowstyle': 'simple', 'facecolor': 'C1',
                                'lw': 0}, fontsize=22)
    
    # subplot decorators
    axs[1].set_ylim(0, max(semblance[:,i_v])*1.1)
    axs[1].text(.99, .93, 'rms "semblance" at {} m/s'.format(vels[i_v]),
                transform=axs[1].transAxes, ha='right', fontsize=22, va='top')
    axs[1].set_ylabel('rms\n"semblance"')
    for spine in axs[1].spines:
        axs[1].spines[spine].set_visible(True)
    xticks      = [((starttime-offset) + h).timestamp \
                   for h in range(0, xdays*DAY, int(interval*HOUR))]
    xticklabels = [UTCDateTime(xt).strftime('%H:%M') for xt in xticks]
    axs[1].set_xticks(xticks, minor=True)
    axs[1].set_xticklabels(xticklabels, minor=True)
    xticks      = [((starttime-offset) + DAY*d).timestamp \
                   for d in range(xdays)]
    xticklabels = [UTCDateTime(xt).strftime('%Y/%m/%d') for xt in xticks]
    axs[1].set_xticks(xticks)
    axs[1].set_xticklabels(xticklabels)
    axs[1].set_xlim(times[0],times[-1])
    axs[1].set_xlabel('time')
    axs[1].set_yticks([])
    
    # main plot decorators
    axs[0].set_yticks(vels[np.where(vels%50 == 0)])
    axs[0].set_ylabel('velocity [m/s]')
    for spine in axs[0].spines:
        axs[0].spines[spine].set_visible(True)
    axs[0].set_title('rms "semblance" analysis', pad=52)
    axs[0].text(1, 1.01, '{}\n{}'.format(target["Comment"],
                         UTCDateTime(target["Date"]).strftime('%d %b %Y')),
            ha='right', transform=axs[0].transAxes, fontsize=22)
    axs[0].text(.8, .71, r'$S_t=\sqrt{\frac{1}{N}\sum_{i=1}^{N}f_{ti}^2}$' + \
                '\n\n',
                fontsize=18, ha='left', va='bottom',
                transform=axs[0].transAxes)
    _text = 'N: traces time corrected\n       with velocity and\n       ' \
            + 'station-event distance'
    axs[0].text(.795, .71, _text, transform=axs[0].transAxes,
                fontsize=11)
    axs[0].add_patch(mpl.patches.Rectangle((.785, .69), .2, .285,
                                    facecolor=mpl.rcParams['figure.facecolor'],
                                    edgecolor='white',
                                    transform=axs[0].transAxes))
    
    # finalize and save
    plt.subplots_adjust(hspace=0)
    plt.savefig('semblance.png', dpi=600)

if __name__ == "__main__":
    main()