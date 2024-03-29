# -*- coding: utf-8 -*-
"""
Section plot for an event.
Data is expected to be available as .mseed files
Use parameters in infrasection.cfg file

@author: felix.eckel@ifg.uni-kiel.de
March 2022
"""
from pandas import read_csv, read_json
from obspy import read, UTCDateTime
from obspy.geodetics.base import locations2degrees
from obspy.signal.trigger import recursive_sta_lta, trigger_onset
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import sys
import utils

mpl.style.use('./style/jupyter-dark.mplstyle')

CIRC_EARTH = 40008000
HOUR       = 3600
DAY        = 86400
MIN        = 60

def main():
    # input parameters
    config     = utils.read_config("infrasection.cfg")
    
    # processing queue and meta data
    csv        = read_csv(config["stationlist"], delimiter=" ",
                        header=None, usecols=[1,2])
    requests   = zip(csv[1].values.tolist(),
                   csv[2].values.tolist())
    n          = len(csv[2].values)
    stations   = read_json(config["metafile"])
    target     = read_json(config["targetfile"], typ='series')
    lat        = target["Latitude"]
    lon        = target["Longitude"]
    rsr        = config["resampling"]
    qc_list    = config["qc_list"]
    if qc_list is not None:
        qc_list = read_csv(qc_list, header=None, index_col=False)[0].values
    
    # load x days of data from 5 h prior to the event
    starttime  = UTCDateTime(target["Date"]) - config["leading_hours"]*HOUR
    starttime -= starttime.microsecond/1000
    offset     = starttime.hour*HOUR + starttime.minute*MIN + starttime.second
    endtime    = starttime + DAY*config["days"]
        
    # plotting: initialize figure
    fig        = plt.figure(figsize=(12, 8))
    ax         = fig.add_axes([0, 0, 1, 1])
    
    # iterate over station queue
    for i, request in enumerate(requests):
        
        if qc_list is not None:
            if not any([request[0] in q and request[1] in q for q in qc_list]):
                continue
            
        # progress info
        s = '\rprocessing file {:04d} from {:d} | {:04.1f} %'.format(i+1, n,
                                                                   100*(i+1)/n)
        sys.stdout.write(s)
        sys.stdout.flush()
        
        # check data availability on disk
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
        
        # sometimes merging fails, then only keep the longest trace in stream
        while len(st) > 1:
            st.remove(st[np.argmin([len(t) for t in st])])
            
        # bandpass filter
        st.filter('bandpass', freqmax=1/config["filter_up"],
                  freqmin=1/config["filter_low"])
        st.resample(rsr)
        
        tr    = st[0]
        
        # get station for current data stream read meta data
        stat  = stations.at[tr.meta.station, tr.meta.network]
        lat_s = stat['latitude']
        lon_s = stat['longitude']
        dist  = locations2degrees(lat, lon, lat_s, lon_s)
        # raw waveform
        if config["plotting_type"] == 0:
            data, times = tr.data, utils.get_times(tr, starttime)
            
        # rms waveform
        elif config["plotting_type"] == 1:
            rms_len = int(config["rms_length"] * rsr)
            times   = utils.get_times(tr, starttime)
            tmp = np.cumsum(abs(tr.data)**2)
            data    = utils.rolling_rms(tr.data, rms_len)
            times   = times[int(rms_len/2):len(times)-int(rms_len/2)]
            if len(data) == 0:
                continue
            
        # get STA/LTA characteristic function
        elif config["plotting_type"] >= 2:
            # get sta/lta parameters
            sta     = config["sta_length"]
            lta     = config["lta_length"]
            
            times   = utils.get_times(tr, starttime)
            
            # init trigger array
            trigger = np.array([])
            
            # calculate characteristic function
            data    = recursive_sta_lta(tr.data, nsta=int(sta*rsr),
                                        nlta=int(lta*rsr))
            
            # the first lta window length has no useful data
            omit    = int(config["lta_length"]*rsr)
            times   = times[omit:]
            data    = data[omit:]
            
            if len(data) == 0:
                continue
            
            # evaluate trigger
            trigger = trigger_onset(data, config["threshold_on"],
                                    config["threshold_off"])
            
            # plot triggers if they exist
            if len(trigger) > 0 and config["plotting_type"] == 3:
                ax.plot(times[trigger[:,0]], np.ones(len(trigger[:,0]))*dist,
                        'o', ms=6, mec='none', mfc='#DC9230', alpha=.25,
                        zorder=181)
            
        # normalize data to 0 - 4 (or -4 - 4) and add distance
        # search for maximum within the first 18 hours (half circumference)
        # max_search = np.where(times > (starttime + \
        #        (18 + config['leading_hours']) * HOUR) - tr.stats.starttime)[0]
        # if len(max_search) == 0:
        #     max_search = None
        # elif max_search[0] == 0:
        #     continue
        # else:
        #     max_search = max_search[0]
        #     # ignore traces that have larger spikes after the first arrival time
        #     if np.argmax(data) > max_search:
        #         continue
        
        data = data / (np.max(data)*.25) + dist
        
        # plot trace as rideline
        ax.fill_between(times, np.ones(len(times))*dist, data,
                        facecolor=mpl.rcParams['figure.facecolor'], 
                        zorder=180-dist)
        ax.plot(times, data, color='C1', lw=.2, zorder=180-dist)
        
    # custom plotting parameters
    xdays       = config["days"]
    ticksperday = config["ticks_per_day"]
    interval    = 24/ticksperday
    vels        = config["velocities"]
    
    # plot theoretical travel times
    if len(vels) > 0:
        origin = config["leading_hours"]*HOUR
        
        for i_vel, vel in enumerate(vels):
            color = mpl.cm.get_cmap('plasma', len(vels)+2)(i_vel+2)
            
            # travel time per degree [s]
            tt_deg = (CIRC_EARTH/360)/vel
            
            # cycles
            c = 0
            while (c*360*tt_deg)/DAY < xdays:
                ax.plot([origin+(c*360+180)*tt_deg, origin+(c*360+360)*tt_deg],
                        [180, 0], lw=1, alpha=.5, color=color, zorder=180)
                ax.plot([origin+(c*360+180)*tt_deg, origin+(c*360)*tt_deg],
                        [180, 0], lw=1, alpha=.5, color=color, zorder=180)
                c += 1
                
            # find correct placement for velocity annotation
            x_text = 24*xdays*HOUR
            y_text = ((24*xdays*HOUR)-origin)/tt_deg
            half_phase = False
            while(y_text > 180):
                y_text -= 180
                half_phase = not half_phase
            if half_phase:
                y_text = 180 - y_text
            ax.text(x_text, y_text,'  {} m/s'.format(vel),
                    ha='left', va='center')
       
    
    
    # indicate midnight by vertical lines
    for xday in range(1, int(np.floor(xdays+1))):
        ax.vlines(x=xday*DAY-offset, ymin=0, ymax=180, lw=1,
                  color='C1', alpha=.4, zorder=-1)
        
    # set x axis tick labels
    ax.tick_params(axis='x', which='both', labelrotation=45)
    ticks = np.arange(-offset, (endtime-starttime), HOUR*interval)
    ticklabels = [(starttime-offset+HOUR*interval*h).strftime("%H:%M") \
                  for h in range(0, len(ticks))]
    ax.set_xticks(ticks, minor=True)
    ax.set_xticklabels(ticklabels, minor=True)
    ticks = np.arange(-offset, (endtime-starttime), DAY)
    ticklabels = [(starttime+DAY*d).strftime("%m/%d") \
                  for d in range(0, len(ticks))]
    ax.set_xticks(ticks)
    ax.set_xticklabels(ticklabels)
    
    # set axes limits
    if config["plotting_type"] > 1:
        ax.set_xlim(config["lta_length"], (24*xdays)*HOUR)
    else:
        ax.set_xlim(0, (24*xdays)*HOUR)
    ax.set_ylim(0, 180)
    
    # set texts
    ax.set_ylabel('distance [deg]', fontsize=22)
    ax.set_xlabel('time', fontsize=22)
    
    # raw waveform
    if config["plotting_type"] == 0:
        title_str = 'waveform'
        
    # rms waveform
    elif config["plotting_type"] == 1:
        title_str = 'rms amplitude'
        
    # STA/LTA characteristic function
    elif config["plotting_type"] > 1:
        title_str = 'characteristic function'
        
        # create a legend
        if config["plotting_type"] == 3:
            ax_legend = fig.add_axes([.84, .02, .14, .07])
            ax_legend.plot(0, 0, 'o', ms=8, color='#DC9230', mec='none')
            ax_legend.text(1.3, 0, 'STA/LTA trigger', ha='right', va='center')
            text_y = -1
        else:
            ax_legend = fig.add_axes([.84, .02, .14, .035])
            text_y = -.8
        ax_legend.set_xticks([])
        ax_legend.set_yticks([])
        sta_text = "STA: {}".format(utils.fmt_time(config["sta_length"]))
        lta_text = "LTA: {}".format(utils.fmt_time(config["lta_length"]))
        ax_legend.text(-.07, text_y, sta_text, ha='left', va='center')
        ax_legend.text(1.3, text_y, lta_text, ha='right', va='center')
        ax_legend.set_xlim([-.3, 1.5])
        ax_legend.set_ylim([-2, 1])
        for spine in ax_legend.spines:
            ax_legend.spines[spine].set_visible(True)
        
    for spine in ax.spines:
        ax.spines[spine].set_zorder(1000)
        
    ax.set_title(title_str, pad=52, loc='right', fontsize=32)
    ax.text(1, 1.01, '{}\n{}'.format(target["Comment"],
                         UTCDateTime(target["Date"]).strftime('%d %b %Y')),
            ha='right', transform=ax.transAxes, fontsize=22)
    
    plt.savefig('section.png', dpi=600)
 
if __name__ == '__main__':
    main()
