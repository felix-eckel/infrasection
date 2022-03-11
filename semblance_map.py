# -*- coding: utf-8 -*-
"""
Created on Tue Jan 18 17:16:39 2022

@author: felix
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
import sys
from pandas import read_json, read_csv
from obspy import read, UTCDateTime
from obspy.geodetics.base import gps2dist_azimuth
from glob import glob
import utils

mpl.style.use('./style/jupyter-dark.mplstyle')

HOUR = 3600
DAY  = 86400
MIN  = 60

def get_requests(csv):
    return zip(csv[1].values.tolist(), csv[2].values.tolist())

def main():
    # input parameters
    config     = utils.read_config("infrasection.cfg")
    
    # processing queue and meta data
    csv        = read_csv(config["stationlist"], delimiter=" ",
                        header=None, usecols=[1,2])
    n          = len(csv[2].values)
    stations   = read_json(config["metafile"])
    target     = read_json(config["targetfile"], typ='series')
    lat        = target["Latitude"]
    lon        = target["Longitude"]
    
    # determine time limits to include one day of data around the event
    starttime  = UTCDateTime(target["Date"]) - config["leading_hours"]*HOUR
    starttime -= starttime.microsecond/1000
    endtime    = starttime + DAY
    
    # resampling rate
    rsr        = config["resampling"]
    rms_len    = int(np.round(config["rms_length"]*rsr))
    
    # semblance parameters
    lons       = np.arange(-180, 180, 1)
    lats       = np.arange(-90, 90, 1)
    vel        = config["velocity"]
    
    semblance = np.zeros((int(DAY*rsr), len(lons), len(lats)))
    n_semb = np.zeros((len(lons), len(lats)))
    
    # iterate over stations
    for j, request in enumerate(get_requests(csv)):
        sys.stdout.write(\
    '\rprocessing file {:04d} from {} | {:04.1f} %'.format(j+1,n,
                                                          100*(j+1)/n))
        sys.stdout.flush()
        
        # read data
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
        st.filter('bandpass', freqmax=1/config["filter_up"],
                  freqmin=1/config["filter_low"])
        st.resample(rsr)
        
        tr    = st[0]
        
        # get station for current data stream read meta data
        stat  = stations.at[tr.meta.station, tr.meta.network]
        lat_s = stat['latitude']
        lon_s = stat['longitude']
        
        times = utils.get_times(tr, starttime)
        times = times[int(rms_len/2):int(len(times)-(rms_len/2))]
        data  = utils.rolling_rms(tr.data, rms_len)
        
        # add shifted data (^2) according to semblance coordinate
        for i, lo in enumerate(lons):
                for j, la in enumerate(lats):
                    dist = np.abs(gps2dist_azimuth(lat_s, lon_s, la, lo)[0])
                    shift = dist / vel
                    try:
                        i_0 = np.where(times-shift >= 0)[0][0]
                    except:
                        continue
                    if len(data[i_0:]) < semblance.shape[0]:
                        semblance[:len(data[i_0:]), i, j] \
                            += np.power(data[i_0:], 2)
                    else:
                        semblance[:, i, j] \
                            += np.power(data[i_0:semblance.shape[0]+i_0], 2)
                    n_semb[i, j] += 1
    
    # rms semblance
    semblance = np.sqrt(semblance / n_semb)
    
    # initialize figure with cartopy
    fig = plt.figure(figsize=(12,12))
    ax = fig.add_axes([0, 0, 1, 1],
                      projection=ccrs.PlateCarree(central_longitude=180))
    
    # create grids and prepare data
    longrid, latgrid = np.meshgrid(lons, lats)
    semblance_norm  = np.max(semblance, axis=0) / np.max(semblance, axis=None)
    
    # semblance map
    pcm = ax.pcolormesh(longrid, latgrid, semblance_norm.T,
                        transform=ccrs.PlateCarree(), shading='nearest')
    
    # plot stations
    for request in get_requests(csv):
            file = '{}/*/{}.{}*.mseed'.format(config["data_directory"],
                                              request[0], request[1])
            search = glob(file)
            if len(search) >= 1:
                lat = stations.at[request[1], request[0]]['latitude']
                lon = stations.at[request[1], request[0]]['longitude']
                ax.plot(lon, lat, '^', ms=14, color='white', mec='none',
                        transform=ccrs.PlateCarree())
                
    # decorators
    ax.set_global()
    ax.coastlines()
    cax = fig.add_axes([1.01, .25, .025, .5])
    cb = plt.colorbar(cax=cax, mappable=pcm)
    cb.set_label('rms semblance', fontsize=22)
    ax.set_title('rms "semblance"', pad=52, loc='right', fontsize=32)
    ax.text(1, 1.01, '{}\n{}'.format(target["Comment"],
                         UTCDateTime(target["Date"]).strftime('%d %b %Y')),
            ha='right', transform=ax.transAxes, fontsize=22)
    plt.savefig('semblance_map.png', dpi=600)

if __name__ == "__main__":
    main()