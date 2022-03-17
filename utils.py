# -*- coding: utf-8 -*-
"""
Created on Fri Feb 18 16:19:50 2022

@author: felix
"""

import configparser
import numpy as np

def rolling_rms(x, N):
    # fancy way to write root-mean-square
    xc = np.cumsum(abs(x)**2)
    return np.sqrt((xc[N:] - xc[:-N]) / N)

def get_times(tr, starttime):
    offset = tr.stats.starttime - starttime
    return tr.times() + offset

def fmt_time(t):
    if t < 60:
        return "{:d} s".format(int(np.round(t)))
    elif t < 5400:
        t = t/60
        if t%1 == 0:
            return "{:d} min".format(int(t))
        else:
            return "{:.1f} min".format(np.round(t, decimals=1))
    else:
        t = t/3600
        if t%1 == 0:
            return "{:d} h".format(int(t))
        else:
            return "{:.1f} h".format(np.round(t, decimals=1))

def read_config(configfile):
    config = configparser.ConfigParser()
    config.read(configfile)
    
    data_directory = config.get("DATA", "data_directory")
    stationlist    = config.get("DATA", "stationlist")
    metafile       = config.get("DATA", "metafile")
    targetfile     = config.get("DATA", "targetfile")

    days       = config.getint("PROCESSING", "days")
    filter_up  = config.getfloat("PROCESSING", "filter_up")
    filter_low = config.getfloat("PROCESSING", "filter_low")
    rms_length = config.getint("PROCESSING", "rms_length")
    resampling = 1/config.getint("PROCESSING", "resampling")
    sta_length = config.getint("PROCESSING", "sta_length")
    lta_length = config.getint("PROCESSING", "lta_length")
    thresh_on  = config.getfloat("PROCESSING", "threshold_on")
    thresh_off = config.getfloat("PROCESSING", "threshold_off")
    velocity   = config.getfloat("PROCESSING", "velocity")

    leading_hours = config.getint("PLOTTING", "leading_hours")
    ticks_per_day = config.getint("PLOTTING", "ticks_per_day")
    plotting_type = config.getint("PLOTTING", "type")
    
    velocities = config.get("PLOTTING", "velocities")
    velocities = velocities.replace(' ', '')
    try:
        velocities = [float(vel) for vel in velocities.split(',')]
    except ValueError:
        velocities = []
    
    return {"data_directory": data_directory,
            "stationlist":    stationlist,
            "metafile":       metafile,
            "targetfile":     targetfile,
            "days":           days,
            "rms_length":     rms_length,
            "resampling":     resampling,
            "filter_up":      filter_up,
            "filter_low":     filter_low,
            "sta_length":     sta_length,
            "lta_length":     lta_length,
            "threshold_on":   thresh_on,
            "threshold_off":  thresh_off,
            "velocity":       velocity,
            "leading_hours":  leading_hours,
            "ticks_per_day":  ticks_per_day,
            "plotting_type":  plotting_type,
            "velocities":     velocities}