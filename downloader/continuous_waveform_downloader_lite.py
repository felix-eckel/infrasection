# -*- coding: utf-8 -*-
"""
Downloads .mseed files from EIDA nodes. Requests *DF channels (infrasound)

Usage: 
    Create an input list (e.g. with collect_metadata_lite) and specify
    parameters in the donwload.cfg file.

@author: felix.eckel@ifg.uni-kiel.de
March 2022

Based on the works of M. Timkó (Institute of Earth Physics and Space Science)
"""

import pandas as pd
import numpy as np
import os
import obspy
import timeit
import time
import configparser
import logging
from obspy import UTCDateTime
from obspy.clients.fdsn.client import FDSNException
from obspy.clients.fdsn import Client

def run():
    # read config
    config = configparser.ConfigParser()
    config.read("./download.cfg")
    
    # various parameters
    input_file      = config.get("DOWNLOAD", "station_list_path")
    processing      = config.getboolean("PROCESSING", "processing")
    name_ext        = "VEL" if config["processing"] else "RAW"
    token           = config.get("DOWNLOAD", "token_path")
    component       = "F"
    channels        = config.get("DOWNLOAD",
                                 "channels").replace(' ','').split(',')
    date_fmt        = "%Y%m%d"
    override        = config.getboolean("DOWNLOAD", "override")
    save_dir        = config.get("DOWNLOAD", "saving_directory")
    dt              = config.getfloat("DOWNLOAD", "dt")
    max_gap         = config.getfloat("DOWNLOAD", "max_gap")
    data_perc       = config.getfloat("DOWNLOAD", "data_percentage")
    sleep_time      = config.getint("DOWNLOAD", "sleep_time")
    anti_alias_filt = config.getfloat("PROCESSING", "anti_aliasing_filter")
    filter_order    = config.getint("PROCESSING", "filter_order")
    zero_phase      = config.getboolean("PROCESSING", "zero_phase")
    apply_bb_filter = config.getboolean("PROCESSING", "apply_broadband_filter")
    bb_filter       = list(map(float, config.get("PROCESSING",
                             "broadband_filter").replace(' ', '').split(",")))
    resample        = config.getboolean("PROCESSING", "resample")
    sampling_rate   = config.getint("PROCESSING", "sampling_freq")
    
    try:
        # read input file
        df = pd.read_csv(input_file, header=None,
                         delimiter=" ", comment='#',
                         names=["client", "network", "station",
                                "start_time", "end_time"])
    except IOError as e:
        print(e)
    
    # initialize clients
    clients = {}
    
    # initialize logger
    dt_string = UTCDateTime.now().strftime("%Y-%m-%d-%H_%M_%S")
    if (not os.path.exists("./logs")):
        os.mkdir("./logs")
    log_format = "%(asctime)s::%(levelname)s::%(filename)s::" \
        + "%(lineno)d::%(message)s"
    logging.basicConfig(filename="./logs/{}.log".format(dt_string),
                        level='INFO',  format=log_format)
    logger = logging.getLogger("data-download")
    
    # iterate over input file requests
    for _, row in df.iterrows():
        print ("{}.{}".format(row["network"],row["station"]))
        
        # add client
        client = row["client"]
        if(client not in clients):
            try:
                clients[client] = Client(client, eida_token=token)
            except obspy.clients.fdsn.client.FDSNException:
                print ("Token is not accepted. " \
                       + "Init {} without token".format(client))
                clients[client] = Client(client)
            except ValueError:
                print ("Token does not exist. " \
                       + "Init {} without token".format(client))
                clients[client] = Client(client)
            
        # get start/end times
        t_end = UTCDateTime(row["end_time"])
        t     = UTCDateTime(row["start_time"])
        
        # loop over intervals from start to end
        while t <= t_end:
            # specify save sub-directory and filename
            subdir = "{}/{}/{}/".format(component, t.year,
                                         t.strftime(date_fmt))
            filename = "{}.{}.{}_{}{}.mseed".format(row["network"],
                                                    row["station"],
                                                    component,
                                                    name_ext,
                                                    t.strftime(date_fmt))
            
            # logging
            start = timeit.default_timer()
            message = "{}::{}::{}::{}::{}".format(row["client"], 
                                                  row["network"],
                                                  row["station"],
                                                  component,
                                                  t.strftime(date_fmt))
            try:
                if (not override \
                    and os.path.exists("{}/{}/{}".format(save_dir, subdir,
                                                         filename))):
                    raise Waveform_exist(filename)
                
                # get waveforms
                
                # iterate over channels
                for channel in channels:
                    con = "{}DF".format(channel)
                    
                    # initialize flags and container
                    quality = True
                    fdsn = True
                    waveform = inventory = None
                    
                    try:
                        waveform = clients[row["client"]].get_waveforms(
                            network = row["network"], 
                            station = row["station"], 
                            location = '*', 
                            channel = con, 
                            starttime = t, 
                            endtime = t + dt, 
                            attach_response = False
                        )
                        inventory = clients[row["client"]].get_stations(
                            network = row["network"], 
                            station = row["station"], 
                            location = '*', 
                            channel = con, 
                            starttime = t, 
                            endtime = t + dt, 
                            level = "response"
                        )
                        
                        # get samples in stream
                        samples = np.sum([tr.count() for tr in waveform])
                        # get max gap in stream
                        gap = max([g[-1] for g in waveform.get_gaps()])
                        
                        # assess waveform quality
                        if not (gap < waveform[0].stats.sampling_rate*max_gap \
                            and samples / (waveform[0].stats.sampling_rate*dt)\
                                > data_perc):
                            raise Quality_error(row["network"],
                                                row["station"], t)
                            
                        waveform.merge(fill_value="interpolate")
                        waveform = waveform.pop()
                        
                    except FDSNException:
                        fdsn = False
                    except Quality_error:
                        quality = False
                
                # handle exceptions and errors
                if(not fdsn):
                    raise FDSN_error(row["network"], row["station"], t)
                elif(waveform is None):
                    raise Waveform_problem(row["network"], row["station"], t)
                elif(inventory is None):
                    raise Inventory_problem(row["network"], row["station"], t)
                elif (not quality):
                    raise Quality_error(row["network"], row["station"], t)
                raise Unknown_problem(row["network"], row["station"], t)
                
                # process waveform
                
                # logging
                start_proc = timeit.default_timer()
                processing = ""
                
                try:
                    if config["processing"]:
                        # detrend
                        processing += "Detrend -> "
                        waveform.detrend(type='linear')
                        processing += "Demean -> "
                        waveform.detrend(type="demean")
                        
                        # taper
                        processing += "Cosine type taper: 0.05 % ->"
                        waveform.taper(type="cosine", max_percentage=0.05)
                        
                        # filter
                        processing += "lowpass filter at: " \
                            + "{} sec. ".format(anti_alias_filt) \
                            + "Filter-order: {}, ".format(filter_order) \
                            + "zerophase: {} -> ".format(zero_phase)
                        waveform.filter(type="lowpass", freq=1./anti_alias_filt,
                                        corners=filter_order,
                                        zerophase=zero_phase)
                        
                        # resample
                        if (resample \
                            and waveform.stats.sampling_rate > sampling_rate):
                            processing += "resample to: {} Hz -> ".format( \
                                                              sampling_rate)
                            waveform.interpolate(sampling_rate=sampling_rate,
                                             method="weighted_average_slopes")
                            
                        # remove respones
                        processing += "remove response -> "
                        waveform.remove_response(inventory=inventory)

                        # filter
                        if (apply_bb_filter):
                            processing += "bandpass filter between: " \
                                + "{} and {} sec. ".format(bb_filter[0],
                                                           bb_filter[1]) \
                                + "Filter-order: {}, ".format(filter_order) \
                                + "zerophase: {} -> ".format(zero_phase)
                            waveform.filter(type="bandpass",
                                            freqmin=1./bb_filter[0],
                                            freqmax=1./bb_filter[1],
                                            corners=filter_order,
                                            zerophase=zero_phase)
                    
                    # only remove sensitivty
                    else:
                        processing += "remove sensitivity ->"
                        waveform.remove_sensitivity(inventory)

                    proc_time = timeit.default_timer() - start_proc
                    processing += "processing time: {} sec".format(proc_time)
                    
                except:
                    processing_time = timeit.default_timer() - start_proc
                
                # handle exceptions and erros
                if ((~np.isfinite(waveform.data)).any()):
                    raise Waveform_problem(row["network"], row["station"], t)
                    
                # save waveform
                directory = "{}/{}".format(save_dir, subdir)
                if not os.path.exists(directory):
                    os.makedirs(directory)
                save = "{}/{}".format(directory, filename)
                print ("save:", save)
                waveform.write("{}/{}".format(directory, filename),
                               format='mseed')
                
                logger.info("{}::{}::{}".format(message,
                                                timeit.default_timer()-start,
                                                processing_time))
            except FDSN_error as e:
                print(e)
                logger.info("{}::{}::{}".format(message,
                                                timeit.default_timer()-start,
                                                -1))
            except Waveform_problem as e:
                print(e)
                logger.info("{}::{}::{}".format(message,
                                                timeit.default_timer()-start,
                                                -2))
            except Inventory_problem as e:
                print(e)
                logger.info("{}::{}::{}".format(message,
                                                timeit.default_timer()-start,
                                                -3))
            except Waveform_exist as e:
                print(e)
                logger.info("{}::{}::{}".format(message,
                                                timeit.default_timer()-start,
                                                -4))
            except Quality_error as e:
                print(e)
                logger.info("{}::{}::{}".format(message,
                                                timeit.default_timer()-start,
                                                -5))
            except Unknown_problem as e:
                print(e)
                logger.info("{}::{}::{}".format(message,
                                                timeit.default_timer()-start,
                                                -6))
            except:
                logger.info("{}::{}::{}".format(message,
                                                timeit.default_timer()-start,
                                                -7))
            # go to next request interval
            t += dt
            
            # wait a bit before starting the next request
            time.sleep(sleep_time)
        
# Exceptions
class Waveform_exist(Exception):
    def __init__(self, path, message="File exist"):
        self.path = path
        self.message = message
        super().__init__(self.message)
        
class Quality_error(Exception):
    def __init__(self, n, s, t, message="Quality error"):
        self.path = "{}.{}.{}".format(n,s,t)
        self.message = message
        super().__init__(self.message)
        
class FDSN_error(Exception):
    def __init__(self, n, s, t, message="FDSN error"):
        self.path = "{}.{}.{}".format(n,s,t)
        self.message = message
        super().__init__(self.message)
        
class Waveform_problem(Exception):
    def __init__(self, n, s, t, message="Waveform error"):
        self.path = "{}.{}.{}".format(n,s,t)
        self.message = message
        super().__init__(self.message)
        
class Inventory_problem(Exception):
    def __init__(self, n, s, t, message="Inventory error"):
        self.path = "{}.{}.{}".format(n,s,t)
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return "{}: {}".format(self.message, self.path)
    
class Unknown_problem(Exception):
    def __init__(self, n, s, t, message="Unknown error"):
        self.path = "{}.{}.{}".format(n,s,t)
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return "{}: {}".format(self.message, self.path)
    
# run
if __name__ == "__main__":
    run()