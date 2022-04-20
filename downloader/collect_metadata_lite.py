# -*- coding: utf-8 -*-
"""
Requests station availability and station meta data from all available
EIDA nodes globally in a given time period and writes the return to files.
Requests *DF (infrasound) channels.

Usage: 
    use the download.cfg input file to specify a start and end date as well
    as paths where to write the meta info .json file and download input file.

@author: felix.eckel@ifg.uni-kiel.de
March 2022

Based on the works of M. Timkó (Institute of Earth Physics and Space Science)
"""
import configparser
from obspy import UTCDateTime
from obspy.clients.fdsn.header import URL_MAPPINGS, FDSNException
from obspy.clients.fdsn import Client
import sys
import json

def run():
    # read config
    config = configparser.ConfigParser()
    config.read("./download_20220115_LD.cfg")
    
    start_date      = config.get("METADATA", "start_date")
    end_date        = config.get("METADATA", "end_date")
    start_date      = UTCDateTime(start_date)
    end_date        = UTCDateTime(end_date)
    json_path       = config.get("METADATA", "json_path")
    input_file_path = config.get("METADATA", "input_file_path")
    station_dictionary = {}
    
    # request ALL EIDA nodes
    EIDA_nodes = URL_MAPPINGS.keys()
    
    print ("Generating metadata...")
    for node in EIDA_nodes:
        try:
            client = Client(node)
            
            # get all globally available stations for channel *DF,
            inv = client.get_stations(network="*", channel="LD?",
                                      starttime=start_date, endtime=end_date,
                                      latitude=0, longitude=0, level="station",
                                      maxradius=180)
            
            # append to station_dictionary
            for network in inv:
                nw = network.code
                if nw not in station_dictionary:
                    station_dictionary[nw] = {}
                for station in network:
                    st = station.code
                    
                    # get station operation time
                    station_start = station.start_date
                    station_end = station.end_date
                    if start_date < station_start:
                        op_start = station_start
                    else:
                        op_start = start_date
                    if station_end == None:
                        op_end = end_date
                    elif end_date < station_end:
                        op_end = end_date
                    else:
                        op_end = station_end
                        
                    op_start = op_start.strftime("%Y-%m-%d")
                    op_end   = op_end.strftime("%Y-%m-%d")
                    
                    if (st not in station_dictionary[nw]):
                        station_dictionary[nw][st] = {
                            "latitude" : station.latitude,
                            "longitude" : station.longitude,
                            "elevation" : station.elevation,
                            "start_date" : [op_start],
                            "end_date" : [op_end],
                            "status" : station.restricted_status,
                            "services" : [node]
                        }
                    elif (node not in station_dictionary[nw][st]["services"]):
                        station_dictionary[nw][st]["services"].append(node)
                    else:
                        station_dictionary[nw][st]["start_date"].append(\
                                                                    op_start)
                        station_dictionary[nw][st]["end_date"].append(op_end)
            
            print("{} -- success".format(node))
            
        # exception handling
        except FDSNException as e:
            print ("{} -- failed: FDSNException".format(node))
            print(e)
        except KeyboardInterrupt:
            print("failed -- killed")
            sys.exit(1)
        except:
            print("{} -- failed: unknown".format(node))
   
    # write meta info to json file
    with open(json_path, 'w') as fid:
        json.dump(station_dictionary, fid, sort_keys=False, indent=2)
    
    # write download input file
    with open(input_file_path,'w') as fid:
        for network in station_dictionary:
            for station in station_dictionary[network]:
                c = len(station_dictionary[network][station]["start_date"])
                for i in range(c):
                    fid.write("%s %s %s %s %s\n" % (
                        station_dictionary[network][station]["services"][0], 
                        network, 
                        station, 
                        station_dictionary[network][station]["start_date"][i],
                        station_dictionary[network][station]["end_date"][i]
                        )  
                )
    print("Metadata collection -- finished")
    
if __name__ == "__main__":
    run()