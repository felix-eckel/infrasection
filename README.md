# infrasection
plot global infrasound sections and perform semblance analysis

## Scripts

### section.py
Plots a global section of waveforms either as raw waveform, root-mean-square amplitudes or STA/LTA characteristic function.

### semblance.py
Plots a semblance analysis for various velocities. Shifts traces by distance and velocity from a given target and calculates a semblance out of the stacked traces. An example output image can be found in the exampels directory.

### semblance_map.py
Similar to semblance.py but instead of searching for velocities from a given origin it searches for the origin itself with a given velocity. An example output image can be found in the exampels directory. An example output image can be found in the exampels directory.

### Download
This repository comes with a lite version of a continuous waveform downloader.

With collect_metadata_lite a global request of infrasound stations (within the EIDA network) is issued in a given time frame (specify start_date and end_date in download.cfg).

With continuous_waveform_downloader_lite the input file created with the metadata request is handled and the data is being downloaded and processed according to the specified parameters in download.cfg


## Requirements
These scripts require
- obspy to read data and perform data processing
- matplotlib for all plotting routines
- pandas to read json and csv files
- cartopy to plot the semblance map and its decorators
- numpy for various essential tasks

## Usage
Download data with the scripts provided in "download".

Define parameters in infrasection.cfg. Run __main__ on the scripts (section, semblance, semblance_map).

### Suggestion for parameters
for infrasound section try these parameters first:
- filter: 100 s - 14400 s
- rms length: 150 s
- resampling: 25 s
- sta/lta length: 600 s/14400 s
- threshold on/off: 6.0/1.0

for semblance analysis try these parameters first:
- filter: 400 s - 900 s
- rms length: 600 s
- resampling: 100 s

## Data
data is expected to be stored in the data directory given in infrasection.cfg - data_directory.

Files are searched as "data_directory"/*/"NW"."STAT"*.mseed with "NW" the station network code and "STAT" the station identifier. Adapt the "file" variable declaration in each script if your data is stored differently.

infrasection.cfg requires a station list. This should be a csv file (" " as delimiter) with the station network identifier in column 2 and station name in column 3. See example directory for an example.

infrasection.cfg requires a station meta file. This should be a .json file. See example directory for an example.

infrasection.cfg requires a json file with info on the event (coordinates, estimated origin time, description). See example directory fo an example.

(see "Download")

## Performance
Performance for the global section is mostly depentend on the amount of data and the time it takes to read data from the disk.

Semblance analysis is quite performance heavy. If it takes too long, one can adapt the velocity sampling in semblance.py in the declaration for "vels" or the coordinate grid sampling in semblance_map.py in the declaration for "lons" and "lats". 