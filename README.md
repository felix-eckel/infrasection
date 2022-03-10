# infrasection
plot global infrasound sections and perform semblance analysis

## Scripts

### section.py
Plots a global section of waveforms either as raw waveform, root-mean-square amplitudes or STA/LTA characteristic function.

### semblance.py
Plots a semblance analysis for various velocities. Shifts traces by distance and velocity from a given target and calculates a semblance out of the stacked traces.

### semblance_map.py
Similar to semblance.py but instead of searching for velocities from a given origin it searches for the origin itself with a given velocity.

## Usage
define parameters in infrasection.cfg. Run __main__ on the scripts.

## Data
data is expected to be stored in the data directory given in infrasection.cfg - data_directory.

Files are searched as "data_directory"/*/"NW"."STAT"*.mseed with "NW" the station network code and "STAT" the station identifier. Adapt the "file" variable declaration in each script if your data is stored differently.

infrasection.cfg requires a station list. This should be a csv file (" " as delimiter) with the station network identifier in column 2 and station name in column 3. See example directory for an example.

infrasection.cfg requires a station meta file. This should be a .json file. See example directory for an example.

infrasection.cfg requires a json file with info on the event (coordinates, estimated origin time, description). See example directory fo an example.


## Requirements
These scripts require
- obspy to read data and perform data processing
- matplotlib for all plotting routines
- pandas to read json and csv files
- cartopy to plot the semblance map and its decorators
- numpy for various essential tasks