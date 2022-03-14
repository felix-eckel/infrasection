Download continuous waveform data from EIDA nodes.

uses the obspy get_waveform and get_station interfaces.

## Usage
- Define search parameters and save directories in the [METADATA] section in "download.cfg"
- Run the collect_metadata script. This scripts requests all available EIDA nodes and searches for all infrasound stations within the time frame specified in the configuration file. It creates an input text file.
- Specify the input file created with collect_metadata in the [DOWNLOAD] section in "download.cfg"
- Set the additional parameters in [DOWNLOAD] and [PROCESSING]
- Run continuous_waveform_data. This script iterates over all lines in the input file.

## Output
Waveforms are stored in MiniSEED format in the specified data directory ([DOWNLOAD] section in "download.cfg"). A log file is created in "./logs". The last entry in each line contains either the processing time or an error code (-1: FDSN error, -2: Waveform issue, -3: Inventory issue, -4: Data already downloaded, -5: Quality criteria not met, -6: Other error)