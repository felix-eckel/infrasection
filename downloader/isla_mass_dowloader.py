"""
[2022-04-20 11:59:28,826] - obspy.clients.fdsn.mass_downloader -
INFO: Initializing FDSN client(s) for AUSPASS, BGR, ETH, EMSC, GEONET, GEOFON, GFZ, ICGC, IESDMC, INGV, IPGP, IRIS,
IRISPH5, ISC, KNMI, KOERI, LMU, NCEDC, NIEP, NOA, ODC, ORFEUS, RESIF, RESIFPH5, RASPISHAKE,
SCEDC, TEXNET, UIB-NORSAR, USGS, USP.
[2022-04-20 11:59:29,099] - obspy.clients.fdsn.mass_downloader -
INFO: Cannot use client 'ISC' as it does not have 'dataselect' and/or 'station' services.
[2022-04-20 11:59:30,219] - obspy.clients.fdsn.mass_downloader -
INFO: Cannot use client 'EMSC' as it does not have 'dataselect' and/or 'station' services.
[2022-04-20 11:59:34,877] - obspy.clients.fdsn.mass_downloader -
INFO: Cannot use client 'USGS' as it does not have 'dataselect' and/or 'station' services.
[2022-04-20 11:59:34,878] - obspy.clients.fdsn.mass_downloader -
INFO: Successfully initialized 27 client(s): AUSPASS, BGR, ETH, GEONET, GEOFON, GFZ, ICGC, IESDMC, INGV, IPGP, IRIS,
IRISPH5, KNMI, KOERI, LMU, NCEDC, NIEP, NOA, ODC, ORFEUS, RESIF, RESIFPH5, RASPISHAKE, SCEDC, TEXNET, UIB-NORSAR, USP.
"""
import obspy
from obspy.clients.fdsn.mass_downloader import GlobalDomain, \
    Restrictions, MassDownloader
from obspy import UTCDateTime
from typing import List

origin_time: UTCDateTime = UTCDateTime(2022, 1, 16, 0, 0, 0)
end_time: UTCDateTime = UTCDateTime(2022, 1, 17, 0, 0, 0)
# If folder does not exist, the folder is created
# Will be downloading multiple days
mseed_storage_directory: str = "/Users/mgarces/Documents/DATA_2022/Tonga/FDSN_20220115_10D/waveforms/20220116"
xml_storage_directory: str = "/Users/mgarces/Documents/DATA_2022/Tonga/FDSN_20220115_10D/stations/20220116"
# Example: [HB]DF are for channels HDF, BDF
# channels: List[str] = ["[HB]DF", "LD[IAO]", "BD[AO]"]
channels: List[str] = ["BDF", "LD[IAO]", "BD[AO]"]

if __name__ == "__main__":

    domain = GlobalDomain()

    restrictions = Restrictions(
        # Temporal bounds of the waveform data 5 minutes before and after the event.
        starttime=origin_time,
        endtime=end_time,
        # If True, any trace with a gap/overlap will be discarded.
        reject_channels_with_gaps=False,
        # Any trace that is shorter than 50 % of the desired total duration will be discarded.
        minimum_length=0.5,
        # No two stations should be closer than 10 km to each other. This is
        # useful to for example filter out stations that are part of different
        # networks but at the same physical station. Settings this option to
        # zero or None will disable that filtering.
        minimum_interstation_distance_in_m=0,
        # Only HH or BH channels. If a station has HH channels, those will be
        # downloaded, otherwise the BH. Nothing will be downloaded if it has
        # neither. You can add more/less patterns if you like.
        channel_priorities=channels,
        # Location codes are arbitrary and there is no rule as to which
        # location is best. Same logic as for the previous setting.
        # location_priorities=["", "00", "10"]
    )

    # No specified providers will result in all known ones being queried.
    mdl = MassDownloader(providers=obspy.clients.fdsn.header.URL_MAPPINGS)
    # The data will be downloaded to the ``./waveforms/`` and ``./stations/``
    # folders with automatically chosen file names.
    mdl.download(domain, restrictions, mseed_storage=mseed_storage_directory,
                 stationxml_storage=xml_storage_directory)